MODULE hairy_probes
 
        USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                                   get_atomic_kind,&
                                                   get_atomic_kind_set
        USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                                   gto_basis_set_type
        USE cp_control_types,                ONLY: hairy_probes_type  
        USE kinds,                           ONLY: dp
        USE orbital_pointers,                ONLY: indco,&
                                                   nco,&                                                     
                                                   nso                                                              
        USE particle_types,                  ONLY: particle_type
        USE qs_kind_types,                   ONLY: get_qs_kind,&
                                                   get_qs_kind_set,&
                                                   qs_kind_type
        USE qs_environment_types,            ONLY: get_qs_env,& 
                                                   qs_environment_type
 
        USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                                   cp_logger_get_default_io_unit,&
                                                   cp_logger_type,&
                                                   cp_to_string
        USE orbital_transformation_matrices, ONLY: orbtramat
        USE qs_mo_types,                     ONLY: mo_set_p_type,&
                                                   mo_set_type, &
                                                   get_mo_set 
        USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                                   cp_fm_get_submatrix,&
                                                   cp_fm_type
        USE cp_units,                        ONLY: cp_unit_from_cp2k
        USE kahan_sum,                       ONLY: accurate_sum

#include "./base/base_uses.f90"

        IMPLICIT NONE 
        PRIVATE
        
        CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hairy_probes'

        INTEGER, PARAMETER, PRIVATE          :: BISECT_MAX_ITER = 400
        PUBLIC :: probe_occupancy, probe_occupancy_kp, AO_boundaries

CONTAINS

!**************************************************************************************
!> \brief subroutine to calculate occupation number and 'Fermi' level using the 
!> \biref HAIR PROBE approach; gamma point calculation.
!> \param occ(:): occupation numbers
!> \param kTS: entropic energy contribution 
!> \param fermi: fermi level
!> \param probe(np): hairy probe 
!> \param atomic_kind_set(:) 
!> \param particle_set(:)
!> \param qs_kind_set(:)
!> \param coeff: MOs coefficient 
!> \param energies(:): MOs eigenvalues 
!> \param maxocc: maximum allowed occupation number of an MO (1 or 2)
!> \param N: number of electrons
!> \param nAO: number of atomic orbitals (# rows in coeff)
!> \param nMO: number of molecular orbitals (# cols in coeff)
!
!> \var smatrix(:,:): spherical MOs matrix 
!> \var nrow_global, ncol_global: # rows and columns in smatrix (~should coincide 
!>                                with nAO and nMO)
!> \var fermi_min, fermi_max, fermi_fit: temporary fermi level needed to apply  
!>                                       bisection method to find 'fermi'
!> \var N_min, N_max, N_now: temoprary number of electrons needed to apply the 
!>                           bisection method to find 'fermi'
!> \var ip, np, iter: cycle indexes 
!**************************************************************************************
SUBROUTINE probe_occupancy ( occ, fermi, kTS, energies, coeff, maxocc, probe, N )

        !i/o variables and arrays
        REAL(KIND=dp), INTENT(out)                          :: occ(:), kTS, fermi
        TYPE(hairy_probes_type), INTENT(INOUT)              :: probe(:) 
        TYPE(cp_fm_type), INTENT(IN), POINTER               :: coeff         
        REAL(KIND=dp), INTENT(IN)                           :: energies(:), maxocc, N  
                                                                                                                          
        !subroutine variables and arrays 
        REAL(KIND=dp), POINTER                              :: smatrix(:,:) !smatrix: Spherical MOs matrix 
        REAL(KIND=dp), ALLOCATABLE                          :: smatrix_squared(:,:) !squared Spherical MOs matrix 
        REAL(KIND=dp)                                       :: fermi_min, fermi_max, fermi_half, fermi_fit, &
                                                               & N_fit, N_min, N_max, N_half, delta_fermi, &
                                                               & h, y0, y1, y2, de, N_now
        INTEGER                                             :: nrow_global, ncol_global, iter 

        REAL(KIND=dp), PARAMETER                            :: epsocc = 1.0e-12_dp

        CALL cp_fm_get_info(coeff, &
                             nrow_global=nrow_global, &
                             ncol_global=ncol_global)
        ALLOCATE (smatrix(nrow_global, ncol_global))
        CALL cp_fm_get_submatrix(coeff, smatrix)

        ALLOCATE(smatrix_squared(nrow_global, ncol_global))
        smatrix_squared(:,:) = smatrix(:,:)**2.0d0 

!**************************************************************************************
!1)calculate Fermi energy using the hairy probe formula
!**************************************************************************************
         de = probe(1)%T *LOG((1.0_dp - epsocc)/epsocc) !FIX probe(1)%T
         de = MAX(de, 0.5_dp)
         fermi_max = MAXVAL(probe%mu) + de 
         fermi_min = MINVAL(probe%mu) - de 

         CALL HP_occupancy(probe=probe, matrix=smatrix_squared, energies=energies, &
                           maxocc=maxocc, fermi=fermi_max, occupancy=occ, kTS=kTS)
         N_max = accurate_sum(occ)

         CALL HP_occupancy(probe=probe, matrix=smatrix_squared, energies=energies, &
                           maxocc=maxocc, fermi=fermi_min, occupancy=occ, kTS=kTS)
         N_min = accurate_sum(occ)

         iter = 0
         DO WHILE (ABS(N_max - N_min) > N*epsocc)
            iter = iter + 1
          
            fermi_half = (fermi_max + fermi_min)/2.0_dp
            CALL HP_occupancy(probe=probe, matrix=smatrix_squared, energies=energies, &
                  maxocc=maxocc, fermi=fermi_half, occupancy=occ, kTS=kTS)
            N_half = accurate_sum(occ)

            h = fermi_half - fermi_min
            IF ( h .gt. N*epsocc*100 ) THEN 
               y0 = N_min - N
               y1 = N_half - N
               y2 = N_max - N

               CALL three_point_zero(y0,y1,y2,h,delta_fermi)
               fermi_fit = fermi_min + delta_fermi 
         
               CALL HP_occupancy(probe=probe, matrix=smatrix_squared, energies=energies, &
                              maxocc=maxocc, fermi=fermi_fit, occupancy=occ, kTS=kTS)
               N_fit = accurate_sum(occ) 
            END IF
            
            !define 1st bracked using fermi_half
            IF (N_half < N) THEN
               fermi_min = fermi_half
               N_min = N_half
            ELSE IF ( N_half > N ) THEN 
               fermi_max = fermi_half
               N_max = N_half
            ELSE 
               fermi_min = fermi_half
               N_min = N_half
               fermi_max = fermi_half
               N_max = N_half
               h = 0.0d0
            ENDIF
  
            !define 2nd bracker using fermi_fit
            IF ( h .gt. N*epsocc*100 ) THEN 
               IF ( fermi_fit .ge. fermi_min .and. fermi_fit .le. fermi_max ) THEN 
                  IF (N_fit < N) THEN
                     fermi_min = fermi_fit
                     N_min = N_fit
                  ELSE IF ( N_fit > N ) THEN 
                     fermi_max = fermi_fit
                     N_max = N_fit
                  ENDIF   
               END IF 
            END IF 

            IF ( ABS(N_max - N) < N*epsocc ) THEN 
                fermi = fermi_max
                EXIT 
            ELSE IF (  ABS(N_min - N ) < N*epsocc ) THEN 
                fermi = fermi_min 
                EXIT 
            END IF  

            IF (iter > BISECT_MAX_ITER) THEN
               CPWARN("Maximum number of iterations reached while finding the Fermi energy")
               EXIT
            ENDIF
         ENDDO

!**************************************************************************************
!2)calculate occupation numbers according to hairy probe formula 
!**************************************************************************************         
         occ(:) = 0.0_dp
         N_now = 0.0_dp
         CALL HP_occupancy(probe=probe, matrix=smatrix_squared, energies=energies, &
                           maxocc=maxocc, fermi=fermi, occupancy=occ, kTS=kTS)
         N_now = accurate_sum(occ)

         IF (ABS(N_now - N) > N*epsocc) CPWARN("Total number of electrons is not accurate - HP")

         DEALLOCATE(smatrix,smatrix_squared)

END SUBROUTINE probe_occupancy 

!**************************************************************************************
!> \brief subroutine to calculate occupation number and 'Fermi' level using the 
!> \biref HAIR PROBE approach; kpoints calculation.
!> \param occ(nMO,nkp): occupation numbers
!> \param kTS: entropic energy contribution 
!> \param fermi: fermi level
!> \param probe(np): hairy probe 
!> \param atomic_kind_set(:) 
!> \param particle_set(:)
!> \param qs_kind_set(:)
!> \param coeff(nkp): MOs coefficient 
!> \param energies(nMO,nkp): MOs eigenvalues 
!> \param wk(nkp): weight of kpoints
!> \param maxocc: maximum allowed occupation number of an MO (1 or 2)
!> \param N: number of electrons
!> \param nAO: number of atomic orbitals (# rows in coeff)
!> \param nMO: number of molecular orbitals (# cols in coeff)
!
!> \var coeff_kp: pointer to elements of coeff(ikp) 
!> \var smatrix(:,:): spherical MOs matrix 
!> \var kTS_kp: entropic contribution per kpoint
!> \var nrow_global, ncol_global: # rows and columns in smatrix (~should coincide 
!>                                with nAO and nMO)
!> \var fermi_min, fermi_max, fermi_now: temporary fermi level needed to apply  
!>                                       bisection method to find 'fermi'
!> \var N_now: temoprary number of electrons needed to apply the bisection method
!>             to find 'fermi'
!> \var ip, np, ikp, nkpo, iter: cycle indexes 
!**************************************************************************************
SUBROUTINE probe_occupancy_kp ( occ, fermi, kTS, energies, rcoeff, icoeff, maxocc, probe, N, wk )

      CHARACTER(LEN=*), PARAMETER :: routineN = 'probe_occupancy_kp'

      REAL(KIND=dp),INTENT(OUT)                          :: occ(:,:,:), fermi, kTS 
      TYPE(hairy_probes_type), INTENT(IN)                :: probe(:) 
      REAL(KIND=dp), INTENT(IN)                          :: energies(:,:,:), wk(:)
      REAL(KIND=dp), INTENT(IN)                          :: rcoeff(:,:,:,:), icoeff(:,:,:,:)
      REAL(KIND=dp), INTENT(IN)                          :: maxocc, N

      REAL(KIND=dp)                                      :: kTS_kp, de, h, fermi_max, fermi_min, fermi_half, fermi_fit, &
                                                            & N_max, N_min, N_half, N_fit, N_now, &
                                                            & y0, y1, y2, delta_fermi
      REAL(KIND=dp), ALLOCATABLE                         :: coeff_squared(:,:,:,:)
      INTEGER                                            :: nMO, nAO, iter, ikp, nkp, ispin, nspin 

      REAL(KIND=dp), PARAMETER                           :: epsocc = 1.0e-12_dp

      INTEGER :: handle 

      CALL timeset(routineN, handle)

!**************************************************************************************
!1)calculate Fermi energy using the hairy probe formula
!**************************************************************************************
      nAO=SIZE(rcoeff, 1)
      nMO=SIZE(rcoeff, 2)
      nkp=SIZE(rcoeff, 3)
      nspin=SIZE(rcoeff, 4)

      ALLOCATE(coeff_squared(nAO,nMO,nkp,nspin))
      coeff_squared(:,:,:,:) = rcoeff(:,:,:,:)**2.0d0 + icoeff(:,:,:,:)**2.0d0        

      occ(:,:,:)=0.0_dp
  
      !define initial brackets 
      de = probe(1)%T *LOG((1.0_dp - epsocc)/epsocc) !FIX probe(1)%T
      de = MAX(de, 0.5_dp)
      fermi_max = MAXVAL(probe%mu) + de 
      fermi_min = MINVAL(probe%mu) - de 

      N_max = 0.0_dp
      kTS = 0.0_dp
      !***HP 
      DO ispin = 1, nspin
         DO ikp = 1, nkp
            CALL HP_occupancy(probe,coeff_squared(:,:,ikp,ispin),energies(:,ikp,ispin), &
                              maxocc,fermi_max,occ(:,ikp,ispin),kTS_kp)

            kTS = kTS + kTS_kp * wk(ikp) !entropic contribution
            N_max = N_max + accurate_sum(occ(1:nMo,ikp,ispin))*wk(ikp) 
         END DO
      END DO 
      !***HP 

      N_min = 0.0_dp
      kTS = 0.0_dp
      !***HP 
      DO ispin = 1, nspin
         DO ikp = 1, nkp
            CALL HP_occupancy(probe,coeff_squared(:,:,ikp,ispin),energies(:,ikp,ispin), &
                              maxocc,fermi_min,occ(:,ikp,ispin),kTS_kp)

            kTS = kTS + kTS_kp * wk(ikp) !entropic contribution
            N_min = N_min + accurate_sum(occ(1:nMo,ikp,ispin))*wk(ikp) 
         END DO
      END DO 
      !***HP 

      iter = 0
      DO WHILE (ABS(N_max - N_min) > N*epsocc)
         iter = iter + 1
         fermi_half = (fermi_max + fermi_min)/2.0_dp
         N_half = 0.0_dp
         kTS = 0.0_dp
         !***HP 
         DO ispin = 1, nspin
            DO ikp = 1, nkp
               CALL HP_occupancy(probe,coeff_squared(:,:,ikp,ispin),energies(:,ikp,ispin), &
                                 maxocc,fermi_half,occ(:,ikp,ispin),kTS_kp)

               kTS = kTS + kTS_kp * wk(ikp) !entropic contribution
               N_half = N_half + accurate_sum(occ(1:nMo,ikp,ispin))*wk(ikp) 
            END DO
         END DO 
         !***HP 

         h = fermi_half - fermi_min
         IF ( h .gt. N*epsocc*100 ) THEN 
            y0 = N_min - N
            y1 = N_half - N
            y2 = N_max - N

            CALL three_point_zero(y0,y1,y2,h,delta_fermi)
            fermi_fit = fermi_min + delta_fermi 
            N_fit = 0.0_dp
            kTS = 0.0_dp
                 
            !***HP
            DO ispin = 1, nspin
               DO ikp = 1, nkp
                  CALL HP_occupancy(probe,coeff_squared(:,:,ikp,ispin),energies(:,ikp,ispin), &
                                 maxocc,fermi_fit,occ(:,ikp,ispin),kTS_kp)

                  kTS = kTS + kTS_kp * wk(ikp) !entropic contribution
                  N_fit = N_fit + accurate_sum(occ(1:nMo,ikp,ispin))*wk(ikp) 
               END DO
            END DO 
            !***HP

         END IF

         !define 1st bracked using fermi_half
         IF (N_half < N) THEN
            fermi_min = fermi_half
            N_min = N_half
         ELSE IF ( N_half > N ) THEN 
            fermi_max = fermi_half
            N_max = N_half
         ELSE 
            fermi_min = fermi_half
            N_min = N_half
            fermi_max = fermi_half
            N_max = N_half
            h = 0.0d0
         ENDIF
  
         !define 2nd bracker using fermi_fit
         IF ( h .gt. N*epsocc*100 ) THEN 
            IF ( fermi_fit .ge. fermi_min .and. fermi_fit .le. fermi_max ) THEN 
               IF (N_fit < N) THEN
                  fermi_min = fermi_fit
                  N_min = N_fit
               ELSE IF ( N_fit > N ) THEN 
                  fermi_max = fermi_fit
                  N_max = N_fit
               ENDIF   
            END IF 
         END IF 

         IF ( ABS(N_max - N) < N*epsocc ) THEN 
             fermi = fermi_max
             EXIT 
         ELSE IF (  ABS(N_min - N ) < N*epsocc ) THEN 
             fermi = fermi_min 
             EXIT 
         END IF  

         IF (iter > BISECT_MAX_ITER) THEN
            CPWARN("Maximum number of iterations reached while finding the Fermi energy")
            EXIT
         ENDIF
      ENDDO

!**************************************************************************************
!3)calculate occupation numbers using the hairy probe formula
!**************************************************************************************
      N_now = 0.0_dp
      kTS = 0.0_dp

      !***HP 
      DO ispin = 1, nspin
         DO ikp = 1, nkp
            CALL HP_occupancy(probe,coeff_squared(:,:,ikp,ispin),energies(:,ikp,ispin), &
                              maxocc,fermi,occ(:,ikp,ispin),kTS_kp)

            kTS = kTS + kTS_kp * wk(ikp) !entropic contribution
            N_now = N_now + accurate_sum(occ(1:nMo,ikp,ispin))*wk(ikp) 
         END DO
      END DO 
      !***HP 
      
      DEALLOCATE(coeff_squared)

      CALL timestop(handle)
END SUBROUTINE probe_occupancy_kp 

!**************************************************************************************
!
!
!
!**************************************************************************************
SUBROUTINE HP_occupancy(probe, matrix, energies, maxocc, fermi, occupancy, kTS)
        IMPLICIT NONE 

        REAL(KIND=dp), INTENT(OUT)                 :: occupancy(:), kTS
        TYPE(hairy_probes_type), INTENT(IN)         :: probe(:)
        REAL(KIND=dp), INTENT(IN)                  :: matrix(:,:) !squared coefficient matrix
        REAL(KIND=dp), INTENT(IN)                  :: energies(:), fermi, maxocc

        INTEGER :: iMO, ip, np, nMO!, nAO
        REAL(KIND=dp) :: mu, fermi_fun, fermi_fun_sol, sum_fermi_fun, sum_coeff, &
                         & sum_fermi_fun_sol, sum_coeff_sol, C, s, f

        nMO = SIZE(matrix, 2)
        np = SIZE(probe)

        kTS = 0.0_dp
        ! kTS is the entropic contribution to the electronic energy

        mos2: do iMO = 1, nMO

           sum_fermi_fun = 0.0_dp
           sum_coeff = 0.0_dp 

           sum_fermi_fun_sol = 0.0_dp
           sum_coeff_sol = 0.0_dp 

           probes2: DO ip = 1, np
              IF ( probe(ip)%alpha .lt. 1.0_dp ) THEN 
                 !fermi distribution, solution probes 
                 CALL fermi_distribution ( energies(iMO), fermi, probe(ip)%T, fermi_fun_sol )                        
                 !sum of coefficients, solution probes 
                 sum_coeff_sol = sum_coeff_sol + SUM(matrix(probe(ip)%first_ao:probe(ip)%last_ao,iMO))

              ELSE  
                 C = SUM(matrix(probe(ip)%first_ao:probe(ip)%last_ao,iMO))  
                 !bias probes 
                 mu = fermi - probe(ip)%mu 
                 !fermi distribution, main probes                   
                 CALL fermi_distribution ( energies(iMO), mu, probe(ip)%T, fermi_fun )                        
                 !sum fermi distribution * coefficients
                 sum_fermi_fun = sum_fermi_fun + ( fermi_fun * C )
                 !sum cofficients 
                 sum_coeff = sum_coeff + C 
              END IF 
           END DO probes2

           sum_fermi_fun_sol = probe(np)%alpha * fermi_fun_sol * sum_coeff_sol
           sum_coeff_sol = probe(np)%alpha * sum_coeff_sol                       
           f = ( sum_fermi_fun_sol + sum_fermi_fun ) / ( sum_coeff_sol + sum_coeff ) 
           occupancy(iMO) = f * maxocc

           !entropy kTS= kT*[f ln f + (1-f) ln (1-f)]
           IF ( f .eq. 0.0d0 .or. f .eq. 1.0d0 ) THEN
              s = 0.0d0
           ELSE 
              s = f * LOG(f) + (1.0d0-f) * LOG(1.0d0-f) 
           END IF 
           kTS = kTS + probe(np)%T*maxocc*s 
        END DO MOs2

END SUBROUTINE HP_occupancy

!**************************************************************************************
!
!
!
!**************************************************************************************
SUBROUTINE AO_boundaries (probe, atomic_kind_set, qs_kind_set, particle_set, nAO)
        IMPLICIT NONE 
        
        TYPE(hairy_probes_type), INTENT(INOUT)           :: probe
        TYPE(atomic_kind_type), INTENT(IN), POINTER      :: atomic_kind_set(:) !from get_atomic_kind_set 
        TYPE(particle_type), INTENT(IN), POINTER         :: particle_set(:)    !from get_atomic_kind
        TYPE(qs_kind_type), INTENT(IN), POINTER          :: qs_kind_set(:)     !from get_qs_kind_set
        INTEGER, INTENT(IN) :: nAO
 
        !subroutine variables and arrays 
        TYPE(gto_basis_set_type), POINTER                :: orb_basis_set   !from get_qs_kind 
        INTEGER, DIMENSION(:), POINTER                   :: nshell          !from get_gto_basis_set 
        INTEGER, DIMENSION(:, :), POINTER                :: l               !from get_gto_basis_set 
        INTEGER                                          :: nset            !from get_gto_basis_set 
                                                           !nset : Number of exponent sets
                                                           !nshell(l) : Number of shells for angular momentum quantum number l
        INTEGER                                          :: lshell      
        INTEGER                                          :: natom, nsgf

        INTEGER :: p_atom, isgf, iatom, ii, ikind, iset, ishell 

        !get sets from different modules: 
        CALL get_atomic_kind_set ( atomic_kind_set=atomic_kind_set, natom=natom )
        CALL get_qs_kind_set ( qs_kind_set=qs_kind_set, nsgf=nsgf ) !indexes for orbital symbols 

        !get iAO boundaries:
        p_atom = SIZE ( probe%atom_ids )
        probe%first_ao = nsgf !nAO
        probe%last_ao = 1
        isgf = 0
        
        atoms: DO iatom = 1, natom 
                NULLIFY ( orb_basis_set ) 
                !1) iatom is used to find the correct atom kind and its index (ikind) is the particle set 
                CALL get_atomic_kind ( particle_set(iatom)%atomic_kind, kind_number=ikind )        

                !2) ikind is used to find the basis set associate to that atomic kind
                CALL get_qs_kind ( qs_kind_set(ikind), basis_set=orb_basis_set )

                !3) orb_basis_set is used to get the gto basis set variables 
                IF ( ASSOCIATED ( orb_basis_set ) ) THEN 
                        CALL get_gto_basis_set ( gto_basis_set=orb_basis_set, &
                                nset=nset, nshell=nshell, l=l ) !? ,cgf_symbol=bcgf_symbol )
                END IF
                 
                !4) get iAO boundaries                 
                sets: DO iset = 1, nset 

                        shells: DO ishell = 1, nshell(iset) 
                                lshell = l ( ishell, iset ) 

                                isgf = isgf + nso(lshell)

                                        boundaries: DO ii = 1, p_atom 
                                                IF ( iatom .ne. probe%atom_ids(ii) ) THEN                                          
                                                        CYCLE boundaries
                                                ELSE 
                                                        !defines iAO boundaries***********
                                                        probe%first_ao = MIN ( probe%first_ao, isgf )
                                                        probe%last_ao = MAX ( probe%last_ao, isgf )
                                                END IF
                                        END DO boundaries 
                                        
                        END DO shells
                END DO sets
        END DO atoms                            

        IF ( isgf .ne. nAO ) CPWARN("row count does not correspond to nAO, number of rows in mo_coeff")

 END SUBROUTINE AO_boundaries 

!*******************************************************************************************************
!*******************************************************************************************************

SUBROUTINE fermi_distribution (E, pot, temp, f_p)
      IMPLICIT NONE 

      REAL(kind=dp),INTENT(OUT) :: f_p
      REAL(kind=dp),INTENT(IN) :: E
      REAL(kind=dp),INTENT(IN) :: pot, temp

      REAL(KIND=dp)            :: arg, exponential, exponential_plus_1, f, one_minus_f

      ! have the result of exp go to zero instead of overflowing  
      IF (E > pot) THEN                                           
         arg = -(E - pot)/temp                                    
         ! exponential is smaller than 1                                  
         exponential = EXP(arg)                                           
         exponential_plus_1 = exponential + 1.0_dp                                      

         one_minus_f = exponential/exponential_plus_1                                                      
         f = 1.0_dp/exponential_plus_1                                       
         f_p = one_minus_f
      ELSE                                                        
         arg = (E - pot)/temp                                     
         ! exponential is smaller than 1                                  
         exponential = EXP(arg)                                           
         exponential_plus_1 = exponential + 1.0_dp                                      

         f = 1.0_dp/exponential_plus_1                                                   
         one_minus_f = exponential/exponential_plus_1                                          
         f_p = f 
      END IF                                                      
                                       
END SUBROUTINE fermi_distribution                                        

!*******************************************************************************************************
!*******************************************************************************************************
SUBROUTINE three_point_zero(y0,y1,y2,h,delta)
   IMPLICIT NONE 

   REAL(kind=dp), INTENT(IN) :: y0, y1, y2, h
   REAL(kind=dp), INTENT(OUT) :: delta 

   REAL(kind=dp) :: a, b, c, d 
   REAL(kind=dp) :: u0

   a = ( y2 - 2.0_dp*y1 + y0 ) / ( 2.0_dp * h * h )
   b = ( 4.0_dp*y1 - 3.0_dp*y0 - y2 ) / ( 2.0_dp * h )
   c = y0

   d = sqrt( b * b - 4.0_dp*a*c )
   
   u0 = ( -b + d ) / ( 2.0_dp * a )

   IF ( u0 .ge. 0.0_dp .and. u0 .le. 2.0_dp*h ) THEN
      delta = u0 
      !RETURN
   ELSE 
      u0 = ( -b - d ) / ( 2.0_dp * a )
      IF ( u0 .ge. 0.0_dp .and. u0 .le. 2.0_dp*h ) THEN
         delta = u0
         !RETURN
      ELSE
         IF ( y1 .lt. 0.0_dp ) delta = 2.0_dp*h
         IF ( y1 .ge. 0.0_dp ) delta = 0.0_dp
      END IF
   END IF 

   
END SUBROUTINE three_point_zero


END MODULE hairy_probes 
