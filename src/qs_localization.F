!   CP2K: A general program to perform molecular dynamics simulations         !

!!****** cp2k/qs_localization [1.0] *
!!
!!   NAME
!!     qs_localization
!!
!!   FUNCTION
!! 
!!     should be able to find properties related to the position operator in PBC
!!     e.g. location and spread of periodic orbitals
!!          Wannier functions
!!
!!   AUTHORS
!!     Joost VandeVondele with some help of Leonardo Bernasconi (01.2003)
!!     CJM (04.2003)
!!     CJM & IKUO (08.2003)
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE qs_localization
! *****************************************************************************
  USE coefficient_types,               ONLY: coeff_type
  USE cp_cfm_basic_linalg,             ONLY: cp_cfm_lu_decompose
  USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                             cp_cfm_release,&
                                             cp_cfm_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_element,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_log_handling,                 ONLY: cp_logger_get_default_unit_nr
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE global_types,                    ONLY: WANNIER_CRAZY,&
                                             WANNIER_DIRECT,&
                                             WANNIER_JACOBI
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: twopi
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_max,&
                                             mp_minloc
  USE molecule_kind_types,             ONLY: atom_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: angstrom
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_give_back_coeff,&
                                             pw_pool_init_coeff,&
                                             pw_pool_type
  USE pw_types,                        ONLY: REALDATA3D,&
                                             REALSPACE
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_integrate_potential,          ONLY: integrate_v_rspace
  USE qs_localization_methods,         ONLY: crazy_rotations,&
                                             direct_mini,&
                                             initialize_weights,&
                                             jacobi_rotations
  USE simulation_cell,                 ONLY: cell_type,&
                                             pbc
  USE sparse_matrix_types,             ONLY: deallocate_matrix,&
                                             real_matrix_p_type,&
                                             replicate_matrix_structure,&
                                             set_matrix
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

! *** Public ***
  PUBLIC :: qs_loc_init
  PUBLIC :: qs_loc_env_type, PIPEK_MENZY, POSITION_OPERATOR
  PUBLIC :: qs_loc_compute_properties
  PUBLIC :: qs_loc_destroy

! *** Localization Data Type ***
  TYPE qs_loc_env_type
       INTEGER :: operator_type
       INTEGER :: localization_method
       LOGICAL :: wannier
       LOGICAL :: molecular_states
       INTEGER :: is_initialized
       TYPE ( real_matrix_p_type ),POINTER :: OPERATOR(:,:)
       TYPE ( cp_para_env_type ), POINTER :: para_env
       TYPE ( cell_type ), POINTER :: cell
       TYPE ( particle_type ), POINTER :: particle_set ( : )
       TYPE(distribution_1d_type), POINTER         :: local_molecules

       REAL (KIND = dp) :: weights ( 3 )
  END TYPE qs_loc_env_type
  
  INTEGER, PARAMETER :: PIPEK_MENZY=301, POSITION_OPERATOR=302

! *****************************************************************************

CONTAINS

!******************************************************************************
!
! allocates the data, and initializes the operators
!
!******************************************************************************
SUBROUTINE qs_loc_init(qs_loc_env,qs_env)

    TYPE(qs_loc_env_type), POINTER           :: qs_loc_env
    TYPE(qs_environment_type)                :: qs_env

    INTEGER                                  :: handle, i, isos, j, natoms
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s

  CALL timeset("qs_loc_init","I",'',handle)
  NULLIFY(qs_loc_env%para_env)
  NULLIFY(qs_loc_env%cell)
  NULLIFY(qs_loc_env%particle_set)
  !
  ! start replicating the s structure, we do have
  ! 6 operators (x,y,z)*(real,imag) for POSITION_OPERATOR
  ! and we have NATOMS operators for PIPEK_MENZY
  ! notice that the operator should have the same structure as S
  ! borrow para env from S
  !
  CALL get_qs_env(qs_env=qs_env,matrix_s=matrix_s,cell=qs_loc_env%cell, &
                   local_molecules=qs_loc_env%local_molecules,&
                   particle_set = qs_loc_env % particle_set, &
                   para_env=qs_loc_env%para_env)
  CALL cp_para_env_retain(qs_loc_env%para_env) 

! dbg
  qs_loc_env % operator_type = POSITION_OPERATOR
! dbg
  IF ( qs_loc_env % operator_type == POSITION_OPERATOR )  THEN
    ALLOCATE ( qs_loc_env % operator ( 2, 3 ), stat = isos )
    IF ( isos /= 0 ) CALL stop_memory ( 'qs_localization', 'operator' )
  ELSEIF ( qs_loc_env % operator_type == PIPEK_MENZY ) THEN
    natoms = SIZE ( qs_loc_env % particle_set )
    ALLOCATE ( qs_loc_env % operator ( 1, natoms ), stat = isos )
    IF ( isos /= 0 ) CALL stop_memory ( 'qs_localization', 'operator' )
  ENDIF

  DO i=1, SIZE ( qs_loc_env % operator, 2 )
     DO j=1, SIZE ( qs_loc_env % operator, 1 )
       NULLIFY(qs_loc_env%operator(j,i)%matrix)
       CALL replicate_matrix_structure(matrix_s(1)%matrix, &
            qs_loc_env%operator(j,i)%matrix,"operator")
       CALL set_matrix(qs_loc_env%operator(j,i)%matrix,0.0_dp)
     ENDDO
  END DO


  IF ( qs_loc_env % operator_type == POSITION_OPERATOR ) THEN 
    CALL initialize_weights ( qs_loc_env % cell, qs_loc_env % weights ) 
    CALL get_boys_foster_operator ( qs_loc_env, qs_env )
  ENDIF

  IF ( qs_loc_env % operator_type == PIPEK_MENZY ) & 
  CALL get_pipek_menzy_operator ( qs_loc_env, qs_env )
  
  qs_loc_env%molecular_states = .FALSE.

! done
  qs_loc_env%is_initialized = 777
  CALL timestop(0.0_dp,handle) 
END SUBROUTINE qs_loc_init

!******************************************************************************
!
! Computes the Boys-Foster operator for periodic systems
!
!******************************************************************************
SUBROUTINE get_boys_foster_operator ( qs_loc_env, qs_env ) 

! globals
    TYPE(qs_loc_env_type), POINTER           :: qs_loc_env
    TYPE(qs_environment_type)                :: qs_env

    INTEGER                                  :: i, ix, iy, iz, j
    REAL(KIND=dp)                            :: thearg, theval, x, y, z
    TYPE(coeff_type)                         :: op_rspace
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool

  NULLIFY ( pw_env, auxbas_pw_pool )
!
! OK, now we construct these operators the easy way, put them on the pw grid and integrate
! analytical construction should be faster.
!
! get a pw grid
  CALL get_qs_env(qs_env=qs_env,pw_env=pw_env)
  CALL pw_env_get(pw_env=pw_env,auxbas_pw_pool=auxbas_pw_pool)
  CALL pw_pool_init_coeff(auxbas_pw_pool,op_rspace,&
                          use_data=REALDATA3D,in_space=REALSPACE)

  DO i=1,3
   DO j=1,2
      ! CALL coeff_zero(op_rspace) 
      ! the definition of the position on the PW grid has to be 
      ! compatible with the one used in the integration procedure / collocation procedure
      ! integrate_v_rspace calculate_rho_elec
      DO iz=op_rspace%pw%pw_grid%bounds_local(1,3),op_rspace%pw%pw_grid%bounds_local(2,3)
      DO iy=op_rspace%pw%pw_grid%bounds_local(1,2),op_rspace%pw%pw_grid%bounds_local(2,2)
      DO ix=op_rspace%pw%pw_grid%bounds_local(1,1),op_rspace%pw%pw_grid%bounds_local(2,1)
         x=REAL(ix-op_rspace%pw%pw_grid%bounds(1,1),dp)*twopi/REAL(op_rspace%pw%pw_grid%npts(1),dp)
         y=REAL(iy-op_rspace%pw%pw_grid%bounds(1,2),dp)*twopi/REAL(op_rspace%pw%pw_grid%npts(2),dp)
         z=REAL(iz-op_rspace%pw%pw_grid%bounds(1,3),dp)*twopi/REAL(op_rspace%pw%pw_grid%npts(3),dp) 
         SELECT CASE (i)
         CASE (1)
            thearg=x
         CASE (2)
            thearg=y
         CASE (3)
            thearg=z
         END SELECT
         SELECT CASE (j)
         CASE (1)
            theval=COS(thearg)
         CASE (2)
            theval=SIN(thearg)
         END SELECT
         ! notice the volume element, slightly ugly I think
         op_rspace%pw%cr3d(ix,iy,iz)=theval*op_rspace%pw%pw_grid%dvol 
      ENDDO
      ENDDO
      ENDDO
! do the actual integration
      CALL set_matrix(qs_loc_env%operator(j,i)%matrix,0.0_dp)
      CALL integrate_v_rspace(op_rspace,h=qs_loc_env%operator(j,i),qs_env=qs_env,calculate_forces=.FALSE.)
   ENDDO
  ENDDO

  CALL pw_pool_give_back_coeff ( auxbas_pw_pool, op_rspace )

END SUBROUTINE get_boys_foster_operator

!******************************************************************************
!
! Computes the Pipek-Menzy operator for periodic systems
!
!******************************************************************************
SUBROUTINE get_pipek_menzy_operator ( qs_loc_env, qs_env )
! globals
    TYPE(qs_loc_env_type), POINTER           :: qs_loc_env
    TYPE(qs_environment_type)                :: qs_env

  CALL stop_program("qs_localization","PIPEK MENZY not implemented")
END SUBROUTINE get_pipek_menzy_operator

!******************************************************************************
!
! applies the operators to obtain the necessary data
!
!******************************************************************************
SUBROUTINE qs_loc_compute_properties(qs_loc_env, molecule_set, vectors, center, ionode, iunit)

    TYPE(qs_loc_env_type), POINTER           :: qs_loc_env
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(cp_fm_type), POINTER                :: vectors
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: center
    LOGICAL                                  :: ionode
    INTEGER                                  :: iunit

    INTEGER                                  :: handle, i, isos, j, k, n, &
                                                natoms, sweeps, scr
    TYPE(cp_fm_p_type), ALLOCATABLE          :: ZIJ( :, : )
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: opvec

  CALL timeset("qs_loc_compute_properties","I",'',handle)
  NULLIFY(tmp_fm_struct,opvec)

  ! get rows and cols of the input
  CALL cp_fm_get_info(vectors,nrow_global=n,ncol_global=k)
  ! replicate the input kind of matrix
  CALL cp_fm_create(opvec,vectors%matrix_struct)
  ! create a new matrix structure k x k in size, borrowing most of info
  CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=k, &
            ncol_global=k,para_env=vectors%matrix_struct%para_env, &
            context=vectors%matrix_struct%context)

  IF (ANY(SHAPE(center).NE.(/k,6/))) CALL stop_program('centers_and_spreads','center')

! Allocating ZIJ
  IF ( qs_loc_env % operator_type == POSITION_OPERATOR ) THEN
    ALLOCATE ( ZIJ ( 2, 3 ), stat = isos )
    IF ( isos /= 0 ) CALL stop_memory ( 'qs_localization', 'ZIJ' )
  ELSEIF ( qs_loc_env % operator_type == PIPEK_MENZY ) THEN
    natoms = SIZE ( qs_loc_env % operator, 2 ) 
    ALLOCATE ( ZIJ ( 1, natoms ), stat = isos )
    IF ( isos /= 0 ) CALL stop_memory ( 'qs_localization', 'ZIJ' )
  ENDIF
! make the actual matrices
  DO i=1,SIZE ( ZIJ, 2 )
    DO j=1, SIZE ( ZIJ, 1 )
        NULLIFY ( ZIJ ( j, i ) % matrix )
        CALL cp_fm_create ( ZIJ ( j, i ) % matrix, tmp_fm_struct )
    ENDDO
  ENDDO

! release one instance of the fm_struct
  CALL cp_fm_struct_release ( tmp_fm_struct ) 

! Compute zij here
  DO i=1,SIZE ( ZIJ, 2 )
    DO j=1,SIZE ( ZIJ, 1 )
       CALL cp_sm_fm_multiply(qs_loc_env%operator(j,i)%matrix, &
                            vectors,opvec,ncol=k,para_env=qs_loc_env%para_env)
       CALL cp_fm_gemm("T","N",k,k,n,1.0_dp,vectors,opvec,0.0_dp,ZIJ(j,i)%matrix)
    ENDDO
  ENDDO

  IF ( qs_loc_env % wannier ) THEN
       SELECT CASE (qs_loc_env % localization_method) 
       CASE ( WANNIER_JACOBI )
          CALL jacobi_rotations ( qs_loc_env % weights, zij, vectors, sweeps ) 
       CASE ( WANNIER_CRAZY )
          CALL crazy_rotations  ( qs_loc_env % weights, zij, vectors, sweeps ) 
       CASE ( WANNIER_DIRECT )
          CALL direct_mini      ( qs_loc_env % weights, zij, vectors, sweeps ) 
       CASE DEFAULT
          CALL stop_program("qs_localization","unknown localization method")
       END SELECT
       IF (ionode) THEN
           scr=cp_logger_get_default_unit_nr()
           WRITE(scr,'(A,I6,A)') " Wannier localization converged in ",sweeps," iterations"
       ENDIF
  ELSE
       IF (ionode) THEN
           scr=cp_logger_get_default_unit_nr()
           WRITE(scr,'(A,I6,A)') " No Wannier localization applied "
       ENDIF
  ENDIF

  IF ( qs_loc_env % operator_type == POSITION_OPERATOR )  THEN
    CALL centers_and_spreads ( ionode, iunit, qs_loc_env, zij, center ) 
    IF (qs_loc_env % molecular_states ) &
    CALL wfc_to_molecule ( qs_loc_env, center, molecule_set )
  ENDIF 


! release results
  DO i=1,SIZE ( ZIJ, 2 )
    DO j=1, SIZE ( ZIJ, 1 )
       CALL cp_fm_release(ZIJ(j,i)%matrix)
    ENDDO
  ENDDO
 
  DEALLOCATE(zij, STAT = isos)
  IF(isos/=0) CALL stop_memory('qs_loc_compute_properties','zij')
  CALL cp_fm_release(opvec)

  CALL timestop(0.0_dp,handle) 

END SUBROUTINE qs_loc_compute_properties

!******************************************************************************
!
! computes spread and centers
!
!******************************************************************************

SUBROUTINE centers_and_spreads ( ionode, iunit, qs_loc_env, zij, center )
! globals
    LOGICAL, INTENT(IN)                      :: ionode
    INTEGER, INTENT(IN)                      :: iunit
    TYPE(qs_loc_env_type), POINTER           :: qs_loc_env
    TYPE(cp_fm_p_type), INTENT(INOUT)        :: zij( :, : )
    REAL(KIND=dp), INTENT(OUT)               :: center( :, : )

    COMPLEX(KIND=dp)                         :: almost_determinant, z
    INTEGER                                  :: i, l, nstates, scr
    REAL(KIND=dp)                            :: imagpart, realpart
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_cfm_type), POINTER               :: complexresult

! get nstates

  CALL cp_fm_get_info ( zij ( 1, 1 ) % matrix, nrow_global = nstates )

! initialize complexresult for boys-foster
  NULLIFY ( complexresult )
  CALL cp_cfm_create ( complexresult, zij ( 1, 1 ) % matrix % matrix_struct )
!
  
! OK print out the centers ( center ( 1:3) )  and spreads ( center ( 4:6 ) ) of these friends
  cell => qs_loc_env % cell
  DO l=1,nstates
    DO i=1,3
       CALL cp_fm_get_element(zij(1,i)%matrix,l,l,realpart)
       CALL cp_fm_get_element(zij(2,i)%matrix,l,l,imagpart)
       z=CMPLX(realpart,imagpart,dp)
       center (l ,i )=( cell % hmat ( i, i ) / twopi ) * AIMAG ( LOG ( z ) )
       center (l, i+3 )=SQRT( -qs_loc_env % weights ( i ) * 2.0_dp * LOG ( ABS ( Z ) ) ) / twopi
    ENDDO
    IF (ionode) WRITE(iunit,'(A,6F16.8)') "He", angstrom*center(l, 1 : 6 )
  ENDDO
  
  ! now, find the spread of the orbitals
  DO i=1,3
     complexresult%local_data ( :, : ) = CMPLX (ZIJ(1,i)%matrix%local_data(:,:), &
                                         ZIJ(2,i)%matrix%local_data(:,:),dp)
     CALL cp_cfm_lu_decompose ( complexresult, almost_determinant )
     IF (ionode) THEN
         scr=cp_logger_get_default_unit_nr()
         WRITE(scr,'(A,I6,F16.10)') " D Eq. 11 ",i,SQRT(-2.0*nstates*LOG(ABS(almost_determinant)))
     ENDIF
  ENDDO

  CALL cp_cfm_release(complexresult)

END SUBROUTINE centers_and_spreads

!******************************************************************************
!
! deallocates all data 
!
!******************************************************************************
SUBROUTINE qs_loc_destroy ( qs_loc_env )

    TYPE(qs_loc_env_type), POINTER           :: qs_loc_env

    INTEGER                                  :: i, j

  IF (qs_loc_env%is_initialized.ne.777) CALL stop_program("qs_localization","problems")

  DO i=1,SIZE ( qs_loc_env % operator, 2 )
   DO j=1,SIZE ( qs_loc_env % operator, 1 )
     CALL deallocate_matrix(qs_loc_env%operator(j,i)%matrix)
   ENDDO
  END DO
  CALL cp_para_env_release(qs_loc_env%para_env) 

END SUBROUTINE qs_loc_destroy

!******************************************************************************
!
! maps wfc's to molecules 
!
!******************************************************************************
SUBROUTINE wfc_to_molecule ( qs_loc_env, center, molecule_set ) 
    TYPE(qs_loc_env_type), INTENT(IN)        :: qs_loc_env
    REAL(KIND=dp), INTENT(IN)                :: center( :, : )
    TYPE(molecule_type), POINTER             :: molecule_set( : )

    CHARACTER(LEN=default_string_length)     :: mol_name
    INTEGER :: counter, first_atom, i, iatom, ikind, imol, imol_now, isos, &
      istate, local_location, natom, natom_loc, natom_max, nkind, nmol, nstate
    INTEGER, POINTER                         :: wfc_to_atom_map( : )
    REAL(KIND=dp)                            :: dr(3), mydist(2), r_pbc( 3 )
    REAL(KIND=dp), POINTER                   :: distance( : ), r( :, : )
    TYPE(atom_type), DIMENSION(:), POINTER   :: atom_list
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(particle_type), POINTER             :: particle_set( : )

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  particle_set => qs_loc_env % particle_set
  para_env => qs_loc_env % para_env
  local_molecules => qs_loc_env % local_molecules
  nstate = SIZE ( center, 1 )

  ALLOCATE(wfc_to_atom_map(nstate),STAT=isos)
  IF(isos/=0) CALL stop_memory('wfc_to_molecule','wfc_to_atom_map')
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  nkind = SIZE (local_molecules%n_el)
  natom = 0
  natom_max = 0
  DO ikind=1,nkind
    nmol = SIZE ( local_molecules%list(ikind)%array )
    DO imol=1,nmol
      i = local_molecules%list(ikind)%array(imol)
!      WRITE(*,*) para_env%mepos,"ABC",ikind,imol,i
      molecule_kind => molecule_set(i)%molecule_kind
      CALL get_molecule_kind(molecule_kind=molecule_kind,&
                             natom=natom)
      natom_max = natom_max + natom
      molecule_set(i)%lmi%nstates=0
      IF(ASSOCIATED(molecule_set(i)%lmi%states)) THEN
        DEALLOCATE(molecule_set(i)%lmi%states,STAT=isos)
        IF(isos/=0) CALL stop_memory('wfc_to_molecule','states')
      END IF
    END DO
  END DO
  natom_loc = natom_max
  natom = natom_max

  CALL mp_max(natom_max,para_env%group)

!  WRITE(*,*) para_env%mepos,"nkind  = ",nkind
!  WRITE(*,*) para_env%mepos,"natom_loc  = ",natom_loc
!  WRITE(*,*) para_env%mepos,"natom_max  = ",natom_max
!  WRITE(*,*) para_env%mepos,"nstate = ",nstate

  ALLOCATE(r(3,natom_max),STAT=isos)
  IF(isos/=0) CALL stop_memory('centers_and_spreads','r')

  ALLOCATE(distance(natom_max),STAT=isos)
  IF(isos/=0) CALL stop_memory('wfc_to_molecule','distance')

  !Zero all the stuff
  r(:,:) = 0.0_dp
  distance(:) = 1.E10_dp

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  counter = 0
  nkind = SIZE (local_molecules%n_el)
  DO ikind=1,nkind
    nmol = SIZE ( local_molecules%list(ikind)%array )
    DO imol=1,nmol
      i = local_molecules%list(ikind)%array(imol)
      molecule_kind => molecule_set(i)%molecule_kind
      first_atom = molecule_set(i)%first_atom
      CALL get_molecule_kind(molecule_kind=molecule_kind,&
                             natom=natom,&
                             name=mol_name,&
                             atom_list=atom_list)
      DO iatom=1,natom
        counter = counter + 1
        r ( :, counter ) = particle_set(first_atom+iatom-1) % r ( : )
      END DO
    END DO
  END DO

  !DO iatom=1,natom_loc
  !  WRITE(*,*) iatom,r(:,iatom)
  !END DO

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  DO istate=1,nstate
    distance(:) = 1.E10_dp
    DO iatom=1,natom_loc
      dr(1) = r(1,iatom) - center(istate,1)
      dr(2) = r(2,iatom) - center(istate,2)
      dr(3) = r(3,iatom) - center(istate,3)
      r_pbc = pbc ( dr, qs_loc_env%cell )
      distance(iatom) = SQRT ( DOT_PRODUCT ( r_pbc, r_pbc ) )
    END DO

    !combine distance() from all procs
    local_location = MINLOC(distance,DIM=1)
    
    mydist(1) = distance(local_location)
    mydist(2) = para_env%mepos

    CALL mp_minloc(mydist,para_env%group)
    
    IF(mydist(2)==para_env%mepos) THEN
      wfc_to_atom_map(istate)  =  local_location
    ELSE
      wfc_to_atom_map(istate)  =  0
    END IF
  END DO
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  IF(natom_loc/=0) THEN
    DO istate = 1,nstate
      iatom = wfc_to_atom_map(istate)
      IF(iatom/=0) THEN
        counter = 0
        nkind = SIZE (local_molecules%n_el)
        DO ikind=1,nkind
          nmol = SIZE ( local_molecules%list(ikind)%array )
          DO imol=1,nmol
            imol_now = local_molecules%list(ikind)%array(imol)
            molecule_kind => molecule_set(imol_now)%molecule_kind
            CALL get_molecule_kind(molecule_kind=molecule_kind,&
                                   natom=natom)
            counter = counter + natom
            IF(counter>=iatom) EXIT
          END DO
          IF(counter>=iatom) EXIT
        END DO
       
!        WRITE(*,*) "PP",para_env%mepos,istate,iatom,imol_now
     
        i = molecule_set(imol_now)%lmi%nstates
        i = i + 1
        molecule_set(imol_now)%lmi%nstates = i
  
        CALL reallocate(molecule_set(imol_now)%lmi%states,1,i)
        molecule_set(imol_now)%lmi%states(i) = istate
      END IF
    END DO
  END IF

!  nkind = SIZE (local_molecules%n_el)
!  DO ikind=1,nkind
!    nmol = SIZE ( local_molecules%list(ikind)%array )
!    DO imol=1,nmol
!      i = local_molecules%list(ikind)%array(imol)
!      WRITE(*,*)"XX",para_env%mepos,i,molecule_set(i)%lmi%states
!    END DO
!  END DO
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  DEALLOCATE(distance,STAT=isos)
  IF(isos/=0) CALL stop_memory('wfc_to_molecule','distance')

  DEALLOCATE(wfc_to_atom_map,STAT=isos)
  IF(isos/=0) CALL stop_memory('wfc_to_molecule','wfc_to_atom_map')
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  
  !CALL stop_program("wfc_to_atom","checkpoint")

END SUBROUTINE wfc_to_molecule
!------------------------------------------------------------------------------

END MODULE qs_localization

