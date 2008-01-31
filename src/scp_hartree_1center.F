!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****** cp2k/scp_hartree_1center
!!
!!   NAME
!!     scp_hartree_1center
!!
!!   FUNCTION
!!     Perform the calculation of the hartree 1 center 3 electron
!!     self-interaction terms for the SCP method
!!
!!   AUTHOR
!!      CJM
!!
!!   SOURCE
!******************************************************************************
MODULE  scp_hartree_1center
  USE ai_coulomb,                      ONLY: coulomb3, coulomb2
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
       get_atomic_kind,&
       get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
       gto_basis_set_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE kinds,                           ONLY: dp,dp_size,&
       int_size
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: ncoset, &
       nco, &
       indco, &
       nsoset, &
       nso
  USE orbital_transformation_matrices, ONLY: orbtramat
  USE qs_environment_types,            ONLY: get_qs_env,&
       qs_environment_type
  USE qs_util,                         ONLY: trace_r_AxB_new, trace_r_AxB
  USE scp_environment_types,           ONLY: scp_environment_type, &
       get_scp_env
  USE scp_coeff_types,                 ONLY: aux_coeff_set_type, & 
       aux_coeff_type, &
       aux_coeff_zero_fc
  USE scp_energy_types,                ONLY: scp_energy_type
  USE sparse_matrix_types,             ONLY: get_block_node,&
       real_matrix_p_type, &
       create_dbm_from_sm, &
       deallocate_dbm_from_sm,&
       allocate_matrix_set
  USE termination,                     ONLY: stop_memory,&
       stop_program
  USE timings,                         ONLY: timeset,&
       timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  ! *** Global parameters (only in this module)

  CHARACTER(LEN=*), PARAMETER :: moduleN = "scp_hartree_1center"

  ! Public Subroutine

  PUBLIC :: integrate_a_vhscp_b, integrate_vhscp_gscp, integrate_vhscp_gcore

  !******************************************************************************
CONTAINS
  !!****f* scp_hartree_1center/integrate_a_vhscp_b [1.0] *
  !!
  !!   NAME
  !!     integrate_a_vscp_b
  !!
  !!   FUNCTION
  !!     Analytic calculation of the 1 center integrals that make up the self-terms
  !!     for SCP.   
  !!
  !!   NOTE
  !!
  !!   INPUTS
  !!
  !!
  !!   AUTHOR
  !!     CJM
  !!
  !!
  !!*** **********************************************************************

  SUBROUTINE integrate_a_vhscp_b (qs_env,h_mat,p_mat, just_energy, error, debug)

    TYPE(qs_environment_type), POINTER         :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
         INTENT(IN)                               :: p_mat
    TYPE(real_matrix_p_type), DIMENSION(:), &
         INTENT(INOUT)                            :: h_mat
    LOGICAL, INTENT ( IN )                     :: just_energy
    LOGICAL, INTENT ( IN ), OPTIONAL           :: debug
    TYPE(cp_error_type), INTENT(inout)  &
         :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'integrate_a_vhscp_b', &
         routineP = moduleN//':'//routineN


    TYPE(atomic_kind_type), POINTER            :: atomic_kind_set ( : )
    TYPE(atomic_kind_type),       POINTER      :: atomic_kind
    TYPE(distribution_1d_type), POINTER        :: local_particles
    TYPE(aux_coeff_type), POINTER              :: local_coeffs
    TYPE(cp_logger_type), POINTER              :: logger
    TYPE(gto_basis_set_type),     POINTER      :: aux_basis, orb_basis
    TYPE(aux_coeff_set_type),     POINTER      :: aux_coeff_set
    TYPE(scp_environment_type ),  POINTER      :: scp_env
    TYPE(scp_energy_type),        POINTER      :: energy
    TYPE(cp_para_env_type), POINTER            :: para_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
         POINTER                                :: diagmat_ks, diagmat_p
    REAL(dp), DIMENSION(:,:),     POINTER      :: gccc
    REAL(dp), DIMENSION(:,:),     POINTER      :: h_block,p_block
    REAL(dp), DIMENSION(:,:),     POINTER      :: h_block2,p_block2
    REAL(dp), DIMENSION(:,:),     POINTER      :: rpgf, sphi
    REAL(dp), DIMENSION(:,:),     POINTER      :: hab, fcoeff, work, work_scp
    REAL(dp)                                   :: e_hartree_1c, fac, zetc, rpgfc, trace
    REAL(dp), DIMENSION ( :, : ), POINTER      :: rpgf_scp
    REAL(dp), DIMENSION ( :, : ), POINTER      :: sphi_scp, zet, zet_scp
    REAL(dp), DIMENSION(:,:,:),    POINTER     :: CPC_ab, Int_abc, int_tmp
    REAL(dp), DIMENSION(:,:),     POINTER      :: Int_ab_sum
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE  :: vv
    REAL(dp), DIMENSION(:),   ALLOCATABLE      :: ff
    REAL(dp), DIMENSION(:), POINTER            :: rpgfa, rpgfb
    REAL(dp), DIMENSION(:),     POINTER        :: zeta,zetb

    INTEGER :: handle, iatom,  i, istat
    INTEGER :: ic1, ic1_min_m1, ic1_max
    INTEGER :: ic2, ic2_min_m1, ic2_max
    INTEGER :: ip1, ip2, m1, m2, ifs1, ifs2, nz1, nz2
    INTEGER :: ns1, ns2, nco1, nco2, n1, n2, nn1, nn2
    INTEGER :: ldintab1, ldintab2, ldhab1, ldhab2, ldcpc1, ldcpc2
    INTEGER :: iparticle_local, nparticle_local, ikind, nkind, nspins, ispin
    INTEGER :: iseta, isetb, ico1, ipgf_scp, iset_scp
    INTEGER :: ncotot, maxco, nset, nsgf
    INTEGER :: ncotot_scp,maxco_scp,nset_scp,nsgf_scp
    INTEGER :: nco_scp, sgf_scp, icoscp, offset
    INTEGER :: maxco_global, maxlgto, maxsgf_set, output_unit
    INTEGER :: lshell, is, iso, ic, ico

    INTEGER, DIMENSION ( : ), POINTER :: l_max, l_min, npgf, nsgf_set
    INTEGER, DIMENSION ( : ), POINTER :: lscp_max, lscp_min, npgf_scp, nsgf_scp_set
    INTEGER, DIMENSION ( :, : ), POINTER :: first_sgf_scp, first_sgf

    LOGICAL                                    :: ionode, my_debug

    !   ---------------------------------------------------------------------------

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) output_unit= cp_logger_get_default_unit_nr(logger)

    NULLIFY(atomic_kind_set, local_particles, aux_basis, aux_coeff_set)
    NULLIFY( para_env, scp_env, local_coeffs, energy )

    my_debug = .FALSE.
    IF ( PRESENT ( debug ) ) my_debug = debug

    ! Zero energy
    e_hartree_1c = 0._dp

    fac = 1.0_dp
    nspins = qs_env%dft_control%nspins

    CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,&
         para_env=para_env, scp_env = scp_env, &
         local_particles=local_particles, error=error)

    CALL get_scp_env(scp_env=scp_env, energy=energy, &
         aux_coeff_set=aux_coeff_set, &
         error=error)

    ! Zero forces on the coefficients
    IF ( .NOT. just_energy ) &
         CALL aux_coeff_zero_fc ( aux_coeff_set, error )

    nkind = SIZE(atomic_kind_set,1)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         maxco=maxco_global,maxlgto=maxlgto,&
         maxsgf_set=maxsgf_set)


    CALL timeset(routineN,"I"," ",handle)

    !   Allocated diagonal matricies
    CALL allocate_matrix_set(diagmat_ks,nspins,error=error)
    CALL allocate_matrix_set(diagmat_p,nspins,error=error)


    !   Allocate the work array for the orbitals
    ALLOCATE(work(maxco_global,maxsgf_set),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "work",maxco_global*maxsgf_set)

    !   Allocate gccc array
    ALLOCATE (gccc(maxco_global,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "gccc",maxco_global*1*dp_size)
    gccc = 0._dp
    !   Allocate work_scp array
    ALLOCATE (work_scp (maxco_global,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "work_scp",maxco_global*1*dp_size)
    !   Allocate fcoeff array
    ALLOCATE (fcoeff (maxco_global,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "fcoeff",maxco_global*1*dp_size)

    NULLIFY(CPC_ab,hab)
    ALLOCATE(CPC_ab(1,1,nspins),STAT=istat)
    ldcpc1 = 1
    ldcpc2 = 1
    ALLOCATE(hab(1,1),STAT=istat)
    ldhab1 = 1
    ldhab2 = 1
    NULLIFY(Int_ab_sum,Int_abc)
    ALLOCATE(Int_ab_sum(1,1))
    ldintab1 = 1
    ldintab2 = 1
    CPC_ab = 0.0_dp
    hab = 0.0_dp

    ALLOCATE(int_tmp(maxco_global,maxco_global,maxco_global),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "int_tmp",(maxco_global*maxco_global*maxco_global)*dp_size)
    ALLOCATE(vv(ncoset(maxlgto),ncoset(maxlgto),ncoset(maxlgto),&
         maxlgto+maxlgto+maxlgto+1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "vv",ncoset(maxlgto)**2*ncoset(maxlgto)*(maxlgto+maxlgto+maxlgto+1)*dp_size)
    ALLOCATE(ff(0:maxlgto+maxlgto+maxlgto),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "ff",(1+maxlgto+maxlgto+maxlgto)*dp_size)

    DO ispin=1,nspins
       ! Allocate diagonal block matrices
       CALL create_dbm_from_sm(diagmat_p(ispin)%matrix,p_mat(ispin)%matrix,&
            "DIAGONAL BLOCK P MATRIX",.TRUE.,para_env,error=error)
       CALL create_dbm_from_sm(diagmat_ks(ispin)%matrix,h_mat(ispin)%matrix,&
            "DIAGONAL BLOCK KS MATRIX",.FALSE.,error=error)
    END DO

    DO ikind=1,nkind

       NULLIFY(atomic_kind,orb_basis)
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,aux_basis_set=aux_basis, &
            orb_basis_set=orb_basis )

       ! Check to see if the atom is an SCP atom
       IF (.NOT.ASSOCIATED(aux_basis)) CYCLE
       local_coeffs=>aux_coeff_set%coeffs_of_kind(ikind)%coeffs
       ! Get SCP basis functions (from AUX_BASIS)  corresponding to atom_A
       NULLIFY(first_sgf_scp,lscp_max,lscp_min,npgf_scp,nsgf_scp_set)
       NULLIFY(rpgf_scp,sphi_scp,zet_scp)
       CALL get_gto_basis_set(gto_basis_set=aux_basis,&
            first_sgf=first_sgf_scp,&
            lmax=lscp_max,&
            lmin=lscp_min,&
            maxco=maxco_scp,&
            npgf=npgf_scp,&
            nset=nset_scp,&
            nsgf=nsgf_scp,&
            nsgf_set=nsgf_scp_set,&
            pgf_radius=rpgf_scp,&
            sphi=sphi_scp,&
            zet=zet_scp)

       ncotot_scp = maxco_scp * nset_scp

       ! Get ORBITAL BASIS functions mu,nu corresponding to atom_A and atom_B
       NULLIFY(first_sgf,l_max,l_min,npgf,nsgf_set,rpgf,sphi,zet)
       CALL get_gto_basis_set(gto_basis_set=orb_basis,&
            first_sgf=first_sgf,&
            lmax=l_max,&
            lmin=l_min,&
            maxco=maxco,&
            npgf=npgf,&
            nset=nset,&
            nsgf=nsgf,&
            nsgf_set=nsgf_set,&
            pgf_radius=rpgf,&
            sphi=sphi,&
            zet=zet)
       ncotot = maxco * nset

       ! Get the number of particles of ikind ( local )
       nparticle_local = local_particles%n_el(ikind)

       ! Loop over the (local) particles
       DO iparticle_local=1,nparticle_local
          iatom = local_particles%list(ikind)%array(iparticle_local)

          ! Retrieve KS and density matrix block for this atom
          NULLIFY ( p_block, p_block2, h_block, h_block2 )
          CALL get_block_node(matrix=diagmat_ks(1)%matrix,&
               block_row=iatom,&
               block_col=iatom,&
               BLOCK=h_block)

          CALL get_block_node(matrix=diagmat_p(1)%matrix,&
               block_row=iatom,&
               block_col=iatom,&
               BLOCK=p_block)
          !         CALL get_block_node(matrix=h_mat(1)%matrix,&
          !                            block_row=iatom,&
          !                            block_col=iatom,&
          !                            BLOCK=h_block)
          !
          !         CALL get_block_node(matrix=p_mat(1)%matrix,&
          !                            block_row=iatom,&
          !                            block_col=iatom,&
          !                            BLOCK=p_block)

          IF(nspins == 2) THEN
             CALL get_block_node(matrix=diagmat_ks(2)%matrix,&
                  block_row=iatom,&
                  block_col=iatom,&
                  BLOCK=h_block2)

             CALL get_block_node(matrix=diagmat_p(2)%matrix,&
                  block_row=iatom,&
                  block_col=iatom,&
                  BLOCK=p_block2)
             !           CALL get_block_node(matrix=h_mat(2)%matrix,&
             !                            block_row=iatom,&
             !                            block_col=iatom,&
             !                            BLOCK=h_block2)
             !
             !           CALL get_block_node(matrix=p_mat(2)%matrix,&
             !                            block_row=iatom,&
             !                            block_col=iatom,&
             !                            BLOCK=p_block2)

          END IF
          ! Initializing arrays for the calculation of the coefficient matrix CPC
          CALL reallocate(CPC_ab,1,ncotot,1,ncotot,1,nspins)
          ldcpc1 = ncotot
          ldcpc2 = ncotot
          CPC_ab = 0.0_dp

          IF(maxco .GT. ldhab1 .OR. maxco .GT. ldhab2) THEN
             CALL reallocate(hab,1,maxco,1,maxco)
             ldhab1 = maxco
             ldhab2 = maxco
          ELSE
             hab = 0.0_dp
          END IF
          !    Build the CPC coefficients that should then multiplied by the
          !    3-center integrals of primitives to get the energy 

          !    Loop over the mu and nu ( iseta and isetb ) of the orbital basis
          m1 = 0
          DO iseta = 1,nset
             ifs1 = first_sgf(1,iseta)
             nco1 = npgf(iseta)*ncoset(l_max(iseta))
             ns1  = nsgf_set(iseta)
             ic1_min_m1 = ncoset(l_min(iseta)-1)
             ic1_max    = ncoset(l_max(iseta))
             m2 = 0
             DO isetb = 1,nset
                ifs2 = first_sgf(1,isetb)
                nco2 = npgf(isetb)*ncoset(l_max(isetb))
                ns2  = nsgf_set(isetb)
                ic2_min_m1 = ncoset(l_min(isetb)-1)
                ic2_max    = ncoset(l_max(isetb))

                CALL dgemm("N","N",nco1,ns2,ns1,&
                     1.0_dp,sphi(1,ifs1),maxco,&
                     p_block(ifs1,ifs2),nsgf,&
                     0.0_dp,work(1,1),maxco_global)
                CALL dgemm("N","T",nco1,nco2,ns2,&
                     1.0_dp,work(1,1),maxco_global,&
                     sphi(1,ifs2), maxco,&
                     0.0_dp,hab(1,1),ldhab1)


                DO ip1 = 1,npgf(iseta)
                   ic1 = ic1_min_m1+1
                   n1 = ic1_max*(ip1-1)
                   nn1 = ic1_max*ip1
                   DO ip2 = 1,npgf(isetb)
                      ic2 = ic2_min_m1+1
                      n2 =  ic2_max*(ip2-1)
                      nn2 = ic2_max*ip2
                      CPC_ab(ic1+n1+m1:nn1+m1,ic2+n2+m2:nn2+m2,1) = &
                           hab(ic1+n1:nn1,ic2+n2:nn2)
                   END DO
                END DO

                IF(nspins == 2) THEN
                   CALL dgemm("N","N",nco1,ns2,ns1,&
                        1.0_dp,sphi(1,ifs1),maxco,&
                        p_block2(ifs1,ifs2),nsgf,&
                        0.0_dp,work(1,1),maxco_global)
                   CALL dgemm("N","T",nco1,nco2,ns2,&
                        1.0_dp,work(1,1),maxco_global,&
                        sphi(1,ifs2), maxco,&
                        0.0_dp,hab(1,1),ldhab1)


                   DO ip1 = 1,npgf(iseta)
                      ic1 = ic1_min_m1+1
                      n1 = ic1_max*(ip1-1)
                      nn1 = ic1_max*ip1
                      DO ip2 = 1,npgf(isetb)
                         ic2 = ic2_min_m1+1
                         n2 =  ic2_max*(ip2-1)
                         nn2 = ic2_max*ip2
                         CPC_ab(ic1+n1+m1:nn1+m1,ic2+n2+m2:nn2+m2,2) = &
                              hab(ic1+n1:nn1,ic2+n2:nn2)
                      END DO
                   END DO

                END IF

                m2 = m2 + maxco
             END DO  ! isetb
             m1 = m1 + maxco
          END DO  ! iseta

          ! Allocate array for sum of integrals
          IF(maxco .GT. ldintab1 .OR. maxco .GT. ldintab2) THEN
             CALL reallocate(Int_ab_sum,1,maxco,1,maxco)
             ldintab1 = maxco
             ldintab2 = maxco
          END IF
          !   Allocate array for the 3-centers integrals
          CALL reallocate(Int_abc,1,ncotot,1,ncotot,1,ncotot_scp)

          !------Now loop over SCP coefficients and compute integrals-----

          offset = 0
          DO iset_scp = 1,nset_scp
             nco_scp = npgf_scp ( iset_scp ) * ncoset ( lscp_max ( iset_scp ) )
             sgf_scp = first_sgf_scp ( 1, iset_scp )

             ! SCP density coefficients of each contracted
             ! gaussian polarization functions (i.e. the dynamical variables associated to
             ! polarization), are decontracted with  sphi and stored in  gccc(:,:).

             DO i=1,nsgf_scp_set ( iset_scp )
                work_scp ( i, 1 )=-local_coeffs%c(iparticle_local,offset+i)
             ENDDO

             CALL dgemm("N","N",nco_scp,1,nsgf_scp_set ( iset_scp ),&
                  1.0_dp,sphi_scp(1,sgf_scp),SIZE(sphi_scp,1),&
                  work_scp(1,1),SIZE(work_scp,1),&
                  0.0_dp,gccc (1,1),SIZE(gccc,1))

             Int_abc = 0.0_dp
             ! Loop over the primitive exponents of the SCP-basis (AUXILLARY BASIS)
             DO ipgf_scp = 1, npgf_scp ( iset_scp )
                zetc = zet_scp ( ipgf_scp, iset_scp ) 
                rpgfc = rpgf_scp ( ipgf_scp, iset_scp )
                ! Loop over the mu and nu ( iseta and isetb ) of the orbital basis
                m1 = 0
                DO iseta = 1, nset
                   ifs1 = first_sgf ( 1, iseta )
                   nz1 = npgf ( iseta ) 
                   nco1 = npgf ( iseta ) * ncoset ( l_max ( iseta ) )
                   ns1 = nsgf_set ( iseta )
                   rpgfa => rpgf ( 1:nz1, iseta )
                   zeta => zet ( 1:nz1, iseta )
                   m2 = 0
                   DO isetb = 1, nset
                      ifs2 = first_sgf ( 1, isetb )
                      nz2 = npgf ( isetb ) 
                      nco2 = npgf ( isetb ) * ncoset ( l_max ( isetb ) )
                      ns2 = nsgf_set ( isetb )
                      rpgfb => rpgf ( 1:nz2, isetb )
                      zetb => zet ( 1:nz2, isetb )
                      ! For energies and KS-matrix update only Int_ab_sum is used on output
                      Int_ab_sum = 0.0_dp
                      CALL coulomb3 ( l_max ( iseta ), nz1, zeta , rpgfa, l_min ( iseta ),       &
                           l_max ( isetb ), nz2, zetb,  rpgfb, l_min ( isetb ),       &
                           lscp_max ( iset_scp ), zetc, rpgfc, 0, &
                           gccc(:,1), (/0.0_dp,0.0_dp,0.0_dp/), 0.0_dp, &   
                           (/0.0_dp,0.0_dp,0.0_dp/), 0.0_dp, &   
                           (/0.0_dp,0.0_dp,0.0_dp/), 0.0_dp, &   
                           Int_ab_sum, int_tmp, vv, ff )
                      ! update KS-matrix
                      IF ( .NOT. just_energy ) THEN
                         CALL dgemm("N","N",nco1,ns2,nco2,&
                              1.0_dp,Int_ab_sum(1,1),ldintab1,&
                              sphi(1,ifs2),SIZE(sphi,1),&
                              0.0_dp,work(1,1),maxco_global)

                         CALL dgemm("T","N",ns1,ns2,nco1,&
                              fac,sphi(1,ifs1),SIZE(sphi,1),&
                              work(1,1),maxco_global,&
                              1.0_dp,h_block(ifs1,ifs2),SIZE(h_block,1))

                         IF(nspins == 2) THEN
                            CALL dgemm("T","N",ns1,ns2,nco1,&
                                 fac,sphi(1,ifs1),SIZE(sphi,1),&
                                 work(1,1),maxco_global,&
                                 1.0_dp,h_block2(ifs1,ifs2),SIZE(h_block,1))
                         END IF
                      END IF ! just_energy
                      ! Compute energy
                      DO ispin = 1,nspins
                         trace = trace_r_AxB_new(CPC_ab(:,:,ispin),ldcpc1,m1,m2,&
                              Int_ab_sum,ldintab1,0,0,nco1,nco2)

                         e_hartree_1c = e_hartree_1c + trace
                      END DO
                      ! For gradients with respect to SCP coefficients only int_tmp is used on output

                      ! Introducing Int_abc so we can trace over nco1, nco2 to obtain a function with
                      ! index of icoscp to get SCP forces, namely:
                      ! Int_abc(1+m1:nco1+m1,1+m2:nco2+m2,ico) = int_tmp(1:nco1,1:nco2,ico)
                      DO icoscp = 1, nco_scp
                         DO ico1 = 1,nco2
                            !  CALL DCOPY(nco1,int_tmp(1,ico1,icoscp),1,Int_abc(m1+1,m2+ico1,icoscp),1)
                            CALL DAXPY(nco1,1.0_dp,int_tmp(1,ico1,icoscp),1,Int_abc(m1+1,m2+ico1,icoscp),1)
                         END DO
                      END DO

                      m2 = m2 + maxco
                   END DO ! isetb
                   m1 = m1 + maxco
                END DO ! iseta
             END DO ! ipgf_scp
             ! update coefficient force
             IF ( .NOT. just_energy ) THEN
                work_scp = 0.0_dp
                DO ispin = 1, nspins
                   DO icoscp = 1, nco_scp
                      work_scp ( icoscp, 1 ) = work_scp ( icoscp, 1 ) + &
                           trace_r_AxB ( CPC_ab ( :, :, ispin ), &
                           ldcpc1, Int_abc ( :, :, icoscp ), &
                           ldcpc1, ldcpc1, ldcpc2 )
                   END DO
                END DO
                ! multiply work_scp by the contraction coeff (sphi_scp) to get the correct fcoeffs 
                CALL dgemm("T","N",nsgf_scp_set ( iset_scp ),1,nco_scp,&
                     1.0_dp,sphi_scp(1,sgf_scp),SIZE(sphi_scp,1),&
                     work_scp(1,1),SIZE(work_scp,1),&
                     0.0_dp,fcoeff(1,1),SIZE(fcoeff,1))
                DO i = 1, nsgf_scp_set ( iset_scp )
                   local_coeffs%fc(iparticle_local,offset+i) = &
                        local_coeffs%fc(iparticle_local,offset+i) + fcoeff(i,1)
                   IF ( my_debug ) WRITE ( *, * ) 'F_ANALYTICAL', iatom, i, fcoeff ( i, 1 )
                END DO
             ENDIF ! just energy
             offset = offset + nsgf_scp_set ( iset_scp )
          END DO ! iset_scp
       END DO ! iparticle
    END DO ! ikind

    DO ispin = 1, nspins
       CALL deallocate_dbm_from_sm(diagmat_p(ispin)%matrix,&
            p_mat(1)%matrix,.FALSE.,error=error)
       CALL deallocate_dbm_from_sm(diagmat_ks(ispin)%matrix,&
            h_mat(ispin)%matrix,.TRUE.,para_env,error=error)
    END DO
    DEALLOCATE(diagmat_ks,diagmat_p)

    !   IF PARALLEL sum up the contributions to e_scp_hartree_1c from
    !   ab pairs handled by different processors
    CALL mp_sum(e_hartree_1c,para_env%group)
    energy % e_scp_ks_self = e_hartree_1c

    IF(ASSOCIATED(Int_abc)) THEN
       DEALLOCATE(Int_abc,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
            "Int_abc")

    END IF
    IF(ASSOCIATED(hab)) THEN
       DEALLOCATE(hab, STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
            "hab")
    END IF

    IF(ASSOCIATED(CPC_ab)) THEN
       DEALLOCATE(CPC_ab, STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
            "CPC_ab")
    END IF

    DEALLOCATE(work,STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,__LINE__,&
            "work")
    END IF

    DEALLOCATE(work_scp, gccc, fcoeff, STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,__LINE__,&
            "work_scp, gccc, fcoeff")
    END IF

    DEALLOCATE(int_tmp,vv,ff,STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,__LINE__,&
            "int_tmp,vv,ff")
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE integrate_a_vhscp_b

  !******************************************************************************
  !!****f* scp_hartree_1center/integrate_vhscp_gscp [1.0] *
  !!
  !!   NAME
  !!     integrate_vhscp_gscp
  !!
  !!   FUNCTION
  !!     Analytic calculation of the 1 center integrals that make up the self-terms
  !!     for SCP.   
  !!
  !!   NOTE
  !!
  !!   INPUTS
  !!
  !!
  !!   AUTHOR
  !!     CJM
  !!
  !!
  !!*** **********************************************************************

  SUBROUTINE integrate_vhscp_gscp ( scp_env, atomic_kind_set, local_particles, just_energy, error, debug )

    TYPE(scp_environment_type), POINTER                 :: scp_env
    TYPE(atomic_kind_type), DIMENSION(:), POINTER       :: atomic_kind_set
    TYPE(distribution_1d_type), POINTER                 :: local_particles
    LOGICAL, INTENT ( IN ), OPTIONAL                    :: just_energy
    TYPE(cp_error_type), INTENT(inout)                  :: error
    LOGICAL, INTENT ( IN ), OPTIONAL                    :: debug

    CHARACTER(len=*), PARAMETER :: routineN = 'integrate_vhscp_gscp', &
         routineP = moduleN//':'//routineN


    TYPE(atomic_kind_type),       POINTER      :: atomic_kind
    TYPE(gto_basis_set_type),     POINTER      :: aux_basis
    TYPE(cp_para_env_type), POINTER            :: para_env
    TYPE ( aux_coeff_set_type ), POINTER       :: aux_coeff_set
    TYPE ( aux_coeff_type ), POINTER           :: local_coeffs
    TYPE ( scp_energy_type ), POINTER          :: energy
    REAL(kind=dp)                              :: e_scp

    REAL(dp), DIMENSION(:,:),       ALLOCATABLE  :: work, gcca, gccb,fcoeff
    REAL(dp), DIMENSION(:,:,:),   ALLOCATABLE  :: vv
    REAL(dp), DIMENSION ( :, : ), POINTER      :: vab

    REAL(dp), DIMENSION ( :, : ), POINTER      :: sphi, zet, rpgf
    REAL(dp), DIMENSION ( : ), POINTER         :: zeta, zetb, rpgfa, rpgfb
    REAL(dp), DIMENSION(:),   ALLOCATABLE      :: ff

    INTEGER :: handle, iatom,  i, istat
    INTEGER :: icoa, icob, ncoa, ncob
    INTEGER :: iparticle_local, nparticle_local, ikind, nkind
    INTEGER :: iseta, isetb, sgfa, sgfb
    INTEGER :: ncotot, maxco, nset
    INTEGER :: offseta, offsetb, nza, nzb, ldvab1, ldvab2
    INTEGER :: maxco_global, maxlgto, maxsgf_set, output_unit
    INTEGER, DIMENSION ( : ), POINTER :: l_max, l_min, npgf, nsgf_set
    INTEGER, DIMENSION ( :, : ), POINTER :: first_sgf

    LOGICAL :: ionode, my_debug, energy_only



    !   ---------------------------------------------------------------------------
    my_debug = .FALSE.
    IF ( PRESENT ( debug ) ) my_debug = debug
    energy_only = .FALSE.
    IF ( PRESENT ( just_energy ) ) energy_only = just_energy

    NULLIFY( atomic_kind, aux_basis, aux_coeff_set)
    NULLIFY( energy, local_coeffs, para_env )

    CALL get_scp_env(scp_env=scp_env, para_env=para_env, & 
         energy=energy, aux_coeff_set=aux_coeff_set, &
         error=error)

    nkind = SIZE(atomic_kind_set,1)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         maxco=maxco_global,maxlgto=maxlgto)


    CALL timeset(routineN,"I"," ",handle)

    ALLOCATE (work (maxco_global,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "work",maxco_global*1*dp_size)
    ALLOCATE (gcca (maxco_global,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "gcca",maxco_global*1*dp_size)
    ALLOCATE (gccb (maxco_global,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "gccb",maxco_global*1*dp_size)
    ALLOCATE (fcoeff (maxco_global,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "fcoeff",maxco_global*1*dp_size)

    ALLOCATE(vab(1,1),STAT = istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "vab",1*1*dp_size)
    ldvab1 = 1
    ldvab2 = 1
    vab = 0.0_dp
    ALLOCATE(vv(ncoset(maxlgto),ncoset(maxlgto),maxlgto+maxlgto+1),STAT = istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "vv",maxco_global*maxco_global*(maxlgto+maxlgto+1)*dp_size)
    ALLOCATE(ff(0:maxlgto+maxlgto),STAT = istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "ff", (maxlgto+maxlgto+1)*dp_size)

    e_scp = 0._dp
    DO ikind = 1,nkind

       NULLIFY(atomic_kind)
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,aux_basis_set=aux_basis )


       ! Check to see if the atom is an SCP atom
       IF (.NOT.ASSOCIATED(aux_basis)) CYCLE
       local_coeffs=>aux_coeff_set%coeffs_of_kind(ikind)%coeffs
       ! Get SCP basis functions (from AUX_BASIS)  corresponding to atom_A
       NULLIFY(first_sgf,l_max,l_min,npgf,nsgf_set)
       NULLIFY(rpgf,sphi,zet)
       CALL get_gto_basis_set(gto_basis_set=aux_basis,&
            first_sgf=first_sgf,&
            lmax=l_max,&
            lmin=l_min,&
            maxco=maxco,&
            npgf=npgf,&
            nset=nset,&
            nsgf_set=nsgf_set,&
            pgf_radius=rpgf,&
            sphi=sphi,&
            zet=zet)
       ncotot = maxco * nset
       ! Get the number of particles of ikind ( local )
       nparticle_local = local_particles%n_el(ikind)
       DO iparticle_local = 1, nparticle_local
          iatom = local_particles%list(ikind)%array(iparticle_local)
          offseta = 0
          DO iseta = 1, nset
             ncoa = npgf ( iseta ) * ncoset ( l_max ( iseta ) )
             WRITE ( *, * ) 'NCOAM', npgf ( iseta ), ncoset ( 0 ), 0
             WRITE ( *, * ) 'NCOA', npgf ( iseta ), ncoset ( l_max ( iseta ) ), l_max ( iseta )
             WRITE ( *, * ) 'NCOAP',npgf ( iseta ), ncoset ( l_max ( iseta ) + 1 ), l_max ( iseta ) + 1
             sgfa = first_sgf ( 1, iseta )
             nza = npgf ( iseta )
             rpgfa => rpgf ( 1:nza, iseta )
             zeta => zet ( 1 : nza, iseta )
             ! Assigning the coefficients
             DO i=1, nsgf_set ( iseta )
                work ( i, 1 ) = local_coeffs % c ( iparticle_local, offseta+i )
             END DO
             CALL dgemm ( "N", "N", ncoa, 1, nsgf_set ( iseta ), &
                  1.0_dp, sphi ( 1, sgfa ), SIZE ( sphi, 1 ), &
                  work ( 1, 1 ), SIZE ( work, 1 ), &
                  0.0_dp, gcca ( 1, 1 ), SIZE ( gcca, 1 ) )
             offsetb = 0
             DO isetb = 1, nset
                ncob = npgf ( isetb ) * ncoset ( l_max ( isetb ) )
                sgfb = first_sgf ( 1, isetb )
                nzb = npgf ( isetb )
                rpgfb => rpgf ( 1 : nzb, isetb )
                zetb => zet ( 1 : nzb, isetb )
                ! Allcating the arrays for the integrals
                IF ( ncoa > ldvab1 .OR.  ncob > ldvab2 ) THEN
                   CALL reallocate ( vab , 1, ncoa, 1, ncob )
                   ldvab1 = ncoa
                   ldvab2 = ncob
                ELSE
                   vab = 0._dp 
                ENDIF
                ! Assigning the coefficients
                DO i = 1, nsgf_set ( isetb )
                   work ( i, 1 ) = local_coeffs % c ( iparticle_local, offsetb+i )
                END DO
                CALL dgemm ( "N", "N", ncob, 1, nsgf_set ( isetb ), &
                     1.0_dp, sphi ( 1, sgfb ), SIZE ( sphi, 1 ), &
                     work ( 1, 1 ), SIZE ( work, 1 ), &
                     0.0_dp, gccb ( 1, 1 ), SIZE ( gccb, 1 ) )
                CALL coulomb2( l_max ( iseta ), npgf ( iseta ), zeta, rpgfa, l_min ( iseta ), &
                     l_max ( isetb ), npgf ( isetb ), zetb, rpgfb, l_min ( isetb ), &
                     (/0.0_dp,0.0_dp,0.0_dp/), 0.0_dp, vab, vv, ff ( 0: ) )
                work = 0.0_dp
                DO icoa = 1,  ncoa ! There might be a stride... Need to check
                   DO icob = 1, ncob ! There might be a stride... Need to check
                      work ( icoa, 1 )  = work ( icoa, 1 ) + gccb ( icob, 1 ) * vab ( icoa, icob )
                   END DO
                   e_scp = e_scp - gcca ( icoa , 1) * work ( icoa, 1 )
                END DO
                offsetb = offsetb + nsgf_set ( isetb )
             END DO ! isetb
             ! Forces on the SCP coefficients
             IF ( .NOT. energy_only ) THEN
                CALL dgemm("T","N",nsgf_set ( iseta ),1,ncoa,&
                     1.0_dp,sphi(1,sgfa),SIZE(sphi,1),&
                     work(1,1),SIZE(work,1),&
                     0.0_dp,fcoeff(1,1),SIZE(fcoeff,1))
                DO i = 1, nsgf_set ( iseta )
                   local_coeffs%fc(iparticle_local,offseta+i) = &
                        local_coeffs%fc(iparticle_local,offseta+i) + fcoeff(i,1)
                   IF ( my_debug ) WRITE ( *, * ) 'F_ANALYTICAL', iatom, i, fcoeff ( i, 1 )
                END DO
             ENDIF ! energy_only
             offseta = offseta + nsgf_set ( iseta )
          END DO ! iseta
       END DO ! iparticle

    END DO  ! ikind
    e_scp = 0.5_dp * e_scp

    CALL mp_sum(e_scp,para_env%group)
    energy % e_scp_self = e_scp

    DEALLOCATE(vab,vv,ff,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "vab.vv.ff")

    DEALLOCATE(work, gcca, gccb, fcoeff,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "work.gcca.gccb,fcoeff")

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE integrate_vhscp_gscp

  !******************************************************************************
  !!****f* scp_hartree_1center/integrate_vhscp_gcore [1.0] *
  !!
  !!   NAME
  !!     integrate_vhscp_gcore
  !!
  !!   FUNCTION
  !!     Analytic calculation of the 1 center integrals that make up the self-terms
  !!     for SCP.   
  !!
  !!   NOTE
  !!
  !!   INPUTS
  !!
  !!
  !!   AUTHOR
  !!     CJM
  !!
  !!
  !!*** **********************************************************************
  SUBROUTINE integrate_vhscp_gcore ( scp_env, atomic_kind_set, local_particles, just_energy, error, debug )

    TYPE(scp_environment_type), POINTER                 :: scp_env
    TYPE(atomic_kind_type), DIMENSION(:), POINTER       :: atomic_kind_set
    TYPE(distribution_1d_type), POINTER                 :: local_particles
    TYPE(cp_error_type), INTENT(inout)                  :: error
    LOGICAL, INTENT ( IN ), OPTIONAL                    :: debug, just_energy

    CHARACTER(len=*), PARAMETER :: routineN = 'integrate_vhscp_gcore', &
         routineP = moduleN//':'//routineN


    TYPE(atomic_kind_type),       POINTER      :: atomic_kind
    TYPE(gto_basis_set_type),     POINTER      :: aux_basis
    TYPE(aux_coeff_set_type),     POINTER      :: aux_coeff_set
    TYPE(aux_coeff_type),     POINTER          :: local_coeffs
    TYPE(scp_energy_type),    POINTER          :: energy
    TYPE(cp_para_env_type), POINTER            :: para_env

    REAL(dp), DIMENSION(:),       ALLOCATABLE  :: ff
    REAL(dp), DIMENSION(:,:),       ALLOCATABLE  :: zetacore
    REAL(dp), DIMENSION(:,:),       ALLOCATABLE  :: work, gcca, fcoeff
    REAL(dp), DIMENSION(:,:),   POINTER          ::  vab
    REAL(dp), DIMENSION(:,:,:),   ALLOCATABLE  :: vv
    REAL(dp) :: e_scp, gccore

    REAL(dp), DIMENSION ( :, : ), POINTER      :: sphi, zet, rpgf
    REAL(dp), DIMENSION ( : ), POINTER         :: zeta, zetcore, rpgfa, rpgfcore

    INTEGER :: handle, iatom,  i, istat
    INTEGER :: icoa, ncoa, ldvab1
    INTEGER :: iparticle_local, nparticle_local, ikind, nkind
    INTEGER :: iseta, sgfa
    INTEGER :: ncotot, maxco, nset
    INTEGER :: offseta, nza
    INTEGER :: maxco_global, maxlgto, maxsgf_set, output_unit
    INTEGER, DIMENSION ( : ), POINTER :: l_max, l_min, npgf, nsgf_set
    INTEGER, DIMENSION ( :, : ), POINTER :: first_sgf

    LOGICAL :: ionode, my_debug, energy_only


    !   ---------------------------------------------------------------------------

    NULLIFY(para_env, atomic_kind, aux_basis)
    NULLIFY(aux_coeff_set,local_coeffs,energy)

    my_debug = .FALSE.
    IF ( PRESENT ( debug ) ) my_debug = debug
    energy_only = .FALSE.
    IF ( PRESENT ( just_energy ) ) energy_only = just_energy

    CALL get_scp_env(scp_env=scp_env, para_env=para_env,  &
         aux_coeff_set = aux_coeff_set, energy=energy, &
         error=error)

    nkind = SIZE(atomic_kind_set,1)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         maxco=maxco_global,maxlgto=maxlgto)


    CALL timeset(routineN,'I',' ',handle)

    ALLOCATE (zetcore ( 1 ),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "zetcore",1*dp_size)
    ALLOCATE (rpgfcore ( 1 ),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "rpgfcore",1*dp_size)
    ALLOCATE (work (maxco_global,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "work",maxco_global*1*dp_size)
    ALLOCATE (gcca (maxco_global,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "gcca",maxco_global*1*dp_size)
    ALLOCATE (fcoeff (maxco_global,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "fcoeff",maxco_global*1*dp_size)

    ALLOCATE(vab(1,1),STAT = istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "vab",maxco_global*1*dp_size)
    ldvab1 = 1
    vab = 0.0_dp
    ALLOCATE(vv(ncoset(maxlgto),1,maxlgto+1+1),STAT = istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "vv",maxco_global*1*(maxlgto+1+1)*dp_size)
    ALLOCATE(ff(0:maxlgto+1),STAT = istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "ff", (maxlgto+1+1)*dp_size)

    e_scp = 0._dp
    DO ikind = 1,nkind

       NULLIFY(atomic_kind)
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,aux_basis_set=aux_basis, &
            alpha_core_charge = zetcore ( 1 ), &
            ccore_charge = gccore )

       ! Check to see if the atom is an SCP atom
       IF (.NOT.ASSOCIATED(aux_basis)) CYCLE
       local_coeffs=>aux_coeff_set%coeffs_of_kind(ikind)%coeffs

       ! Get SCP basis functions (from AUX_BASIS)  corresponding to atom_A
       NULLIFY(first_sgf,l_max,l_min,npgf,nsgf_set)
       NULLIFY(rpgf,sphi,zet)
       CALL get_gto_basis_set(gto_basis_set=aux_basis,&
            first_sgf=first_sgf,&
            lmax=l_max,&
            lmin=l_min,&
            maxco=maxco,&
            npgf=npgf,&
            nset=nset,&
            nsgf_set=nsgf_set,&
            pgf_radius=rpgf,&
            sphi=sphi,&
            zet=zet)
       ncotot = maxco * nset

       ! Get the number of particles of ikind ( local )
       nparticle_local = local_particles%n_el(ikind)
       DO iparticle_local = 1, nparticle_local
          iatom = local_particles%list(ikind)%array(iparticle_local)
          offseta = 0._dp
          DO iseta = 1, nset
             ncoa = npgf ( iseta ) * ncoset ( l_max ( iseta ) )
             sgfa = first_sgf ( 1, iseta )
             nza = npgf ( iseta )
             rpgfa => rpgf ( 1:nza, iseta )
             zeta => zet ( 1 : nza, iseta )
             ! Allcating the arrays for the integrals
             IF ( ncoa > ldvab1 ) THEN
                CALL reallocate ( vab, 1, ncoa, 1, 1 )
                ldvab1 = ncoa
             ELSE
                vab = 0._dp 
             ENDIF
             DO i=1, nsgf_set ( iseta )
                work ( i, 1 ) = local_coeffs % c ( iparticle_local, offseta+i )
             END DO

             CALL dgemm ( "N", "N", ncoa, 1, nsgf_set ( iseta ), &
                  1.0_dp, sphi ( 1, sgfa ), SIZE ( sphi, 1 ), &
                  work ( 1, 1 ), SIZE ( work, 1 ), &
                  0.0_dp, gcca ( 1, 1 ), SIZE ( gcca, 1 ) )

             rpgfcore ( 1 ) = 0._dp ! rpgf is only used for screening. Not applicable for self-terms

             CALL coulomb2( l_max(iseta),npgf(iseta),zeta,rpgfa,l_min(iseta), &
                  0, 1, zetcore, rpgfa, 0, & 
                  (/0.0_dp,0.0_dp,0.0_dp/),0.0_dp,vab,vv,ff (0:))

             work = 0._dp
             DO icoa = 1,  ncoa ! There might be a stride... Need to check
                e_scp = e_scp - gcca (icoa, 1 ) * (-1._dp ) * gccore * vab ( icoa, 1 )
                work ( icoa, 1 )  =  work ( icoa, 1 ) + ( -1._dp ) * gccore * vab ( icoa, 1 )
             END DO
             ! Forces on the SCP coefficients
             IF ( .NOT. energy_only ) THEN
                CALL dgemm("T","N",nsgf_set ( iseta ),1,ncoa,&
                     1.0_dp,sphi(1,sgfa),SIZE(sphi,1),&
                     work(1,1),SIZE(work,1),&
                     0.0_dp,fcoeff(1,1),SIZE(fcoeff,1))
                DO i = 1, nsgf_set ( iseta )
                   local_coeffs%fc(iparticle_local,offseta+i) = &
                        local_coeffs%fc(iparticle_local,offseta+i) + fcoeff(i,1)
                   IF ( my_debug ) WRITE ( *, * ) 'FCOEFF', i, fcoeff ( i, 1 )
                END DO
             END IF ! energy  only
             offseta = offseta + nsgf_set ( iseta )
          END DO ! iseta
       END DO ! iparticle

    END DO  ! ikind

    !   *** sum up in the energy of hartree local from 2 centers terms
    CALL mp_sum(e_scp,para_env%group)
    energy % e_scp_core = e_scp

    DEALLOCATE (zetcore,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "zetcore")
    DEALLOCATE (rpgfcore,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "rpgfcore")
    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "work")
    DEALLOCATE (gcca,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "gcca")
    DEALLOCATE (fcoeff,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "fcoeff")
    DEALLOCATE(vab,vv,ff,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "vab.vv.ff")
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE integrate_vhscp_gcore

  !******************************************************************************
END MODULE scp_hartree_1center
!******************************************************************************
