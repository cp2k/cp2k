!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
!> \brief Perform the calculation of the dispersion contributions
!>        for the SCP method
!> \author CJM
!******************************************************************************

MODULE  scp_dispersion
  USE ai_coulomb,                      ONLY: coulomb2_new, coulomb2
  USE ai_derivatives,                  ONLY: dabdr
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp,&
                                             dp_size,int_size
  USE machine,                         ONLY: m_flush
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: ncoset
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_set_type, &
       neighbor_list_type, neighbor_node_type, next
  USE scp_coeff_types,                 ONLY: aux_coeff_set_type,&
                                             aux_coeff_type,&
                                             aux_coeff_zero_fc
  USE scp_energy_types,                ONLY: scp_energy_type
  USE scp_force_types,                 ONLY: scp_force_type
  USE scp_environment_types,           ONLY: get_scp_env,&
                                             scp_environment_type
  USE sparse_matrix_types,             ONLY: allocate_matrix_set,&
                                             create_dbm_from_sm,&
                                             deallocate_dbm_from_sm,&
                                             get_block_node,&
                                             real_matrix_p_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  ! *** Global parameters (only in this module)

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scp_dispersion'

  ! Public Subroutine

  PUBLIC :: scp_qs_dispersion

CONTAINS

! *****************************************************************************
!> \brief Analytic calculation of the dispersion energy and force
!>        for SCP using the QS environment
!> \author CJM
! *****************************************************************************

  SUBROUTINE scp_qs_dispersion ( qs_env, calculate_forces, error, debug )

    TYPE(qs_environment_type), POINTER                 :: qs_env
    LOGICAL, INTENT ( IN ), OPTIONAL                   :: calculate_forces
    TYPE(cp_error_type), INTENT(inout)                 :: error
    LOGICAL, INTENT ( IN ), OPTIONAL                   :: debug

    CHARACTER(LEN=*), PARAMETER :: routineN = "scp_qs_dispersion"

    TYPE(atomic_kind_type), DIMENSION(:), POINTER       :: atomic_kind_set
    
    TYPE(cell_type), POINTER                 :: cell
    TYPE(atomic_kind_type),       POINTER    :: atomic_kind
    TYPE(gto_basis_set_type),     POINTER    :: aux_basis
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(scp_energy_type), POINTER           :: energy
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(scp_force_type), POINTER            :: force(:)
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_scp
    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    TYPE(neighbor_list_type), POINTER        :: neighbor_list
    TYPE(neighbor_node_type), POINTER        :: neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(kind=dp)                            :: e_scp

    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE  :: vv
    REAL(dp), DIMENSION ( :, : ), POINTER    :: vab, vab_plus

    REAL(dp), DIMENSION ( :, : ), POINTER    :: sphia, sphib, dabdx, dabdy, dabdz
    REAL(dp), DIMENSION ( :, : ), POINTER    :: I_ab, I_dabdx, I_dabdy, I_dabdz, I_tmp
    REAL(dp), DIMENSION ( :, : ), POINTER    :: zeta, zetb, rpgfa, rpgfb
    REAL(dp), DIMENSION(:), POINTER          :: r1, r2, z1, z2
    REAL(dp), DIMENSION(:),   ALLOCATABLE    :: ff
    REAL(dp), DIMENSION(3)                   :: force_a 
    REAL(dp)                                 :: ra(3), rab(3), rab2, rb(3), dab, I_a, I_b, alpha_a, alpha_b
    REAL(dp)                                 :: disp_factor

    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, kind_of
    INTEGER :: ab, ilist, nlist, inode, nnode, maxder, natom
    INTEGER :: maxcoa, maxcob, nseta, nsetb
    INTEGER :: handle, iatom, jatom, i, istat, idir, i_atom, j_atom
    INTEGER :: ia, ib, ncoa, ncob, ncoap, ncobp
    INTEGER :: iparticle_local, nparticle_local, ikind, nkind, jkind
    INTEGER :: iseta, isetb, sgfa, sgfb
    INTEGER :: nz1, nz2, ldvab1, ldvab2
    INTEGER :: maxco_global, maxlgto, maxlgto1, output_unit
    INTEGER, DIMENSION ( : ), POINTER :: la_max, la_min, npgfa, nsgf_seta
    INTEGER, DIMENSION ( : ), POINTER :: lb_max, lb_min, npgfb, nsgf_setb
    INTEGER, DIMENSION ( :, : ), POINTER :: first_sgfa, first_sgfb

    LOGICAL :: ionode, my_debug, forces
    


!   ---------------------------------------------------------------------------
    CALL timeset(routineN,"I"," ",handle)
    output_unit = -1
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) output_unit= cp_logger_get_default_unit_nr(logger)
    my_debug = .FALSE.
    IF ( PRESENT ( debug ) ) my_debug = debug
    forces = .FALSE.
    IF ( PRESENT ( calculate_forces ) ) forces = calculate_forces

    NULLIFY( atomic_kind, aux_basis, scp_env)
    NULLIFY( energy, force, para_env, atomic_kind_set )
    NULLIFY( cell, sab_scp )

    CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,&
                    cell=cell, scp_env = scp_env, &
                    para_env=para_env,&
                    particle_set=particle_set,&
                    sab_scp=sab_scp,error=error)

    CALL get_scp_env(scp_env=scp_env, & 
                     energy=energy, force = force, error=error)
   
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco_global,maxlgto=maxlgto)

    nkind = SIZE(atomic_kind_set,1)
    natom = SIZE(particle_set,1)

    maxder = 0
    IF ( forces ) THEN
      maxder = 1

      ALLOCATE (atom_of_kind(natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
           "atom_of_kind",natom*int_size)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
            atom_of_kind=atom_of_kind)

      IF ( my_debug ) THEN
        ALLOCATE (kind_of(natom),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
             "kind_of",natom*int_size)
        CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                                 kind_of=kind_of)
      END IF

      CALL get_scp_env(scp_env=scp_env, force=force,error=error)

      maxlgto1 = maxlgto + maxder
    ELSE
      maxlgto1 = maxlgto
    END IF

    ALLOCATE( I_tmp ( maxco_global, maxco_global ), STAT=istat )
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "I_tmp",maxco_global*maxco_global*dp_size)
    ALLOCATE(vab(1,1),STAT = istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "vab",1*1*dp_size)
    ldvab1 = 1
    ldvab2 = 1
    vab = 0.0_dp

    ALLOCATE(vv(ncoset(maxlgto1),ncoset(maxlgto1),maxlgto1+maxlgto1+1),STAT = istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "vv",maxco_global*maxco_global*(maxlgto1+maxlgto+1)*dp_size)
    ALLOCATE(ff(0:maxlgto1+maxlgto1),STAT = istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "ff", (maxlgto1+maxlgto1+1)*dp_size)

    e_scp = 0._dp
    DO ikind = 1,nkind

      NULLIFY(atomic_kind)
      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,aux_basis_set=aux_basis, &
                           alpha_scp = alpha_a, I_scp = I_a )
                           

! Check to see if the atom is an SCP atom
      IF (.NOT.ASSOCIATED(aux_basis)) CYCLE
      NULLIFY(first_sgfa,la_max,la_min,npgfa,nsgf_seta)
      NULLIFY(rpgfa,sphia,zeta)
      CALL get_gto_basis_set(gto_basis_set=aux_basis,&
                            first_sgf=first_sgfa,&
                            lmax=la_max,&
                            lmin=la_min,&
                            maxco=maxcoa,&
                            npgf=npgfa,&
                            nset=nseta,&
                            nsgf_set=nsgf_seta,&
                            pgf_radius=rpgfa,&
                            sphi=sphia,&
                            zet=zeta)
      DO jkind = 1, nkind 
        NULLIFY(atomic_kind)
        atomic_kind => atomic_kind_set(jkind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,aux_basis_set=aux_basis, & 
                             alpha_scp = alpha_b, I_scp = I_b )
                           
! Computing the constants for the dispersion interaction
        disp_factor = -1._dp/8._dp * alpha_a * alpha_b * I_a * I_b / ( I_a + I_b )

! Check to see if the atom is an SCP atom
        IF (.NOT.ASSOCIATED(aux_basis)) CYCLE
        NULLIFY(first_sgfb,lb_max,lb_min,npgfb,nsgf_setb)
        NULLIFY(rpgfb,sphib,zetb)
        CALL get_gto_basis_set(gto_basis_set=aux_basis,&
                              first_sgf=first_sgfb,&
                              lmax=lb_max,&
                              lmin=lb_min,&
                              maxco=maxcob,&
                              npgf=npgfb,&
                              nset=nsetb,&
                              nsgf_set=nsgf_setb,&
                              pgf_radius=rpgfb,&
                              sphi=sphib,&
                              zet=zetb)

        ab = ikind + nkind*(jkind-1)

        IF( .NOT.ASSOCIATED(sab_scp(ab)%neighbor_list_set)) CYCLE

        neighbor_list_set => sab_scp(ab)%neighbor_list_set

        CALL get_neighbor_list_set(neighbor_list_set=neighbor_list_set,&
                                   nlist=nlist)

        neighbor_list => first_list(neighbor_list_set)

        DO ilist = 1,nlist

          CALL get_neighbor_list(neighbor_list=neighbor_list,&
                                 atom=iatom,&
                                 nnode=nnode)

          ra(:) = pbc(particle_set(iatom)%r(:),cell)

          neighbor_node => first_node(neighbor_list)

          DO inode=1,nnode

            CALL get_neighbor_node(neighbor_node=neighbor_node,&
                                   neighbor=jatom,&
                                   r=rab(:))

            rab2  = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)

            IF ( forces ) dab = SQRT ( rab2 )

            rb(:) = pbc(particle_set(jatom)%r,cell)

            IF ( .NOT. my_debug ) THEN
             rab ( : ) = ra ( : ) - rb ( : )
             rab2  = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
            END IF 

            IF ( iatom == jatom ) THEN
              neighbor_node => next(neighbor_node)
              CYCLE
            END IF
  

            DO iseta = 1, nseta
              ncoa = npgfa ( iseta ) * ncoset ( la_max ( iseta ) )
! Redifine ncoap if there are forces. This is needed to get the derivative
              IF ( forces ) ncoap = npgfa ( iseta ) * ncoset ( la_max ( iseta ) + 1 )
              sgfa = first_sgfa ( 1, iseta )
              nz1 = npgfa ( iseta )
              r1 => rpgfa ( 1:nz1, iseta )
              z1 => zeta ( 1 : nz1, iseta )
              DO isetb = 1, nsetb
                 ncob = npgfb ( isetb ) * ncoset ( lb_max ( isetb ) )
! Redifine ncobp if there are forces. This is needed to get the derivative
                 IF ( forces ) ncobp = npgfb ( isetb ) * ncoset ( lb_max ( isetb ) + 1 )
                 sgfb = first_sgfb ( 1, isetb )
                 nz2 = npgfb ( isetb )
                 r2 => rpgfb ( 1 : nz2, isetb )
                 z2 => zetb ( 1 : nz2, isetb )
! Allcating the arrays for the integrals
                 IF ( ncoa > ldvab1 .OR.  ncob > ldvab2 ) THEN
                   CALL reallocate ( vab , 1, ncoa, 1, ncob )
                   ldvab1 = ncoa
                   ldvab2 = ncob
                 ELSE
                   vab = 0._dp 
                 ENDIF 

                 ALLOCATE( I_ab ( nsgf_seta ( iseta ), nsgf_setb ( isetb ) ), STAT=istat )
                 IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                 "I_ab",nsgf_seta(iseta)*nsgf_setb(isetb)*dp_size)

                 IF ( forces ) THEN

                   ALLOCATE( I_dabdx ( nsgf_seta ( iseta ), nsgf_setb ( isetb ) ), STAT=istat )
                   IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                   "I_dabdx",nsgf_seta(iseta)*nsgf_setb(isetb)*dp_size)
 
                   ALLOCATE( I_dabdy ( nsgf_seta ( iseta ), nsgf_setb ( isetb ) ), STAT=istat )
                   IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                   "I_dabdy",nsgf_seta(iseta)*nsgf_setb(isetb)*dp_size)

                   ALLOCATE( I_dabdz ( nsgf_seta ( iseta ), nsgf_setb ( isetb ) ), STAT=istat )
                   IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                   "I_dabdz",nsgf_seta(iseta)*nsgf_setb(isetb)*dp_size)

                   ALLOCATE( vab_plus ( ncoap, ncobp ), STAT=istat )
                   IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                   "vab_plus",ncoap*ncobp*dp_size)

                   ALLOCATE(dabdx(ncoa,ncob),dabdy(ncoa,ncob),dabdz(ncoa,ncob),STAT=istat)
                   IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                   "dbdx.dbdy.dbdz",3*ncoa*ncob*dp_size)

                   CALL coulomb2_new( la_max ( iseta )+1, npgfa ( iseta ), z1, r1, 0, &
                                      lb_max ( isetb )+1, npgfb ( isetb ), z2, r2, 0, &
                                      rab, rab2, vab, vv, ff ( 0: ), maxder, vab_plus )

! Use vab_plus integral block to compute derivative
                   CALL dabdr ( la_max ( iseta ), npgfa ( iseta ), z1, r1, 0, &
                                lb_max ( isetb ), npgfb ( isetb ), z2, r2, 0, &
                                dab, vab_plus, dabdx, dabdy, dabdz )

! Contract over the cartesian orbital index to obtain integral
! block labeled by molecular orbitals
! vab->I_ab
                   CALL dgemm ("N","N", ncoa,nsgf_setb(isetb),ncob, &
                         1.0_dp,vab (1,1),SIZE (vab,1), &
                         sphib(1,sgfb),SIZE (sphib,1), &
                         0.0_dp, I_tmp(1,1), SIZE(I_tmp,1))
                   CALL dgemm("T","N",nsgf_seta (iseta),nsgf_setb(isetb),ncoa,&
                         1.0_dp,sphia(1,sgfa),SIZE(sphia,1),&
                         I_tmp(1,1),SIZE(I_tmp,1),&
                         0.0_dp,I_ab(1,1),SIZE(I_ab,1))
! dabdx->I_dabdx
                   CALL dgemm ("N","N", ncoa,nsgf_setb(isetb),ncob, &
                         1.0_dp,dabdx (1,1),SIZE (dabdx,1), &
                         sphib(1,sgfb),SIZE (sphib,1), &
                         0.0_dp, I_tmp(1,1), SIZE(I_tmp,1))
                   CALL dgemm("T","N",nsgf_seta (iseta),nsgf_setb(isetb),ncoa,&
                         1.0_dp,sphia(1,sgfa),SIZE(sphia,1),&
                         I_tmp(1,1),SIZE(I_tmp,1),&
                         0.0_dp,I_dabdx(1,1),SIZE(I_dabdx,1))
! dabdy->I_dabdy
                   CALL dgemm ("N","N", ncoa,nsgf_setb(isetb),ncob, &
                         1.0_dp,dabdy (1,1),SIZE (dabdy,1), &
                         sphib(1,sgfb),SIZE (sphib,1), &
                         0.0_dp, I_tmp(1,1), SIZE(I_tmp,1))
                   CALL dgemm("T","N",nsgf_seta (iseta),nsgf_setb(isetb),ncoa,&
                         1.0_dp,sphia(1,sgfa),SIZE(sphia,1),&
                         I_tmp(1,1),SIZE(I_tmp,1),&
                         0.0_dp,I_dabdy(1,1),SIZE(I_dabdy,1))
! dabdz->I_dabdz
                   CALL dgemm ("N","N", ncoa,nsgf_setb(isetb),ncob, &
                         1.0_dp,dabdz (1,1),SIZE (dabdz,1), &
                         sphib(1,sgfb),SIZE (sphib,1), &
                         0.0_dp, I_tmp(1,1), SIZE(I_tmp,1))
                   CALL dgemm("T","N",nsgf_seta (iseta),nsgf_setb(isetb),ncoa,&
                         1.0_dp,sphia(1,sgfa),SIZE(sphia,1),&
                         I_tmp(1,1),SIZE(I_tmp,1),&
                         0.0_dp,I_dabdz(1,1),SIZE(I_dabdz,1))
                   force_a = 0.0_dp
                   DO ia = 1, nsgf_seta ( iseta )
                     DO ib = 1, nsgf_setb ( isetb )
             ! **** need to omit (s|s) ****
                       e_scp = e_scp + disp_factor * I_ab ( ia, ib )**2 
                       force_a(1)=force_a(1)+2.0_dp*disp_factor*I_ab(ia,ib)*I_dabdx(ia,ib)
                       force_a(2)=force_a(2)+2.0_dp*disp_factor*I_ab(ia,ib)*I_dabdy(ia,ib)
                       force_a(3)=force_a(3)+2.0_dp*disp_factor*I_ab(ia,ib)*I_dabdz(ia,ib)
                     END DO
                   END DO
!             ***  sum up in the force arrey  ***
              ! NOTE: This is the actual force, e.g. -d/dR
                   i_atom = atom_of_kind(iatom)
                   j_atom = atom_of_kind(jatom)
                   DO idir = 1,3
                     force(ikind)%f_scp_disp(idir,i_atom) = &
                       force(ikind)%f_scp_disp(idir,i_atom) - force_a(idir)
                     force(jkind)%f_scp_disp(idir,j_atom) = &
                        force(jkind)%f_scp_disp(idir,j_atom) + force_a(idir)
                   END DO
                   
                   DEALLOCATE(dabdx,dabdy,dabdz,STAT=istat)
                   IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                   "dbdx.dbdy.dbdz")
                   DEALLOCATE(vab_plus, STAT=istat)
                   IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                   "vab_plus")
                   DEALLOCATE( I_dabdx, I_dabdy, I_dabdz, STAT=istat)
                   IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                   "I_dabdz.I_dabdy.I_dabdz")
                 ELSE
                   CALL coulomb2( la_max ( iseta ), npgfa ( iseta ), z1, r1, la_min ( iseta ), &
                                  lb_max ( isetb ), npgfb ( isetb ), z2, r2, lb_min ( isetb ), &
                                  rab, rab2, vab, vv, ff ( 0: ) )
! Contract over the cartesian orbital index to obtain integral
! block labeled by molecular orbitals
! vab->I_ab
 
                   CALL dgemm ("N","N", ncoa,nsgf_setb(isetb),ncob, &
                         1.0_dp,vab (1,1),SIZE (vab,1), &
                         sphib(1,sgfb),SIZE (sphib,1), &
                         0.0_dp, I_tmp(1,1), SIZE(I_tmp,1))
                   CALL dgemm("T","N",nsgf_seta (iseta),nsgf_setb(isetb),ncoa,&
                         1.0_dp,sphia(1,sgfa),SIZE(sphia,1),&
                         I_tmp(1,1),SIZE(I_tmp,1),&
                         0.0_dp,I_ab(1,1),SIZE(I_ab,1))

                   DO ia = 1, nsgf_seta ( iseta )
                     DO ib = 1, nsgf_setb ( isetb )
             ! **** need to omit (s|s) ****
                       e_scp= e_scp +  disp_factor * I_ab ( ia, ib )**2 
                     END DO
                   END DO

                 ENDIF

                 DEALLOCATE( I_ab, STAT=istat)
                 IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                 "I_ab")

               END DO ! isetb
             END DO ! iseta
             neighbor_node => next(neighbor_node)
           END DO ! inode
           neighbor_list => next(neighbor_list)
         END DO ! ilist 
       END DO ! jkind
     END DO  ! ikind

     IF ( my_debug .AND. ionode ) THEN
        DO iatom = 1, natom
          ikind = kind_of ( iatom )
          i_atom = atom_of_kind ( iatom )
          WRITE ( *, * ) 'F_DISPERSION', iatom,  force ( ikind ) % f_scp_disp ( :, i_atom )
        END DO
     END IF

     CALL mp_sum(e_scp,para_env%group)
     energy % e_scp_disp =  e_scp 

    ! print
    IF ( .NOT. my_debug ) THEN
      IF (output_unit>0) THEN
         WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
              "Post Processed SCP Dispersion                          ",energy % e_scp_disp

         CALL m_flush(output_unit)
      END IF
    END IF

   
    DEALLOCATE(vab,vv,ff,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
          "vab.vv.ff")

    DEALLOCATE(I_tmp,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
          "I_tmp")
    IF (forces) THEN
       DEALLOCATE(atom_of_kind,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
            "atom_of_kind")
       IF ( my_debug ) THEN
         DEALLOCATE(kind_of,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
            "kind_of")
       END IF
    END IF

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE scp_qs_dispersion

!******************************************************************************
END MODULE  scp_dispersion

