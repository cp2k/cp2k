!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Space Group Symmetry Module  (version 1.0, January 16, 2020)
!> \par History
!>      Pierre-André Cazade [pcazade] 01.2020 - University of Limerick
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
MODULE space_groups

   USE bibliography,                    ONLY: Togo2018,&
                                              cite_reference
   USE kinds,                           ONLY: dp, sp
   USE spglib_f08,                      ONLY: &
        spg_get_international, spg_get_ir_reciprocal_mesh, spg_get_major_version, &
        spg_get_micro_version, spg_get_minor_version, spg_get_multiplicity, spg_get_pointgroup, &
        spg_get_schoenflies, spg_get_symmetry
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE gopt_f_types,                    ONLY: gopt_f_type
   USE cell_types,                      ONLY: cell_type, cell_copy,&
                                              cell_create, init_cell,&
                                              cell_release,&
                                              real_to_scaled, scaled_to_real
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE input_constants,                 ONLY: default_cell_direct_id,&
                                              default_cell_geo_opt_id,&
                                              default_cell_md_id,&
                                              default_cell_method_id,&
                                              default_minimization_method_id,&
                                              default_shellcore_method_id,&
                                              default_ts_method_id
   USE cp_fm_types,                     ONLY: cp_fm_type, cp_fm_get_info
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE message_passing,                 ONLY: mp_sum
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE mathlib,                         ONLY: det_3x3,&
                                              inv_3x3
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'space_groups'
   
   TYPE spgr_type
      LOGICAL                                     :: symlib = .FALSE.
      INTEGER                                     :: plevel = 0
      INTEGER                                     :: punit = -1
      INTEGER                                     :: istriz = -1
      REAL(KIND=dp)                               :: delta = 1.0e-8_dp
      INTEGER                                     :: nat
      INTEGER, DIMENSION(:), ALLOCATABLE          :: atype
      INTEGER, DIMENSION(:,:), ALLOCATABLE        :: eqatom
      !SPGLIB
      CHARACTER(len=11)                           :: international_symbol
      CHARACTER(len=6)                            :: pointgroup_symbol
      CHARACTER(len=7)                            :: schoenflies
      INTEGER                                     :: n_operations
      INTEGER, DIMENSION(:, :, :), ALLOCATABLE    :: rotations
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: translations
      TYPE(cell_type),POINTER                     :: cell_ref
   END TYPE spgr_type
   
   PUBLIC :: release_spgr_type, spgr_gen, identify_space_group, spgr_find_equivalent_atoms
   PUBLIC :: spgr_apply_rotations_coord, spgr_apply_rotations_force, print_spgr, spgr_type
   PUBLIC :: spgr_print_cell, spgr_print_xg, spgr_print_hess, spgr_print_hess_sp
   PUBLIC :: spgr_apply_rotations_hess

CONTAINS

! **************************************************************************************************
!> \brief Release the SPGR type
!> \param spgr  The SPGR type
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
   SUBROUTINE release_spgr_type(spgr)
      
      IMPLICIT NONE
      
      TYPE(spgr_type)                                    :: spgr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'release_spgr_type', &
         routineP = moduleN//':'//routineN

      IF (ALLOCATED(spgr%rotations)) THEN
         DEALLOCATE (spgr%rotations)
      END IF
      IF (ALLOCATED(spgr%translations)) THEN
         DEALLOCATE (spgr%translations)
      END IF
      IF (ALLOCATED(spgr%atype)) THEN
         DEALLOCATE (spgr%atype)
      END IF
      IF (ALLOCATED(spgr%eqatom)) THEN
         DEALLOCATE (spgr%eqatom)
      END IF

      CALL cell_release(spgr%cell_ref)

   END SUBROUTINE release_spgr_type

! **************************************************************************************************
!> \brief routine creates the space group structure
!> \param spgr ...
!> \param scoor ...
!> \param types ...
!> \param cell ...
!> \param gopt_env ...
!> \param delta ...
!> \param iounit ...
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
   SUBROUTINE spgr_gen(spgr, scoor, types, cell, gopt_env, delta, eps, iounit)
      
      IMPLICIT NONE
      
      TYPE(spgr_type)                                    :: spgr
      TYPE(gopt_f_type), POINTER                         :: gopt_env
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: scoor
      INTEGER, DIMENSION(:), INTENT(IN)                  :: types
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: delta, eps
      INTEGER, INTENT(IN), OPTIONAL                      :: iounit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'spgr_gen', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ierr, major, micro, minor, nat, &
                                                            nop, tra_mat(3, 3)
      LOGICAL                                            :: spglib

      CALL timeset(routineN, handle)

      !..total number of atoms
      nat = SIZE(scoor, 2)
      spgr%nat = nat
      !write(6,*)"SPGR| NUMBER OF ATOMS:",routineN,nat
      ! output unit
      IF (PRESENT(iounit)) THEN
         spgr%punit = iounit
      ELSE
         spgr%punit = 6
      END IF

      ! accuracy for symmetry
      IF (PRESENT(delta)) THEN
         spgr%delta = delta
      ELSE
         spgr%delta = 1.e-4_dp
      END IF

      !..set cell values
      NULLIFY(spgr%cell_ref)
      CALL cell_create(spgr%cell_ref)
      CALL cell_copy(cell, spgr%cell_ref)
      SELECT CASE (gopt_env%type_id)
      CASE (default_minimization_method_id, default_ts_method_id)
         CALL init_cell(spgr%cell_ref, hmat=cell%hmat)
      CASE (default_cell_method_id)
         SELECT CASE (gopt_env%cell_method_id)
         CASE (default_cell_direct_id)
            CALL init_cell(spgr%cell_ref, hmat=gopt_env%h_ref)
         CASE (default_cell_geo_opt_id, default_cell_md_id)
            CPABORT("SPACE_GROUP_SYMMETRY SHOULD NOT BE INVOKED DURING THE CELL STEP.")
         CASE DEFAULT
            CPABORT("")
         END SELECT
      CASE DEFAULT
         CPABORT("")
      END SELECT

      ! atom types
      ALLOCATE (spgr%atype(nat))
      spgr%atype(1:nat) = types(1:nat)

      spgr%n_operations = 0

      !..try spglib
      major = spg_get_major_version()
      minor = spg_get_minor_version()
      micro = spg_get_micro_version()
      IF (major == 0) THEN
         CALL cp_warn(__LOCATION__, "Symmetry library SPGLIB not available")
         spglib = .FALSE.
      ELSE
         !IF (major /= 1 .OR. minor /= 9 .OR. micro /= 9) THEN
         !   CALL cp_warn(__LOCATION__, "Version of Symmetry Library SPGLIB not tested")
         !END IF
         spglib = .TRUE.
         CALL cite_reference(Togo2018)
         spgr%international_symbol = '           '
         ierr = spg_get_international(spgr%international_symbol, TRANSPOSE(cell%hmat), scoor, types, nat, eps)
         IF (ierr == 0) THEN
            CALL cp_warn(__LOCATION__, "Symmetry Library SPGLIB failed")
            spglib = .FALSE.
         ELSE
            nop = spg_get_multiplicity(TRANSPOSE(cell%hmat), scoor, types, nat, eps)
            ALLOCATE (spgr%rotations(3, 3, nop), spgr%translations(3, nop))
            ALLOCATE (spgr%eqatom(nop,nat))
            spgr%n_operations = nop
            ierr = spg_get_symmetry(spgr%rotations, spgr%translations, nop, &
                                    TRANSPOSE(cell%hmat), scoor, types, nat, eps)
            ! Schoenflies Symbol
            spgr%schoenflies = '       '
            ierr = spg_get_schoenflies(spgr%schoenflies, TRANSPOSE(cell%hmat), scoor, types, nat, eps)
            ! Point Group
            spgr%pointgroup_symbol = '      '
            tra_mat = 0
            ierr = spg_get_pointgroup(spgr%pointgroup_symbol, tra_mat, &
                                      spgr%rotations, spgr%n_operations)
         END IF
      END IF
      spgr%symlib = spglib

      CALL timestop(handle)

   END SUBROUTINE spgr_gen
   
! **************************************************************************************************
!> \brief routine indentifies the space group and finds rotation matrices.
!> \param spgr ...
!> \param particle_set ...
!> \param cell ...
!> \param geo_section ...
!> \param gopt_env ...
!> \param iounit ...
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
!> \note  rotation matrices innclude translations and translation symmetry:
!>        it works with supercells as well.
! **************************************************************************************************
   SUBROUTINE identify_space_group(spgr,particle_set, cell, geo_section, gopt_env, iounit)
    
      IMPLICIT NONE
      
      TYPE(section_vals_type), POINTER                   :: geo_section
      TYPE(gopt_f_type), POINTER                         :: gopt_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell

      CHARACTER(LEN=*), PARAMETER :: routineN = 'identify_space_group', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, iounit, natom
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atype
      REAL(KIND=dp)                                      :: eps_geo, eps
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: scoord
      TYPE(spgr_type)                                    :: spgr

      CALL timeset(routineN, handle)
      
      natom = SIZE(particle_set)
      ALLOCATE (scoord(3, natom), atype(natom))
      DO i = 1, natom
         CALL get_atomic_kind(atomic_kind=particle_set(i)%atomic_kind, kind_number=atype(i))
         CALL real_to_scaled(scoord(1:3, i), particle_set(i)%r(1:3), cell)
      END DO
      CALL section_vals_val_get(geo_section, "EPS_GEO", r_val=eps_geo)
      CALL section_vals_val_get(geo_section, "EPS_SPGR", r_val=eps)
      CALL spgr_gen(spgr, scoord, atype, cell, gopt_env, delta=eps_geo, eps=eps, iounit=iounit)
      CALL spgr_find_equivalent_atoms(spgr, scoord)
      DEALLOCATE(scoord, atype)
      
      CALL timestop(handle)
      
   END SUBROUTINE identify_space_group
   
! **************************************************************************************************
!> \brief routine indentifies the equivalent atoms for each rotation matrix.
!> \param spgr ...
!> \param scoord ...
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
   SUBROUTINE spgr_find_equivalent_atoms(spgr, scoord)
      
      IMPLICIT NONE
      
      TYPE(spgr_type)                                    :: spgr
      
      CHARACTER(LEN=*), PARAMETER :: routineN = 'spgr_find_equivalent_atoms', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ia, ib, ir, natom, nop, handle !, i
      REAL(KIND=dp)                                      :: diff
      REAL(KIND=dp), DIMENSION(3)                        :: rb, ri, ro, tr
      REAL(KIND=dp), DIMENSION(3, 3)                     :: rot
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: scoord
      
      CALL timeset(routineN, handle)
      
      nop=spgr%n_operations
      natom = spgr%nat

      DO ia = 1, natom
         spgr%eqatom(1,ia) = ia
         spgr%eqatom(2,ia) = 1
      END DO

      !$OMP PARALLEL DO PRIVATE (ia,ib,ir,ri,rb,ro,rot,tr,diff) SHARED (spgr,scoord,natom,nop)
      DO ia = 1, natom
         ri(1:3) = scoord(1:3, ia)
         !ri(1) = ri(1) - FLOOR(ri(1))
         !ri(2) = ri(2) - FLOOR(ri(2))
         !ri(3) = ri(3) - FLOOR(ri(3))
         DO ir = 1, nop
            rot(1:3, 1:3) = spgr%rotations(1:3, 1:3, ir)
            tr(1:3) = spgr%translations(1:3, ir)
            DO ib = 1, natom
               rb(1:3) = scoord(1:3, ib)
               ro(1) = REAL(rot(1, 1), dp)*rb(1)+REAL(rot(2, 1), dp)*rb(2)+REAL(rot(3, 1), dp)*rb(3)+tr(1)
               ro(2) = REAL(rot(1, 2), dp)*rb(1)+REAL(rot(2, 2), dp)*rb(2)+REAL(rot(3, 2), dp)*rb(3)+tr(2)
               ro(3) = REAL(rot(1, 3), dp)*rb(1)+REAL(rot(2, 3), dp)*rb(2)+REAL(rot(3, 3), dp)*rb(3)+tr(3)
               ro(1) = ro(1) - REAL(NINT(ro(1)-ri(1)))
               ro(2) = ro(2) - REAL(NINT(ro(2)-ri(2)))
               ro(3) = ro(3) - REAL(NINT(ro(3)-ri(3)))
               !ro(1) = ro(1) - FLOOR(ro(1))
               !ro(2) = ro(2) - FLOOR(ro(2))
               !ro(3) = ro(3) - FLOOR(ro(3))
               diff = SQRT(SUM((ri(:)-ro(:))**2))
               IF ( (diff < spgr%delta) .AND. (spgr%atype(ia)==spgr%atype(ib)) ) THEN
                  !IF (ib < spgr%eqatom(1,ia)) THEN
                     spgr%eqatom(ir,ia) = ib
                     !spgr%eqatom(2,ia) = ir
                     EXIT
                  !END IF
               END IF
            END DO
         END DO
      END DO
      !$OMP END PARALLEL DO

      CALL timestop(handle)
      
   END SUBROUTINE spgr_find_equivalent_atoms

! **************************************************************************************************
!> \brief routine applies the rotation matrices to the coordinates.
!> \param spgr ...
!> \param coord ...
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
   SUBROUTINE spgr_apply_rotations_coord(spgr, coord)
      
      IMPLICIT NONE
      
      TYPE(spgr_type)                                    :: spgr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'spgr_apply_rotations_coord', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ia, ib, ir, natom, nop, handle, &
                                                            ja, jb
      REAL(KIND=dp), DIMENSION(3)                        :: ri, ro, tr, rf, rn
      REAL(KIND=dp), DIMENSION(3, 3)                     :: rot
      REAL(KIND=dp), DIMENSION(:)                        :: coord
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: cold

      INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
      
      CALL timeset(routineN, handle)

      allocate(cold(size(coord)))
      cold=coord
      nop=spgr%n_operations
      natom = spgr%nat
      
      !$OMP PARALLEL DO PRIVATE (ia,ib,ja,jb,ir,ri,ro,rf,rn,rot,tr) SHARED (spgr,coord,natom,nop)
      DO ia = 1, natom
         ja = 3*(ia-1)
         CALL real_to_scaled(rf(1:3), coord(ja+1:ja+3), spgr%cell_ref)
         rn(1:3) = 0.d0
         DO ir = 1, nop
            ib = spgr%eqatom(ir, ia)
            rot(1:3, 1:3) = spgr%rotations(1:3, 1:3, ir)
            tr(1:3) = spgr%translations(1:3, ir)
            jb = 3*(ib-1)
            CALL real_to_scaled(ri(1:3), coord(jb+1:jb+3), spgr%cell_ref)
            ro(1) = REAL(rot(1, 1), dp)*ri(1)+REAL(rot(2, 1), dp)*ri(2)+REAL(rot(3, 1), dp)*ri(3)+tr(1)
            ro(2) = REAL(rot(1, 2), dp)*ri(1)+REAL(rot(2, 2), dp)*ri(2)+REAL(rot(3, 2), dp)*ri(3)+tr(2)
            ro(3) = REAL(rot(1, 3), dp)*ri(1)+REAL(rot(2, 3), dp)*ri(2)+REAL(rot(3, 3), dp)*ri(3)+tr(3)
            ro(1) = ro(1) - REAL(NINT(ro(1)-rf(1)))
            ro(2) = ro(2) - REAL(NINT(ro(2)-rf(2)))
            ro(3) = ro(3) - REAL(NINT(ro(3)-rf(3)))
            rn(1:3) = rn(1:3) + ro(1:3)
         END DO
         rn = rn / REAL(nop)
         CALL scaled_to_real(coord(ja+1:ja+3), rn(1:3), spgr%cell_ref)
      END DO
      !$OMP END PARALLEL DO

      deallocate(cold)
      CALL timestop(handle)
   
   END SUBROUTINE spgr_apply_rotations_coord

! **************************************************************************************************
!> \brief routine applies the rotation matrices to the forces.
!> \param spgr ...
!> \param force ...
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
   SUBROUTINE spgr_apply_rotations_force(spgr, force)
      
      IMPLICIT NONE
      
      TYPE(spgr_type)                                    :: spgr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'spgr_apply_rotations_force', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ia, ib, ir, natom, nop, handle, &
                                                            ja, jb
      REAL(KIND=dp), DIMENSION(3)                        :: ri, ro, rn
      REAL(KIND=dp), DIMENSION(3, 3)                     :: rot
      REAL(KIND=dp), DIMENSION(:)                        :: force
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: fold
      
      CALL timeset(routineN, handle)

      allocate(fold(size(force)))
      fold=force
      nop=spgr%n_operations
      natom = spgr%nat

      !$OMP PARALLEL DO PRIVATE (ia,ib,ja,jb,ir,ri,ro,rn,rot) SHARED (spgr,force,natom,nop)
      DO ia = 1, natom
         ja = 3*(ia-1)
         rn(1:3) = 0.d0
         DO ir=1,nop
            ib = spgr%eqatom(ir, ia)
            rot(1:3, 1:3) = spgr%rotations(1:3, 1:3, ir)
            jb = 3*(ib-1)
            CALL real_to_scaled(ri(1:3), force(jb+1:jb+3), spgr%cell_ref)
            ro(1) = REAL(rot(1, 1), dp)*ri(1)+REAL(rot(2, 1), dp)*ri(2)+REAL(rot(3, 1), dp)*ri(3)
            ro(2) = REAL(rot(1, 2), dp)*ri(1)+REAL(rot(2, 2), dp)*ri(2)+REAL(rot(3, 2), dp)*ri(3)
            ro(3) = REAL(rot(1, 3), dp)*ri(1)+REAL(rot(2, 3), dp)*ri(2)+REAL(rot(3, 3), dp)*ri(3)
            rn(1:3) = rn(1:3) + ro(1:3)
         END DO
         rn = rn / REAL(nop)
         CALL scaled_to_real(force(ja+1:ja+3), rn(1:3), spgr%cell_ref)
      END DO
      !$OMP END PARALLEL DO

      deallocate(fold)
      CALL timestop(handle)
   
   END SUBROUTINE spgr_apply_rotations_force

! **************************************************************************************************
!> \brief routine applies the rotation matrices to the hessian.
!> \param ndf ...
!> \param spgr ...
!> \param hess_mat ...
!> \param para_env ...
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
!> \note  Not currently in use and most likely unnecessary. Will be removed in later version.
! **************************************************************************************************
   SUBROUTINE spgr_apply_rotations_hess(ndf,spgr,hess_mat,para_env)

      IMPLICIT NONE

      INTEGER, INTENT(IN)                                :: ndf
      TYPE(spgr_type)                                    :: spgr
      TYPE(cp_fm_type), POINTER                          :: hess_mat
      TYPE(cp_para_env_type), OPTIONAL, POINTER          :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'pgr_apply_rotations_hess',&
              routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, j, k, l, ncol_local, &
                                                            nrow_local, ia, ib, natom, &
                                                            nop, ja, jb, ka, kb, la, lb, &
                                                            iri, irk
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp), DIMENSION(3, 3)                     :: rot, roti, rotk, r, s, t
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:)         :: work
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: local_hes


      CALL timeset(routineN,handle)

      ALLOCATE(work(ndf,ndf))

      CALL cp_fm_get_info(hess_mat,row_indices=row_indices,col_indices=col_indices, &
           local_data=local_hes,nrow_local=nrow_local,ncol_local=ncol_local)

      work=0.0_dp
      DO i=1,nrow_local
         j=row_indices(i)
         DO k=1,ncol_local
            l=col_indices(k)
            work(j,l)=local_hes(i,k)
         END DO
      END DO

      CALL mp_sum(work,para_env%group)

      nop=spgr%n_operations
      natom = spgr%nat
      DO ia = 1, natom
         ib = spgr%eqatom(1, ia)
         iri = spgr%eqatom(2, ia)
         roti(1:3, 1:3) = spgr%rotations(1:3, 1:3, iri)
         ja = 3*(ia-1)
         jb = 3*(ib-1)
         DO ka = ia, natom
            r=0.0_dp
            kb = spgr%eqatom(1, ka)
            IF (ib==ia .AND. kb==ka) THEN
               CYCLE
            END IF
            irk = spgr%eqatom(2, ka)
            rotk(1:3, 1:3) = spgr%rotations(1:3, 1:3, irk)
            la = 3*(ka-1)
            lb = 3*(kb-1)
            rot=TRANSPOSE(MATMUL(rotk,roti))
            r(:,:) = work(jb+1:jb+3,lb+1:lb+3)
            s = MATMUL(spgr%cell_ref%h_inv,r)
            t = MATMUL(rot,s)
            r = MATMUL(spgr%cell_ref%hmat,t)
            work(ja+1:ja+3,la+1:la+3) = r(:,:)
            work(la+1:la+3,ja+1:ja+3) = TRANSPOSE(r)
         END DO
      END DO
      
      DO i=1,nrow_local
         j=row_indices(i)
         DO k=1,ncol_local
            l=col_indices(k)
            local_hes(i,k) = work(j,l)
         END DO
      END DO

      DEALLOCATE(work)

      CALL timestop(handle)
      
   END SUBROUTINE spgr_apply_rotations_hess

! **************************************************************************************************
!> \brief routine prints Space Group Information.
!> \param spgr ...
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
   SUBROUTINE print_spgr(spgr)
      
      IMPLICIT NONE
      
      TYPE(spgr_type)                                    :: spgr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_spgr', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, iunit, j, plevel

      iunit = spgr%punit
      IF (iunit >= 0) THEN
         plevel = spgr%plevel+5
         write(iunit,*)"output unit:",iunit
         WRITE (iunit, '(/,T2,A,A)') "----------------------------------------", &
            "---------------------------------------"
         WRITE (iunit, "(T2,A,T25,A,T77,A)") "----", "SPACE GROUP SYMMETRY INFORMATION", "----"
         WRITE (iunit, '(T2,A,A)') "----------------------------------------", &
            "---------------------------------------"
         IF (spgr%symlib) THEN
            CALL spgr_print_cell(spgr%cell_ref,0,iunit)
            WRITE (iunit, *) "SPGR| INTERNATIONAL SYMBOL: ", TRIM(ADJUSTR(spgr%international_symbol))
            WRITE (iunit, *) "SPGR| POINT GROUP SYMBOL  : ", TRIM(ADJUSTR(spgr%pointgroup_symbol))
            WRITE (iunit, *) "SPGR| SCHOENFLIES SYMBOL  : ", TRIM(ADJUSTR(spgr%schoenflies))
            WRITE (iunit, *) "SPGR| NUMBER OF SYMMETRY OPERATIONS:", spgr%n_operations
            IF (plevel > 10) THEN
               WRITE (iunit, '(/,T2,A,A)') "----------------------------------------", &
               "---------------------------------------"
               WRITE (iunit, '(T2,A,T34,A,T77,A)') "----", "EQUIVALENT ATOMS", "----"
               WRITE (iunit, '(T2,A,A)') "----------------------------------------", &
               "---------------------------------------"
               DO i = 1, spgr%nat
                  DO j = 1, spgr%n_operations
                     WRITE (iunit, '(T36,I8,I8,I8)') i, j, spgr%eqatom(j,i)
                  END DO
               END DO
               WRITE (iunit, '(T2,A,A)') "----------------------------------------", &
               "---------------------------------------"
               DO i = 1, spgr%n_operations
                  WRITE (iunit, '(A,i4,T51,3I10,/,T51,3I10,/,T51,3I10)') &
                     "           Rotation #: ", i, (spgr%rotations(j, :, i), j=1, 3)
                  WRITE (iunit, '(T36,3F15.7)') spgr%translations(:, i)
               END DO
            END IF
         ELSE
            WRITE (iunit, "(T2,A)") "SPGLIB for Crystal Symmetry Information determination is not availale"
         END IF
      END IF

   END SUBROUTINE print_spgr

! **************************************************************************************************
!> \brief routine prints cell details. Only for debugging purpose.
!> \param cell ...
!> \param iter ...
!> \param io output unit
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
   SUBROUTINE spgr_print_cell(cell,iter,io)
      IMPLICIT NONE
      
      INTEGER                                            :: iter,io
      TYPE(cell_type)                                    :: cell

      WRITE (io, *) "SPGR| ITERATION: ",iter
      WRITE (io, *) "SPGR| CELL_REF| A: ", cell%hmat(:,1)*0.52917721067
      WRITE (io, *) "SPGR| CELL_REF| B: ", cell%hmat(:,2)*0.52917721067
      WRITE (io, *) "SPGR| CELL_REF| C: ", cell%hmat(:,3)*0.52917721067

   END SUBROUTINE spgr_print_cell

! **************************************************************************************************
!> \brief routine prints coordinates or forces. Only for debugging purpose.
!> \param x ...
!> \param iter ...
!> \param io output unit
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
   SUBROUTINE spgr_print_xg(x,iter,io)
      IMPLICIT NONE

      INTEGER                                            :: iter,nat,i,ia,io
      REAL(KIND=dp), DIMENSION(:)                        :: x

      nat=size(x)/3
      write(io,*) nat
      write(io,*) "Iteration: ",iter
      do i = 1,nat
         ia = (i-1)*3
         write(io,*) x(ia+1)*0.52917721067,x(ia+2)*0.52917721067,x(ia+3)*0.52917721067
      end do

   END SUBROUTINE spgr_print_xg

! **************************************************************************************************
!> \brief routine prints hessian matrix. Only for debugging purpose.
!> \param hess ...
!> \param iter ...
!> \param io output unit
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
   SUBROUTINE spgr_print_hess(hess,iter,io)
      IMPLICIT NONE

      INTEGER                                            :: iter,nl,i,j,io
      REAL(KIND=dp), DIMENSION(:,:)                      :: hess

      nl=size(hess(1,:))
      write(io,*) nl/3
      write(io,*) "Iteration: ",iter
      do i = 1,nl/2
         write(io,*) (hess(j,i),j=1,nl/2)
      end do

      do i = nl/2+1,nl
         write(io+1,*) (hess(j,i),j=nl/2+1,nl)
      end do

      do i = 1,nl/2
         write(io+2,*) (hess(j,i),j=nl/2+1,nl)
      end do

      do i = nl/2+1,nl
         write(io+3,*) (hess(j,i),j=1,nl/2)
      end do

   END SUBROUTINE spgr_print_hess

! **************************************************************************************************
!> \brief routine prints single precision hessian matrix. Only for debugging purpose.
!> \param hess ...
!> \param iter ...
!> \param io output unit
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
   SUBROUTINE spgr_print_hess_sp(hess,iter,io)
      IMPLICIT NONE

      INTEGER                                            :: iter,nl,i,j,io
      REAL(KIND=sp), DIMENSION(:,:)                      :: hess

      nl=size(hess(1,:))
      write(io,*) nl/3
      write(io,*) "Iteration: ",iter
      do i = 1,nl/2
         write(io,*) (hess(j,i),j=1,nl/2)
      end do

      do i = nl/2+1,nl
         write(io+1,*) (hess(j,i),j=nl/2+1,nl)
      end do

      do i = 1,nl/2
         write(io+2,*) (hess(j,i),j=nl/2+1,nl)
      end do

      do i = nl/2+1,nl
         write(io+3,*) (hess(j,i),j=1,nl/2)
      end do

   END SUBROUTINE spgr_print_hess_sp
  
END MODULE space_groups
