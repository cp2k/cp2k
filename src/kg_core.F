!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/KG_CORE[1.0] *
!!
!!   NAME
!!     kg_core
!!
!!   FUNCTION
!!     Calculation of the energies concerning the core charge distribution
!!   AUTHOR
!!     gt (14.10.2002)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE kg_core

! *****************************************************************************

  USE kinds, ONLY: int_size,&
                   wp => dp,&
                   wp_size => dp_size

  USE atomic_kind_types,        ONLY: atomic_kind_type,&
                                      get_atomic_kind,&
                                      get_atomic_kind_set
  USE global_types,             ONLY: LOW,global_environment_type
  USE kg_environment_types,     ONLY: kg_environment_type, &
                                      get_kg_env
  USE kg_force_types,           ONLY: kg_force_type 
  USE mathconstants,            ONLY: oorootpi,twopi
  USE message_passing,          ONLY: mp_sum
  USE particle_types,           ONLY: particle_type
  USE simulation_cell,          ONLY: cell_type,&
                                      get_cell,&
                                      pbc,&
                                      real_to_scaled,&
                                      scaled_to_real
  USE termination,              ONLY: stop_memory
  USE timings,                  ONLY: timeset,&
                                      timestop

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: calculate_ecore_rspace,&
            calculate_eself

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE calculate_ecore_rspace(kg_env,globenv,e_rspace_core,force)

!   Purpose: Calculate the overlap energy of the core charge distribution.
!   preliminar version:
!   adapted to KG from the corresponding qs routine - no linked lists!

!   History: - Creation (14.10.2002)

!   ***************************************************************************

    IMPLICIT NONE

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(kg_environment_type), INTENT(INOUT)  :: kg_env
    REAL(wp), intent (out) :: e_rspace_core
    TYPE(kg_force_type), DIMENSION(:), POINTER, OPTIONAL  :: force

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE calculate_ecore_rspace (MODULE kg_core)"

!   *** Local data types ***

    TYPE pbc_coord_type
      REAL(wp), DIMENSION(:,:), POINTER :: r,s
    END TYPE pbc_coord_type

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER :: atomic_kind
    TYPE(cell_type), POINTER        :: cell
    CHARACTER(LEN=8)                :: unit_of_length_name
    REAL(wp)                        :: a_max,a_min,aab,&
                                       b_max,b_min,&
                                       c_max,c_min,dab,eab,f0,&
                                       rab_max,rab2,rab2_max,rootaab,subcells,&
                                       unit_of_length,zab
    INTEGER                         :: atom_a,handle,iatom,icell,igrid,iijk,&
                                       ikind,istat,jatom,jcell,jgrid,jkind,&
                                       kcell,kgrid,maxatom,natom_a,nkind,&
                                       output_unit
    LOGICAL                         :: equal_kinds,ionode

    REAL(wp), DIMENSION(3) :: deab,r,ra_pbc,rab,rb,sa_pbc,sab,sab_max,sb,&
                              sb_max,sb_min,sb_pbc
    INTEGER, DIMENSION(3)  :: cell_b,ncell,ngrid,periodic

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), POINTER    :: part
    INTEGER, DIMENSION(:), POINTER                :: atom_list

    TYPE(pbc_coord_type), DIMENSION(:), ALLOCATABLE :: pbc_coord
    REAL(wp), DIMENSION(:), ALLOCATABLE             :: alpha,radius,qeff
    INTEGER, DIMENSION(:), ALLOCATABLE              :: natom
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE       :: grid_max,grid_min
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE          :: nijk
    INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE        :: ijk

    REAL(wp), EXTERNAL :: erfc

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_ecore_rspace","I","",handle)

    ionode = globenv%ionode
    output_unit = globenv%scr

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    box=cell,&
                    part=part) 

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxatom=maxatom)

    CALL get_cell(cell=cell,&
                  periodic=periodic,&
                  subcells=subcells,&
                  unit_of_length=unit_of_length,&
                  unit_of_length_name=unit_of_length_name)

!   *** Allocate work storage ***

    nkind = SIZE(atomic_kind_set)

    ALLOCATE (alpha(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"alpha",nkind*wp_size)
    alpha(:) = 0.0_wp

    ALLOCATE (radius(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"radius",nkind*wp_size)
    radius(:) = 0.0_wp

    ALLOCATE (qeff(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"qeff",nkind*int_size)
    qeff(:) = 0

    ALLOCATE (natom(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"natom",nkind*int_size)

    ALLOCATE (pbc_coord(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pbc_coord",nkind*int_size)

!   *** Calculate PBC coordinates ***

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      NULLIFY (atom_list)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           natom=natom(ikind),&
                           atom_list=atom_list,&
                           alpha_core_charge=alpha(ikind),&
                           core_charge_radius=radius(ikind),&
                           qeff=qeff(ikind))

      IF (alpha(ikind) == 0.0_wp) THEN

        NULLIFY (pbc_coord(ikind)%r)
        NULLIFY (pbc_coord(ikind)%s)

      ELSE

        ALLOCATE (pbc_coord(ikind)%r(3,natom(ikind)),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"pbc_coord(ikind)%r",&
                                         3*natom(ikind)*wp_size)

        ALLOCATE (pbc_coord(ikind)%s(3,natom(ikind)),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"pbc_coord(ikind)%s",&
                                         3*natom(ikind)*wp_size)


        DO iatom=1,natom(ikind)
          atom_a = atom_list(iatom)
          ra_pbc(:) = pbc(part(atom_a)%r(:),cell)
          pbc_coord(ikind)%r(:,iatom) = ra_pbc(:)
          pbc_coord(ikind)%s(:,iatom) = real_to_scaled(ra_pbc(:),cell)
        END DO

      END IF

    END DO

    IF ((ionode.AND.globenv%print%cell_parameters).AND.&
        (globenv%print%level > LOW)) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/,/,T3,A,T29,A,T54,A)")&
        "SUBCELL GRID FOR THE CALCULATION OF THE CORE CHARGE OVERLAP",&
        "Atomic kind pair","Grid size",&
        "Subcell size in "//unit_of_length_name
    END IF

    e_rspace_core = 0.0_wp

    DO ikind=1,nkind

      IF (alpha(ikind) == 0.0_wp) CYCLE

      DO jkind=ikind,nkind

        IF (alpha(jkind) == 0.0_wp) CYCLE

        equal_kinds = (ikind == jkind)

        zab = REAL(qeff(ikind)*qeff(jkind),wp)
        aab = alpha(ikind)*alpha(jkind)/(alpha(ikind) + alpha(jkind))
        rootaab = SQRT(aab)
        f0 = 2.0_wp*oorootpi*zab*rootaab

!       *** Calculate the square of the maximum interaction distance ***

        rab_max = radius(ikind) + radius(jkind)
        rab2_max = rab_max*rab_max

        r(:) = rab_max
        sab_max(:) = real_to_scaled(r(:),cell)

        ncell(:) = (INT(sab_max(:)) + 1)*periodic(:)
        ngrid(:) = MAX(1,NINT(0.5_wp*subcells/sab_max(:)))

!       *** Print subcell information for the current atomic kind pair ***

        IF ((ionode.AND.globenv%print%cell_parameters).AND.&
            (globenv%print%level > LOW)) THEN
          WRITE (UNIT=output_unit,FMT="(T3,2I8,4X,3I5,6X,3F12.6)")&
            ikind,jkind,ngrid,&
            scaled_to_real(1.0_wp/REAL(ngrid(:),wp),cell)/unit_of_length
        END IF

!       *** Check, if we have to consider a subcell grid ***

        IF (SUM(ngrid) == 3) THEN

          DO jatom=1,natom(jkind)

            sb_pbc(:) = pbc_coord(jkind)%s(:,jatom)

            DO icell=-ncell(1),ncell(1)
              cell_b(1) = icell
              DO jcell=-ncell(2),ncell(2)
                cell_b(2) = jcell
                DO kcell=-ncell(3),ncell(3)
                  cell_b(3) = kcell

                  sb(:) = sb_pbc(:) + REAL(cell_b(:),wp)
                  rb(:) = scaled_to_real(sb(:),cell)

                  IF (equal_kinds) THEN
                    natom_a = jatom - 1
                  ELSE
                    natom_a = natom(ikind)
                  END IF

                  DO iatom=1,natom_a
                    rab(:) = rb(:) - pbc_coord(ikind)%r(:,iatom)
                    rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                    IF (rab2 < rab2_max) THEN
                      dab = SQRT(rab2)
                      eab = zab*erfc(rootaab*dab)/dab
                      e_rspace_core = e_rspace_core + eab
                      IF (PRESENT (force)) THEN
                        deab(:) = rab(:)*(eab + f0*exp(-aab*rab2))/rab2
                        force(ikind)%f_rspace_core(:,iatom) =&
                          force(ikind)%f_rspace_core(:,iatom) - deab(:)
                        force(jkind)%f_rspace_core(:,jatom) =&
                          force(jkind)%f_rspace_core(:,jatom) + deab(:)
                      END IF
                    END IF
                  END DO

                END DO
              END DO
            END DO

          END DO

        ELSE

          ALLOCATE (grid_min(3,ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"grid_min",3*PRODUCT(ngrid)*wp_size)
          END IF

          ALLOCATE (grid_max(3,ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"grid_max",3*PRODUCT(ngrid)*wp_size)
          END IF

          ALLOCATE (nijk(ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"nijk",PRODUCT(ngrid)*int_size)
          END IF
          nijk(:,:,:) = 0

          ALLOCATE (ijk(natom(ikind),ngrid(1),ngrid(2),ngrid(3)),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"ijk",&
                             natom(ikind)*PRODUCT(ngrid)*int_size)
          END IF

          DO igrid=1,ngrid(1)
            a_min = REAL(igrid-1,wp)/REAL(ngrid(1),wp) - 0.5_wp
            a_max = REAL(igrid,wp)/REAL(ngrid(1),wp) - 0.5_wp
            DO jgrid=1,ngrid(2)
              b_min = REAL(jgrid-1,wp)/REAL(ngrid(2),wp) - 0.5_wp
               b_max = REAL(jgrid,wp)/REAL(ngrid(2),wp) - 0.5_wp
              DO kgrid=1,ngrid(3)
                c_min = REAL(kgrid-1,wp)/REAL(ngrid(3),wp) - 0.5_wp
                c_max = REAL(kgrid,wp)/REAL(ngrid(3),wp) - 0.5_wp
                grid_min(:,igrid,jgrid,kgrid) = (/a_min,b_min,c_min/)
                grid_max(:,igrid,jgrid,kgrid) = (/a_max,b_max,c_max/)
              END DO
            END DO
          END DO

          DO iatom=1,natom(ikind)
            sa_pbc(:) = pbc_coord(ikind)%s(:,iatom)
            igrid = MAX(1,CEILING((sa_pbc(1) + 0.5_wp)*ngrid(1)))
            jgrid = MAX(1,CEILING((sa_pbc(2) + 0.5_wp)*ngrid(2)))
            kgrid = MAX(1,CEILING((sa_pbc(3) + 0.5_wp)*ngrid(3)))
            nijk(igrid,jgrid,kgrid) = nijk(igrid,jgrid,kgrid) + 1
            ijk(nijk(igrid,jgrid,kgrid),igrid,jgrid,kgrid) = iatom
          END DO

          DO jatom=1,natom(jkind)

            sb_pbc(:) = pbc_coord(jkind)%s(:,jatom)

            DO icell=-ncell(1),ncell(1)
              cell_b(1) = icell
              DO jcell=-ncell(2),ncell(2)
                cell_b(2) = jcell
                DO kcell=-ncell(3),ncell(3)
                  cell_b(3) = kcell

                  sb(:) = sb_pbc(:) + REAL(cell_b(:),wp)
                  rb(:) = scaled_to_real(sb(:),cell)
                  sb_min(:) = sb(:) - sab_max(:)
                  sb_max(:) = sb(:) + sab_max(:)

                  IF (sb_max(1) < grid_min(1,1,1,1)) CYCLE
                  IF (sb_max(2) < grid_min(2,1,1,1)) CYCLE
                  IF (sb_max(3) < grid_min(3,1,1,1)) CYCLE

                  igrid = ngrid(1)
                  jgrid = ngrid(2)
                  kgrid = ngrid(3)

                  IF (sb_min(1) >= grid_max(1,igrid,jgrid,kgrid)) CYCLE
                  IF (sb_min(2) >= grid_max(2,igrid,jgrid,kgrid)) CYCLE
                  IF (sb_min(3) >= grid_max(3,igrid,jgrid,kgrid)) CYCLE

                  DO igrid=1,ngrid(1)
                    DO jgrid=1,ngrid(2)
                      DO kgrid=1,ngrid(3)

                        IF (nijk(igrid,jgrid,kgrid) == 0) CYCLE

                        IF (sb_max(1) < grid_min(1,igrid,jgrid,kgrid)) CYCLE
                        IF (sb_max(2) < grid_min(2,igrid,jgrid,kgrid)) CYCLE
                        IF (sb_max(3) < grid_min(3,igrid,jgrid,kgrid)) CYCLE

                        IF (sb_min(1) >= grid_max(1,igrid,jgrid,kgrid)) CYCLE
                        IF (sb_min(2) >= grid_max(2,igrid,jgrid,kgrid)) CYCLE
                        IF (sb_min(3) >= grid_max(3,igrid,jgrid,kgrid)) CYCLE

                        DO iijk=1,nijk(igrid,jgrid,kgrid)
                          iatom = ijk(iijk,igrid,jgrid,kgrid)
                          IF (equal_kinds) THEN
                            IF (jatom <= iatom) CYCLE
                          END IF
                          rab(:) = rb(:) - pbc_coord(ikind)%r(:,iatom)
                          rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                          IF (rab2 < rab2_max) THEN
                            dab = SQRT(rab2)
                            e_rspace_core = e_rspace_core +&
                                            zab*erfc(aab*dab)/dab
                            eab = zab*erfc(rootaab*dab)/dab
                            e_rspace_core = e_rspace_core + eab
                            IF (PRESENT (force)) THEN
                              deab(:) = rab(:)*(eab + f0*exp(-aab*rab2))/rab2
                              force(ikind)%f_rspace_core(:,iatom) =&
                                force(ikind)%f_rspace_core(:,iatom) - deab(:)
                              force(jkind)%f_rspace_core(:,jatom) =&
                                force(jkind)%f_rspace_core(:,jatom) + deab(:)
                            END IF
                          END IF
                        END DO

                      END DO
                    END DO
                  END DO

                END DO
              END DO
            END DO

          END DO

          DEALLOCATE (grid_min,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"grid_min")

          DEALLOCATE (grid_max,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"grid_max")

          DEALLOCATE (nijk,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"nijk")

          DEALLOCATE (ijk,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"ijk")

        END IF

      END DO
    END DO

!   *** Release work storage ***

    DEALLOCATE (alpha,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"alpha")

    DEALLOCATE (radius,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"radius")

    DEALLOCATE (qeff,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"qeff")

    DEALLOCATE (natom,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"natom")

    DO ikind=1,nkind
      IF (ASSOCIATED(pbc_coord(ikind)%r)) THEN
        DEALLOCATE (pbc_coord(ikind)%r,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"pbc_coord(ikind)%r")
      END IF
      IF (ASSOCIATED(pbc_coord(ikind)%s)) THEN
        DEALLOCATE (pbc_coord(ikind)%s,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"pbc_coord(ikind)%s")
      END IF
    END DO

    DEALLOCATE (pbc_coord,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pbc_coord")

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE calculate_ecore_rspace

! *****************************************************************************

  SUBROUTINE calculate_eself(kg_env,ecore_self)

!   Purpose: Calculate the self energy of the core charge distribution.

!   History: - Creation (14.10.2002) 

!   ***************************************************************************

    TYPE(kg_environment_type), INTENT(IN) :: kg_env
    REAL (wp), intent(OUT) :: ecore_self

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER :: atomic_kind
    REAL(wp)                        :: alpha_core_charge,qeff
    INTEGER                         :: ikind,natom

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set

!   ---------------------------------------------------------------------------

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set) 

    ecore_self = 0.0_wp

    DO ikind=1,SIZE(atomic_kind_set)

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           natom=natom,&
                           qeff=qeff,&
                           alpha_core_charge=alpha_core_charge)

      ecore_self = ecore_self - REAL(natom*qeff**2,wp)*SQRT(alpha_core_charge)

    END DO

    ecore_self = ecore_self/SQRT(twopi)

  END SUBROUTINE calculate_eself

! *****************************************************************************

END MODULE kg_core
