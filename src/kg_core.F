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
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE global_types,                    ONLY: LOW,&
                                             global_environment_type
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type
  USE kg_force_types,                  ONLY: kg_force_type
  USE kinds,                           ONLY: dbl,&
                                             int_size,&
                                             dp_size
  USE linklist_utilities,              ONLY: get_match
  USE mathconstants,                   ONLY: oorootpi,&
                                             twopi, pi
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell,&
                                             pbc,&
                                             real_to_scaled,&
                                             scaled_to_real
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
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
    REAL(dbl), INTENT (out) :: e_rspace_core
    TYPE(kg_force_type), DIMENSION(:), POINTER, OPTIONAL  :: force

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE calculate_ecore_rspace (MODULE kg_core)"

!   *** Local data types ***

    TYPE pbc_coord_type
      REAL(dbl), DIMENSION(:,:), POINTER :: r,s
    END TYPE pbc_coord_type
    TYPE alist_type
      INTEGER , DIMENSION(:), POINTER :: ipart
    END TYPE alist_type

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER :: atomic_kind
    TYPE(cell_type), POINTER        :: cell
    CHARACTER(LEN=8)                :: unit_of_length_name
    INTEGER, PARAMETER              :: els = 50
    REAL(dbl)                        :: a_max,a_min,aab,&
                                       b_max,b_min,&
                                       c_max,c_min,dab,eab,f0,&
                                       rab_max,rab2,rab2_max,rootaab,subcells,&
                                       unit_of_length,zab
    INTEGER                         :: atom_a,atom_b,handle,iatom,icell,igrid,iijk,&
                                       ikind,istat,jatom,jcell,jgrid,jkind,&
                                       kcell,kgrid,maxatom,natom_a,nkind,&
                                       nexcl,iexclude,output_unit
    LOGICAL                         :: equal_kinds,ionode

    REAL(dbl), DIMENSION(3) :: deab,r,ra_pbc,rab,rb,sa_pbc,sab,sab_max,sb,&
                              sb_max,sb_min,sb_pbc
    INTEGER, DIMENSION(3)  :: cell_b,ncell,ngrid,periodic

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set
    INTEGER, DIMENSION(:), POINTER                :: atom_list
    TYPE(pbc_coord_type), DIMENSION(:), ALLOCATABLE :: pbc_coord
    TYPE(alist_type), DIMENSION(:), ALLOCATABLE :: alist
    REAL(dbl), DIMENSION(:), ALLOCATABLE             :: alpha,radius,ccore_charge
    INTEGER, DIMENSION(:), ALLOCATABLE              :: natom
!    INTEGER, DIMENSION(:), ALLOCATABLE              :: elist
    INTEGER, DIMENSION(:), ALLOCATABLE              :: ind
    REAL(dbl), DIMENSION(:,:,:,:), ALLOCATABLE       :: grid_max,grid_min
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE          :: nijk
    INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE        :: ijk

    REAL(dbl), EXTERNAL :: erfc

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_ecore_rspace","I","",handle)

    ionode = globenv%ionode
    output_unit = globenv%scr

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    particle_set=particle_set) 

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
    IF (istat /= 0) CALL stop_memory(routine,"alpha",nkind*dp_size)
    alpha(:) = 0.0_dbl

    ALLOCATE (radius(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"radius",nkind*dp_size)
    radius(:) = 0.0_dbl

    ALLOCATE (ccore_charge(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"ccore_charge",nkind*int_size)
    ccore_charge(:) = 0

    ALLOCATE (natom(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"natom",nkind*int_size)

    ALLOCATE (alist(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"alist",nkind*int_size)

    ALLOCATE (pbc_coord(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pbc_coord",nkind*int_size)
! gt excl to be done
!    ALLOCATE (elist(els),STAT=istat)
!    IF ( istat /= 0 ) CALL stop_memory (routine,"elist",els*int_size )
! gt excl to be done

    ALLOCATE (ind(els),STAT=istat)
    IF ( istat /= 0 ) CALL stop_memory (routine,"ind",els*int_size )


!   *** Calculate PBC coordinates ***

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      NULLIFY (atom_list)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           natom=natom(ikind),&
                           atom_list=atom_list,&
                           alpha_core_charge=alpha(ikind),&
                           core_charge_radius=radius(ikind),&
                           ccore_charge=ccore_charge(ikind))

      IF (alpha(ikind) == 0.0_dbl) THEN

        NULLIFY (pbc_coord(ikind)%r)
        NULLIFY (pbc_coord(ikind)%s)
        NULLIFY (alist(ikind)%ipart)

      ELSE

        ALLOCATE (pbc_coord(ikind)%r(3,natom(ikind)),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"pbc_coord(ikind)%r",&
                                         3*natom(ikind)*dp_size)

        ALLOCATE (pbc_coord(ikind)%s(3,natom(ikind)),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"pbc_coord(ikind)%s",&
                                         3*natom(ikind)*dp_size)

        ALLOCATE (alist(ikind)%ipart(natom(ikind)),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"alist(ikind)%ipart",&
                                         natom(ikind)*int_size)


        DO iatom=1,natom(ikind)
          atom_a = atom_list(iatom)
          ra_pbc(:) = pbc(particle_set(atom_a)%r(:),cell)
          alist(ikind)%ipart(iatom) = atom_a
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

    e_rspace_core = 0.0_dbl

    DO ikind=1,nkind

      IF (alpha(ikind) == 0.0_dbl) CYCLE

      DO jkind=ikind,nkind

        IF (alpha(jkind) == 0.0_dbl) CYCLE

        equal_kinds = (ikind == jkind)

        zab = (ccore_charge(ikind)/SQRT((alpha(ikind)/pi)**3))*  &
              (ccore_charge(jkind)/SQRT((alpha(jkind)/pi)**3))
        aab = alpha(ikind)*alpha(jkind)/(alpha(ikind) + alpha(jkind))
        rootaab = SQRT(aab)
        f0 = 2.0_dbl*oorootpi*zab*rootaab

!       *** Calculate the square of the maximum interaction distance ***

        rab_max = radius(ikind) + radius(jkind)
        rab2_max = rab_max*rab_max

        r(:) = rab_max
        sab_max(:) = real_to_scaled(r(:),cell)

        ncell(:) = (INT(sab_max(:)) + 1)*periodic(:)
        ngrid(:) = MAX(1,NINT(0.5_dbl*subcells/sab_max(:)))

!       *** Print subcell information for the current atomic kind pair ***

        IF ((ionode.AND.globenv%print%cell_parameters).AND.&
            (globenv%print%level > LOW)) THEN
          WRITE (UNIT=output_unit,FMT="(T3,2I8,4X,3I5,6X,3F12.6)")&
            ikind,jkind,ngrid,&
            scaled_to_real(1.0_dbl/REAL(ngrid(:),dbl),cell)/unit_of_length
        END IF

        IF (SUM(ngrid) == 3) THEN

          DO jatom=1,natom(jkind)

            atom_a = alist(jkind)%ipart(jatom) 
!gt exclusion to be done
! By now it doesnt work for molecules:
! wait for the new exclusion lists
!            nexcl = pnode(atom_a) %nexcl
!            ll_exclude => pnode(atom_a) %ex
! list number of exclusion
!              DO iexclude = 1, nexcl
!               elist(iexclude) = ll_exclude%p%iatom
!               ll_exclude => ll_exclude%next
!            END DO
! sort elist
!            CALL sort(elist,nexcl,ind)
!gt exclusion to be done

            sb_pbc(:) = pbc_coord(jkind)%s(:,jatom)

            DO icell=-ncell(1),ncell(1)
              cell_b(1) = icell
              DO jcell=-ncell(2),ncell(2)
                cell_b(2) = jcell
                DO kcell=-ncell(3),ncell(3)
                  cell_b(3) = kcell

                  sb(:) = sb_pbc(:) + REAL(cell_b(:),dbl)
                  rb(:) = scaled_to_real(sb(:),cell)

                  IF (equal_kinds) THEN
                    natom_a = jatom - 1
                  ELSE
                    natom_a = natom(ikind)
                  END IF

                  DO iatom=1,natom_a
                    atom_b = alist(ikind)%ipart(iatom) 
!gt exclusion to be done
!                    IF ( get_match(atom_b,elist,nexcl)) CYCLE
!gt exclusion to be done
                    rab(:) = rb(:) - pbc_coord(ikind)%r(:,iatom)
                    rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                    IF (rab2 < rab2_max) THEN
                      dab = SQRT(rab2)
                      eab = zab*erfc(rootaab*dab)/dab
                      e_rspace_core = e_rspace_core + eab
                      IF (PRESENT (force)) THEN
                        deab(:) = rab(:)*(eab + f0*EXP(-aab*rab2))/rab2
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

          CALL stop_program(routine,"subcell grids not imlemented")

        END IF

      END DO
    END DO

!   *** Release work storage ***

    DEALLOCATE (alpha,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"alpha")

    DEALLOCATE (radius,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"radius")

    DEALLOCATE (ccore_charge,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"ccore_charge")

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
      IF (ASSOCIATED(alist(ikind)%ipart)) THEN
        DEALLOCATE (alist(ikind)%ipart,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"alist(ikind)%ipart")
      END IF
    END DO

    DEALLOCATE (pbc_coord,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pbc_coord")

    DEALLOCATE (alist,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"alist")

    DEALLOCATE (ind,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"ind")

! gt excl to be done
!    DEALLOCATE (elist,STAT=istat)
!    IF (istat /= 0) CALL stop_memory(routine,"elist")
! gt excl to be done

    CALL timestop(0.0_dbl,handle)

  END SUBROUTINE calculate_ecore_rspace

! *****************************************************************************

  SUBROUTINE calculate_eself(kg_env,ecore_self)

!   Purpose: Calculate the self energy of the core charge distribution.

!   History: - Creation (14.10.2002) 

!   ***************************************************************************

    TYPE(kg_environment_type), INTENT(IN) :: kg_env
    REAL (dbl), INTENT(OUT) :: ecore_self

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER :: atomic_kind
    REAL(dbl)                        :: alpha_core_charge,qeff,ccore_charge
    INTEGER                         :: ikind,natom

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set

!   ---------------------------------------------------------------------------

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set) 

    ecore_self = 0.0_dbl

    DO ikind=1,SIZE(atomic_kind_set)

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           natom=natom,&
                           ccore_charge=ccore_charge,&
                           alpha_core_charge=alpha_core_charge)

      qeff=ccore_charge/SQRT(alpha_core_charge/pi)**3
      ecore_self = ecore_self - REAL(natom*qeff**2,dbl)*SQRT(alpha_core_charge)

    END DO

    ecore_self = ecore_self/SQRT(twopi)

  END SUBROUTINE calculate_eself

! *****************************************************************************

END MODULE kg_core
