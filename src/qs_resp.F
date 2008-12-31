!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief provides a resp fit for gas phase systems
!> \par History
!>      created
!> \author Joost VandeVondele (02.2007)
! *****************************************************************************
MODULE qs_resp
  USE atomic_charges,                  ONLY: print_atomic_charges
  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE cell_types,                      ONLY: cell_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE pw_types,                        ONLY: pw_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_resp'

  PUBLIC :: resp_fit

CONTAINS

! *****************************************************************************
  SUBROUTINE resp_fit(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'resp_fit', &
      routineP = moduleN//':'//routineN

    INTEGER :: bo(2,3), gbo(2,3), handle, i, info, jx, jy, jz, k, m, natom, &
      ncons, ncons_i, ncons_sec, ncons_v, np(3), npoints, nres, nrest_sec, &
      nvar, output_unit, stride(3), z
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ipiv
    INTEGER, DIMENSION(:), POINTER           :: atom_list, my_stride
    LOGICAL                                  :: explicit, failure, has_resp, &
                                                itc, rheavies
    REAL(KIND=dp) :: dr(3), dvol, mycharge, r(3), rheavies_strength, rmax, &
      rmin, shift(3), strength, TARGET, vec(3), vec_pbc(3), vj
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: dist, rhs
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: matrix
    REAL(KIND=dp), DIMENSION(:), POINTER     :: atom_coef
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_type), POINTER                   :: v_hartree_pw
    TYPE(section_vals_type), POINTER         :: cons_section, input, &
                                                resp_section, rest_section

    CALL timeset(routineN,handle)

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
       CALL get_qs_env(qs_env, input=input, subsys=subsys, particle_set=particle_set, &
                               cell=cell, error=error)
       resp_section => section_vals_get_subs_vals(input,"PROPERTIES%RESP",error=error)
       CALL section_vals_get(resp_section, explicit=has_resp, error=error)
    END IF

    IF (.NOT. failure .AND. has_resp) THEN
       logger => cp_error_get_logger(error)
       output_unit=cp_print_key_unit_nr(logger,resp_section,"PROGRAM_RUN_INFO",&
                                        extension=".resp",error=error)

       cons_section=>section_vals_get_subs_vals(resp_section,"CONSTRAINT",error=error)
       rest_section=>section_vals_get_subs_vals(resp_section,"RESTRAINT",error=error)

       CALL section_vals_get(cons_section, explicit=explicit, error=error)
       IF (explicit) THEN
          CALL section_vals_get(cons_section,n_repetition=ncons_sec,error=error) 
       ELSE
          ncons_sec=0
       ENDIF
       CALL section_vals_get(rest_section, explicit=explicit, error=error)
       IF (explicit) THEN
          CALL section_vals_get(rest_section,n_repetition=nrest_sec,error=error) 
       ELSE
          nrest_sec=0
       ENDIF
       ncons=ncons_sec
       nres=nrest_sec

       CALL section_vals_val_get(resp_section,"INTEGER_TOTAL_CHARGE",l_val=itc,error=error)
       IF (itc) ncons=ncons+1

       CALL section_vals_val_get(resp_section,"RESTRAIN_HEAVIES_TO_ZERO",l_val=rheavies,error=error)
       IF (rheavies) THEN
          CALL section_vals_val_get(resp_section,"RESTRAIN_HEAVIES_STRENGTH",r_val=rheavies_strength,error=error)
       ENDIF

       CALL section_vals_val_get(resp_section,"RMIN",r_val=rmin,error=error) 
       CALL section_vals_val_get(resp_section,"RMAX",r_val=rmax,error=error) 
       CALL section_vals_val_get(resp_section,"STRIDE",i_vals=my_stride,error=error) 
       CALL cp_assert(SIZE(my_stride)==1.OR.SIZE(my_stride)==3,cp_fatal_level,cp_assertion_failed,routineP,&
            "STRIDE keyword can accept only 1 (the same for X,Y,Z) or 3 values. Correct your input file."//&
CPSourceFileRef,&
            only_ionode=.TRUE.)       
       IF (SIZE(my_stride)==1) THEN
          DO i = 1,3
             stride(i) = my_stride(1)
          END DO
       ELSE
          stride = my_stride(1:3)
       END IF
       IF (output_unit>0) THEN
          WRITE(output_unit,'(2X,A)')        "Starting RESP fit"
          WRITE(output_unit,'(2X,A,I0)')     "Number of explicit constraints: ",ncons_sec
          WRITE(output_unit,'(2X,A,I0)')     "Number of explicit restraints: ",nrest_sec
          WRITE(output_unit,'(2X,A,A)')      "Constrain total charge ",MERGE("T","F",itc) 
          WRITE(output_unit,'(2X,A,A)')      "Restrain heavy atoms ",MERGE("T","F",rheavies) 
          IF (rheavies) THEN
             WRITE(output_unit,'(2X,A,F10.5)') "Heavy atom restraint strength: ",rheavies_strength
          ENDIF
          WRITE(output_unit,'(2X,A,F10.5)')  "Rmin [a.u.]: ",rmin
          WRITE(output_unit,'(2X,A,F10.5)')  "Rmax [a.u.]: ",rmax
          WRITE(output_unit,'(2X,A,3I5)')     "Stride: ",stride
       ENDIF

       CALL cp_subsys_get(subsys,particles=particles,error=error)
       natom=particles%n_els
       nvar=natom+ncons

       ALLOCATE(matrix(nvar,nvar))
       matrix = 0.0_dp 
       ALLOCATE(rhs(nvar))
       rhs = 0.0_dp
       ALLOCATE(ipiv(nvar))

       ! loop over all grid points that need to be fitted
       v_hartree_pw => qs_env%ks_env%v_hartree_rspace%pw
       bo=v_hartree_pw%pw_grid%bounds_local
       gbo=v_hartree_pw%pw_grid%bounds
       np=v_hartree_pw%pw_grid%npts
       dr=v_hartree_pw%pw_grid%dr
       dvol=v_hartree_pw%pw_grid%dvol
       shift=-REAL(MODULO(np,2),dp)*dr/2.0_dp
       ALLOCATE(dist(natom))

       mycharge=0.0_dp
       npoints = 0
       DO jz=bo(1,3),bo(2,3)
       DO jy=bo(1,2),bo(2,2)
       DO jx=bo(1,1),bo(2,1)  
          IF (.NOT.(MODULO(jz,stride(3))==0)) CYCLE
          IF (.NOT.(MODULO(jy,stride(2))==0)) CYCLE
          IF (.NOT.(MODULO(jx,stride(1))==0)) CYCLE

          ! compute distance from the grid point to all atoms
          r(3)=jz*dr(3)+shift(3)
          r(2)=jy*dr(2)+shift(2)
          r(1)=jx*dr(1)+shift(1)
          DO i=1,natom
             vec=r-particles%els(i)%r
             ! always wrap this 
             vec_pbc(3)=MODULO(vec(3),cell%hmat(3,3))-cell%hmat(3,3)/2._dp
             vec_pbc(2)=MODULO(vec(2),cell%hmat(2,2))-cell%hmat(2,2)/2._dp
             vec_pbc(1)=MODULO(vec(1),cell%hmat(1,1))-cell%hmat(1,1)/2._dp
             dist(i)=SQRT(SUM(vec_pbc**2))
          ENDDO
          ! check the point is sufficiently close and sufficiently far
          ! if OK, we can use the point for our fitting procedure
          IF (ALL(dist>RMAX).OR.ANY(dist<RMIN)) CYCLE
          npoints = npoints + 1

          ! correct for the fact that v_hartree is scaled by dvol, and has the opposite sign
          vj=-v_hartree_pw%cr3d(jx,jy,jz)/dvol
          dist=1.0_dp/dist

          DO i=1,natom
           DO m=1,natom
              matrix(m,i)=matrix(m,i)+2.0_dp*dist(i)*dist(m)
           ENDDO
           rhs(i)     =rhs(i)+2.0_dp*vj*dist(i)
          ENDDO
       ENDDO
       ENDDO
       ENDDO

       CALL mp_sum(npoints,v_hartree_pw%pw_grid%para%group)
       CALL mp_sum(matrix,v_hartree_pw%pw_grid%para%group)
       CALL mp_sum(rhs,v_hartree_pw%pw_grid%para%group)
       matrix=matrix/npoints
       rhs=rhs/npoints

       IF (output_unit>0) THEN
          WRITE(output_unit,'(2X,A,I0)') "Number of potential fitting points found: ",npoints
          WRITE(output_unit,'()') 
       ENDIF

       ! add the restraints
       DO i=1,nrest_sec
          CALL section_vals_val_get(rest_section,"ATOM_LIST",i_rep_section=i,i_vals=atom_list,error=error)
          CALL section_vals_val_get(rest_section,"ATOM_COEF",i_rep_section=i,r_vals=atom_coef,error=error)
          CALL section_vals_val_get(rest_section,"TARGET",i_rep_section=i,r_val=TARGET,error=error)
          CALL section_vals_val_get(rest_section,"STRENGTH",i_rep_section=i,r_val=strength,error=error)
          CPPrecondition(SIZE(atom_list)==SIZE(atom_coef),cp_failure_level,routineP,error,failure)
          DO m=1,SIZE(atom_list)
             DO k=1,SIZE(atom_list)
                matrix(atom_list(m),atom_list(k))=matrix(atom_list(m),atom_list(k))+ &
                                                   atom_coef(m)*atom_coef(k)*2.0_dp*strength
             ENDDO
             rhs(atom_list(m))=rhs(atom_list(m))+2.0_dp*TARGET*strength*atom_coef(m)
          ENDDO
       ENDDO

       ! if heavies are restrained to zero, add these as well
       IF (rheavies) THEN
          DO i=1,natom
             CALL get_atomic_kind(atomic_kind=particle_set(i)%atomic_kind,z=z)
             IF (z.NE.1) THEN
                matrix(i,i)=matrix(i,i)+2.0_dp*rheavies_strength
             ENDIF
          ENDDO
       ENDIF

       ! add the constraints
       ncons_i=0
       IF (itc) THEN
          ncons_i=ncons_i+1        
          ncons_v=natom+ncons_i
          matrix(1:natom,ncons_v)=1.0_dp
          matrix(ncons_v,1:natom)=1.0_dp
          rhs(ncons_v)=qs_env%dft_control%charge
       ENDIF

       DO i=1,ncons_sec
          ncons_i=ncons_i+1
          ncons_v=ncons_v+1
          CALL section_vals_val_get(cons_section,"ATOM_LIST",i_rep_section=i,i_vals=atom_list,error=error)
          CALL section_vals_val_get(cons_section,"ATOM_COEF",i_rep_section=i,r_vals=atom_coef,error=error)
          CALL section_vals_val_get(cons_section,"TARGET",i_rep_section=i,r_val=TARGET,error=error)
          CPPrecondition(SIZE(atom_list)==SIZE(atom_coef),cp_failure_level,routineP,error,failure)
          DO m=1,SIZE(atom_list)
             matrix(atom_list(m),ncons_v)=atom_coef(m)
             matrix(ncons_v,atom_list(m))=atom_coef(m)
          ENDDO
          rhs(ncons_v)=TARGET
       ENDDO

       ! solve system for the values of the charges and the lagrangian multipliers
       CALL DGETRF(nvar,nvar,matrix,nvar,ipiv,info)
       CPPrecondition(info==0,cp_failure_level,routineP,error,failure)
       
       CALL DGETRS('N',nvar,1,matrix,nvar,ipiv,rhs,nvar,info)
       CPPrecondition(info==0,cp_failure_level,routineP,error,failure)

       CALL print_atomic_charges(particle_set,output_unit,title="RESP charges:",atomic_charges=rhs(1:natom))

       DEALLOCATE(matrix)
       DEALLOCATE(rhs)
       DEALLOCATE(ipiv)

       CALL cp_print_key_finished_output(output_unit,logger,resp_section,&
            "PROGRAM_RUN_INFO", error=error)

    END IF

    CALL timestop(handle)

  END SUBROUTINE resp_fit

END MODULE qs_resp
