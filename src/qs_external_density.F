!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
!> \brief Routines to handle an external density 
!>        The external density can be generic and is provided by user input
!> \author D. Varsano
! *****************************************************************************
MODULE qs_external_density
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_files,                        ONLY: close_file,&
                                             get_unit_number,&
                                             open_file
  USE gaussian_gridlevels,             ONLY: gridlevel_info_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE message_passing,                 ONLY: mp_recv,&
                                             mp_send,&
                                             mp_sync
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_integrate_function
  USE pw_types,                        ONLY: pw_p_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
  USE realspace_grid_types,            ONLY: realspace_grid_desc_p_type,&
                                             realspace_grid_p_type,&
                                             rs_grid_create,&
                                             rs_grid_zero
  USE rs_pw_interface,                 ONLY: density_rs2pw
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_external_density'

  PUBLIC :: external_read_density

CONTAINS

! *****************************************************************************
!> \brief  Computes the external density on the grid
!> \param qs_env ...
!> \date   03.2011
!> \author D. Varsano 
! *****************************************************************************
  SUBROUTINE external_read_density(qs_env)

    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(len=*), PARAMETER :: routineN = 'external_read_density', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: filename
    INTEGER                                  :: extunit, gid, handle, i, &
                                                igrid_level, ip, j, k, &
                                                master, my_rank, nat, ndum, &
                                                num_pe, stat, tag
    INTEGER, DIMENSION(3)                    :: lbounds, lbounds_local, &
                                                npoints, npoints_local, &
                                                ubounds, ubounds_local
    LOGICAL                                  :: failure
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: buffer
    REAL(kind=dp), DIMENSION(3)              :: dr, rdum
    REAL(KIND=dp), DIMENSION(:), POINTER     :: tot_rho_r_ext
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_ext_g, rho_ext_r
    TYPE(qs_rho_type), POINTER               :: rho_external
    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), POINTER                  :: rs_descs
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_rho_ext
    TYPE(section_vals_type), POINTER         :: ext_den_section, input

    failure=.FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(cell, rs_rho_ext, input, ext_den_section, rs_descs, dft_control)
    NULLIFY(rho_ext_r, rho_ext_g, tot_rho_r_ext)

    CALL get_qs_env(qs_env,&
                    cell=cell,&
                    rho_external=rho_external,&
                    input=input,&
                    pw_env=pw_env,&
                    dft_control=dft_control)

    IF (dft_control%apply_external_density) THEN
      CALL qs_rho_get(rho_external,&
                    rho_r=rho_ext_r,&
                    rho_g=rho_ext_g,&
                    tot_rho_r=tot_rho_r_ext)

      gridlevel_info=>pw_env%gridlevel_info

      CALL pw_env_get(pw_env, rs_descs=rs_descs)

      ALLOCATE (rs_rho_ext(gridlevel_info%ngrid_levels),STAT=stat)

      DO igrid_level=1,gridlevel_info%ngrid_levels
        CALL rs_grid_create(rs_rho_ext(igrid_level)%rs_grid,&
                            rs_descs(igrid_level)%rs_desc)
        CALL rs_grid_zero(rs_rho_ext(igrid_level)%rs_grid)
      ENDDO

      igrid_level=igrid_level-1

      ext_den_section => section_vals_get_subs_vals(input,"DFT%EXTERNAL_DENSITY")
      CALL section_vals_val_get(ext_den_section,"FILE_DENSITY",c_val=filename)

      gid     = rho_ext_r(1)%pw%pw_grid%para%group
      my_rank = rho_ext_r(1)%pw%pw_grid%para%my_pos
      num_pe  = rho_ext_r(1)%pw%pw_grid%para%group_size
      tag     = 1

      IF (dft_control%read_external_density) THEN

        DO i =1,3
          dr(i) = rs_descs(igrid_level)%rs_desc%dh(i,i)  
        END DO
        npoints = rs_descs(igrid_level)%rs_desc%npts
        lbounds = rs_descs(igrid_level)%rs_desc%lb
        ubounds = rs_descs(igrid_level)%rs_desc%ub

        npoints_local = rho_ext_r(1)%pw%pw_grid%npts_local  
        lbounds_local = rho_ext_r(1)%pw%pw_grid%bounds_local(1,:)
        ubounds_local = rho_ext_r(1)%pw%pw_grid%bounds_local(2,:)
        
!        IF ( my_rank .eq. 0 ) THEN
!          WRITE(*,*) my_rank,npoints, lbounds, ubounds
!          WRITE(*,*) my_rank,npoints_local
!          WRITE(*,*) my_rank,lbounds_local
!          WRITE(*,*) my_rank,ubounds_local
!        END IF

        ALLOCATE(buffer(lbounds_local(3):ubounds_local(3)))
 
        IF ( my_rank == 0 ) THEN
          WRITE (*,FMT="(/,/,T2,A)") "INITIALIZING ZMP CONSTRAINED DENSITY METHOD"
          WRITE(*,FMT="(/,(T3,A,T51,A30))")"ZMP| Reading the target density:     ",filename

          CALL open_file(file_name=filename,&
                         file_status="OLD",&
                         file_form="FORMATTED",&
                         file_action="READ", &
                         unit_number=extunit)

          DO i =1,2
            READ(extunit,*) 
          ENDDO
          READ(extunit,*) nat,rdum
          DO i=1,3
            READ(extunit,*) ndum,rdum
            IF (ndum/=npoints(i).OR.(ABS(rdum(i)-dr(i)) > 1e-4 )) THEN   
              WRITE(*,*) "ZMP | ERROR! | CUBE FILE NOT COINCIDENT WITH INTERNAL GRID ", i
              WRITE(*,*) "ZMP | ",ndum," DIFFERS FROM ", npoints(i) 
              WRITE(*,*) "ZMP | ",rdum," DIFFERS FROM ", dr(i)
            ENDIF
          ENDDO
          DO i =1,nat
            READ(extunit,*)
          ENDDO
        ENDIF


        DO i=lbounds(1),ubounds(1)
          DO j=lbounds(2),ubounds(2) 
            IF ( my_rank .EQ. 0 ) THEN 
              READ(extunit,*) (buffer(k),k=lbounds(3),ubounds(3))   
              IF ( num_pe .GT. 1 ) THEN 
                DO ip=1,num_pe-1
                  CALL mp_send(buffer(lbounds(3):ubounds(3)),ip,tag,gid)
!                  WRITE(*,*) my_rank," sending to : ",ip
                ENDDO
!                CALL mp_sync(gid)
              ENDIF
            ELSE
               master=0  
               CALL mp_recv(buffer(lbounds(3):ubounds(3)),master,tag,gid)
!               WRITE(*,*) my_rank," recieving from master"
            ENDIF
            CALL mp_sync(gid)

            IF ( (lbounds_local(1).le.i).AND.(i.le.ubounds_local(1)).AND.(lbounds_local(2).le.j)&
                      .AND.(j.le.ubounds_local(2)) ) THEN
               rs_rho_ext(igrid_level)%rs_grid%r(i,j,lbounds(3):ubounds(3))=buffer(lbounds(3):ubounds(3))
            ENDIF

          ENDDO
        ENDDO
        IF (my_rank == 0) CALL close_file(unit_number=extunit)
      ENDIF


      CALL density_rs2pw(pw_env,rs_rho_ext,rho=rho_ext_r(1),rho_gspace=rho_ext_g(1))
      tot_rho_r_ext(1) = pw_integrate_function(rho_ext_r(1)%pw,isign=1)
      IF (my_rank == 0) THEN
        WRITE(*,FMT="(T3,A,T61,F20.10)") "ZMP| Total external charge:                    ",&
               tot_rho_r_ext(1)
      ENDIF
      DEALLOCATE(buffer,rs_rho_ext)
      CALL mp_sync(gid)
    ENDIF

    CALL timestop(handle)
    
  END SUBROUTINE external_read_density

! *****************************************************************************
!> \brief  Read external vxc potential
!> \param qs_env ...
!> \date   03.2011
!> \author DV 
! *****************************************************************************
  SUBROUTINE external_read_vxc(qs_env)

    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(len=*), PARAMETER :: routineN = 'external_read_vxc', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: filename
    INTEGER                                  :: extunit, handle, i, j, k, &
                                                nat, ndum
    INTEGER, DIMENSION(2, 3)                 :: bo_global, bo_local
    INTEGER, DIMENSION(3)                    :: lbounds, npoints, ubounds
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: dvol
    REAL(kind=dp), DIMENSION(3)              :: dr, rdum
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_p_type), POINTER                 :: vxc_ext
    TYPE(section_vals_type), POINTER         :: ext_vxc_section, input

    failure=.FALSE.
    NULLIFY(cell, vxc_ext, input, ext_vxc_section, dft_control)
    CALL timeset(routineN,handle)
    CALL get_qs_env(qs_env,&
                    cell=cell,&
                    external_vxc=vxc_ext,&
                    input=input,&
                    dft_control=dft_control)
    IF (dft_control%apply_external_vxc) THEN
      ext_vxc_section => section_vals_get_subs_vals(input,"DFT%EXTERNAL_VXC")
      IF (dft_control%read_external_vxc) THEN
        dr=vxc_ext%pw%pw_grid%dr
        dvol=vxc_ext%pw%pw_grid%dvol
        vxc_ext%pw%cr3d=0.0_dp

        bo_local=vxc_ext%pw%pw_grid%bounds_local
        bo_global=vxc_ext%pw%pw_grid%bounds
        DO i=1,3
          npoints(i) = bo_local(2,i)-bo_local(1,i)+1
          lbounds(i) = bo_local(1,i)
          ubounds(i) = bo_local(2,i)
        ENDDO
        CALL section_vals_val_get(ext_vxc_section,"FILE_VXC",c_val=filename)
        WRITE(*,*) "NOW READING VXC FILE:   ", filename
        extunit=get_unit_number()
        CALL open_file(file_name=filename,file_status="OLD",&
               file_form="FORMATTED", file_action="READ", unit_number=extunit)
        DO i =1,2
          READ(extunit,*)
        ENDDO
        READ(extunit,*) nat,rdum
        DO i=1,3
           READ(extunit,*) ndum,rdum
           IF (ndum/=npoints(i).OR.(ABS(rdum(i)-dr(i)) > 1e-4 )) THEN
!Here valid only for cubic cells
             WRITE(*,*) "check dimensions ", i
             WRITE(*,*) ndum,"differs from", npoints(i)
             WRITE(*,*) rdum(i),"differs from", dr(i)
             CPABORT("")
           ENDIF
        ENDDO
        DO i =1,nat
          READ(extunit,*)
        ENDDO 
        DO i=lbounds(1),ubounds(1)
          DO j=lbounds(2),ubounds(2)
            READ(extunit,'(6E13.5)') (vxc_ext%pw%cr3d(i,j,k),k=lbounds(3),ubounds(3))
          ENDDO
        ENDDO
        CALL close_file(unit_number=extunit)
        dft_control%read_external_vxc = .FALSE.
      END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE external_read_vxc

END MODULE qs_external_density

