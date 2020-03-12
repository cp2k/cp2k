! collocate driver, reads relevant data from file and sets up the needed arrays
use kinds
use l_utils
use cube_utils
use orbital_pointers
use qs_collocate_density
use qs_integrate_potential
use basic_types
IMPLICIT NONE

REAL(kind=dp) :: eps_rho_rspace,dr(3),dh(3,3),error,max_error,ele
REAL :: time,time_min
INTEGER :: Ntest
INTEGER :: lmax,cmax,ngrid_levels,igrid_level,cmax_i,i,nprimitives,iostat
INTEGER, ALLOCATABLE :: max_radius(:)
TYPE(realspace_grid_type), ALLOCATABLE :: rs_grids(:)
TYPE(l_info_type) :: l_info
TYPE(cube_info_type), ALLOCATABLE :: cube_info(:)
TYPE(primitive_type), ALLOCATABLE :: p(:)
TYPE(primitive_type) :: ptest
REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: rtest

LOGICAL :: autogen_one, do_collocate
INTEGER :: one_lmax_a, one_lmax_b, one_lmin_a, one_lmin_b,input_collo_inte

write(6,*) "Check collocate (0) or integrate (1) "
read(5,*) input_collo_inte 
do_collocate=.FALSE.
IF (input_collo_inte==0) do_collocate=.TRUE.
IF (do_collocate) THEN
   write(6,*) "Checking collocate routines"
ELSE
   write(6,*) "Checking integrate routines"
ENDIF

write(6,*) "Make only one special case (T/F)"
read(5,*) autogen_one
write(6,*) autogen_one
IF (autogen_one) THEN
  write(6,*) "lmax_a, max_b, lmin_a, lmin_b"
  read(5,*) one_lmax_a, one_lmax_b, one_lmin_a, one_lmin_b
  write(6,*) one_lmax_a, one_lmax_b, one_lmin_a, one_lmin_b
ENDIF
write(6,*) "Number of times to run the test"
read(5,*) Ntest
write(6,*) "Ntest =",Ntest


! read data about the reference run
    OPEN(UNIT=1000,FILE='setup.dat')
    read(1000,*) eps_rho_rspace
    read(1000,*) lmax,cmax
    IF (autogen_one) THEN
       lmax=MAX(lmax,one_lmax_a,one_lmax_b)
    ENDIF
    read(1000,*) ngrid_levels
    ALLOCATE(max_radius(ngrid_levels))
    ALLOCATE(rs_grids(ngrid_levels))
    ALLOCATE(cube_info(ngrid_levels))
    DO igrid_level=1,ngrid_levels
        read(1000,*) rs_grids(igrid_level)%npts, &
                     rs_grids(igrid_level)%lb, &
                     rs_grids(igrid_level)%ub, &
                     rs_grids(igrid_level)%dr, &
                     rs_grids(igrid_level)%dh, &
                     rs_grids(igrid_level)%perd, &
                     rs_grids(igrid_level)%border, &
                     rs_grids(igrid_level)%ngpts_local, &
                     rs_grids(igrid_level)%npts_local, &
                     rs_grids(igrid_level)%lb_local, &
                     rs_grids(igrid_level)%ub_local

        ! allocate grids, assume the lb_local stup is actually right..
        ALLOCATE(rs_grids(igrid_level)%r(rs_grids(igrid_level)%lb_local(1):rs_grids(igrid_level)%ub_local(1), &
                                         rs_grids(igrid_level)%lb_local(2):rs_grids(igrid_level)%ub_local(2), &
                                         rs_grids(igrid_level)%lb_local(3):rs_grids(igrid_level)%ub_local(3)))

        read(1000,*) dr, &
                     dh, &
                     max_radius(igrid_level) ! notice not the same as the argument to init
        IF (ANY(dr.NE.rs_grids(igrid_level)%dr)) write(6,*) "Hi ... dr?"
        IF (ANY(dh.NE.rs_grids(igrid_level)%dh)) write(6,*) "Hi ... dh?"
    END DO
    CLOSE(1000)
! call init routines
    CALL init_orbital_pointers(2*lmax+1)
    CALL init_l_info(l_info,lmax,cmax)
    DO igrid_level=1,ngrid_levels
       CALL init_cube_info(cube_info(igrid_level),rs_grids(igrid_level)%dr,rs_grids(igrid_level)%dh, &
                            .TRUE.,max_radius(igrid_level),cmax_i)
       IF (cmax_i>cmax) write(6,*) "Hi ... cmax?"
    ENDDO

! read primitives data 
    OPEN(UNIT=1000,FILE='primitives.dat')
    REWIND(1000)
    nprimitives=0
    DO 
      READ(1000,*,IOSTAT=iostat) ptest%la_max,ptest%zeta,ptest%la_min,ptest%lb_max,ptest%zetb, &
                                 ptest%lb_min,ptest%ra,ptest%rab,ptest%rab2,ptest%igrid
      IF (IOSTAT==0) THEN
          nprimitives=nprimitives+1
      ELSE
          exit
      ENDIF
    ENDDO
    REWIND(1000)
    ALLOCATE(p(nprimitives))
    DO i=1,nprimitives
       READ(1000,*) p(i)%la_max,p(i)%zeta,p(i)%la_min,p(i)%lb_max,p(i)%zetb,p(i)%lb_min,p(i)%ra,p(i)%rab,p(i)%rab2,p(i)%igrid
       IF (autogen_one) THEN
          p(i)%la_max=one_lmax_a
          p(i)%lb_max=one_lmax_b
          p(i)%la_min=one_lmin_a
          p(i)%lb_min=one_lmin_b
       ENDIF
       p(i)%intsum=0.0_dp
    ENDDO
    write(6,*) '# of primitives ',nprimitives
! we're setup

    write(6,*) '# computational kernel timings ',Ntest


    IF (do_collocate) THEN

            time_min=HUGE(time_min)
            DO i=1,Ntest
               CALL collocate_all_type_1(p,rs_grids,cube_info,l_info,eps_rho_rspace,time)
               write(6,*) "Kernel time",time
               IF (time<time_min) time_min=time
            ENDDO
            write(6,*) '# best time ',time_min
            
            ! read/write for checking
            IF (.NOT. autogen_one) THEN
             max_error=0.0
             DO i=1,SIZE(rs_grids)
                ALLOCATE(rtest(rs_grids(i)%lb_local(1):rs_grids(i)%ub_local(1), &
                               rs_grids(i)%lb_local(2):rs_grids(i)%ub_local(2), &
                               rs_grids(i)%lb_local(3):rs_grids(i)%ub_local(3)))
                OPEN(300+i)
                read(300+i,*) rtest
                CLOSE(300+i)
                error=MAXVAL(ABS(rtest-rs_grids(i)%r))
                write(6,'(1X,A,E10.2)') "# Max difference between computed and reference ",error
                if (error>max_error) max_error=error
                DEALLOCATE(rtest)
             ENDDO
             write(6,*) "# largest error ",max_error
            ENDIF

    ELSE

            time_min=HUGE(time_min)
            DO i=1,Ntest
               CALL integrate_all(p,rs_grids,cube_info,l_info,eps_rho_rspace,time)
               write(6,*) "Kernel time",time
               IF (time<time_min) time_min=time
            ENDDO
            write(6,*) '# best time ',time_min

            ! read/write for checking
            IF (.NOT. autogen_one) THEN
               max_error=0.0_dp
               OPEN(200)
               DO i=1,SIZE(p)
                  read(200,*) ele
                  error=ABS(p(i)%intsum-ele)/MAX(1.0_dp,ABS(ele))
                  IF (error>max_error) max_error=error
               ENDDO
               CLOSE(200)
               write(6,*) "# largest error ",max_error
            ENDIF

    ENDIF

! cleanup

    CALL destroy_l_info(l_info)
    DO igrid_level=1,ngrid_levels
       CALL destroy_cube_info(cube_info(igrid_level))
    ENDDO
    CALL deallocate_orbital_pointers()
    DO igrid_level=1,ngrid_levels
       DEALLOCATE(rs_grids(igrid_level)%r)
    ENDDO
END
