!***********************************************************************************************************************************
!
!
!
!       COLLOCATE COLLOCATE COLLOCATE COLLOCATE COLLOCATE COLLOCATE COLLOCATE COLLOCATE COLLOCATE COLLOCATE COLLOCATE COLLOCATE 
!
!
!
!***********************************************************************************************************************************
MODULE collocate_generate
  USE option_module

  IMPLICIT NONE


CONTAINS

SUBROUTINE collocate_generate_stub(iunit,lp)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iunit,lp

write(iunit,'(A,I0,A)')"  SUBROUTINE collocate_core_",lp,"(grid,coef_xyz,pol_x,pol_y,pol_z,map,sphere_bounds,cmax,gridbounds)"
write(iunit,'(A)')""
write(iunit,'(A)')"    IMPLICIT NONE"
write(iunit,'(A)')""
write(iunit,'(A)')"  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND ( 14, 200 )"
write(iunit,'(A,I0)')" INTEGER, PARAMETER :: lp=",lp
write(iunit,'(A)')"    INTEGER, INTENT(IN)     :: cmax"
write(iunit,'(A)')"    INTEGER, INTENT(IN)     :: gridbounds(2,3)"
write(iunit,'(A)')"    INTEGER, INTENT(IN)     :: sphere_bounds(*)"
write(iunit,'(A)')"    INTEGER, INTENT(IN)     :: map(-cmax:cmax,1:3)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(IN)    :: pol_x(0:lp,-cmax:cmax)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(IN)    :: pol_y(1:2,0:lp,-cmax:0)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(IN)    :: pol_z(1:2,0:lp,-cmax:0)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(IN)    :: coef_xyz(((lp+1)*(lp+2)*(lp+3))/6)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(INOUT) :: grid(gridbounds(1,1):gridbounds(2,1), & "
write(iunit,'(A)')"                                    gridbounds(1,2):gridbounds(2,2), &"
write(iunit,'(A)')"                                    gridbounds(1,3):gridbounds(2,3))"
write(iunit,'(A)')" CALL collocate_core_default(grid,coef_xyz,pol_x,pol_y,pol_z,&"
write(iunit,'(A,I0,A)')  "map,sphere_bounds,",lp,",cmax,gridbounds)" 
write(iunit,'(A,I0)')" END SUBROUTINE collocate_core_",lp

END SUBROUTINE collocate_generate_stub

SUBROUTINE collocate_generate_specialized(iunit,lp,options)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iunit,lp
  TYPE (option_type)  :: options

write(iunit,'(A,I0,A)')"  SUBROUTINE collocate_core_",lp,"(grid,coef_xyz,pol_x,pol_y,pol_z,map,sphere_bounds,cmax,gridbounds)"
write(iunit,'(A)')""
write(iunit,'(A)')"    IMPLICIT NONE"
write(iunit,'(A)')""
write(iunit,'(A)')"  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND ( 14, 200 )"
write(iunit,'(A,I0)')" INTEGER, PARAMETER :: lp=",lp
write(iunit,'(A)')"    INTEGER, INTENT(IN)     :: cmax"
write(iunit,'(A)')"    INTEGER, INTENT(IN)     :: gridbounds(2,3)"
write(iunit,'(A)')"    INTEGER, INTENT(IN)     :: sphere_bounds(*)"
write(iunit,'(A)')"    INTEGER, INTENT(IN)     :: map(-cmax:cmax,1:3)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(IN)    :: pol_x(0:lp,-cmax:cmax)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(IN)    :: pol_y(1:2,0:lp,-cmax:0)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(IN)    :: pol_z(1:2,0:lp,-cmax:0)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(IN)    :: coef_xyz(((lp+1)*(lp+2)*(lp+3))/6)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(INOUT) :: grid(gridbounds(1,1):gridbounds(2,1), & "
write(iunit,'(A)')"                                    gridbounds(1,2):gridbounds(2,2), &"
write(iunit,'(A)')"                                    gridbounds(1,3):gridbounds(2,3))"
write(iunit,'(A)')""
write(iunit,'(A)')"    REAL(KIND=dp) ::  coef_xy(2,(lp+1)*(lp+2)/2)"
write(iunit,'(A)')"    REAL(KIND=dp) ::  coef_x(4,0:lp)"
write(iunit,'(A)')""
write(iunit,'(A)')"    INTEGER kg,k,jgmin,jg,j,j2,igmax,ig,i,kgmin,igmin,k2,i2,jg2,kg2"
write(iunit,'(A)')"    INTEGER sci,lxp,lyp,lzp"
write(iunit,'(A)')"    REAL(KIND=dp) s01,s02,s03,s04,s(4)"
write(iunit,'(A)')"    INTEGER lxyz,lxy,lx"
write(iunit,'(A)')""
write(iunit,'(A)')"    sci=1"
write(iunit,'(A)')""
write(iunit,'(A)')"    kgmin=sphere_bounds(sci)"
write(iunit,'(A)')"    sci=sci+1"
write(iunit,'(A)')"    DO kg=kgmin,0"
write(iunit,'(A)')"       kg2=1-kg"
write(iunit,'(A)')"       k=map(kg,3)"
write(iunit,'(A)')"       k2=map(kg2,3)"

CALL write_kg_body()

write(iunit,'(A)')"       jgmin=sphere_bounds(sci)"
write(iunit,'(A)')"       sci=sci+1"
write(iunit,'(A)')"       DO jg=jgmin,0"
write(iunit,'(A)')"          jg2=1-jg"
write(iunit,'(A)')"          j=map(jg,2)"
write(iunit,'(A)')"          j2=map(jg2,2)"
write(iunit,'(A)')"          igmin=sphere_bounds(sci)"
write(iunit,'(A)')"          sci=sci+1"
write(iunit,'(A)')"          igmax=1-igmin"

CALL write_jg_body()

CALL write_ig_loop()

write(iunit,'(A)')"       END DO"
write(iunit,'(A)')"    END DO"
write(iunit,'(A)')"  "
write(iunit,'(A,I0)')"  END SUBROUTINE collocate_core_",lp

CONTAINS

!***********************************************************************************************************************************
!
! write_kg_body
!
!***********************************************************************************************************************************

SUBROUTINE write_kg_body()
  INTEGER :: lxp,lyp,lzp,lxy,lxyz  

  write(iunit,'(A)')"       coef_xy=0.0_dp"
  SELECT CASE (options%kg_loop_unroll_lxp)
  CASE (0)
    write(iunit,'(A)')"       lxyz = 0"
    write(iunit,'(A)')"       DO lzp=0,lp"
    write(iunit,'(A)')"          lxy=0"
    write(iunit,'(A)')"          DO lyp=0,lp-lzp"
    write(iunit,'(A)')"             DO lxp=0,lp-lzp-lyp"
    write(iunit,'(A)')"                lxyz=lxyz+1 ; lxy=lxy+1"
    CALL write_kg_lxy_body(lzp,lxy,lxyz)
    write(iunit,'(A)')"             ENDDO"
    write(iunit,'(A)')"             lxy=lxy+lzp"
    write(iunit,'(A)')"          ENDDO"
    write(iunit,'(A)')"       ENDDO"
  CASE (1)
    lxyz = 0
    DO lzp=0,lp
       lxy=0
       DO lyp=0,lp-lzp
          DO lxp=0,lp-lzp-lyp
             lxyz=lxyz+1 ; lxy=lxy+1
             CALL write_kg_lxy_body(lzp,lxy,lxyz)
          ENDDO
          lxy=lxy+lzp
       ENDDO
    ENDDO
  END SELECT
END SUBROUTINE

SUBROUTINE write_kg_lxy_body(lzp,lxy,lxyz)
  INTEGER :: lzp,lxy,lxyz
  SELECT CASE (options%kg_loop_unroll_lxp)
  CASE (0)
    SELECT CASE (options%kg_loop_vector_notation)
    CASE(0)
      write(iunit,'(A)')"                coef_xy(1,lxy)=coef_xy(1,lxy)+coef_xyz(lxyz)*pol_z(1,lzp,kg)"
      write(iunit,'(A)')"                coef_xy(2,lxy)=coef_xy(2,lxy)+coef_xyz(lxyz)*pol_z(2,lzp,kg)"
    CASE(1)
      write(iunit,'(A)')"                coef_xy(:,lxy)=coef_xy(:,lxy)+coef_xyz(lxyz)*pol_z(:,lzp,kg)"
    END SELECT
  CASE(1)
    SELECT CASE (options%kg_loop_vector_notation)
    CASE(0)
      write(iunit,'(A,I0,A,I0,A,I0,A,I0,A)')"                coef_xy(1,",lxy,")=coef_xy(1,",lxy,")+"//&
                                                                                "coef_xyz(",lxyz,")*pol_z(1,",lzp,",kg)"
      write(iunit,'(A,I0,A,I0,A,I0,A,I0,A)')"                coef_xy(2,",lxy,")=coef_xy(2,",lxy,")+"//&
                                                                                "coef_xyz(",lxyz,")*pol_z(2,",lzp,",kg)"
    CASE(1)
      write(iunit,'(A,I0,A,I0,A,I0,A,I0,A)')"                coef_xy(:,",lxy,")=coef_xy(:,",lxy,")+"//&
                                                                                "coef_xyz(",lxyz,")*pol_z(:,",lzp,",kg)"
    END SELECT
  END SELECT
END SUBROUTINE

!***********************************************************************************************************************************
!
! write_jg_body
!
!***********************************************************************************************************************************

SUBROUTINE write_jg_body()
  INTEGER :: lxp,lyp,lxy

  write(iunit,'(A)')"          coef_x=0.0_dp"
  SELECT CASE (options%jg_loop_unroll_lxp)
  CASE (0)
    write(iunit,'(A)')"          lxy=0"
    write(iunit,'(A)')"          DO lyp=0,lp"
    write(iunit,'(A)')"          DO lxp=0,lp-lyp"
    write(iunit,'(A)')"             lxy=lxy+1"
    CALL write_jg_lxy_body(lyp,lxp,lxy)
    write(iunit,'(A)')"          ENDDO"
    write(iunit,'(A)')"          ENDDO"
  CASE (1)
    lxy=0
    DO lyp=0,lp
    DO lxp=0,lp-lyp
       lxy=lxy+1
       CALL write_jg_lxy_body(lyp,lxp,lxy)
    ENDDO
    ENDDO
  END SELECT
END SUBROUTINE write_jg_body

SUBROUTINE write_jg_lxy_body(lyp,lxp,lxy)
  INTEGER :: lxp,lyp,lxy

  SELECT CASE (options%jg_loop_unroll_lxp)
  CASE (0)
    SELECT CASE (options%jg_loop_vector_notation)
    CASE(0)
      write(iunit,'(A)')"             coef_x(1,lxp)=coef_x(1,lxp)+coef_xy(1,lxy)*pol_y(1,lyp,jg)"
      write(iunit,'(A)')"             coef_x(2,lxp)=coef_x(2,lxp)+coef_xy(2,lxy)*pol_y(1,lyp,jg)"
      write(iunit,'(A)')"             coef_x(3,lxp)=coef_x(3,lxp)+coef_xy(1,lxy)*pol_y(2,lyp,jg)"
      write(iunit,'(A)')"             coef_x(4,lxp)=coef_x(4,lxp)+coef_xy(2,lxy)*pol_y(2,lyp,jg)"
    CASE(1)
      write(iunit,'(A)')"             coef_x(1:2,lxp)=coef_x(1:2,lxp)+coef_xy(1:2,lxy)*pol_y(1,lyp,jg)"
      write(iunit,'(A)')"             coef_x(3:4,lxp)=coef_x(3:4,lxp)+coef_xy(1:2,lxy)*pol_y(2,lyp,jg)"
    END SELECT
  CASE (1)
    SELECT CASE (options%jg_loop_vector_notation)
    CASE(0)
      write(iunit,'(A,I0,A,I0,A,I0,A,I0,A)')"             coef_x(1,",lxp,")=coef_x(1,",lxp,&
                                                                            ")+coef_xy(1,",lxy,")*pol_y(1,",lyp,",jg)"
      write(iunit,'(A,I0,A,I0,A,I0,A,I0,A)')"             coef_x(2,",lxp,")=coef_x(2,",lxp,&
                                                                            ")+coef_xy(2,",lxy,")*pol_y(1,",lyp,",jg)"
      write(iunit,'(A,I0,A,I0,A,I0,A,I0,A)')"             coef_x(3,",lxp,")=coef_x(3,",lxp,&
                                                                            ")+coef_xy(1,",lxy,")*pol_y(2,",lyp,",jg)"
      write(iunit,'(A,I0,A,I0,A,I0,A,I0,A)')"             coef_x(4,",lxp,")=coef_x(4,",lxp,&
                                                                            ")+coef_xy(2,",lxy,")*pol_y(2,",lyp,",jg)"
    CASE(1)
      write(iunit,'(A,I0,A,I0,A,I0,A,I0,A)')"             coef_x(1:2,",lxp,")=coef_x(1:2,",lxp,&
                                                                            ")+coef_xy(1:2,",lxy,")*pol_y(1,",lyp,",jg)"
      write(iunit,'(A,I0,A,I0,A,I0,A,I0,A)')"             coef_x(3:4,",lxp,")=coef_x(3:4,",lxp,&
                                                                            ")+coef_xy(1:2,",lxy,")*pol_y(2,",lyp,",jg)"
    END SELECT
  END SELECT
END SUBROUTINE


!***********************************************************************************************************************************
!
! write_ig_loop
!
!***********************************************************************************************************************************

SUBROUTINE write_ig_loop()

  INTEGER :: lxp

  write(iunit,'(A)')"          DO ig=igmin,igmax"
  write(iunit,'(A)')"             i=map(ig,1)"

  SELECT CASE(options%ig_loop_vector_notation)
  CASE(0)
     write(iunit,'(A)')"             s01=0.0_dp"
     write(iunit,'(A)')"             s02=0.0_dp"
     write(iunit,'(A)')"             s03=0.0_dp"
     write(iunit,'(A)')"             s04=0.0_dp"
  CASE(1)
     write(iunit,'(A)')"             s(:)=0.0_dp"
  END SELECT

  SELECT CASE (options%ig_loop_unroll_lxp)
  CASE (0)
    write(iunit,'(A)')"             DO lxp=0,lp"
    CALL write_ig_loop_lxp_body(lxp)
    write(iunit,'(A)')"             ENDDO"
  CASE (1)
    DO lxp=0,lp
    CALL write_ig_loop_lxp_body(lxp)
    ENDDO
  END SELECT


  SELECT CASE(options%ig_loop_vector_notation)
  CASE(0)
    write(iunit,'(A)')"             grid(i,j,k) = grid(i,j,k)     + s01 "
    write(iunit,'(A)')"             grid(i,j2,k) = grid(i,j2,k)   + s03 "
    write(iunit,'(A)')"             grid(i,j,k2) = grid(i,j,k2)   + s02 "
    write(iunit,'(A)')"             grid(i,j2,k2) = grid(i,j2,k2) + s04 "
  CASE(1)
    write(iunit,'(A)')"             grid(i,j,k) = grid(i,j,k)     + s(1) "
    write(iunit,'(A)')"             grid(i,j2,k) = grid(i,j2,k)   + s(3) "
    write(iunit,'(A)')"             grid(i,j,k2) = grid(i,j,k2)   + s(2) "
    write(iunit,'(A)')"             grid(i,j2,k2) = grid(i,j2,k2) + s(4) "
  END SELECT

  write(iunit,'(A)')"          END DO"

END SUBROUTINE write_ig_loop

SUBROUTINE write_ig_loop_lxp_body(lxp)
  INTEGER :: lxp

  SELECT CASE (options%ig_loop_unroll_lxp)
  CASE (0)
    SELECT CASE (options%ig_loop_vector_notation)
    CASE(0)
      write(iunit,'(A)')"                s01=s01+coef_x(1,lxp)*pol_x(lxp,ig)"
      write(iunit,'(A)')"                s02=s02+coef_x(2,lxp)*pol_x(lxp,ig)"
      write(iunit,'(A)')"                s03=s03+coef_x(3,lxp)*pol_x(lxp,ig)"
      write(iunit,'(A)')"                s04=s04+coef_x(4,lxp)*pol_x(lxp,ig)"
    CASE(1)
      write(iunit,'(A)')"                s(:)=s(:)+coef_x(:,lxp)*pol_x(lxp,ig)"
    END SELECT
  CASE (1)
    SELECT CASE (options%ig_loop_vector_notation)
    CASE(0)
      write(iunit,'(A,I0,A,I0,A)')"                s01=s01+coef_x(1,",lxp,")*pol_x(",lxp,",ig)"
      write(iunit,'(A,I0,A,I0,A)')"                s02=s02+coef_x(2,",lxp,")*pol_x(",lxp,",ig)"
      write(iunit,'(A,I0,A,I0,A)')"                s03=s03+coef_x(3,",lxp,")*pol_x(",lxp,",ig)"
      write(iunit,'(A,I0,A,I0,A)')"                s04=s04+coef_x(4,",lxp,")*pol_x(",lxp,",ig)"
    CASE(1)
      write(iunit,'(A,I0,A,I0,A)')"                s(:)=s(:)+coef_x(:,",lxp,")*pol_x(",lxp,",ig)"
    END SELECT
  END SELECT
END SUBROUTINE


END SUBROUTINE collocate_generate_specialized



!***********************************************************************************************************************************
!
!
!
! this writes the default code for a given l quantum number
!
!
!
!  2012.Jun  RR: if lp_stub > 0 , the name of the routine is changed to collocate_core_ + lp_stub
!***********************************************************************************************************************************

SUBROUTINE collocate_generate_default(iunit)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iunit
!collocate_core_,lp

write(iunit,'(A)')"  SUBROUTINE collocate_core_default(grid,coef_xyz,pol_x,pol_y,pol_z,map,sphere_bounds,lp,cmax,gridbounds)"
write(iunit,'(A)')""
write(iunit,'(A)')"    IMPLICIT NONE"
write(iunit,'(A)')""
write(iunit,'(A)')"  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND ( 14, 200 )"
write(iunit,'(A)')"    INTEGER, INTENT(IN)     :: lp"
write(iunit,'(A)')"    INTEGER, INTENT(IN)     :: cmax"
write(iunit,'(A)')"    INTEGER, INTENT(IN)     :: gridbounds(2,3)"
write(iunit,'(A)')"    INTEGER, INTENT(IN)     :: sphere_bounds(*)"
write(iunit,'(A)')"    INTEGER, INTENT(IN)     :: map(-cmax:cmax,1:3)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(IN)    :: pol_x(0:lp,-cmax:cmax)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(IN)    :: pol_y(1:2,0:lp,-cmax:0)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(IN)    :: pol_z(1:2,0:lp,-cmax:0)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(IN)    :: coef_xyz(((lp+1)*(lp+2)*(lp+3))/6)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(INOUT) :: grid(gridbounds(1,1):gridbounds(2,1), & "
write(iunit,'(A)')"                                    gridbounds(1,2):gridbounds(2,2), &"
write(iunit,'(A)')"                                    gridbounds(1,3):gridbounds(2,3))"
write(iunit,'(A)')""
write(iunit,'(A)')"    REAL(KIND=dp) ::  coef_xy(2,(lp+1)*(lp+2)/2)"
write(iunit,'(A)')"    REAL(KIND=dp) ::  coef_x(4,0:lp)"
write(iunit,'(A)')""
write(iunit,'(A)')"    INTEGER kg,k,jgmin,jg,j,j2,igmax,ig,i,kgmin,igmin,k2,i2,jg2,kg2"
write(iunit,'(A)')"    INTEGER sci,lxp,lyp,lzp"
write(iunit,'(A)')"    REAL(KIND=dp) s01,s02,s03,s04,s(4)"
write(iunit,'(A)')"    INTEGER lxyz,lxy,lx"
write(iunit,'(A)')""
write(iunit,'(A)')"    sci=1"
write(iunit,'(A)')""
write(iunit,'(A)')"    kgmin=sphere_bounds(sci)"
write(iunit,'(A)')"    sci=sci+1"
write(iunit,'(A)')"    DO kg=kgmin,0"
write(iunit,'(A)')"       kg2=1-kg"
write(iunit,'(A)')"       k=map(kg,3)"
write(iunit,'(A)')"       k2=map(kg2,3)"
write(iunit,'(A)')""
write(iunit,'(A)')"       coef_xy=0.0_dp"
write(iunit,'(A)')"       lxyz = 0"
write(iunit,'(A)')"       DO lzp=0,lp"
write(iunit,'(A)')"          lxy=0"
write(iunit,'(A)')"          DO lyp=0,lp-lzp"
write(iunit,'(A)')"             DO lxp=0,lp-lzp-lyp"
write(iunit,'(A)')"                lxyz=lxyz+1 ; lxy=lxy+1"
write(iunit,'(A)')"                coef_xy(1,lxy)=coef_xy(1,lxy)+coef_xyz(lxyz)*pol_z(1,lzp,kg)"
write(iunit,'(A)')"                coef_xy(2,lxy)=coef_xy(2,lxy)+coef_xyz(lxyz)*pol_z(2,lzp,kg)"
write(iunit,'(A)')"             ENDDO"
write(iunit,'(A)')"             lxy=lxy+lzp"
write(iunit,'(A)')"          ENDDO"
write(iunit,'(A)')"       ENDDO"
write(iunit,'(A)')""
write(iunit,'(A)')"       jgmin=sphere_bounds(sci)"
write(iunit,'(A)')"       sci=sci+1"
write(iunit,'(A)')"       DO jg=jgmin,0"
write(iunit,'(A)')"          jg2=1-jg"
write(iunit,'(A)')"          j=map(jg,2)"
write(iunit,'(A)')"          j2=map(jg2,2)"
write(iunit,'(A)')"          igmin=sphere_bounds(sci)"
write(iunit,'(A)')"          sci=sci+1"
write(iunit,'(A)')"          igmax=1-igmin"
write(iunit,'(A)')""
write(iunit,'(A)')"          coef_x=0.0_dp"
write(iunit,'(A)')"          lxy=0"
write(iunit,'(A)')"          DO lyp=0,lp"
write(iunit,'(A)')"          DO lxp=0,lp-lyp"
write(iunit,'(A)')"             lxy=lxy+1"
write(iunit,'(A)')"             coef_x(1,lxp)=coef_x(1,lxp)+coef_xy(1,lxy)*pol_y(1,lyp,jg)"
write(iunit,'(A)')"             coef_x(2,lxp)=coef_x(2,lxp)+coef_xy(2,lxy)*pol_y(1,lyp,jg)"
write(iunit,'(A)')"             coef_x(3,lxp)=coef_x(3,lxp)+coef_xy(1,lxy)*pol_y(2,lyp,jg)"
write(iunit,'(A)')"             coef_x(4,lxp)=coef_x(4,lxp)+coef_xy(2,lxy)*pol_y(2,lyp,jg)"
write(iunit,'(A)')"          ENDDO"
write(iunit,'(A)')"          ENDDO"
write(iunit,'(A)')""
write(iunit,'(A)')"          DO ig=igmin,igmax"
write(iunit,'(A)')"             i=map(ig,1)"
write(iunit,'(A)')"             s01=0.0_dp"
write(iunit,'(A)')"             s02=0.0_dp"
write(iunit,'(A)')"             s03=0.0_dp"
write(iunit,'(A)')"             s04=0.0_dp"
write(iunit,'(A)')"             DO lxp=0,lp"
write(iunit,'(A)')"                s01=s01+coef_x(1,lxp)*pol_x(lxp,ig)"
write(iunit,'(A)')"                s02=s02+coef_x(2,lxp)*pol_x(lxp,ig)"
write(iunit,'(A)')"                s03=s03+coef_x(3,lxp)*pol_x(lxp,ig)"
write(iunit,'(A)')"                s04=s04+coef_x(4,lxp)*pol_x(lxp,ig)"
write(iunit,'(A)')"             ENDDO"
write(iunit,'(A)')"             grid(i,j,k) = grid(i,j,k)     + s01 "
write(iunit,'(A)')"             grid(i,j2,k) = grid(i,j2,k)   + s03 "
write(iunit,'(A)')"             grid(i,j,k2) = grid(i,j,k2)   + s02 "
write(iunit,'(A)')"             grid(i,j2,k2) = grid(i,j2,k2) + s04 "
write(iunit,'(A)')"          END DO"
write(iunit,'(A)')""
write(iunit,'(A)')"       END DO"
write(iunit,'(A)')"    END DO"
write(iunit,'(A)')"  "
write(iunit,'(A)')"  END SUBROUTINE collocate_core_default"

END SUBROUTINE collocate_generate_default

SUBROUTINE collocate_generate_case_statement(iunit,lmax)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: iunit,lmax

   INTEGER :: lp
   
   write(iunit,'(A)')"  SELECT CASE(lp)"
   DO lp=0,lmax
      write(iunit,'(A,I0,A)')" CASE(",lp,")"
      write(iunit,'(A,I0,A)')" CALL collocate_core_",lp,"(grid(1,1,1),coef_xyz(1),pol_x(0,-cmax),"//&
                                                 " pol_y(1,0,-cmax),pol_z(1,0,-cmax), &"
      write(iunit,'(A)')"                        map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))"
   ENDDO
   write(iunit,'(A)')"  CASE DEFAULT"
   write(iunit,'(A)')" CALL collocate_core_default(grid(1,1,1),coef_xyz(1),pol_x(0,-cmax),pol_y(1,0,-cmax),pol_z(1,0,-cmax), &"
   write(iunit,'(A)')"                        map(-cmax,1),sphere_bounds(1),lp,cmax,gridbounds(1,1))"
   write(iunit,'(A)')"  END SELECT"
END SUBROUTINE collocate_generate_case_statement

END MODULE collocate_generate

!***********************************************************************************************************************************
!
!
!
!       INTEGRATE INTEGRATE INTEGRATE INTEGRATE INTEGRATE INTEGRATE INTEGRATE INTEGRATE INTEGRATE INTEGRATE INTEGRATE INTEGRATE
!
!
!
!***********************************************************************************************************************************

MODULE integrate_generate
  USE option_module
  IMPLICIT NONE

CONTAINS

SUBROUTINE integrate_generate_stub(iunit,lp)
  INTEGER, INTENT(IN) :: iunit,lp

write(iunit,'(A,I0,A)')"  SUBROUTINE integrate_core_",lp,"(grid,coef_xyz,pol_x,pol_y,pol_z,map,sphere_bounds,cmax,gridbounds)"
write(iunit,'(A)')""
write(iunit,'(A)')"    IMPLICIT NONE"
write(iunit,'(A)')""
write(iunit,'(A)')"  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND ( 14, 200 )"
write(iunit,'(A,I0)')"    INTEGER, PARAMETER      :: lp=",lp
write(iunit,'(A)')"    INTEGER, INTENT(IN)     :: cmax"
write(iunit,'(A)')"    INTEGER, INTENT(IN)     :: gridbounds(2,3)"
write(iunit,'(A)')"    INTEGER, INTENT(IN)     :: sphere_bounds(*)"
write(iunit,'(A)')"    INTEGER, INTENT(IN)     :: map(-cmax:cmax,1:3)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(IN)    :: pol_x(0:lp,-cmax:cmax)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(IN)    :: pol_y(1:2,0:lp,-cmax:0)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(IN)    :: pol_z(1:2,0:lp,-cmax:0)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(OUT)   :: coef_xyz(((lp+1)*(lp+2)*(lp+3))/6)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(IN)    :: grid(gridbounds(1,1):gridbounds(2,1), & "
write(iunit,'(A)')"                                    gridbounds(1,2):gridbounds(2,2), &"
write(iunit,'(A)')"                                    gridbounds(1,3):gridbounds(2,3))"
write(iunit,'(A)')" CALL integrate_core_default(grid,coef_xyz,pol_x,pol_y, &"
write(iunit,'(A)')                      "pol_z,map,sphere_bounds, &"
write(iunit,'(I0,A)') lp,",cmax,gridbounds)"
write(iunit,'(A,I0,A)')"END SUBROUTINE integrate_core_",lp," "

END SUBROUTINE integrate_generate_stub

SUBROUTINE integrate_generate_specialized(iunit,lp,options)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iunit,lp
  TYPE (option_type)  :: options

write(iunit,'(A,I0,A)')"  SUBROUTINE integrate_core_",lp,"(grid,coef_xyz,pol_x,pol_y,pol_z,map,sphere_bounds,cmax,gridbounds)"
write(iunit,'(A)')""
write(iunit,'(A)')"    IMPLICIT NONE"
write(iunit,'(A)')""
write(iunit,'(A)')"  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND ( 14, 200 )"
write(iunit,'(A,I0)')"    INTEGER, PARAMETER      :: lp=",lp
write(iunit,'(A)')"    INTEGER, INTENT(IN)     :: cmax"
write(iunit,'(A)')"    INTEGER, INTENT(IN)     :: gridbounds(2,3)"
write(iunit,'(A)')"    INTEGER, INTENT(IN)     :: sphere_bounds(*)"
write(iunit,'(A)')"    INTEGER, INTENT(IN)     :: map(-cmax:cmax,1:3)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(IN)    :: pol_x(0:lp,-cmax:cmax)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(IN)    :: pol_y(1:2,0:lp,-cmax:0)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(IN)    :: pol_z(1:2,0:lp,-cmax:0)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(OUT)   :: coef_xyz(((lp+1)*(lp+2)*(lp+3))/6)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(IN)    :: grid(gridbounds(1,1):gridbounds(2,1), & "
write(iunit,'(A)')"                                    gridbounds(1,2):gridbounds(2,2), &"
write(iunit,'(A)')"                                    gridbounds(1,3):gridbounds(2,3))"
write(iunit,'(A)')""
write(iunit,'(A)')"    REAL(KIND=dp) ::  coef_xy(2,((lp+1)*(lp+2))/2)"
write(iunit,'(A)')"    REAL(KIND=dp) ::  coef_x(4,0:lp)"
write(iunit,'(A)')""
write(iunit,'(A)')"    INTEGER kg,k,jgmin,jg,j,j2,igmax,ig,i,kgmin,igmin,k2,i2,jg2,kg2"
write(iunit,'(A)')"    INTEGER sci,lxp,lyp,lzp"
write(iunit,'(A)')"    REAL(KIND=dp) s01,s02,s03,s04,s(4)"
write(iunit,'(A)')"    INTEGER lxyz,lxy,lx"
write(iunit,'(A)')""
write(iunit,'(A)')"    sci=1"
write(iunit,'(A)')""
write(iunit,'(A)')"    coef_xyz=0.0_dp"
write(iunit,'(A)')""
write(iunit,'(A)')"    kgmin=sphere_bounds(sci)"
write(iunit,'(A)')"    sci=sci+1"
write(iunit,'(A)')"    DO kg=kgmin,0"
write(iunit,'(A)')"       kg2=1-kg"
write(iunit,'(A)')"       k=map(kg,3)"
write(iunit,'(A)')"       k2=map(kg2,3)"
write(iunit,'(A)')""
write(iunit,'(A)')"       coef_xy=0.0_dp"
write(iunit,'(A)')""
write(iunit,'(A)')"       jgmin=sphere_bounds(sci)"
write(iunit,'(A)')"       sci=sci+1"
write(iunit,'(A)')"       DO jg=jgmin,0"
write(iunit,'(A)')"          jg2=1-jg"
write(iunit,'(A)')"          j=map(jg,2)"
write(iunit,'(A)')"          j2=map(jg2,2)"
write(iunit,'(A)')"          igmin=sphere_bounds(sci)"
write(iunit,'(A)')"          sci=sci+1"
write(iunit,'(A)')"          igmax=1-igmin"
write(iunit,'(A)')""
write(iunit,'(A)')"          coef_x=0.0_dp"

CALL write_ig_loop()

CALL write_jg_body()

write(iunit,'(A)')"       END DO"

CALL write_kg_body

write(iunit,'(A)')"    END DO"
write(iunit,'(A)')"  "
write(iunit,'(A,I0)')"  END SUBROUTINE integrate_core_",lp

CONTAINS

SUBROUTINE write_ig_loop()
  INTEGER :: lxp

  write(iunit,'(A)')"          DO ig=igmin,igmax"
  write(iunit,'(A)')"             i=map(ig,1)"

  SELECT CASE(options%ig_loop_vector_notation)
  CASE(0)
    write(iunit,'(A)')"             s01=grid(i,j,k)"
    write(iunit,'(A)')"             s02=grid(i,j,k2)"
    write(iunit,'(A)')"             s03=grid(i,j2,k)"
    write(iunit,'(A)')"             s04=grid(i,j2,k2)"
  CASE(1)
    write(iunit,'(A)')"             s(1)=grid(i,j,k)"
    write(iunit,'(A)')"             s(2)=grid(i,j,k2)"
    write(iunit,'(A)')"             s(3)=grid(i,j2,k)"
    write(iunit,'(A)')"             s(4)=grid(i,j2,k2)"
  END SELECT

  SELECT CASE (options%ig_loop_unroll_lxp)
  CASE (0)
    write(iunit,'(A)')"             DO lxp=0,lp"
    CALL write_ig_loop_lxp_body(lxp)
    write(iunit,'(A)')"             ENDDO"
  CASE (1)
    DO lxp=0,lp
    CALL write_ig_loop_lxp_body(lxp)
    ENDDO
  END SELECT

  write(iunit,'(A)')"          END DO"
END SUBROUTINE write_ig_loop 

SUBROUTINE write_ig_loop_lxp_body(lxp)
  INTEGER :: lxp

  SELECT CASE (options%ig_loop_unroll_lxp)
  CASE (0)
    SELECT CASE (options%ig_loop_vector_notation)
    CASE(0)
      write(iunit,'(A)')"                coef_x(1,lxp)=coef_x(1,lxp)+s01*pol_x(lxp,ig)"
      write(iunit,'(A)')"                coef_x(2,lxp)=coef_x(2,lxp)+s02*pol_x(lxp,ig)"
      write(iunit,'(A)')"                coef_x(3,lxp)=coef_x(3,lxp)+s03*pol_x(lxp,ig)"
      write(iunit,'(A)')"                coef_x(4,lxp)=coef_x(4,lxp)+s04*pol_x(lxp,ig)"
    CASE(1)
      write(iunit,'(A)')"                coef_x(:,lxp)=coef_x(:,lxp)+s(:)*pol_x(lxp,ig)"
    END SELECT
  CASE (1)
    SELECT CASE (options%ig_loop_vector_notation)
    CASE(0)
      write(iunit,'(A,I0,A,I0,A,I0,A)')"                coef_x(1,",lxp,")=coef_x(1,",lxp,")+s01*pol_x(",lxp,",ig)"
      write(iunit,'(A,I0,A,I0,A,I0,A)')"                coef_x(2,",lxp,")=coef_x(2,",lxp,")+s02*pol_x(",lxp,",ig)"
      write(iunit,'(A,I0,A,I0,A,I0,A)')"                coef_x(3,",lxp,")=coef_x(3,",lxp,")+s03*pol_x(",lxp,",ig)"
      write(iunit,'(A,I0,A,I0,A,I0,A)')"                coef_x(4,",lxp,")=coef_x(4,",lxp,")+s04*pol_x(",lxp,",ig)"
    CASE(1)
      write(iunit,'(A,I0,A,I0,A,I0,A)')"                coef_x(:,",lxp,")=coef_x(:,",lxp,")+s(:)*pol_x(",lxp,",ig)"
    END SELECT
  END SELECT
END SUBROUTINE write_ig_loop_lxp_body

SUBROUTINE write_jg_body()
  INTEGER :: lxp,lyp,lxy

  SELECT CASE (options%jg_loop_unroll_lxp)
  CASE (0)
    write(iunit,'(A)')"          lxy=0"
    write(iunit,'(A)')"          DO lyp=0,lp"
    write(iunit,'(A)')"          DO lxp=0,lp-lyp"
    write(iunit,'(A)')"             lxy=lxy+1"
    CALL write_jg_lxy_body(lyp,lxp,lxy)
    write(iunit,'(A)')"          ENDDO"
    write(iunit,'(A)')"          ENDDO"
  CASE (1)
    lxy=0
    DO lyp=0,lp
    DO lxp=0,lp-lyp
       lxy=lxy+1
       CALL write_jg_lxy_body(lyp,lxp,lxy)
    ENDDO
    ENDDO
  END SELECT


END SUBROUTINE write_jg_body

SUBROUTINE write_jg_lxy_body(lyp,lxp,lxy)
  INTEGER :: lxp,lyp,lxy

  SELECT CASE (options%jg_loop_unroll_lxp)
  CASE (0)
    SELECT CASE (options%jg_loop_vector_notation)
    CASE(0)
      write(iunit,'(A)')"             coef_xy(1,lxy)=coef_xy(1,lxy)+coef_x(1,lxp)*pol_y(1,lyp,jg)"
      write(iunit,'(A)')"             coef_xy(2,lxy)=coef_xy(2,lxy)+coef_x(2,lxp)*pol_y(1,lyp,jg)"
      write(iunit,'(A)')"             coef_xy(1,lxy)=coef_xy(1,lxy)+coef_x(3,lxp)*pol_y(2,lyp,jg)"
      write(iunit,'(A)')"             coef_xy(2,lxy)=coef_xy(2,lxy)+coef_x(4,lxp)*pol_y(2,lyp,jg)"
    CASE(1)
      write(iunit,'(A)')"             coef_xy(:,lxy)=coef_xy(:,lxy)+coef_x(1:2,lxp)*pol_y(1,lyp,jg)"
      write(iunit,'(A)')"             coef_xy(:,lxy)=coef_xy(:,lxy)+coef_x(3:4,lxp)*pol_y(2,lyp,jg)"
    END SELECT
  CASE (1)
    SELECT CASE (options%jg_loop_vector_notation)
    CASE(0)
      write(iunit,'(A,I0,A,I0,A,I0,A,I0,A)')"             coef_xy(1,",lxy,")=coef_xy(1,",lxy,")+"//&
                                                                  "coef_x(1,",lxp,")*pol_y(1,",lyp,",jg)"
      write(iunit,'(A,I0,A,I0,A,I0,A,I0,A)')"             coef_xy(2,",lxy,")=coef_xy(2,",lxy,")+"//&
                                                                  "coef_x(2,",lxp,")*pol_y(1,",lyp,",jg)"
      write(iunit,'(A,I0,A,I0,A,I0,A,I0,A)')"             coef_xy(1,",lxy,")=coef_xy(1,",lxy,")+"//&
                                                                  "coef_x(3,",lxp,")*pol_y(2,",lyp,",jg)"
      write(iunit,'(A,I0,A,I0,A,I0,A,I0,A)')"             coef_xy(2,",lxy,")=coef_xy(2,",lxy,")+"//&
                                                                  "coef_x(4,",lxp,")*pol_y(2,",lyp,",jg)"
    CASE(1)
      write(iunit,'(A,I0,A,I0,A,I0,A,I0,A)')"             coef_xy(:,",lxy,")=coef_xy(:,",lxy,")+"//&
                                                                  "coef_x(1:2,",lxp,")*pol_y(1,",lyp,",jg)"
      write(iunit,'(A,I0,A,I0,A,I0,A,I0,A)')"             coef_xy(:,",lxy,")=coef_xy(:,",lxy,")+"//&
                                                                  "coef_x(3:4,",lxp,")*pol_y(2,",lyp,",jg)"
    END SELECT
  END SELECT
END SUBROUTINE

SUBROUTINE write_kg_body()
  INTEGER :: lxp,lyp,lzp,lxy,lxyz

  SELECT CASE (options%kg_loop_unroll_lxp)
  CASE (0)
    write(iunit,'(A)')"       lxyz = 0"
    write(iunit,'(A)')"       DO lzp=0,lp"
    write(iunit,'(A)')"          lxy=0"
    write(iunit,'(A)')"          DO lyp=0,lp-lzp"
    write(iunit,'(A)')"             DO lxp=0,lp-lzp-lyp"
    write(iunit,'(A)')"                lxyz=lxyz+1 ; lxy=lxy+1"
    CALL write_kg_lxy_body(lzp,lxy,lxyz)
    write(iunit,'(A)')"             ENDDO"
    write(iunit,'(A)')"             lxy=lxy+lzp"
    write(iunit,'(A)')"          ENDDO"
    write(iunit,'(A)')"       ENDDO"
  CASE (1)
    lxyz = 0
    DO lzp=0,lp
       lxy=0
       DO lyp=0,lp-lzp
          DO lxp=0,lp-lzp-lyp
             lxyz=lxyz+1 ; lxy=lxy+1
             CALL write_kg_lxy_body(lzp,lxy,lxyz)
          ENDDO
          lxy=lxy+lzp
       ENDDO
    ENDDO
  END SELECT
END SUBROUTINE

SUBROUTINE write_kg_lxy_body(lzp,lxy,lxyz)
  INTEGER :: lzp,lxy,lxyz
  SELECT CASE (options%kg_loop_unroll_lxp)
  CASE (0)
    SELECT CASE (options%kg_loop_vector_notation)
    CASE(0)
      write(iunit,'(A)')"                coef_xyz(lxyz)=coef_xyz(lxyz)+coef_xy(1,lxy)*pol_z(1,lzp,kg)"
      write(iunit,'(A)')"                coef_xyz(lxyz)=coef_xyz(lxyz)+coef_xy(2,lxy)*pol_z(2,lzp,kg)"
    CASE(1)
      write(iunit,'(A)')"                coef_xyz(lxyz)=coef_xyz(lxyz)+SUM(coef_xy(:,lxy)*pol_z(:,lzp,kg))"
    END SELECT
  CASE(1)
    SELECT CASE (options%kg_loop_vector_notation)
    CASE(0)
      write(iunit,'(A,I0,A,I0,A,I0,A,I0,A)')"                coef_xyz(",lxyz,")=coef_xyz(",lxyz,")+"//&
                                                               "coef_xy(1,",lxy,")*pol_z(1,",lzp,",kg)"
      write(iunit,'(A,I0,A,I0,A,I0,A,I0,A)')"                coef_xyz(",lxyz,")=coef_xyz(",lxyz,")+"//&
                                                               "coef_xy(2,",lxy,")*pol_z(2,",lzp,",kg)"
    CASE(1)
      write(iunit,'(A,I0,A,I0,A,I0,A,I0,A)')"                coef_xyz(",lxyz,")=coef_xyz(",lxyz,")+"//&
                                                               "SUM(coef_xy(:,",lxy,")*pol_z(:,",lzp,",kg))"
    END SELECT
  END SELECT
END SUBROUTINE

END SUBROUTINE integrate_generate_specialized


!***********************************************************************************************************************************
!
!
!
! this writes the default code for a given l quantum number
!
!
! 2012.Jun RR If lp_stub > 0 then name of routine is integrate_core_ + lp value
!***********************************************************************************************************************************

SUBROUTINE integrate_generate_default(iunit)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iunit 

write(iunit,'(A)')"  SUBROUTINE integrate_core_default(grid,coef_xyz,pol_x,pol_y,pol_z,map,sphere_bounds,lp,cmax,gridbounds)"
write(iunit,'(A)')""
write(iunit,'(A)')"    IMPLICIT NONE"
write(iunit,'(A)')""
write(iunit,'(A)')"  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND ( 14, 200 )"
write(iunit,'(A)')"    INTEGER, INTENT(IN)     :: lp"
write(iunit,'(A)')"    INTEGER, INTENT(IN)     :: cmax"
write(iunit,'(A)')"    INTEGER, INTENT(IN)     :: gridbounds(2,3)"
write(iunit,'(A)')"    INTEGER, INTENT(IN)     :: sphere_bounds(*)"
write(iunit,'(A)')"    INTEGER, INTENT(IN)     :: map(-cmax:cmax,1:3)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(IN)    :: pol_x(0:lp,-cmax:cmax)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(IN)    :: pol_y(1:2,0:lp,-cmax:0)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(IN)    :: pol_z(1:2,0:lp,-cmax:0)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(OUT)   :: coef_xyz(((lp+1)*(lp+2)*(lp+3))/6)"
write(iunit,'(A)')"    REAL(KIND=dp), INTENT(IN)    :: grid(gridbounds(1,1):gridbounds(2,1), & "
write(iunit,'(A)')"                                    gridbounds(1,2):gridbounds(2,2), &"
write(iunit,'(A)')"                                    gridbounds(1,3):gridbounds(2,3))"
write(iunit,'(A)')""
write(iunit,'(A)')"    REAL(KIND=dp) ::  coef_xy(2,((lp+1)*(lp+2))/2)"
write(iunit,'(A)')"    REAL(KIND=dp) ::  coef_x(4,0:lp)"
write(iunit,'(A)')""
write(iunit,'(A)')"    INTEGER kg,k,jgmin,jg,j,j2,igmax,ig,i,kgmin,igmin,k2,i2,jg2,kg2"
write(iunit,'(A)')"    INTEGER sci,lxp,lyp,lzp"
write(iunit,'(A)')"    REAL(KIND=dp) s01,s02,s03,s04,s(4)"
write(iunit,'(A)')"    INTEGER lxyz,lxy,lx"
write(iunit,'(A)')""
write(iunit,'(A)')"    sci=1"
write(iunit,'(A)')""
write(iunit,'(A)')"    coef_xyz=0.0_dp"
write(iunit,'(A)')""
write(iunit,'(A)')"    kgmin=sphere_bounds(sci)"
write(iunit,'(A)')"    sci=sci+1"
write(iunit,'(A)')"    DO kg=kgmin,0"
write(iunit,'(A)')"       kg2=1-kg"
write(iunit,'(A)')"       k=map(kg,3)"
write(iunit,'(A)')"       k2=map(kg2,3)"
write(iunit,'(A)')""
write(iunit,'(A)')"       coef_xy=0.0_dp"
write(iunit,'(A)')""
write(iunit,'(A)')"       jgmin=sphere_bounds(sci)"
write(iunit,'(A)')"       sci=sci+1"
write(iunit,'(A)')"       DO jg=jgmin,0"
write(iunit,'(A)')"          jg2=1-jg"
write(iunit,'(A)')"          j=map(jg,2)"
write(iunit,'(A)')"          j2=map(jg2,2)"
write(iunit,'(A)')"          igmin=sphere_bounds(sci)"
write(iunit,'(A)')"          sci=sci+1"
write(iunit,'(A)')"          igmax=1-igmin"
write(iunit,'(A)')""
write(iunit,'(A)')"          coef_x=0.0_dp"
write(iunit,'(A)')""
write(iunit,'(A)')"          DO ig=igmin,igmax"
write(iunit,'(A)')"             i=map(ig,1)"
write(iunit,'(A)')"             s01=grid(i,j,k)"
write(iunit,'(A)')"             s02=grid(i,j,k2)"
write(iunit,'(A)')"             s03=grid(i,j2,k)"
write(iunit,'(A)')"             s04=grid(i,j2,k2)"
write(iunit,'(A)')"             DO lxp=0,lp"
write(iunit,'(A)')"                coef_x(1,lxp)=coef_x(1,lxp)+s01*pol_x(lxp,ig)"
write(iunit,'(A)')"                coef_x(2,lxp)=coef_x(2,lxp)+s02*pol_x(lxp,ig)"
write(iunit,'(A)')"                coef_x(3,lxp)=coef_x(3,lxp)+s03*pol_x(lxp,ig)"
write(iunit,'(A)')"                coef_x(4,lxp)=coef_x(4,lxp)+s04*pol_x(lxp,ig)"
write(iunit,'(A)')"             ENDDO"
write(iunit,'(A)')"          END DO"
write(iunit,'(A)')""
write(iunit,'(A)')"          lxy=0"
write(iunit,'(A)')"          DO lyp=0,lp"
write(iunit,'(A)')"          DO lxp=0,lp-lyp"
write(iunit,'(A)')"             lxy=lxy+1"
write(iunit,'(A)')"             coef_xy(1,lxy)=coef_xy(1,lxy)+coef_x(1,lxp)*pol_y(1,lyp,jg)"
write(iunit,'(A)')"             coef_xy(2,lxy)=coef_xy(2,lxy)+coef_x(2,lxp)*pol_y(1,lyp,jg)"
write(iunit,'(A)')"             coef_xy(1,lxy)=coef_xy(1,lxy)+coef_x(3,lxp)*pol_y(2,lyp,jg)"
write(iunit,'(A)')"             coef_xy(2,lxy)=coef_xy(2,lxy)+coef_x(4,lxp)*pol_y(2,lyp,jg)"
write(iunit,'(A)')"          ENDDO"
write(iunit,'(A)')"          ENDDO"
write(iunit,'(A)')""
write(iunit,'(A)')"       END DO"
write(iunit,'(A)')""
write(iunit,'(A)')"       lxyz = 0"
write(iunit,'(A)')"       DO lzp=0,lp"
write(iunit,'(A)')"          lxy=0"
write(iunit,'(A)')"          DO lyp=0,lp-lzp"
write(iunit,'(A)')"             DO lxp=0,lp-lzp-lyp"
write(iunit,'(A)')"                lxyz=lxyz+1 ; lxy=lxy+1"
write(iunit,'(A)')"                coef_xyz(lxyz)=coef_xyz(lxyz)+coef_xy(1,lxy)*pol_z(1,lzp,kg)"
write(iunit,'(A)')"                coef_xyz(lxyz)=coef_xyz(lxyz)+coef_xy(2,lxy)*pol_z(2,lzp,kg)"
write(iunit,'(A)')"             ENDDO"
write(iunit,'(A)')"             lxy=lxy+lzp"
write(iunit,'(A)')"          ENDDO"
write(iunit,'(A)')"       ENDDO"
write(iunit,'(A)')""
write(iunit,'(A)')"    END DO"
write(iunit,'(A)')"  "
write(iunit,'(A)')"  END SUBROUTINE integrate_core_default"

END SUBROUTINE integrate_generate_default

SUBROUTINE integrate_generate_case_statement(iunit,lmax)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: iunit,lmax

   INTEGER :: lp
   
   write(iunit,'(A)')"  SELECT CASE(lp)"
   DO lp=0,lmax
      write(iunit,'(A,I0,A)')" CASE(",lp,")"
      write(iunit,'(A,I0,A)')" CALL integrate_core_",lp,"(grid(1,1,1),coef_xyz(1),pol_x(0,-cmax),"//&
                                                 " pol_y(1,0,-cmax),pol_z(1,0,-cmax), &"
      write(iunit,'(A)')"                        map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))"
   ENDDO
   write(iunit,'(A)')"  CASE DEFAULT"
   write(iunit,'(A)')" CALL integrate_core_default(grid(1,1,1),coef_xyz(1),pol_x(0,-cmax),pol_y(1,0,-cmax),pol_z(1,0,-cmax), &"
   write(iunit,'(A)')"                        map(-cmax,1),sphere_bounds(1),lp,cmax,gridbounds(1,1))"
   write(iunit,'(A)')"  END SELECT"
END SUBROUTINE integrate_generate_case_statement

END MODULE integrate_generate

!===================================================================================================================================
!
!
! CODE GENERATION  CODE GENERATION  CODE GENERATION  CODE GENERATION  CODE GENERATION  CODE GENERATION  CODE GENERATION  
!
!
!===================================================================================================================================

SUBROUTINE  generate_collocate(do_best)
  USE option_module
  USE collocate_generate
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: do_best

  INTEGER, PARAMETER :: grid_fast_unit=10, call_collocate_unit=11
  CHARACTER(LEN=200) :: directory_name, file_name
  INTEGER :: lmax,lp,Nopts,iopt
  TYPE(option_type), DIMENSION(:), ALLOCATABLE :: options


  CALL generate_all_options()
  Nopts=SIZE(all_options)

  WRITE(6,*) "Collocate: generate up to lmax (= lmax_a + lmax_b)"
  READ(5,*)  lmax
  WRITE(6,*) "Using lmax = ",lmax
  ALLOCATE(options(0:lmax))
  DO lp=0,lmax
     WRITE(6,*) "Which option to be used for code generation [1-",Nopts,"] for l=",lp
     READ(5,*)  iopt
     WRITE(6,*) "Using opt = ",iopt," for l=",lp
     options(lp)=all_options(iopt)
     WRITE(6,*) all_options(iopt)%ig_loop_unroll_lxp

  ENDDO

  IF (lmax .gt. -1) THEN
   IF (do_best .eq. 1) THEN
     write(file_name,"(A)") "out_best/call_collocate.f90"
   ELSE 
     write(file_name,"(A,I0,A,I0,A)") "out_",lmax,"_",iopt,"/call_collocate.f90"
   ENDIF
   write(6,*) "Writting collocate to file ",file_name

  OPEN(UNIT=call_collocate_unit,FILE=file_name)
  CALL collocate_generate_case_statement(call_collocate_unit,lmax)
  CLOSE(call_collocate_unit)

  IF (do_best .eq. 1) THEN
     write(file_name,"(A)") "out_best/collocate_fast.F"
  ELSE
     write(file_name,"(A,I0,A,I0,A)") "out_",lmax,"_",iopt,"/collocate_fast.F"
  ENDIF


  OPEN(UNIT=grid_fast_unit,FILE=file_name)

  CALL collocate_generate_default(grid_fast_unit)
  DO lp=0,lmax
     write(6,*) "Generate l=",lp,"using ",options(lp)
     CALL collocate_generate_specialized(grid_fast_unit,lp,options(lp))
  ENDDO
  IF (do_best .eq. 1) THEN
     ! There should be at least 9 routines
     DO lp=(lmax+1),DEFAULT_LMAX
         CALL collocate_generate_stub(grid_fast_unit,lp)
     END DO
  ENDIF

  CLOSE(grid_fast_unit)
  ENDIF

  CALL deallocate_all_options()

END SUBROUTINE generate_collocate

SUBROUTINE  generate_integrate(do_best)
  USE option_module
  USE integrate_generate
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: do_best

  INTEGER, PARAMETER :: grid_fast_unit=10, call_integrate_unit=11
  CHARACTER(LEN=200) :: directory_name, file_name

  INTEGER :: lmax,lp,Nopts,iopt
  TYPE(option_type), DIMENSION(:), ALLOCATABLE :: options

  CALL generate_all_options()
  Nopts=SIZE(all_options)

  WRITE(6,*) "Integrate: generate up to lmax (= lmax_a + lmax_b)"
  READ(5,*)  lmax
  WRITE(6,*) "Using lmax = ",lmax
  ALLOCATE(options(0:lmax))
  DO lp=0,lmax
     WRITE(6,*) "Which option to be used for code generation [1-",Nopts,"] for l=",lp
     READ(5,*)  iopt
     WRITE(6,*) "Using opt = ",iopt," for l=",lp
     options(lp)=all_options(iopt)
  ENDDO


  IF (lmax .gt. -1) THEN 

   IF (do_best .eq. 1) THEN
    write(file_name,"(A)") "out_best/call_integrate.f90"
   ELSE
    write(file_name,"(A,I0,A,I0,A)") "out_",lmax,"_",iopt,"/call_integrate.f90"
   ENDIF

   write(6,*) "Writting integrate to file ",file_name

  OPEN(UNIT=call_integrate_unit,FILE=file_name)
  CALL integrate_generate_case_statement(call_integrate_unit,lmax)
  CLOSE(call_integrate_unit)
 
  IF (do_best .eq. 1) THEN
    write(file_name,"(A)") "out_best/integrate_fast.F"
  ELSE
    write(file_name,"(A,I0,A,I0,A)") "out_",lmax,"_",iopt,"/integrate_fast.F"
  ENDIF

  OPEN(UNIT=grid_fast_unit,FILE=file_name)

  CALL integrate_generate_default(grid_fast_unit)
  DO lp=0,lmax
     write(6,*) "Generate l=",lp,"using ",options(lp)
     CALL integrate_generate_specialized(grid_fast_unit,lp,options(lp))
  ENDDO

  IF (do_best .eq. 1) THEN
     ! There should be at least 9 routines
     DO lp=(lmax+1),DEFAULT_LMAX
         CALL integrate_generate_stub(grid_fast_unit,lp)
     END DO
  ENDIF


  CLOSE(grid_fast_unit)

  ENDIF

  CALL deallocate_all_options()

END SUBROUTINE generate_integrate



PROGRAM generate_fast
   INTEGER :: do_best
          
   ! If program receives an argument it will generate the output in out_best directory
   IF (iargc() > 0) THEN
      write(6,*) "Generate best combination" 
      do_best = 1
   ELSE
      do_best = 0
   ENDIF

   CALL generate_collocate(do_best)

   CALL generate_integrate(do_best)

END PROGRAM generate_fast
