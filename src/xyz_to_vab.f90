
! *****************************************************************************
!> \brief provides specific implementations for common cases of xyz_to_vab routine.
!> 
!> \note
!>     ____              _ _     __  __           _ _  __         _____ _     _       _____ _ _      _
!>    |  _ \  ___  _ __ ( ) |_  |  \/  | ___   __| (_)/ _|_   _  |_   _| |__ (_)___  |  ___(_) | ___| |
!>    | | | |/ _ \| '_ \|/| __| | |\/| |/ _ \ / _` | | |_| | | |   | | | '_ \| / __| | |_  | | |/ _ \ |
!>    | |_| | (_) | | | | | |_  | |  | | (_) | (_| | |  _| |_| |   | | | | | | \__ \ |  _| | | |  __/_|
!>    |____/ \___/|_| |_|  \__| |_|  |_|\___/ \__,_|_|_|  \__, |   |_| |_| |_|_|___/ |_|   |_|_|\___(_)
!>                                                        |___/
!>      ____ _                  ___                              _ _       _       _
!>     / ___| | ___  ___  ___  |_ _|_ __ ___  _ __ ___   ___  __| (_) __ _| |_ ___| |_   _
!>    | |   | |/ _ \/ __|/ _ \  | || '_ ` _ \| '_ ` _ \ / _ \/ _` | |/ _` | __/ _ \ | | | |
!>    | |___| | (_) \__ \  __/  | || | | | | | | | | | |  __/ (_| | | (_| | ||  __/ | |_| |
!>     \____|_|\___/|___/\___| |___|_| |_| |_|_| |_| |_|\___|\__,_|_|\__,_|\__\___|_|\__, |
!>                                                                                   |___/
!>     _____ _     _       _____ _ _      _
!>    |_   _| |__ (_)___  |  ___(_) | ___| |
!>      | | | '_ \| / __| | |_  | | |/ _ \ |
!>      | | | | | | \__ \ |  _| | | |  __/_|
!>      |_| |_| |_|_|___/ |_|   |_|_|\___(_)
!>
!>      This is a template
!>
!>      **** DO NOT MODIFY THIS .f90 FILE ****
!>      modify the .template instead
!>      The script to regenerate this file can be found at cp2k/tools/generate_xyz_to_vab.py
!> \par History
!>      05.2012 created
!>
!> \author Ruyman Reyes
! *****************************************************************************
! *****************************************************************************
    SUBROUTINE xyz_to_vab_1_0

    !INTEGER :: my_rank,ierror
    !REAL(kind=dp) :: alpha(0:(1+0),0:1,0:0,3)

    INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp
    REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb

    REAL(kind=dp) :: alpha_1(0:(1+0),0:1,0:0)
    REAL(kind=dp) :: alpha_2(0:(1+0),0:1,0:0)
    REAL(kind=dp) :: alpha_3(0:(1+0),0:1,0:0)
    REAL(kind=dp) :: coef_ttz(0:1,0:0)
    REAL(kind=dp) :: coef_tyz(0:1,0:0,0:1,0:0)
    coef_xyz=coef_xyz*prefactor
!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max_local + lb_max_local
!   (current implementation is l**7)
!
!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***
  !RR    DO iaxis=1,3
  alpha_1(:,:,:)=0.0_dp
  alpha_2(:,:,:)=0.0_dp
  alpha_3(:,:,:)=0.0_dp

         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,0,0)=alpha_1(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,1,0)=alpha_1(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,1,0)=alpha_1(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,0,0)=alpha_2(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,1,0)=alpha_2(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,1,0)=alpha_2(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,0,0)=alpha_3(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,1,0)=alpha_3(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,1,0)=alpha_3(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
  !RR    ENDDO
    !
    !   compute v_{lxa,lya,lza,lxb,lyb,lzb} given v_{lxp,lyp,lzp} and alpha(ls,lxa,lxb,1)
    !   use a three step procedure
    !

    lxyz=0
!   DO lzp=0,1+0
       coef_tyz=0.0_dp
       DO lyp=0,1+0-0
          coef_ttz=0.0_dp
          DO lxp=0,1+0-0-lyp
             lxyz=lxyz+1
                coef_ttz(0,0)=coef_ttz(0,0)+coef_xyz(lxyz)*alpha_1(lxp,0,0)
                coef_ttz(1,0)=coef_ttz(1,0)+coef_xyz(lxyz)*alpha_1(lxp,1,0)
          ENDDO

     ! Improve locality
     !temp(:,:) = alpha(lyp,0:1,0:0,2)
     ! RR do lyb = 0, lb_max_local
     ! RR do lya = 0, la_max_local
             !DO lxb=0,0-lyb
             !DO lxa=0,1-lya
                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
                coef_tyz(1,0,0,0)=coef_tyz(1,0,0,0)+coef_ttz(1,0)*alpha_2(lyp,0,0)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,0-lyb
             !DO lxa=0,1-lya
                coef_tyz(0,0,1,0)=coef_tyz(0,0,1,0)+coef_ttz(0,0)*alpha_2(lyp,1,0)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
       ENDDO


          !DO lyb=0,0-lzb
          !DO lya=0,1-lza
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=la_min_local,1
                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(0,0,0)

             ENDDO

          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(0,0,0)


          !ENDDO
          !ENDDO
          !DO lyb=0,0-lzb
          !DO lya=0,1-lza
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(0,1,0)


          !ENDDO
          !ENDDO

!  ENDDO
       coef_tyz=0.0_dp
       DO lyp=0,1+0-1
          coef_ttz=0.0_dp
          DO lxp=0,1+0-1-lyp
             lxyz=lxyz+1
                coef_ttz(0,0)=coef_ttz(0,0)+coef_xyz(lxyz)*alpha_1(lxp,0,0)
                coef_ttz(1,0)=coef_ttz(1,0)+coef_xyz(lxyz)*alpha_1(lxp,1,0)
          ENDDO

     ! Improve locality
     !temp(:,:) = alpha(lyp,0:1,0:0,2)
     ! RR do lyb = 0, lb_max_local
     ! RR do lya = 0, la_max_local
             !DO lxb=0,0-lyb
             !DO lxa=0,1-lya
                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
                coef_tyz(1,0,0,0)=coef_tyz(1,0,0,0)+coef_ttz(1,0)*alpha_2(lyp,0,0)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,0-lyb
             !DO lxa=0,1-lya
                coef_tyz(0,0,1,0)=coef_tyz(0,0,1,0)+coef_ttz(0,0)*alpha_2(lyp,1,0)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
       ENDDO


          !DO lyb=0,0-lzb
          !DO lya=0,1-lza
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=la_min_local,1
                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(1,0,0)

             ENDDO

          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(1,0,0)


          !ENDDO
          !ENDDO
          !DO lyb=0,0-lzb
          !DO lya=0,1-lza
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(1,1,0)


          !ENDDO
          !ENDDO

!  ENDDO

    END SUBROUTINE xyz_to_vab_1_0

! *****************************************************************************
    SUBROUTINE xyz_to_vab_1_1

    !INTEGER :: my_rank,ierror
    !REAL(kind=dp) :: alpha(0:(1+1),0:1,0:1,3)

    INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp
    REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb

    REAL(kind=dp) :: alpha_1(0:(1+1),0:1,0:1)
    REAL(kind=dp) :: alpha_2(0:(1+1),0:1,0:1)
    REAL(kind=dp) :: alpha_3(0:(1+1),0:1,0:1)
    REAL(kind=dp) :: coef_ttz(0:1,0:1)
    REAL(kind=dp) :: coef_tyz(0:1,0:1,0:1,0:1)
    coef_xyz=coef_xyz*prefactor
!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max_local + lb_max_local
!   (current implementation is l**7)
!
!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***
  !RR    DO iaxis=1,3
  alpha_1(:,:,:)=0.0_dp
  alpha_2(:,:,:)=0.0_dp
  alpha_3(:,:,:)=0.0_dp

         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,0,0)=alpha_1(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,0,1)=alpha_1(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,0,1)=alpha_1(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,1,0)=alpha_1(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,1,0)=alpha_1(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,1,1)=alpha_1(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,1,1)=alpha_1(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,1,1)=alpha_1(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,1,1)=alpha_1(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,0,0)=alpha_2(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,0,1)=alpha_2(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,0,1)=alpha_2(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,1,0)=alpha_2(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,1,0)=alpha_2(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,1,1)=alpha_2(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,1,1)=alpha_2(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,1,1)=alpha_2(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,1,1)=alpha_2(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,0,0)=alpha_3(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,0,1)=alpha_3(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,0,1)=alpha_3(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,1,0)=alpha_3(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,1,0)=alpha_3(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,1,1)=alpha_3(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,1,1)=alpha_3(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,1,1)=alpha_3(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,1,1)=alpha_3(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
  !RR    ENDDO
    !
    !   compute v_{lxa,lya,lza,lxb,lyb,lzb} given v_{lxp,lyp,lzp} and alpha(ls,lxa,lxb,1)
    !   use a three step procedure
    !

    lxyz=0
    DO lzp=0,1+1
       coef_tyz=0.0_dp
       DO lyp=0,1+1-lzp
          coef_ttz=0.0_dp
          DO lxp=0,1+1-lzp-lyp
             lxyz=lxyz+1
                coef_ttz(0,0)=coef_ttz(0,0)+coef_xyz(lxyz)*alpha_1(lxp,0,0)
                coef_ttz(1,0)=coef_ttz(1,0)+coef_xyz(lxyz)*alpha_1(lxp,1,0)
                coef_ttz(0,1)=coef_ttz(0,1)+coef_xyz(lxyz)*alpha_1(lxp,0,1)
                coef_ttz(1,1)=coef_ttz(1,1)+coef_xyz(lxyz)*alpha_1(lxp,1,1)
          ENDDO

     ! Improve locality
     !temp(:,:) = alpha(lyp,0:1,0:1,2)
     ! RR do lyb = 0, lb_max_local
     ! RR do lya = 0, la_max_local
             !DO lxb=0,1-lyb
             !DO lxa=0,1-lya
                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
                coef_tyz(1,0,0,0)=coef_tyz(1,0,0,0)+coef_ttz(1,0)*alpha_2(lyp,0,0)
             !ENDDO
             !DO lxa=0,1-lya
                coef_tyz(0,1,0,0)=coef_tyz(0,1,0,0)+coef_ttz(0,1)*alpha_2(lyp,0,0)
                coef_tyz(1,1,0,0)=coef_tyz(1,1,0,0)+coef_ttz(1,1)*alpha_2(lyp,0,0)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,1-lyb
             !DO lxa=0,1-lya
                coef_tyz(0,0,1,0)=coef_tyz(0,0,1,0)+coef_ttz(0,0)*alpha_2(lyp,1,0)
             !ENDDO
             !DO lxa=0,1-lya
                coef_tyz(0,1,1,0)=coef_tyz(0,1,1,0)+coef_ttz(0,1)*alpha_2(lyp,1,0)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
     ! RR do lya = 0, la_max_local
             !DO lxb=0,1-lyb
             !DO lxa=0,1-lya
                coef_tyz(0,0,0,1)=coef_tyz(0,0,0,1)+coef_ttz(0,0)*alpha_2(lyp,0,1)
                coef_tyz(1,0,0,1)=coef_tyz(1,0,0,1)+coef_ttz(1,0)*alpha_2(lyp,0,1)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,1-lyb
             !DO lxa=0,1-lya
                coef_tyz(0,0,1,1)=coef_tyz(0,0,1,1)+coef_ttz(0,0)*alpha_2(lyp,1,1)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
       ENDDO


          !DO lyb=0,1-lzb
          !DO lya=0,1-lza
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO

          !ENDDO
          !ENDDO
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,1-lza
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,0,0)



          !ENDDO
          !ENDDO
          !DO lyb=0,1-lzb
          !DO lya=0,1-lza
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,1-lza
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,1,0)



          !ENDDO
          !ENDDO
          !DO lyb=0,1-lzb
          !DO lya=0,1-lza
             lxb = 0
             !jco=coset(0,0,1)
             jco=4

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,1)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,1)
             jco=4

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,1)



          !ENDDO
          !ENDDO
          !DO lyb=0,1-lzb
          !DO lya=0,1-lza
             lxb = 0
             !jco=coset(0,0,1)
             jco=4

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,1)



          !ENDDO
          !ENDDO

   ENDDO

    END SUBROUTINE xyz_to_vab_1_1

! *****************************************************************************
    SUBROUTINE xyz_to_vab_0_0

    !INTEGER :: my_rank,ierror
    !REAL(kind=dp) :: alpha(0:(0+0),0:0,0:0,3)

    INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp
    REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb

    REAL(kind=dp) :: alpha_1(0:(0+0),0:0,0:0)
    REAL(kind=dp) :: alpha_2(0:(0+0),0:0,0:0)
    REAL(kind=dp) :: alpha_3(0:(0+0),0:0,0:0)
!    REAL(kind=dp) :: coef_ttz, coef_tyz
    coef_xyz=coef_xyz*prefactor
!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max_local + lb_max_local
!   (current implementation is l**7)
!
!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***
  !RR    DO iaxis=1,3
    !alpha_1(0,0,0)=1.0_dp
    !alpha_2(0,0,0)=1.0_dp
    !alpha_3(0,0,0)=1.0_dp
    !
    !   compute v_{lxa,lya,lza,lxb,lyb,lzb} given v_{lxp,lyp,lzp} and alpha(ls,lxa,lxb,1)
    !   use a three step procedure
    !

       !coef_tyz=0.0_dp
          !coef_ttz=0.0_dp
             !   coef_ttz=0.0_dp+coef_xyz(1)*1.0_dp
             !   coef_tyz=0.0_dp+coef_ttz*1.0_dp
             !DO lxb=MAX(lb_min_local,0),0
             jco=1 ! coset(0,0,0)
             !DO lxa=MAX(la_min_local,0),0
                ico=1 ! coset(0,0,0)
                vab(ico,jco)=vab(ico,jco)+coef_xyz(1)
             !ENDDO
             !ENDDO

    END SUBROUTINE xyz_to_vab_0_0

! *****************************************************************************
    SUBROUTINE xyz_to_vab_0_1

    !INTEGER :: my_rank,ierror
    !REAL(kind=dp) :: alpha(0:(0+1),0:0,0:1,3)

    INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp
    REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb

    REAL(kind=dp) :: alpha_1(0:(0+1),0:0,0:1)
    REAL(kind=dp) :: alpha_2(0:(0+1),0:0,0:1)
    REAL(kind=dp) :: alpha_3(0:(0+1),0:0,0:1)
    REAL(kind=dp) :: coef_ttz(0:0,0:1)
    REAL(kind=dp) :: coef_tyz(0:0,0:1,0:0,0:1)
    coef_xyz=coef_xyz*prefactor
!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max_local + lb_max_local
!   (current implementation is l**7)
!
!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***
  !RR    DO iaxis=1,3
  alpha_1(:,:,:)=0.0_dp
  alpha_2(:,:,:)=0.0_dp
  alpha_3(:,:,:)=0.0_dp

         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,0,0)=alpha_1(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,0,1)=alpha_1(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,0,1)=alpha_1(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,0,0)=alpha_2(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,0,1)=alpha_2(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,0,1)=alpha_2(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,0,0)=alpha_3(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,0,1)=alpha_3(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,0,1)=alpha_3(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
  !RR    ENDDO
    !
    !   compute v_{lxa,lya,lza,lxb,lyb,lzb} given v_{lxp,lyp,lzp} and alpha(ls,lxa,lxb,1)
    !   use a three step procedure
    !

    lxyz=0
!   DO lzp=0,0+1
       coef_tyz=0.0_dp
       DO lyp=0,0+1-0
          coef_ttz=0.0_dp
          DO lxp=0,0+1-0-lyp
             lxyz=lxyz+1
                coef_ttz(0,0)=coef_ttz(0,0)+coef_xyz(lxyz)*alpha_1(lxp,0,0)
                coef_ttz(0,1)=coef_ttz(0,1)+coef_xyz(lxyz)*alpha_1(lxp,0,1)
          ENDDO

     ! Improve locality
     !temp(:,:) = alpha(lyp,0:0,0:1,2)
     ! RR do lyb = 0, lb_max_local
     ! RR do lya = 0, la_max_local
             !DO lxb=0,1-lyb
             !DO lxa=0,0-lya
                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
             !ENDDO
             !DO lxa=0,0-lya
                coef_tyz(0,1,0,0)=coef_tyz(0,1,0,0)+coef_ttz(0,1)*alpha_2(lyp,0,0)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
     ! RR do lya = 0, la_max_local
             !DO lxb=0,1-lyb
             !DO lxa=0,0-lya
                coef_tyz(0,0,0,1)=coef_tyz(0,0,0,1)+coef_ttz(0,0)*alpha_2(lyp,0,1)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
       ENDDO


          !DO lyb=0,1-lzb
          !DO lya=0,0-lza
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(0,0,0)


             ENDDO
          !ENDDO
          !ENDDO
          !DO lya=0,0-lza
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(0,0,0)


          !ENDDO
          !ENDDO
          !DO lyb=0,1-lzb
          !DO lya=0,0-lza
             lxb = 0
             !jco=coset(0,0,1)
             jco=4

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(0,0,1)


          !ENDDO
          !ENDDO

!  ENDDO
       coef_tyz=0.0_dp
       DO lyp=0,0+1-1
          coef_ttz=0.0_dp
          DO lxp=0,0+1-1-lyp
             lxyz=lxyz+1
                coef_ttz(0,0)=coef_ttz(0,0)+coef_xyz(lxyz)*alpha_1(lxp,0,0)
                coef_ttz(0,1)=coef_ttz(0,1)+coef_xyz(lxyz)*alpha_1(lxp,0,1)
          ENDDO

     ! Improve locality
     !temp(:,:) = alpha(lyp,0:0,0:1,2)
     ! RR do lyb = 0, lb_max_local
     ! RR do lya = 0, la_max_local
             !DO lxb=0,1-lyb
             !DO lxa=0,0-lya
                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
             !ENDDO
             !DO lxa=0,0-lya
                coef_tyz(0,1,0,0)=coef_tyz(0,1,0,0)+coef_ttz(0,1)*alpha_2(lyp,0,0)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
     ! RR do lya = 0, la_max_local
             !DO lxb=0,1-lyb
             !DO lxa=0,0-lya
                coef_tyz(0,0,0,1)=coef_tyz(0,0,0,1)+coef_ttz(0,0)*alpha_2(lyp,0,1)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
       ENDDO


          !DO lyb=0,1-lzb
          !DO lya=0,0-lza
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(1,0,0)


             ENDDO
          !ENDDO
          !ENDDO
          !DO lya=0,0-lza
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(1,0,0)


          !ENDDO
          !ENDDO
          !DO lyb=0,1-lzb
          !DO lya=0,0-lza
             lxb = 0
             !jco=coset(0,0,1)
             jco=4

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(1,0,1)


          !ENDDO
          !ENDDO

!  ENDDO

    END SUBROUTINE xyz_to_vab_0_1

! *****************************************************************************
    SUBROUTINE xyz_to_vab_0_2

    !INTEGER :: my_rank,ierror
    !REAL(kind=dp) :: alpha(0:(0+2),0:0,0:2,3)

    INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp
    REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb

    REAL(kind=dp) :: alpha_1(0:(0+2),0:0,0:2)
    REAL(kind=dp) :: alpha_2(0:(0+2),0:0,0:2)
    REAL(kind=dp) :: alpha_3(0:(0+2),0:0,0:2)
    REAL(kind=dp) :: coef_ttz(0:0,0:2)
    REAL(kind=dp) :: coef_tyz(0:0,0:2,0:0,0:2)
    coef_xyz=coef_xyz*prefactor
!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max_local + lb_max_local
!   (current implementation is l**7)
!
!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***
  !RR    DO iaxis=1,3
  alpha_1(:,:,:)=0.0_dp
  alpha_2(:,:,:)=0.0_dp
  alpha_3(:,:,:)=0.0_dp

         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,0,0)=alpha_1(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,0,1)=alpha_1(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,0,1)=alpha_1(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,0,2)=alpha_1(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,0,2)=alpha_1(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,0,2)=alpha_1(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,0,0)=alpha_2(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,0,1)=alpha_2(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,0,1)=alpha_2(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,0,2)=alpha_2(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,0,2)=alpha_2(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,0,2)=alpha_2(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,0,0)=alpha_3(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,0,1)=alpha_3(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,0,1)=alpha_3(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,0,2)=alpha_3(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,0,2)=alpha_3(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,0,2)=alpha_3(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
  !RR    ENDDO
    !
    !   compute v_{lxa,lya,lza,lxb,lyb,lzb} given v_{lxp,lyp,lzp} and alpha(ls,lxa,lxb,1)
    !   use a three step procedure
    !

    lxyz=0
    DO lzp=0,0+2
       coef_tyz=0.0_dp
       DO lyp=0,0+2-lzp
          coef_ttz=0.0_dp
          DO lxp=0,0+2-lzp-lyp
             lxyz=lxyz+1
                coef_ttz(0,0)=coef_ttz(0,0)+coef_xyz(lxyz)*alpha_1(lxp,0,0)
                coef_ttz(0,1)=coef_ttz(0,1)+coef_xyz(lxyz)*alpha_1(lxp,0,1)
                coef_ttz(0,2)=coef_ttz(0,2)+coef_xyz(lxyz)*alpha_1(lxp,0,2)
          ENDDO

     ! Improve locality
     !temp(:,:) = alpha(lyp,0:0,0:2,2)
     ! RR do lyb = 0, lb_max_local
     ! RR do lya = 0, la_max_local
             !DO lxb=0,2-lyb
             !DO lxa=0,0-lya
                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
             !ENDDO
             !DO lxa=0,0-lya
                coef_tyz(0,1,0,0)=coef_tyz(0,1,0,0)+coef_ttz(0,1)*alpha_2(lyp,0,0)
             !ENDDO
             !DO lxa=0,0-lya
                coef_tyz(0,2,0,0)=coef_tyz(0,2,0,0)+coef_ttz(0,2)*alpha_2(lyp,0,0)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
     ! RR do lya = 0, la_max_local
             !DO lxb=0,2-lyb
             !DO lxa=0,0-lya
                coef_tyz(0,0,0,1)=coef_tyz(0,0,0,1)+coef_ttz(0,0)*alpha_2(lyp,0,1)
             !ENDDO
             !DO lxa=0,0-lya
                coef_tyz(0,1,0,1)=coef_tyz(0,1,0,1)+coef_ttz(0,1)*alpha_2(lyp,0,1)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
     ! RR do lya = 0, la_max_local
             !DO lxb=0,2-lyb
             !DO lxa=0,0-lya
                coef_tyz(0,0,0,2)=coef_tyz(0,0,0,2)+coef_ttz(0,0)*alpha_2(lyp,0,2)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
       ENDDO


          !DO lyb=0,2-lzb
          !DO lya=0,0-lza
              DO lxb=lb_min_local,2
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,0-lza
              DO lxb=MAX(lb_min_local-0-1,0),1
             jco=coset(lxb,1,0)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,0-lza
             lxb = 0
             !jco=coset(0,2,0)
             jco=8

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,0,0)



          !ENDDO
          !ENDDO
          !DO lyb=0,2-lzb
          !DO lya=0,0-lza
              DO lxb=MAX(lb_min_local-1-0,0),1
             jco=coset(lxb,0,1)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,1)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,0-lza
             lxb = 0
             !jco=coset(0,1,1)
             jco=9

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,1)



          !ENDDO
          !ENDDO
          !DO lyb=0,2-lzb
          !DO lya=0,0-lza
             lxb = 0
             !jco=coset(0,0,2)
             jco=10

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,2)



          !ENDDO
          !ENDDO

   ENDDO

    END SUBROUTINE xyz_to_vab_0_2

! *****************************************************************************
    SUBROUTINE xyz_to_vab_0_3

    !INTEGER :: my_rank,ierror
    !REAL(kind=dp) :: alpha(0:(0+3),0:0,0:3,3)

    INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp
    REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb

    REAL(kind=dp) :: alpha_1(0:(0+3),0:0,0:3)
    REAL(kind=dp) :: alpha_2(0:(0+3),0:0,0:3)
    REAL(kind=dp) :: alpha_3(0:(0+3),0:0,0:3)
    REAL(kind=dp) :: coef_ttz(0:0,0:3)
    REAL(kind=dp) :: coef_tyz(0:0,0:3,0:0,0:3)
    coef_xyz=coef_xyz*prefactor
!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max_local + lb_max_local
!   (current implementation is l**7)
!
!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***
  !RR    DO iaxis=1,3
  alpha_1(:,:,:)=0.0_dp
  alpha_2(:,:,:)=0.0_dp
  alpha_3(:,:,:)=0.0_dp

         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,0,0)=alpha_1(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,0,1)=alpha_1(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,0,1)=alpha_1(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,0,2)=alpha_1(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,0,2)=alpha_1(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,0,2)=alpha_1(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(3,0,3)=alpha_1(3,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(2,0,3)=alpha_1(2,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,0,3)=alpha_1(1,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,0,3)=alpha_1(0,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,0,0)=alpha_2(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,0,1)=alpha_2(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,0,1)=alpha_2(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,0,2)=alpha_2(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,0,2)=alpha_2(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,0,2)=alpha_2(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(3,0,3)=alpha_2(3,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(2,0,3)=alpha_2(2,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,0,3)=alpha_2(1,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,0,3)=alpha_2(0,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,0,0)=alpha_3(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,0,1)=alpha_3(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,0,1)=alpha_3(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,0,2)=alpha_3(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,0,2)=alpha_3(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,0,2)=alpha_3(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(3,0,3)=alpha_3(3,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(2,0,3)=alpha_3(2,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,0,3)=alpha_3(1,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,0,3)=alpha_3(0,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
  !RR    ENDDO
    !
    !   compute v_{lxa,lya,lza,lxb,lyb,lzb} given v_{lxp,lyp,lzp} and alpha(ls,lxa,lxb,1)
    !   use a three step procedure
    !

    lxyz=0
    DO lzp=0,0+3
       coef_tyz=0.0_dp
       DO lyp=0,0+3-lzp
          coef_ttz=0.0_dp
          DO lxp=0,0+3-lzp-lyp
             lxyz=lxyz+1
                coef_ttz(0,0)=coef_ttz(0,0)+coef_xyz(lxyz)*alpha_1(lxp,0,0)
                coef_ttz(0,1)=coef_ttz(0,1)+coef_xyz(lxyz)*alpha_1(lxp,0,1)
                coef_ttz(0,2)=coef_ttz(0,2)+coef_xyz(lxyz)*alpha_1(lxp,0,2)
                coef_ttz(0,3)=coef_ttz(0,3)+coef_xyz(lxyz)*alpha_1(lxp,0,3)
          ENDDO

     ! Improve locality
     !temp(:,:) = alpha(lyp,0:0,0:3,2)
     ! RR do lyb = 0, lb_max_local
     ! RR do lya = 0, la_max_local
             !DO lxb=0,3-lyb
             !DO lxa=0,0-lya
                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
             !ENDDO
             !DO lxa=0,0-lya
                coef_tyz(0,1,0,0)=coef_tyz(0,1,0,0)+coef_ttz(0,1)*alpha_2(lyp,0,0)
             !ENDDO
             !DO lxa=0,0-lya
                coef_tyz(0,2,0,0)=coef_tyz(0,2,0,0)+coef_ttz(0,2)*alpha_2(lyp,0,0)
             !ENDDO
             !DO lxa=0,0-lya
                coef_tyz(0,3,0,0)=coef_tyz(0,3,0,0)+coef_ttz(0,3)*alpha_2(lyp,0,0)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
     ! RR do lya = 0, la_max_local
             !DO lxb=0,3-lyb
             !DO lxa=0,0-lya
                coef_tyz(0,0,0,1)=coef_tyz(0,0,0,1)+coef_ttz(0,0)*alpha_2(lyp,0,1)
             !ENDDO
             !DO lxa=0,0-lya
                coef_tyz(0,1,0,1)=coef_tyz(0,1,0,1)+coef_ttz(0,1)*alpha_2(lyp,0,1)
             !ENDDO
             !DO lxa=0,0-lya
                coef_tyz(0,2,0,1)=coef_tyz(0,2,0,1)+coef_ttz(0,2)*alpha_2(lyp,0,1)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
     ! RR do lya = 0, la_max_local
             !DO lxb=0,3-lyb
             !DO lxa=0,0-lya
                coef_tyz(0,0,0,2)=coef_tyz(0,0,0,2)+coef_ttz(0,0)*alpha_2(lyp,0,2)
             !ENDDO
             !DO lxa=0,0-lya
                coef_tyz(0,1,0,2)=coef_tyz(0,1,0,2)+coef_ttz(0,1)*alpha_2(lyp,0,2)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
     ! RR do lya = 0, la_max_local
             !DO lxb=0,3-lyb
             !DO lxa=0,0-lya
                coef_tyz(0,0,0,3)=coef_tyz(0,0,0,3)+coef_ttz(0,0)*alpha_2(lyp,0,3)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
       ENDDO


          !DO lyb=0,3-lzb
          !DO lya=0,0-lza
              DO lxb=lb_min_local,3
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,0-lza
              DO lxb=MAX(lb_min_local-0-1,0),2
             jco=coset(lxb,1,0)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,0-lza
              DO lxb=MAX(lb_min_local-0-2,0),1
             jco=coset(lxb,2,0)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,0,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,0-lza
             lxb = 0
             !jco=coset(0,3,0)
             jco=17

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,3)*alpha_3(lzp,0,0)



          !ENDDO
          !ENDDO
          !DO lyb=0,3-lzb
          !DO lya=0,0-lza
              DO lxb=MAX(lb_min_local-1-0,0),2
             jco=coset(lxb,0,1)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,1)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,0-lza
              DO lxb=MAX(lb_min_local-1-1,0),1
             jco=coset(lxb,1,1)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,1)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,0-lza
             lxb = 0
             !jco=coset(0,2,1)
             jco=18

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,0,1)



          !ENDDO
          !ENDDO
          !DO lyb=0,3-lzb
          !DO lya=0,0-lza
              DO lxb=MAX(lb_min_local-2-0,0),1
             jco=coset(lxb,0,2)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,2)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,0-lza
             lxb = 0
             !jco=coset(0,1,2)
             jco=19

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,2)



          !ENDDO
          !ENDDO
          !DO lyb=0,3-lzb
          !DO lya=0,0-lza
             lxb = 0
             !jco=coset(0,0,3)
             jco=20

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,3)



          !ENDDO
          !ENDDO

   ENDDO

    END SUBROUTINE xyz_to_vab_0_3

! *****************************************************************************
    SUBROUTINE xyz_to_vab_0_4

    !INTEGER :: my_rank,ierror
    !REAL(kind=dp) :: alpha(0:(0+4),0:0,0:4,3)

    INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp
    REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb

    REAL(kind=dp) :: alpha_1(0:(0+4),0:0,0:4)
    REAL(kind=dp) :: alpha_2(0:(0+4),0:0,0:4)
    REAL(kind=dp) :: alpha_3(0:(0+4),0:0,0:4)
    REAL(kind=dp) :: coef_ttz(0:0,0:4)
    REAL(kind=dp) :: coef_tyz(0:0,0:4,0:0,0:4)
    coef_xyz=coef_xyz*prefactor
!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max_local + lb_max_local
!   (current implementation is l**7)
!
!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***
  !RR    DO iaxis=1,3
  alpha_1(:,:,:)=0.0_dp
  alpha_2(:,:,:)=0.0_dp
  alpha_3(:,:,:)=0.0_dp

         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,0,0)=alpha_1(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,0,1)=alpha_1(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,0,1)=alpha_1(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,0,2)=alpha_1(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,0,2)=alpha_1(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,0,2)=alpha_1(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(3,0,3)=alpha_1(3,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(2,0,3)=alpha_1(2,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,0,3)=alpha_1(1,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,0,3)=alpha_1(0,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(4,0,4)=alpha_1(4,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(4,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(3,0,4)=alpha_1(3,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)/REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(2,0,4)=alpha_1(2,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)/REAL(3,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,0,4)=alpha_1(1,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(4,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,0,4)=alpha_1(0,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,0,0)=alpha_2(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,0,1)=alpha_2(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,0,1)=alpha_2(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,0,2)=alpha_2(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,0,2)=alpha_2(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,0,2)=alpha_2(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(3,0,3)=alpha_2(3,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(2,0,3)=alpha_2(2,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,0,3)=alpha_2(1,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,0,3)=alpha_2(0,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(4,0,4)=alpha_2(4,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(4,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(3,0,4)=alpha_2(3,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)/REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(2,0,4)=alpha_2(2,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)/REAL(3,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,0,4)=alpha_2(1,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(4,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,0,4)=alpha_2(0,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,0,0)=alpha_3(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,0,1)=alpha_3(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,0,1)=alpha_3(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,0,2)=alpha_3(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,0,2)=alpha_3(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,0,2)=alpha_3(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(3,0,3)=alpha_3(3,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(2,0,3)=alpha_3(2,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,0,3)=alpha_3(1,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,0,3)=alpha_3(0,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(4,0,4)=alpha_3(4,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(4,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(3,0,4)=alpha_3(3,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)/REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(2,0,4)=alpha_3(2,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)/REAL(3,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,0,4)=alpha_3(1,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(4,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,0,4)=alpha_3(0,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
  !RR    ENDDO
    !
    !   compute v_{lxa,lya,lza,lxb,lyb,lzb} given v_{lxp,lyp,lzp} and alpha(ls,lxa,lxb,1)
    !   use a three step procedure
    !

    lxyz=0
    DO lzp=0,0+4
       coef_tyz=0.0_dp
       DO lyp=0,0+4-lzp
          coef_ttz=0.0_dp
          DO lxp=0,0+4-lzp-lyp
             lxyz=lxyz+1
                coef_ttz(0,0)=coef_ttz(0,0)+coef_xyz(lxyz)*alpha_1(lxp,0,0)
                coef_ttz(0,1)=coef_ttz(0,1)+coef_xyz(lxyz)*alpha_1(lxp,0,1)
                coef_ttz(0,2)=coef_ttz(0,2)+coef_xyz(lxyz)*alpha_1(lxp,0,2)
                coef_ttz(0,3)=coef_ttz(0,3)+coef_xyz(lxyz)*alpha_1(lxp,0,3)
                coef_ttz(0,4)=coef_ttz(0,4)+coef_xyz(lxyz)*alpha_1(lxp,0,4)
          ENDDO

     ! Improve locality
     !temp(:,:) = alpha(lyp,0:0,0:4,2)
     ! RR do lyb = 0, lb_max_local
     ! RR do lya = 0, la_max_local
             !DO lxb=0,4-lyb
             !DO lxa=0,0-lya
                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
             !ENDDO
             !DO lxa=0,0-lya
                coef_tyz(0,1,0,0)=coef_tyz(0,1,0,0)+coef_ttz(0,1)*alpha_2(lyp,0,0)
             !ENDDO
             !DO lxa=0,0-lya
                coef_tyz(0,2,0,0)=coef_tyz(0,2,0,0)+coef_ttz(0,2)*alpha_2(lyp,0,0)
             !ENDDO
             !DO lxa=0,0-lya
                coef_tyz(0,3,0,0)=coef_tyz(0,3,0,0)+coef_ttz(0,3)*alpha_2(lyp,0,0)
             !ENDDO
             !DO lxa=0,0-lya
                coef_tyz(0,4,0,0)=coef_tyz(0,4,0,0)+coef_ttz(0,4)*alpha_2(lyp,0,0)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
     ! RR do lya = 0, la_max_local
             !DO lxb=0,4-lyb
             !DO lxa=0,0-lya
                coef_tyz(0,0,0,1)=coef_tyz(0,0,0,1)+coef_ttz(0,0)*alpha_2(lyp,0,1)
             !ENDDO
             !DO lxa=0,0-lya
                coef_tyz(0,1,0,1)=coef_tyz(0,1,0,1)+coef_ttz(0,1)*alpha_2(lyp,0,1)
             !ENDDO
             !DO lxa=0,0-lya
                coef_tyz(0,2,0,1)=coef_tyz(0,2,0,1)+coef_ttz(0,2)*alpha_2(lyp,0,1)
             !ENDDO
             !DO lxa=0,0-lya
                coef_tyz(0,3,0,1)=coef_tyz(0,3,0,1)+coef_ttz(0,3)*alpha_2(lyp,0,1)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
     ! RR do lya = 0, la_max_local
             !DO lxb=0,4-lyb
             !DO lxa=0,0-lya
                coef_tyz(0,0,0,2)=coef_tyz(0,0,0,2)+coef_ttz(0,0)*alpha_2(lyp,0,2)
             !ENDDO
             !DO lxa=0,0-lya
                coef_tyz(0,1,0,2)=coef_tyz(0,1,0,2)+coef_ttz(0,1)*alpha_2(lyp,0,2)
             !ENDDO
             !DO lxa=0,0-lya
                coef_tyz(0,2,0,2)=coef_tyz(0,2,0,2)+coef_ttz(0,2)*alpha_2(lyp,0,2)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
     ! RR do lya = 0, la_max_local
             !DO lxb=0,4-lyb
             !DO lxa=0,0-lya
                coef_tyz(0,0,0,3)=coef_tyz(0,0,0,3)+coef_ttz(0,0)*alpha_2(lyp,0,3)
             !ENDDO
             !DO lxa=0,0-lya
                coef_tyz(0,1,0,3)=coef_tyz(0,1,0,3)+coef_ttz(0,1)*alpha_2(lyp,0,3)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
     ! RR do lya = 0, la_max_local
             !DO lxb=0,4-lyb
             !DO lxa=0,0-lya
                coef_tyz(0,0,0,4)=coef_tyz(0,0,0,4)+coef_ttz(0,0)*alpha_2(lyp,0,4)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
       ENDDO


          !DO lyb=0,4-lzb
          !DO lya=0,0-lza
              DO lxb=lb_min_local,4
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,0-lza
              DO lxb=MAX(lb_min_local-0-1,0),3
             jco=coset(lxb,1,0)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,0-lza
              DO lxb=MAX(lb_min_local-0-2,0),2
             jco=coset(lxb,2,0)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,0,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,0-lza
              DO lxb=MAX(lb_min_local-0-3,0),1
             jco=coset(lxb,3,0)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,3)*alpha_3(lzp,0,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,0-lza
             lxb = 0
             !jco=coset(0,4,0)
             jco=31

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,4)*alpha_3(lzp,0,0)



          !ENDDO
          !ENDDO
          !DO lyb=0,4-lzb
          !DO lya=0,0-lza
              DO lxb=MAX(lb_min_local-1-0,0),3
             jco=coset(lxb,0,1)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,1)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,0-lza
              DO lxb=MAX(lb_min_local-1-1,0),2
             jco=coset(lxb,1,1)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,1)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,0-lza
              DO lxb=MAX(lb_min_local-1-2,0),1
             jco=coset(lxb,2,1)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,0,1)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,0-lza
             lxb = 0
             !jco=coset(0,3,1)
             jco=32

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,3)*alpha_3(lzp,0,1)



          !ENDDO
          !ENDDO
          !DO lyb=0,4-lzb
          !DO lya=0,0-lza
              DO lxb=MAX(lb_min_local-2-0,0),2
             jco=coset(lxb,0,2)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,2)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,0-lza
              DO lxb=MAX(lb_min_local-2-1,0),1
             jco=coset(lxb,1,2)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,2)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,0-lza
             lxb = 0
             !jco=coset(0,2,2)
             jco=33

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,0,2)



          !ENDDO
          !ENDDO
          !DO lyb=0,4-lzb
          !DO lya=0,0-lza
              DO lxb=MAX(lb_min_local-3-0,0),1
             jco=coset(lxb,0,3)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,3)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,0-lza
             lxb = 0
             !jco=coset(0,1,3)
             jco=34

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,3)



          !ENDDO
          !ENDDO
          !DO lyb=0,4-lzb
          !DO lya=0,0-lza
             lxb = 0
             !jco=coset(0,0,4)
             jco=35

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,4)



          !ENDDO
          !ENDDO

   ENDDO

    END SUBROUTINE xyz_to_vab_0_4

! *****************************************************************************
    SUBROUTINE xyz_to_vab_1_2

    !INTEGER :: my_rank,ierror
    !REAL(kind=dp) :: alpha(0:(1+2),0:1,0:2,3)

    INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp
    REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb

    REAL(kind=dp) :: alpha_1(0:(1+2),0:1,0:2)
    REAL(kind=dp) :: alpha_2(0:(1+2),0:1,0:2)
    REAL(kind=dp) :: alpha_3(0:(1+2),0:1,0:2)
    REAL(kind=dp) :: coef_ttz(0:1,0:2)
    REAL(kind=dp) :: coef_tyz(0:1,0:2,0:1,0:2)
    coef_xyz=coef_xyz*prefactor
!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max_local + lb_max_local
!   (current implementation is l**7)
!
!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***
  !RR    DO iaxis=1,3
  alpha_1(:,:,:)=0.0_dp
  alpha_2(:,:,:)=0.0_dp
  alpha_3(:,:,:)=0.0_dp

         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,0,0)=alpha_1(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,0,1)=alpha_1(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,0,1)=alpha_1(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,0,2)=alpha_1(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,0,2)=alpha_1(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,0,2)=alpha_1(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,1,0)=alpha_1(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,1,0)=alpha_1(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,1,1)=alpha_1(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,1,1)=alpha_1(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,1,1)=alpha_1(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,1,1)=alpha_1(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(3,1,2)=alpha_1(3,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(2,1,2)=alpha_1(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,1,2)=alpha_1(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,1,2)=alpha_1(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,1,2)=alpha_1(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,1,2)=alpha_1(0,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,0,0)=alpha_2(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,0,1)=alpha_2(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,0,1)=alpha_2(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,0,2)=alpha_2(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,0,2)=alpha_2(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,0,2)=alpha_2(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,1,0)=alpha_2(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,1,0)=alpha_2(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,1,1)=alpha_2(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,1,1)=alpha_2(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,1,1)=alpha_2(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,1,1)=alpha_2(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(3,1,2)=alpha_2(3,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(2,1,2)=alpha_2(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,1,2)=alpha_2(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,1,2)=alpha_2(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,1,2)=alpha_2(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,1,2)=alpha_2(0,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,0,0)=alpha_3(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,0,1)=alpha_3(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,0,1)=alpha_3(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,0,2)=alpha_3(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,0,2)=alpha_3(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,0,2)=alpha_3(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,1,0)=alpha_3(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,1,0)=alpha_3(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,1,1)=alpha_3(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,1,1)=alpha_3(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,1,1)=alpha_3(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,1,1)=alpha_3(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(3,1,2)=alpha_3(3,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(2,1,2)=alpha_3(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,1,2)=alpha_3(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,1,2)=alpha_3(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,1,2)=alpha_3(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,1,2)=alpha_3(0,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
  !RR    ENDDO
    !
    !   compute v_{lxa,lya,lza,lxb,lyb,lzb} given v_{lxp,lyp,lzp} and alpha(ls,lxa,lxb,1)
    !   use a three step procedure
    !

    lxyz=0
    DO lzp=0,1+2
       coef_tyz=0.0_dp
       DO lyp=0,1+2-lzp
          coef_ttz=0.0_dp
          DO lxp=0,1+2-lzp-lyp
             lxyz=lxyz+1
                coef_ttz(0,0)=coef_ttz(0,0)+coef_xyz(lxyz)*alpha_1(lxp,0,0)
                coef_ttz(1,0)=coef_ttz(1,0)+coef_xyz(lxyz)*alpha_1(lxp,1,0)
                coef_ttz(0,1)=coef_ttz(0,1)+coef_xyz(lxyz)*alpha_1(lxp,0,1)
                coef_ttz(1,1)=coef_ttz(1,1)+coef_xyz(lxyz)*alpha_1(lxp,1,1)
                coef_ttz(0,2)=coef_ttz(0,2)+coef_xyz(lxyz)*alpha_1(lxp,0,2)
                coef_ttz(1,2)=coef_ttz(1,2)+coef_xyz(lxyz)*alpha_1(lxp,1,2)
          ENDDO

     ! Improve locality
     !temp(:,:) = alpha(lyp,0:1,0:2,2)
     ! RR do lyb = 0, lb_max_local
     ! RR do lya = 0, la_max_local
             !DO lxb=0,2-lyb
             !DO lxa=0,1-lya
                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
                coef_tyz(1,0,0,0)=coef_tyz(1,0,0,0)+coef_ttz(1,0)*alpha_2(lyp,0,0)
             !ENDDO
             !DO lxa=0,1-lya
                coef_tyz(0,1,0,0)=coef_tyz(0,1,0,0)+coef_ttz(0,1)*alpha_2(lyp,0,0)
                coef_tyz(1,1,0,0)=coef_tyz(1,1,0,0)+coef_ttz(1,1)*alpha_2(lyp,0,0)
             !ENDDO
             !DO lxa=0,1-lya
                coef_tyz(0,2,0,0)=coef_tyz(0,2,0,0)+coef_ttz(0,2)*alpha_2(lyp,0,0)
                coef_tyz(1,2,0,0)=coef_tyz(1,2,0,0)+coef_ttz(1,2)*alpha_2(lyp,0,0)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,2-lyb
             !DO lxa=0,1-lya
                coef_tyz(0,0,1,0)=coef_tyz(0,0,1,0)+coef_ttz(0,0)*alpha_2(lyp,1,0)
             !ENDDO
             !DO lxa=0,1-lya
                coef_tyz(0,1,1,0)=coef_tyz(0,1,1,0)+coef_ttz(0,1)*alpha_2(lyp,1,0)
             !ENDDO
             !DO lxa=0,1-lya
                coef_tyz(0,2,1,0)=coef_tyz(0,2,1,0)+coef_ttz(0,2)*alpha_2(lyp,1,0)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
     ! RR do lya = 0, la_max_local
             !DO lxb=0,2-lyb
             !DO lxa=0,1-lya
                coef_tyz(0,0,0,1)=coef_tyz(0,0,0,1)+coef_ttz(0,0)*alpha_2(lyp,0,1)
                coef_tyz(1,0,0,1)=coef_tyz(1,0,0,1)+coef_ttz(1,0)*alpha_2(lyp,0,1)
             !ENDDO
             !DO lxa=0,1-lya
                coef_tyz(0,1,0,1)=coef_tyz(0,1,0,1)+coef_ttz(0,1)*alpha_2(lyp,0,1)
                coef_tyz(1,1,0,1)=coef_tyz(1,1,0,1)+coef_ttz(1,1)*alpha_2(lyp,0,1)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,2-lyb
             !DO lxa=0,1-lya
                coef_tyz(0,0,1,1)=coef_tyz(0,0,1,1)+coef_ttz(0,0)*alpha_2(lyp,1,1)
             !ENDDO
             !DO lxa=0,1-lya
                coef_tyz(0,1,1,1)=coef_tyz(0,1,1,1)+coef_ttz(0,1)*alpha_2(lyp,1,1)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
     ! RR do lya = 0, la_max_local
             !DO lxb=0,2-lyb
             !DO lxa=0,1-lya
                coef_tyz(0,0,0,2)=coef_tyz(0,0,0,2)+coef_ttz(0,0)*alpha_2(lyp,0,2)
                coef_tyz(1,0,0,2)=coef_tyz(1,0,0,2)+coef_ttz(1,0)*alpha_2(lyp,0,2)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,2-lyb
             !DO lxa=0,1-lya
                coef_tyz(0,0,1,2)=coef_tyz(0,0,1,2)+coef_ttz(0,0)*alpha_2(lyp,1,2)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
       ENDDO


          !DO lyb=0,2-lzb
          !DO lya=0,1-lza
              DO lxb=lb_min_local,2
             jco=coset(lxb,0,0)

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO

          !ENDDO
          !ENDDO
              DO lxb=lb_min_local,2
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,1-lza
              DO lxb=MAX(lb_min_local-0-1,0),1
             jco=coset(lxb,1,0)

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO

          !ENDDO
          !ENDDO
              DO lxb=MAX(lb_min_local-0-1,0),1
             jco=coset(lxb,1,0)

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,0,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,1-lza
             lxb = 0
             !jco=coset(0,2,0)
             jco=8

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,0,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,2,0)
             jco=8

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,2)*alpha_3(lzp,0,0)



          !ENDDO
          !ENDDO
          !DO lyb=0,2-lzb
          !DO lya=0,1-lza
              DO lxb=lb_min_local,2
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,1-lza
              DO lxb=MAX(lb_min_local-0-1,0),1
             jco=coset(lxb,1,0)

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,1,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,1-lza
             lxb = 0
             !jco=coset(0,2,0)
             jco=8

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,1,0)



          !ENDDO
          !ENDDO
          !DO lyb=0,2-lzb
          !DO lya=0,1-lza
              DO lxb=MAX(lb_min_local-1-0,0),1
             jco=coset(lxb,0,1)

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,1)

             ENDDO

             ENDDO

          !ENDDO
          !ENDDO
              DO lxb=MAX(lb_min_local-1-0,0),1
             jco=coset(lxb,0,1)

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,1)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,1-lza
             lxb = 0
             !jco=coset(0,1,1)
             jco=9

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,1)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,1,1)
             jco=9

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,0,1)



          !ENDDO
          !ENDDO
          !DO lyb=0,2-lzb
          !DO lya=0,1-lza
              DO lxb=MAX(lb_min_local-1-0,0),1
             jco=coset(lxb,0,1)

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,1)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,1-lza
             lxb = 0
             !jco=coset(0,1,1)
             jco=9

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,1,1)



          !ENDDO
          !ENDDO
          !DO lyb=0,2-lzb
          !DO lya=0,1-lza
             lxb = 0
             !jco=coset(0,0,2)
             jco=10

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,2)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,2)
             jco=10

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,2)



          !ENDDO
          !ENDDO
          !DO lyb=0,2-lzb
          !DO lya=0,1-lza
             lxb = 0
             !jco=coset(0,0,2)
             jco=10

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,2)



          !ENDDO
          !ENDDO

   ENDDO

    END SUBROUTINE xyz_to_vab_1_2

! *****************************************************************************
    SUBROUTINE xyz_to_vab_1_3

    !INTEGER :: my_rank,ierror
    !REAL(kind=dp) :: alpha(0:(1+3),0:1,0:3,3)

    INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp
    REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb

    REAL(kind=dp) :: alpha_1(0:(1+3),0:1,0:3)
    REAL(kind=dp) :: alpha_2(0:(1+3),0:1,0:3)
    REAL(kind=dp) :: alpha_3(0:(1+3),0:1,0:3)
    REAL(kind=dp) :: coef_ttz(0:1,0:3)
    REAL(kind=dp) :: coef_tyz(0:1,0:3,0:1,0:3)
    coef_xyz=coef_xyz*prefactor
!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max_local + lb_max_local
!   (current implementation is l**7)
!
!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***
  !RR    DO iaxis=1,3
  alpha_1(:,:,:)=0.0_dp
  alpha_2(:,:,:)=0.0_dp
  alpha_3(:,:,:)=0.0_dp

         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,0,0)=alpha_1(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,0,1)=alpha_1(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,0,1)=alpha_1(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,0,2)=alpha_1(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,0,2)=alpha_1(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,0,2)=alpha_1(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(3,0,3)=alpha_1(3,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(2,0,3)=alpha_1(2,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,0,3)=alpha_1(1,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,0,3)=alpha_1(0,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,1,0)=alpha_1(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,1,0)=alpha_1(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,1,1)=alpha_1(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,1,1)=alpha_1(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,1,1)=alpha_1(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,1,1)=alpha_1(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(3,1,2)=alpha_1(3,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(2,1,2)=alpha_1(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,1,2)=alpha_1(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,1,2)=alpha_1(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,1,2)=alpha_1(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,1,2)=alpha_1(0,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(4,1,3)=alpha_1(4,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(3,1,3)=alpha_1(3,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(2,1,3)=alpha_1(2,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,1,3)=alpha_1(1,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(3,1,3)=alpha_1(3,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(2,1,3)=alpha_1(2,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,1,3)=alpha_1(1,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,1,3)=alpha_1(0,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,0,0)=alpha_2(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,0,1)=alpha_2(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,0,1)=alpha_2(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,0,2)=alpha_2(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,0,2)=alpha_2(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,0,2)=alpha_2(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(3,0,3)=alpha_2(3,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(2,0,3)=alpha_2(2,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,0,3)=alpha_2(1,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,0,3)=alpha_2(0,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,1,0)=alpha_2(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,1,0)=alpha_2(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,1,1)=alpha_2(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,1,1)=alpha_2(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,1,1)=alpha_2(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,1,1)=alpha_2(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(3,1,2)=alpha_2(3,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(2,1,2)=alpha_2(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,1,2)=alpha_2(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,1,2)=alpha_2(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,1,2)=alpha_2(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,1,2)=alpha_2(0,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(4,1,3)=alpha_2(4,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(3,1,3)=alpha_2(3,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(2,1,3)=alpha_2(2,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,1,3)=alpha_2(1,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(3,1,3)=alpha_2(3,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(2,1,3)=alpha_2(2,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,1,3)=alpha_2(1,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,1,3)=alpha_2(0,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,0,0)=alpha_3(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,0,1)=alpha_3(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,0,1)=alpha_3(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,0,2)=alpha_3(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,0,2)=alpha_3(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,0,2)=alpha_3(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(3,0,3)=alpha_3(3,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(2,0,3)=alpha_3(2,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,0,3)=alpha_3(1,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,0,3)=alpha_3(0,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,1,0)=alpha_3(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,1,0)=alpha_3(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,1,1)=alpha_3(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,1,1)=alpha_3(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,1,1)=alpha_3(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,1,1)=alpha_3(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(3,1,2)=alpha_3(3,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(2,1,2)=alpha_3(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,1,2)=alpha_3(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,1,2)=alpha_3(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,1,2)=alpha_3(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,1,2)=alpha_3(0,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(4,1,3)=alpha_3(4,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(3,1,3)=alpha_3(3,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(2,1,3)=alpha_3(2,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,1,3)=alpha_3(1,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(3,1,3)=alpha_3(3,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(2,1,3)=alpha_3(2,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,1,3)=alpha_3(1,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,1,3)=alpha_3(0,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
  !RR    ENDDO
    !
    !   compute v_{lxa,lya,lza,lxb,lyb,lzb} given v_{lxp,lyp,lzp} and alpha(ls,lxa,lxb,1)
    !   use a three step procedure
    !

    lxyz=0
    DO lzp=0,1+3
       coef_tyz=0.0_dp
       DO lyp=0,1+3-lzp
          coef_ttz=0.0_dp
          DO lxp=0,1+3-lzp-lyp
             lxyz=lxyz+1
                coef_ttz(0,0)=coef_ttz(0,0)+coef_xyz(lxyz)*alpha_1(lxp,0,0)
                coef_ttz(1,0)=coef_ttz(1,0)+coef_xyz(lxyz)*alpha_1(lxp,1,0)
                coef_ttz(0,1)=coef_ttz(0,1)+coef_xyz(lxyz)*alpha_1(lxp,0,1)
                coef_ttz(1,1)=coef_ttz(1,1)+coef_xyz(lxyz)*alpha_1(lxp,1,1)
                coef_ttz(0,2)=coef_ttz(0,2)+coef_xyz(lxyz)*alpha_1(lxp,0,2)
                coef_ttz(1,2)=coef_ttz(1,2)+coef_xyz(lxyz)*alpha_1(lxp,1,2)
                coef_ttz(0,3)=coef_ttz(0,3)+coef_xyz(lxyz)*alpha_1(lxp,0,3)
                coef_ttz(1,3)=coef_ttz(1,3)+coef_xyz(lxyz)*alpha_1(lxp,1,3)
          ENDDO

     ! Improve locality
     !temp(:,:) = alpha(lyp,0:1,0:3,2)
     ! RR do lyb = 0, lb_max_local
     ! RR do lya = 0, la_max_local
             !DO lxb=0,3-lyb
             !DO lxa=0,1-lya
                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
                coef_tyz(1,0,0,0)=coef_tyz(1,0,0,0)+coef_ttz(1,0)*alpha_2(lyp,0,0)
             !ENDDO
             !DO lxa=0,1-lya
                coef_tyz(0,1,0,0)=coef_tyz(0,1,0,0)+coef_ttz(0,1)*alpha_2(lyp,0,0)
                coef_tyz(1,1,0,0)=coef_tyz(1,1,0,0)+coef_ttz(1,1)*alpha_2(lyp,0,0)
             !ENDDO
             !DO lxa=0,1-lya
                coef_tyz(0,2,0,0)=coef_tyz(0,2,0,0)+coef_ttz(0,2)*alpha_2(lyp,0,0)
                coef_tyz(1,2,0,0)=coef_tyz(1,2,0,0)+coef_ttz(1,2)*alpha_2(lyp,0,0)
             !ENDDO
             !DO lxa=0,1-lya
                coef_tyz(0,3,0,0)=coef_tyz(0,3,0,0)+coef_ttz(0,3)*alpha_2(lyp,0,0)
                coef_tyz(1,3,0,0)=coef_tyz(1,3,0,0)+coef_ttz(1,3)*alpha_2(lyp,0,0)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,3-lyb
             !DO lxa=0,1-lya
                coef_tyz(0,0,1,0)=coef_tyz(0,0,1,0)+coef_ttz(0,0)*alpha_2(lyp,1,0)
             !ENDDO
             !DO lxa=0,1-lya
                coef_tyz(0,1,1,0)=coef_tyz(0,1,1,0)+coef_ttz(0,1)*alpha_2(lyp,1,0)
             !ENDDO
             !DO lxa=0,1-lya
                coef_tyz(0,2,1,0)=coef_tyz(0,2,1,0)+coef_ttz(0,2)*alpha_2(lyp,1,0)
             !ENDDO
             !DO lxa=0,1-lya
                coef_tyz(0,3,1,0)=coef_tyz(0,3,1,0)+coef_ttz(0,3)*alpha_2(lyp,1,0)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
     ! RR do lya = 0, la_max_local
             !DO lxb=0,3-lyb
             !DO lxa=0,1-lya
                coef_tyz(0,0,0,1)=coef_tyz(0,0,0,1)+coef_ttz(0,0)*alpha_2(lyp,0,1)
                coef_tyz(1,0,0,1)=coef_tyz(1,0,0,1)+coef_ttz(1,0)*alpha_2(lyp,0,1)
             !ENDDO
             !DO lxa=0,1-lya
                coef_tyz(0,1,0,1)=coef_tyz(0,1,0,1)+coef_ttz(0,1)*alpha_2(lyp,0,1)
                coef_tyz(1,1,0,1)=coef_tyz(1,1,0,1)+coef_ttz(1,1)*alpha_2(lyp,0,1)
             !ENDDO
             !DO lxa=0,1-lya
                coef_tyz(0,2,0,1)=coef_tyz(0,2,0,1)+coef_ttz(0,2)*alpha_2(lyp,0,1)
                coef_tyz(1,2,0,1)=coef_tyz(1,2,0,1)+coef_ttz(1,2)*alpha_2(lyp,0,1)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,3-lyb
             !DO lxa=0,1-lya
                coef_tyz(0,0,1,1)=coef_tyz(0,0,1,1)+coef_ttz(0,0)*alpha_2(lyp,1,1)
             !ENDDO
             !DO lxa=0,1-lya
                coef_tyz(0,1,1,1)=coef_tyz(0,1,1,1)+coef_ttz(0,1)*alpha_2(lyp,1,1)
             !ENDDO
             !DO lxa=0,1-lya
                coef_tyz(0,2,1,1)=coef_tyz(0,2,1,1)+coef_ttz(0,2)*alpha_2(lyp,1,1)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
     ! RR do lya = 0, la_max_local
             !DO lxb=0,3-lyb
             !DO lxa=0,1-lya
                coef_tyz(0,0,0,2)=coef_tyz(0,0,0,2)+coef_ttz(0,0)*alpha_2(lyp,0,2)
                coef_tyz(1,0,0,2)=coef_tyz(1,0,0,2)+coef_ttz(1,0)*alpha_2(lyp,0,2)
             !ENDDO
             !DO lxa=0,1-lya
                coef_tyz(0,1,0,2)=coef_tyz(0,1,0,2)+coef_ttz(0,1)*alpha_2(lyp,0,2)
                coef_tyz(1,1,0,2)=coef_tyz(1,1,0,2)+coef_ttz(1,1)*alpha_2(lyp,0,2)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,3-lyb
             !DO lxa=0,1-lya
                coef_tyz(0,0,1,2)=coef_tyz(0,0,1,2)+coef_ttz(0,0)*alpha_2(lyp,1,2)
             !ENDDO
             !DO lxa=0,1-lya
                coef_tyz(0,1,1,2)=coef_tyz(0,1,1,2)+coef_ttz(0,1)*alpha_2(lyp,1,2)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
     ! RR do lya = 0, la_max_local
             !DO lxb=0,3-lyb
             !DO lxa=0,1-lya
                coef_tyz(0,0,0,3)=coef_tyz(0,0,0,3)+coef_ttz(0,0)*alpha_2(lyp,0,3)
                coef_tyz(1,0,0,3)=coef_tyz(1,0,0,3)+coef_ttz(1,0)*alpha_2(lyp,0,3)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,3-lyb
             !DO lxa=0,1-lya
                coef_tyz(0,0,1,3)=coef_tyz(0,0,1,3)+coef_ttz(0,0)*alpha_2(lyp,1,3)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
       ENDDO


          !DO lyb=0,3-lzb
          !DO lya=0,1-lza
              DO lxb=lb_min_local,3
             jco=coset(lxb,0,0)

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO

          !ENDDO
          !ENDDO
              DO lxb=lb_min_local,3
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,1-lza
              DO lxb=MAX(lb_min_local-0-1,0),2
             jco=coset(lxb,1,0)

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO

          !ENDDO
          !ENDDO
              DO lxb=MAX(lb_min_local-0-1,0),2
             jco=coset(lxb,1,0)

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,0,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,1-lza
              DO lxb=MAX(lb_min_local-0-2,0),1
             jco=coset(lxb,2,0)

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO

          !ENDDO
          !ENDDO
              DO lxb=MAX(lb_min_local-0-2,0),1
             jco=coset(lxb,2,0)

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,2)*alpha_3(lzp,0,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,1-lza
             lxb = 0
             !jco=coset(0,3,0)
             jco=17

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,3)*alpha_3(lzp,0,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,3,0)
             jco=17

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,3)*alpha_3(lzp,0,0)



          !ENDDO
          !ENDDO
          !DO lyb=0,3-lzb
          !DO lya=0,1-lza
              DO lxb=lb_min_local,3
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,1-lza
              DO lxb=MAX(lb_min_local-0-1,0),2
             jco=coset(lxb,1,0)

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,1,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,1-lza
              DO lxb=MAX(lb_min_local-0-2,0),1
             jco=coset(lxb,2,0)

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,1,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,1-lza
             lxb = 0
             !jco=coset(0,3,0)
             jco=17

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,3)*alpha_3(lzp,1,0)



          !ENDDO
          !ENDDO
          !DO lyb=0,3-lzb
          !DO lya=0,1-lza
              DO lxb=MAX(lb_min_local-1-0,0),2
             jco=coset(lxb,0,1)

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,1)

             ENDDO

             ENDDO

          !ENDDO
          !ENDDO
              DO lxb=MAX(lb_min_local-1-0,0),2
             jco=coset(lxb,0,1)

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,1)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,1-lza
              DO lxb=MAX(lb_min_local-1-1,0),1
             jco=coset(lxb,1,1)

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,1)

             ENDDO

             ENDDO

          !ENDDO
          !ENDDO
              DO lxb=MAX(lb_min_local-1-1,0),1
             jco=coset(lxb,1,1)

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,0,1)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,1-lza
             lxb = 0
             !jco=coset(0,2,1)
             jco=18

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,0,1)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,2,1)
             jco=18

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,2)*alpha_3(lzp,0,1)



          !ENDDO
          !ENDDO
          !DO lyb=0,3-lzb
          !DO lya=0,1-lza
              DO lxb=MAX(lb_min_local-1-0,0),2
             jco=coset(lxb,0,1)

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,1)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,1-lza
              DO lxb=MAX(lb_min_local-1-1,0),1
             jco=coset(lxb,1,1)

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,1,1)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,1-lza
             lxb = 0
             !jco=coset(0,2,1)
             jco=18

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,1,1)



          !ENDDO
          !ENDDO
          !DO lyb=0,3-lzb
          !DO lya=0,1-lza
              DO lxb=MAX(lb_min_local-2-0,0),1
             jco=coset(lxb,0,2)

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,2)

             ENDDO

             ENDDO

          !ENDDO
          !ENDDO
              DO lxb=MAX(lb_min_local-2-0,0),1
             jco=coset(lxb,0,2)

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,2)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,1-lza
             lxb = 0
             !jco=coset(0,1,2)
             jco=19

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,2)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,1,2)
             jco=19

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,0,2)



          !ENDDO
          !ENDDO
          !DO lyb=0,3-lzb
          !DO lya=0,1-lza
              DO lxb=MAX(lb_min_local-2-0,0),1
             jco=coset(lxb,0,2)

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,2)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,1-lza
             lxb = 0
             !jco=coset(0,1,2)
             jco=19

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,1,2)



          !ENDDO
          !ENDDO
          !DO lyb=0,3-lzb
          !DO lya=0,1-lza
             lxb = 0
             !jco=coset(0,0,3)
             jco=20

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,3)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,3)
             jco=20

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,3)



          !ENDDO
          !ENDDO
          !DO lyb=0,3-lzb
          !DO lya=0,1-lza
             lxb = 0
             !jco=coset(0,0,3)
             jco=20

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,3)



          !ENDDO
          !ENDDO

   ENDDO

    END SUBROUTINE xyz_to_vab_1_3

! *****************************************************************************
    SUBROUTINE xyz_to_vab_1_4

    INTEGER :: my_rank,ierror

    REAL(kind=dp) :: alpha(0:(1+4),0:1,0:4,3)
    REAL(kind=dp) :: coef_ttz(0:1,0:4)
    REAL(kind=dp) :: coef_tyz(0:1,0:4,0:1,0:4)

    coef_xyz=coef_xyz*prefactor

!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max + lb_max
!   (current implementation is l**7)
!
!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***

!**** RR : Values of la_max and lb_max
!    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)
!    WRITE(100+my_rank,*) "la=",la_max_local,"lb=",lb_max_local
!****
    alpha(:,:,:,:)=0.0_dp

    DO lxa=0,1
    DO lxb=0,4
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,1)=alpha(lxa-l+lxb-k,lxa,lxb,1)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(1)-(ra(1)+rab(1)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(1)+rp(1))
       ENDDO
    ENDDO
    ENDDO
    DO lxa=0,1
    DO lxb=0,4
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,2)=alpha(lxa-l+lxb-k,lxa,lxb,2)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(2)-(ra(2)+rab(2)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(2)+rp(2))
       ENDDO
    ENDDO
    ENDDO
    DO lxa=0,1
    DO lxb=0,4
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,3)=alpha(lxa-l+lxb-k,lxa,lxb,3)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(3)-(ra(3)+rab(3)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(3)+rp(3))
       ENDDO
    ENDDO
    ENDDO

    !
    !   compute v_{lxa,lya,lza,lxb,lyb,lzb} given v_{lxp,lyp,lzp} and alpha(ls,lxa,lxb,1)
    !   use a three step procedure
    !

    lxyz=0
    DO lzp=0,5
       coef_tyz=0.0_dp
       DO lyp=0,5-lzp
          coef_ttz=0.0_dp
          DO lxp=0,5-lzp-lyp
             lxyz=lxyz+1
             DO lxb=0,4
             DO lxa=0,1
                coef_ttz(lxa,lxb)=coef_ttz(lxa,lxb)+coef_xyz(lxyz)*alpha(lxp,lxa,lxb,1)
             ENDDO
             ENDDO

          ENDDO ! lxp

          DO lyb=0,4
          DO lya=0,1
             DO lxb=0,4-lyb
             DO lxa=0,1-lya
                coef_tyz(lxa,lxb,lya,lyb)=coef_tyz(lxa,lxb,lya,lyb)+coef_ttz(lxa,lxb)*alpha(lyp,lya,lyb,2)
             ENDDO
             ENDDO
          ENDDO
          ENDDO

       ENDDO !lyp

       DO lzb=0,4
       DO lza=0,1
          DO lyb=0,4-lzb
          DO lya=0,1-lza
             DO lxb=MAX(lb_min_local-lzb-lyb,0),4-lzb-lyb
             jco=coset(lxb,lyb,lzb)
             DO lxa=MAX(la_min_local-lza-lya,0),1-lza-lya
                ico=coset(lxa,lya,lza)
                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,lya,lyb)*alpha(lzp,lza,lzb,3)
             ENDDO
             ENDDO
          ENDDO
          ENDDO
       ENDDO
       ENDDO

    ENDDO

    END SUBROUTINE xyz_to_vab_1_4! *****************************************************************************
    SUBROUTINE xyz_to_vab_2_0

    !INTEGER :: my_rank,ierror
    !REAL(kind=dp) :: alpha(0:(2+0),0:2,0:0,3)

    INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp
    REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb

    REAL(kind=dp) :: alpha_1(0:(2+0),0:2,0:0)
    REAL(kind=dp) :: alpha_2(0:(2+0),0:2,0:0)
    REAL(kind=dp) :: alpha_3(0:(2+0),0:2,0:0)
    REAL(kind=dp) :: coef_ttz(0:2,0:0)
    REAL(kind=dp) :: coef_tyz(0:2,0:0,0:2,0:0)
    coef_xyz=coef_xyz*prefactor
!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max_local + lb_max_local
!   (current implementation is l**7)
!
!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***
  !RR    DO iaxis=1,3
  alpha_1(:,:,:)=0.0_dp
  alpha_2(:,:,:)=0.0_dp
  alpha_3(:,:,:)=0.0_dp

         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,0,0)=alpha_1(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,1,0)=alpha_1(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,1,0)=alpha_1(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,2,0)=alpha_1(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,2,0)=alpha_1(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,2,0)=alpha_1(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,0,0)=alpha_2(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,1,0)=alpha_2(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,1,0)=alpha_2(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,2,0)=alpha_2(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,2,0)=alpha_2(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,2,0)=alpha_2(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,0,0)=alpha_3(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,1,0)=alpha_3(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,1,0)=alpha_3(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,2,0)=alpha_3(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,2,0)=alpha_3(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,2,0)=alpha_3(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
  !RR    ENDDO
    !
    !   compute v_{lxa,lya,lza,lxb,lyb,lzb} given v_{lxp,lyp,lzp} and alpha(ls,lxa,lxb,1)
    !   use a three step procedure
    !

    lxyz=0
    DO lzp=0,2+0
       coef_tyz=0.0_dp
       DO lyp=0,2+0-lzp
          coef_ttz=0.0_dp
          DO lxp=0,2+0-lzp-lyp
             lxyz=lxyz+1
                coef_ttz(0,0)=coef_ttz(0,0)+coef_xyz(lxyz)*alpha_1(lxp,0,0)
                coef_ttz(1,0)=coef_ttz(1,0)+coef_xyz(lxyz)*alpha_1(lxp,1,0)
                coef_ttz(2,0)=coef_ttz(2,0)+coef_xyz(lxyz)*alpha_1(lxp,2,0)
          ENDDO

     ! Improve locality
     !temp(:,:) = alpha(lyp,0:2,0:0,2)
     ! RR do lyb = 0, lb_max_local
     ! RR do lya = 0, la_max_local
             !DO lxb=0,0-lyb
             !DO lxa=0,2-lya
                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
                coef_tyz(1,0,0,0)=coef_tyz(1,0,0,0)+coef_ttz(1,0)*alpha_2(lyp,0,0)
                coef_tyz(2,0,0,0)=coef_tyz(2,0,0,0)+coef_ttz(2,0)*alpha_2(lyp,0,0)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,0-lyb
             !DO lxa=0,2-lya
                coef_tyz(0,0,1,0)=coef_tyz(0,0,1,0)+coef_ttz(0,0)*alpha_2(lyp,1,0)
                coef_tyz(1,0,1,0)=coef_tyz(1,0,1,0)+coef_ttz(1,0)*alpha_2(lyp,1,0)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,0-lyb
             !DO lxa=0,2-lya
                coef_tyz(0,0,2,0)=coef_tyz(0,0,2,0)+coef_ttz(0,0)*alpha_2(lyp,2,0)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
       ENDDO


          !DO lyb=0,0-lzb
          !DO lya=0,2-lza
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=la_min_local,2

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-0-1,0),1

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,2,0)
                ico=8

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,0,0)



          !ENDDO
          !ENDDO
          !DO lyb=0,0-lzb
          !DO lya=0,2-lza
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-1-0,0),1

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,1,1)
                ico=9

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,1,0)



          !ENDDO
          !ENDDO
          !DO lyb=0,0-lzb
          !DO lya=0,2-lza
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,0,2)
                ico=10

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,2,0)



          !ENDDO
          !ENDDO

   ENDDO

    END SUBROUTINE xyz_to_vab_2_0

! *****************************************************************************
    SUBROUTINE xyz_to_vab_2_1

    !INTEGER :: my_rank,ierror
    !REAL(kind=dp) :: alpha(0:(2+1),0:2,0:1,3)

    INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp
    REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb

    REAL(kind=dp) :: alpha_1(0:(2+1),0:2,0:1)
    REAL(kind=dp) :: alpha_2(0:(2+1),0:2,0:1)
    REAL(kind=dp) :: alpha_3(0:(2+1),0:2,0:1)
    REAL(kind=dp) :: coef_ttz(0:2,0:1)
    REAL(kind=dp) :: coef_tyz(0:2,0:1,0:2,0:1)
    coef_xyz=coef_xyz*prefactor
!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max_local + lb_max_local
!   (current implementation is l**7)
!
!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***
  !RR    DO iaxis=1,3
  alpha_1(:,:,:)=0.0_dp
  alpha_2(:,:,:)=0.0_dp
  alpha_3(:,:,:)=0.0_dp

         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,0,0)=alpha_1(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,0,1)=alpha_1(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,0,1)=alpha_1(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,1,0)=alpha_1(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,1,0)=alpha_1(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,1,1)=alpha_1(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,1,1)=alpha_1(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,1,1)=alpha_1(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,1,1)=alpha_1(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,2,0)=alpha_1(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,2,0)=alpha_1(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,2,0)=alpha_1(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(3,2,1)=alpha_1(3,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(2,2,1)=alpha_1(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,2,1)=alpha_1(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,2,1)=alpha_1(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,2,1)=alpha_1(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,2,1)=alpha_1(0,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,0,0)=alpha_2(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,0,1)=alpha_2(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,0,1)=alpha_2(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,1,0)=alpha_2(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,1,0)=alpha_2(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,1,1)=alpha_2(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,1,1)=alpha_2(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,1,1)=alpha_2(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,1,1)=alpha_2(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,2,0)=alpha_2(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,2,0)=alpha_2(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,2,0)=alpha_2(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(3,2,1)=alpha_2(3,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(2,2,1)=alpha_2(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,2,1)=alpha_2(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,2,1)=alpha_2(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,2,1)=alpha_2(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,2,1)=alpha_2(0,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,0,0)=alpha_3(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,0,1)=alpha_3(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,0,1)=alpha_3(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,1,0)=alpha_3(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,1,0)=alpha_3(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,1,1)=alpha_3(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,1,1)=alpha_3(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,1,1)=alpha_3(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,1,1)=alpha_3(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,2,0)=alpha_3(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,2,0)=alpha_3(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,2,0)=alpha_3(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(3,2,1)=alpha_3(3,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(2,2,1)=alpha_3(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,2,1)=alpha_3(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,2,1)=alpha_3(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,2,1)=alpha_3(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,2,1)=alpha_3(0,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
  !RR    ENDDO
    !
    !   compute v_{lxa,lya,lza,lxb,lyb,lzb} given v_{lxp,lyp,lzp} and alpha(ls,lxa,lxb,1)
    !   use a three step procedure
    !

    lxyz=0
    DO lzp=0,2+1
       coef_tyz=0.0_dp
       DO lyp=0,2+1-lzp
          coef_ttz=0.0_dp
          DO lxp=0,2+1-lzp-lyp
             lxyz=lxyz+1
                coef_ttz(0,0)=coef_ttz(0,0)+coef_xyz(lxyz)*alpha_1(lxp,0,0)
                coef_ttz(1,0)=coef_ttz(1,0)+coef_xyz(lxyz)*alpha_1(lxp,1,0)
                coef_ttz(2,0)=coef_ttz(2,0)+coef_xyz(lxyz)*alpha_1(lxp,2,0)
                coef_ttz(0,1)=coef_ttz(0,1)+coef_xyz(lxyz)*alpha_1(lxp,0,1)
                coef_ttz(1,1)=coef_ttz(1,1)+coef_xyz(lxyz)*alpha_1(lxp,1,1)
                coef_ttz(2,1)=coef_ttz(2,1)+coef_xyz(lxyz)*alpha_1(lxp,2,1)
          ENDDO

     ! Improve locality
     !temp(:,:) = alpha(lyp,0:2,0:1,2)
     ! RR do lyb = 0, lb_max_local
     ! RR do lya = 0, la_max_local
             !DO lxb=0,1-lyb
             !DO lxa=0,2-lya
                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
                coef_tyz(1,0,0,0)=coef_tyz(1,0,0,0)+coef_ttz(1,0)*alpha_2(lyp,0,0)
                coef_tyz(2,0,0,0)=coef_tyz(2,0,0,0)+coef_ttz(2,0)*alpha_2(lyp,0,0)
             !ENDDO
             !DO lxa=0,2-lya
                coef_tyz(0,1,0,0)=coef_tyz(0,1,0,0)+coef_ttz(0,1)*alpha_2(lyp,0,0)
                coef_tyz(1,1,0,0)=coef_tyz(1,1,0,0)+coef_ttz(1,1)*alpha_2(lyp,0,0)
                coef_tyz(2,1,0,0)=coef_tyz(2,1,0,0)+coef_ttz(2,1)*alpha_2(lyp,0,0)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,1-lyb
             !DO lxa=0,2-lya
                coef_tyz(0,0,1,0)=coef_tyz(0,0,1,0)+coef_ttz(0,0)*alpha_2(lyp,1,0)
                coef_tyz(1,0,1,0)=coef_tyz(1,0,1,0)+coef_ttz(1,0)*alpha_2(lyp,1,0)
             !ENDDO
             !DO lxa=0,2-lya
                coef_tyz(0,1,1,0)=coef_tyz(0,1,1,0)+coef_ttz(0,1)*alpha_2(lyp,1,0)
                coef_tyz(1,1,1,0)=coef_tyz(1,1,1,0)+coef_ttz(1,1)*alpha_2(lyp,1,0)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,1-lyb
             !DO lxa=0,2-lya
                coef_tyz(0,0,2,0)=coef_tyz(0,0,2,0)+coef_ttz(0,0)*alpha_2(lyp,2,0)
             !ENDDO
             !DO lxa=0,2-lya
                coef_tyz(0,1,2,0)=coef_tyz(0,1,2,0)+coef_ttz(0,1)*alpha_2(lyp,2,0)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
     ! RR do lya = 0, la_max_local
             !DO lxb=0,1-lyb
             !DO lxa=0,2-lya
                coef_tyz(0,0,0,1)=coef_tyz(0,0,0,1)+coef_ttz(0,0)*alpha_2(lyp,0,1)
                coef_tyz(1,0,0,1)=coef_tyz(1,0,0,1)+coef_ttz(1,0)*alpha_2(lyp,0,1)
                coef_tyz(2,0,0,1)=coef_tyz(2,0,0,1)+coef_ttz(2,0)*alpha_2(lyp,0,1)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,1-lyb
             !DO lxa=0,2-lya
                coef_tyz(0,0,1,1)=coef_tyz(0,0,1,1)+coef_ttz(0,0)*alpha_2(lyp,1,1)
                coef_tyz(1,0,1,1)=coef_tyz(1,0,1,1)+coef_ttz(1,0)*alpha_2(lyp,1,1)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,1-lyb
             !DO lxa=0,2-lya
                coef_tyz(0,0,2,1)=coef_tyz(0,0,2,1)+coef_ttz(0,0)*alpha_2(lyp,2,1)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
       ENDDO


          !DO lyb=0,1-lzb
          !DO lya=0,2-lza
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

              DO lxa=la_min_local,2

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO

          !ENDDO
          !ENDDO
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

              DO lxa=MAX(la_min_local-0-1,0),1

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO

          !ENDDO
          !ENDDO
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,2,0)
                ico=8

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,0,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,2-lza
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

              DO lxa=la_min_local,2

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

              DO lxa=MAX(la_min_local-0-1,0),1

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,0,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

                lxa = 0
                !ico=coset(0,2,0)
                ico=8

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,1)*alpha_3(lzp,0,0)



          !ENDDO
          !ENDDO
          !DO lyb=0,1-lzb
          !DO lya=0,2-lza
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

              DO lxa=MAX(la_min_local-1-0,0),1

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,0)

             ENDDO

             ENDDO

          !ENDDO
          !ENDDO
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,1,1)
                ico=9

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,1,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,2-lza
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

              DO lxa=MAX(la_min_local-1-0,0),1

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,1,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

                lxa = 0
                !ico=coset(0,1,1)
                ico=9

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,1,0)



          !ENDDO
          !ENDDO
          !DO lyb=0,1-lzb
          !DO lya=0,2-lza
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,0,2)
                ico=10

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,2,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,2-lza
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

                lxa = 0
                !ico=coset(0,0,2)
                ico=10

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,2,0)



          !ENDDO
          !ENDDO
          !DO lyb=0,1-lzb
          !DO lya=0,2-lza
             lxb = 0
             !jco=coset(0,0,1)
             jco=4

              DO lxa=la_min_local,2

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,1)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,1)
             jco=4

              DO lxa=MAX(la_min_local-0-1,0),1

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,1)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,1)
             jco=4

                lxa = 0
                !ico=coset(0,2,0)
                ico=8

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,0,1)



          !ENDDO
          !ENDDO
          !DO lyb=0,1-lzb
          !DO lya=0,2-lza
             lxb = 0
             !jco=coset(0,0,1)
             jco=4

              DO lxa=MAX(la_min_local-1-0,0),1

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,1)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,1)
             jco=4

                lxa = 0
                !ico=coset(0,1,1)
                ico=9

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,1,1)



          !ENDDO
          !ENDDO
          !DO lyb=0,1-lzb
          !DO lya=0,2-lza
             lxb = 0
             !jco=coset(0,0,1)
             jco=4

                lxa = 0
                !ico=coset(0,0,2)
                ico=10

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,2,1)



          !ENDDO
          !ENDDO

   ENDDO

    END SUBROUTINE xyz_to_vab_2_1

! *****************************************************************************
    SUBROUTINE xyz_to_vab_2_2

    !INTEGER :: my_rank,ierror
    !REAL(kind=dp) :: alpha(0:(2+2),0:2,0:2,3)

    INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp
    REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb

    REAL(kind=dp) :: alpha_1(0:(2+2),0:2,0:2)
    REAL(kind=dp) :: alpha_2(0:(2+2),0:2,0:2)
    REAL(kind=dp) :: alpha_3(0:(2+2),0:2,0:2)
    REAL(kind=dp) :: coef_ttz(0:2,0:2)
    REAL(kind=dp) :: coef_tyz(0:2,0:2,0:2,0:2)
    coef_xyz=coef_xyz*prefactor
!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max_local + lb_max_local
!   (current implementation is l**7)
!
!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***
  !RR    DO iaxis=1,3
  alpha_1(:,:,:)=0.0_dp
  alpha_2(:,:,:)=0.0_dp
  alpha_3(:,:,:)=0.0_dp

         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,0,0)=alpha_1(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,0,1)=alpha_1(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,0,1)=alpha_1(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,0,2)=alpha_1(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,0,2)=alpha_1(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,0,2)=alpha_1(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,1,0)=alpha_1(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,1,0)=alpha_1(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,1,1)=alpha_1(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,1,1)=alpha_1(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,1,1)=alpha_1(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,1,1)=alpha_1(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(3,1,2)=alpha_1(3,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(2,1,2)=alpha_1(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,1,2)=alpha_1(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,1,2)=alpha_1(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,1,2)=alpha_1(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,1,2)=alpha_1(0,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,2,0)=alpha_1(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,2,0)=alpha_1(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,2,0)=alpha_1(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(3,2,1)=alpha_1(3,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(2,2,1)=alpha_1(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,2,1)=alpha_1(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,2,1)=alpha_1(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,2,1)=alpha_1(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,2,1)=alpha_1(0,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(4,2,2)=alpha_1(4,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(3,2,2)=alpha_1(3,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(2,2,2)=alpha_1(2,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(3,2,2)=alpha_1(3,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(2,2,2)=alpha_1(2,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,2,2)=alpha_1(1,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,2,2)=alpha_1(2,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,2,2)=alpha_1(1,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,2,2)=alpha_1(0,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,0,0)=alpha_2(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,0,1)=alpha_2(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,0,1)=alpha_2(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,0,2)=alpha_2(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,0,2)=alpha_2(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,0,2)=alpha_2(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,1,0)=alpha_2(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,1,0)=alpha_2(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,1,1)=alpha_2(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,1,1)=alpha_2(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,1,1)=alpha_2(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,1,1)=alpha_2(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(3,1,2)=alpha_2(3,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(2,1,2)=alpha_2(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,1,2)=alpha_2(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,1,2)=alpha_2(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,1,2)=alpha_2(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,1,2)=alpha_2(0,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,2,0)=alpha_2(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,2,0)=alpha_2(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,2,0)=alpha_2(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(3,2,1)=alpha_2(3,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(2,2,1)=alpha_2(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,2,1)=alpha_2(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,2,1)=alpha_2(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,2,1)=alpha_2(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,2,1)=alpha_2(0,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(4,2,2)=alpha_2(4,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(3,2,2)=alpha_2(3,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(2,2,2)=alpha_2(2,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(3,2,2)=alpha_2(3,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(2,2,2)=alpha_2(2,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,2,2)=alpha_2(1,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,2,2)=alpha_2(2,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,2,2)=alpha_2(1,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,2,2)=alpha_2(0,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,0,0)=alpha_3(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,0,1)=alpha_3(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,0,1)=alpha_3(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,0,2)=alpha_3(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,0,2)=alpha_3(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,0,2)=alpha_3(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,1,0)=alpha_3(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,1,0)=alpha_3(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,1,1)=alpha_3(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,1,1)=alpha_3(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,1,1)=alpha_3(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,1,1)=alpha_3(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(3,1,2)=alpha_3(3,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(2,1,2)=alpha_3(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,1,2)=alpha_3(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,1,2)=alpha_3(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,1,2)=alpha_3(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,1,2)=alpha_3(0,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,2,0)=alpha_3(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,2,0)=alpha_3(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,2,0)=alpha_3(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(3,2,1)=alpha_3(3,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(2,2,1)=alpha_3(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,2,1)=alpha_3(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,2,1)=alpha_3(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,2,1)=alpha_3(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,2,1)=alpha_3(0,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(4,2,2)=alpha_3(4,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(3,2,2)=alpha_3(3,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(2,2,2)=alpha_3(2,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(3,2,2)=alpha_3(3,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(2,2,2)=alpha_3(2,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,2,2)=alpha_3(1,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,2,2)=alpha_3(2,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,2,2)=alpha_3(1,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,2,2)=alpha_3(0,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
  !RR    ENDDO
    !
    !   compute v_{lxa,lya,lza,lxb,lyb,lzb} given v_{lxp,lyp,lzp} and alpha(ls,lxa,lxb,1)
    !   use a three step procedure
    !

    lxyz=0
    DO lzp=0,2+2
       coef_tyz=0.0_dp
       DO lyp=0,2+2-lzp
          coef_ttz=0.0_dp
          DO lxp=0,2+2-lzp-lyp
             lxyz=lxyz+1
                coef_ttz(0,0)=coef_ttz(0,0)+coef_xyz(lxyz)*alpha_1(lxp,0,0)
                coef_ttz(1,0)=coef_ttz(1,0)+coef_xyz(lxyz)*alpha_1(lxp,1,0)
                coef_ttz(2,0)=coef_ttz(2,0)+coef_xyz(lxyz)*alpha_1(lxp,2,0)
                coef_ttz(0,1)=coef_ttz(0,1)+coef_xyz(lxyz)*alpha_1(lxp,0,1)
                coef_ttz(1,1)=coef_ttz(1,1)+coef_xyz(lxyz)*alpha_1(lxp,1,1)
                coef_ttz(2,1)=coef_ttz(2,1)+coef_xyz(lxyz)*alpha_1(lxp,2,1)
                coef_ttz(0,2)=coef_ttz(0,2)+coef_xyz(lxyz)*alpha_1(lxp,0,2)
                coef_ttz(1,2)=coef_ttz(1,2)+coef_xyz(lxyz)*alpha_1(lxp,1,2)
                coef_ttz(2,2)=coef_ttz(2,2)+coef_xyz(lxyz)*alpha_1(lxp,2,2)
          ENDDO

     ! Improve locality
     !temp(:,:) = alpha(lyp,0:2,0:2,2)
     ! RR do lyb = 0, lb_max_local
     ! RR do lya = 0, la_max_local
             !DO lxb=0,2-lyb
             !DO lxa=0,2-lya
                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
                coef_tyz(1,0,0,0)=coef_tyz(1,0,0,0)+coef_ttz(1,0)*alpha_2(lyp,0,0)
                coef_tyz(2,0,0,0)=coef_tyz(2,0,0,0)+coef_ttz(2,0)*alpha_2(lyp,0,0)
             !ENDDO
             !DO lxa=0,2-lya
                coef_tyz(0,1,0,0)=coef_tyz(0,1,0,0)+coef_ttz(0,1)*alpha_2(lyp,0,0)
                coef_tyz(1,1,0,0)=coef_tyz(1,1,0,0)+coef_ttz(1,1)*alpha_2(lyp,0,0)
                coef_tyz(2,1,0,0)=coef_tyz(2,1,0,0)+coef_ttz(2,1)*alpha_2(lyp,0,0)
             !ENDDO
             !DO lxa=0,2-lya
                coef_tyz(0,2,0,0)=coef_tyz(0,2,0,0)+coef_ttz(0,2)*alpha_2(lyp,0,0)
                coef_tyz(1,2,0,0)=coef_tyz(1,2,0,0)+coef_ttz(1,2)*alpha_2(lyp,0,0)
                coef_tyz(2,2,0,0)=coef_tyz(2,2,0,0)+coef_ttz(2,2)*alpha_2(lyp,0,0)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,2-lyb
             !DO lxa=0,2-lya
                coef_tyz(0,0,1,0)=coef_tyz(0,0,1,0)+coef_ttz(0,0)*alpha_2(lyp,1,0)
                coef_tyz(1,0,1,0)=coef_tyz(1,0,1,0)+coef_ttz(1,0)*alpha_2(lyp,1,0)
             !ENDDO
             !DO lxa=0,2-lya
                coef_tyz(0,1,1,0)=coef_tyz(0,1,1,0)+coef_ttz(0,1)*alpha_2(lyp,1,0)
                coef_tyz(1,1,1,0)=coef_tyz(1,1,1,0)+coef_ttz(1,1)*alpha_2(lyp,1,0)
             !ENDDO
             !DO lxa=0,2-lya
                coef_tyz(0,2,1,0)=coef_tyz(0,2,1,0)+coef_ttz(0,2)*alpha_2(lyp,1,0)
                coef_tyz(1,2,1,0)=coef_tyz(1,2,1,0)+coef_ttz(1,2)*alpha_2(lyp,1,0)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,2-lyb
             !DO lxa=0,2-lya
                coef_tyz(0,0,2,0)=coef_tyz(0,0,2,0)+coef_ttz(0,0)*alpha_2(lyp,2,0)
             !ENDDO
             !DO lxa=0,2-lya
                coef_tyz(0,1,2,0)=coef_tyz(0,1,2,0)+coef_ttz(0,1)*alpha_2(lyp,2,0)
             !ENDDO
             !DO lxa=0,2-lya
                coef_tyz(0,2,2,0)=coef_tyz(0,2,2,0)+coef_ttz(0,2)*alpha_2(lyp,2,0)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
     ! RR do lya = 0, la_max_local
             !DO lxb=0,2-lyb
             !DO lxa=0,2-lya
                coef_tyz(0,0,0,1)=coef_tyz(0,0,0,1)+coef_ttz(0,0)*alpha_2(lyp,0,1)
                coef_tyz(1,0,0,1)=coef_tyz(1,0,0,1)+coef_ttz(1,0)*alpha_2(lyp,0,1)
                coef_tyz(2,0,0,1)=coef_tyz(2,0,0,1)+coef_ttz(2,0)*alpha_2(lyp,0,1)
             !ENDDO
             !DO lxa=0,2-lya
                coef_tyz(0,1,0,1)=coef_tyz(0,1,0,1)+coef_ttz(0,1)*alpha_2(lyp,0,1)
                coef_tyz(1,1,0,1)=coef_tyz(1,1,0,1)+coef_ttz(1,1)*alpha_2(lyp,0,1)
                coef_tyz(2,1,0,1)=coef_tyz(2,1,0,1)+coef_ttz(2,1)*alpha_2(lyp,0,1)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,2-lyb
             !DO lxa=0,2-lya
                coef_tyz(0,0,1,1)=coef_tyz(0,0,1,1)+coef_ttz(0,0)*alpha_2(lyp,1,1)
                coef_tyz(1,0,1,1)=coef_tyz(1,0,1,1)+coef_ttz(1,0)*alpha_2(lyp,1,1)
             !ENDDO
             !DO lxa=0,2-lya
                coef_tyz(0,1,1,1)=coef_tyz(0,1,1,1)+coef_ttz(0,1)*alpha_2(lyp,1,1)
                coef_tyz(1,1,1,1)=coef_tyz(1,1,1,1)+coef_ttz(1,1)*alpha_2(lyp,1,1)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,2-lyb
             !DO lxa=0,2-lya
                coef_tyz(0,0,2,1)=coef_tyz(0,0,2,1)+coef_ttz(0,0)*alpha_2(lyp,2,1)
             !ENDDO
             !DO lxa=0,2-lya
                coef_tyz(0,1,2,1)=coef_tyz(0,1,2,1)+coef_ttz(0,1)*alpha_2(lyp,2,1)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
     ! RR do lya = 0, la_max_local
             !DO lxb=0,2-lyb
             !DO lxa=0,2-lya
                coef_tyz(0,0,0,2)=coef_tyz(0,0,0,2)+coef_ttz(0,0)*alpha_2(lyp,0,2)
                coef_tyz(1,0,0,2)=coef_tyz(1,0,0,2)+coef_ttz(1,0)*alpha_2(lyp,0,2)
                coef_tyz(2,0,0,2)=coef_tyz(2,0,0,2)+coef_ttz(2,0)*alpha_2(lyp,0,2)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,2-lyb
             !DO lxa=0,2-lya
                coef_tyz(0,0,1,2)=coef_tyz(0,0,1,2)+coef_ttz(0,0)*alpha_2(lyp,1,2)
                coef_tyz(1,0,1,2)=coef_tyz(1,0,1,2)+coef_ttz(1,0)*alpha_2(lyp,1,2)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,2-lyb
             !DO lxa=0,2-lya
                coef_tyz(0,0,2,2)=coef_tyz(0,0,2,2)+coef_ttz(0,0)*alpha_2(lyp,2,2)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
       ENDDO


          !DO lyb=0,2-lzb
          !DO lya=0,2-lza
              DO lxb=lb_min_local,2
             jco=coset(lxb,0,0)

              DO lxa=la_min_local,2

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO

          !ENDDO
          !ENDDO
              DO lxb=lb_min_local,2
             jco=coset(lxb,0,0)

              DO lxa=MAX(la_min_local-0-1,0),1

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO

          !ENDDO
          !ENDDO
              DO lxb=lb_min_local,2
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,2,0)
                ico=8

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,0,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,2-lza
              DO lxb=MAX(lb_min_local-0-1,0),1
             jco=coset(lxb,1,0)

              DO lxa=la_min_local,2

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO

          !ENDDO
          !ENDDO
              DO lxb=MAX(lb_min_local-0-1,0),1
             jco=coset(lxb,1,0)

              DO lxa=MAX(la_min_local-0-1,0),1

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO

          !ENDDO
          !ENDDO
              DO lxb=MAX(lb_min_local-0-1,0),1
             jco=coset(lxb,1,0)

                lxa = 0
                !ico=coset(0,2,0)
                ico=8

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,1)*alpha_3(lzp,0,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,2-lza
             lxb = 0
             !jco=coset(0,2,0)
             jco=8

              DO lxa=la_min_local,2

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,0,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,2,0)
             jco=8

              DO lxa=MAX(la_min_local-0-1,0),1

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,2)*alpha_3(lzp,0,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,2,0)
             jco=8

                lxa = 0
                !ico=coset(0,2,0)
                ico=8

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,2)*alpha_3(lzp,0,0)



          !ENDDO
          !ENDDO
          !DO lyb=0,2-lzb
          !DO lya=0,2-lza
              DO lxb=lb_min_local,2
             jco=coset(lxb,0,0)

              DO lxa=MAX(la_min_local-1-0,0),1

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,0)

             ENDDO

             ENDDO

          !ENDDO
          !ENDDO
              DO lxb=lb_min_local,2
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,1,1)
                ico=9

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,1,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,2-lza
              DO lxb=MAX(lb_min_local-0-1,0),1
             jco=coset(lxb,1,0)

              DO lxa=MAX(la_min_local-1-0,0),1

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,1,0)

             ENDDO

             ENDDO

          !ENDDO
          !ENDDO
              DO lxb=MAX(lb_min_local-0-1,0),1
             jco=coset(lxb,1,0)

                lxa = 0
                !ico=coset(0,1,1)
                ico=9

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,1,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,2-lza
             lxb = 0
             !jco=coset(0,2,0)
             jco=8

              DO lxa=MAX(la_min_local-1-0,0),1

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,1,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,2,0)
             jco=8

                lxa = 0
                !ico=coset(0,1,1)
                ico=9

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,2)*alpha_3(lzp,1,0)



          !ENDDO
          !ENDDO
          !DO lyb=0,2-lzb
          !DO lya=0,2-lza
              DO lxb=lb_min_local,2
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,0,2)
                ico=10

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,2,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,2-lza
              DO lxb=MAX(lb_min_local-0-1,0),1
             jco=coset(lxb,1,0)

                lxa = 0
                !ico=coset(0,0,2)
                ico=10

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,2,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,2-lza
             lxb = 0
             !jco=coset(0,2,0)
             jco=8

                lxa = 0
                !ico=coset(0,0,2)
                ico=10

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,2,0)



          !ENDDO
          !ENDDO
          !DO lyb=0,2-lzb
          !DO lya=0,2-lza
              DO lxb=MAX(lb_min_local-1-0,0),1
             jco=coset(lxb,0,1)

              DO lxa=la_min_local,2

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,1)

             ENDDO

             ENDDO

          !ENDDO
          !ENDDO
              DO lxb=MAX(lb_min_local-1-0,0),1
             jco=coset(lxb,0,1)

              DO lxa=MAX(la_min_local-0-1,0),1

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,1)

             ENDDO

             ENDDO

          !ENDDO
          !ENDDO
              DO lxb=MAX(lb_min_local-1-0,0),1
             jco=coset(lxb,0,1)

                lxa = 0
                !ico=coset(0,2,0)
                ico=8

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,0,1)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,2-lza
             lxb = 0
             !jco=coset(0,1,1)
             jco=9

              DO lxa=la_min_local,2

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,1)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,1,1)
             jco=9

              DO lxa=MAX(la_min_local-0-1,0),1

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,0,1)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,1,1)
             jco=9

                lxa = 0
                !ico=coset(0,2,0)
                ico=8

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,1)*alpha_3(lzp,0,1)



          !ENDDO
          !ENDDO
          !DO lyb=0,2-lzb
          !DO lya=0,2-lza
              DO lxb=MAX(lb_min_local-1-0,0),1
             jco=coset(lxb,0,1)

              DO lxa=MAX(la_min_local-1-0,0),1

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,1)

             ENDDO

             ENDDO

          !ENDDO
          !ENDDO
              DO lxb=MAX(lb_min_local-1-0,0),1
             jco=coset(lxb,0,1)

                lxa = 0
                !ico=coset(0,1,1)
                ico=9

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,1,1)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,2-lza
             lxb = 0
             !jco=coset(0,1,1)
             jco=9

              DO lxa=MAX(la_min_local-1-0,0),1

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,1,1)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,1,1)
             jco=9

                lxa = 0
                !ico=coset(0,1,1)
                ico=9

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,1,1)



          !ENDDO
          !ENDDO
          !DO lyb=0,2-lzb
          !DO lya=0,2-lza
              DO lxb=MAX(lb_min_local-1-0,0),1
             jco=coset(lxb,0,1)

                lxa = 0
                !ico=coset(0,0,2)
                ico=10

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,2,1)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,2-lza
             lxb = 0
             !jco=coset(0,1,1)
             jco=9

                lxa = 0
                !ico=coset(0,0,2)
                ico=10

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,2,1)



          !ENDDO
          !ENDDO
          !DO lyb=0,2-lzb
          !DO lya=0,2-lza
             lxb = 0
             !jco=coset(0,0,2)
             jco=10

              DO lxa=la_min_local,2

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,2)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,2)
             jco=10

              DO lxa=MAX(la_min_local-0-1,0),1

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,2)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,2)
             jco=10

                lxa = 0
                !ico=coset(0,2,0)
                ico=8

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,0,2)



          !ENDDO
          !ENDDO
          !DO lyb=0,2-lzb
          !DO lya=0,2-lza
             lxb = 0
             !jco=coset(0,0,2)
             jco=10

              DO lxa=MAX(la_min_local-1-0,0),1

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,2)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,2)
             jco=10

                lxa = 0
                !ico=coset(0,1,1)
                ico=9

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,1,2)



          !ENDDO
          !ENDDO
          !DO lyb=0,2-lzb
          !DO lya=0,2-lza
             lxb = 0
             !jco=coset(0,0,2)
             jco=10

                lxa = 0
                !ico=coset(0,0,2)
                ico=10

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,2,2)



          !ENDDO
          !ENDDO

   ENDDO

    END SUBROUTINE xyz_to_vab_2_2

! *****************************************************************************
    SUBROUTINE xyz_to_vab_2_3

    INTEGER :: my_rank,ierror

    REAL(kind=dp) :: alpha(0:(2+3),0:2,0:3,3)
    REAL(kind=dp) :: coef_ttz(0:2,0:3)
    REAL(kind=dp) :: coef_tyz(0:2,0:3,0:2,0:3)

    coef_xyz=coef_xyz*prefactor

!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max + lb_max
!   (current implementation is l**7)
!
!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***

!**** RR : Values of la_max and lb_max
!    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)
!    WRITE(100+my_rank,*) "la=",la_max_local,"lb=",lb_max_local
!****
    alpha(:,:,:,:)=0.0_dp

    DO lxa=0,2
    DO lxb=0,3
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,1)=alpha(lxa-l+lxb-k,lxa,lxb,1)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(1)-(ra(1)+rab(1)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(1)+rp(1))
       ENDDO
    ENDDO
    ENDDO
    DO lxa=0,2
    DO lxb=0,3
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,2)=alpha(lxa-l+lxb-k,lxa,lxb,2)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(2)-(ra(2)+rab(2)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(2)+rp(2))
       ENDDO
    ENDDO
    ENDDO
    DO lxa=0,2
    DO lxb=0,3
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,3)=alpha(lxa-l+lxb-k,lxa,lxb,3)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(3)-(ra(3)+rab(3)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(3)+rp(3))
       ENDDO
    ENDDO
    ENDDO

    !
    !   compute v_{lxa,lya,lza,lxb,lyb,lzb} given v_{lxp,lyp,lzp} and alpha(ls,lxa,lxb,1)
    !   use a three step procedure
    !

    lxyz=0
    DO lzp=0,5
       coef_tyz=0.0_dp
       DO lyp=0,5-lzp
          coef_ttz=0.0_dp
          DO lxp=0,5-lzp-lyp
             lxyz=lxyz+1
             DO lxb=0,3
             DO lxa=0,2
                coef_ttz(lxa,lxb)=coef_ttz(lxa,lxb)+coef_xyz(lxyz)*alpha(lxp,lxa,lxb,1)
             ENDDO
             ENDDO

          ENDDO ! lxp

          DO lyb=0,3
          DO lya=0,2
             DO lxb=0,3-lyb
             DO lxa=0,2-lya
                coef_tyz(lxa,lxb,lya,lyb)=coef_tyz(lxa,lxb,lya,lyb)+coef_ttz(lxa,lxb)*alpha(lyp,lya,lyb,2)
             ENDDO
             ENDDO
          ENDDO
          ENDDO

       ENDDO !lyp

       DO lzb=0,3
       DO lza=0,2
          DO lyb=0,3-lzb
          DO lya=0,2-lza
             DO lxb=MAX(lb_min_local-lzb-lyb,0),3-lzb-lyb
             jco=coset(lxb,lyb,lzb)
             DO lxa=MAX(la_min_local-lza-lya,0),2-lza-lya
                ico=coset(lxa,lya,lza)
                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,lya,lyb)*alpha(lzp,lza,lzb,3)
             ENDDO
             ENDDO
          ENDDO
          ENDDO
       ENDDO
       ENDDO

    ENDDO

    END SUBROUTINE xyz_to_vab_2_3! *****************************************************************************
    SUBROUTINE xyz_to_vab_2_4

    INTEGER :: my_rank,ierror

    REAL(kind=dp) :: alpha(0:(2+4),0:2,0:4,3)
    REAL(kind=dp) :: coef_ttz(0:2,0:4)
    REAL(kind=dp) :: coef_tyz(0:2,0:4,0:2,0:4)

    coef_xyz=coef_xyz*prefactor

!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max + lb_max
!   (current implementation is l**7)
!
!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***

!**** RR : Values of la_max and lb_max
!    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)
!    WRITE(100+my_rank,*) "la=",la_max_local,"lb=",lb_max_local
!****
    alpha(:,:,:,:)=0.0_dp

    DO lxa=0,2
    DO lxb=0,4
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,1)=alpha(lxa-l+lxb-k,lxa,lxb,1)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(1)-(ra(1)+rab(1)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(1)+rp(1))
       ENDDO
    ENDDO
    ENDDO
    DO lxa=0,2
    DO lxb=0,4
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,2)=alpha(lxa-l+lxb-k,lxa,lxb,2)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(2)-(ra(2)+rab(2)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(2)+rp(2))
       ENDDO
    ENDDO
    ENDDO
    DO lxa=0,2
    DO lxb=0,4
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,3)=alpha(lxa-l+lxb-k,lxa,lxb,3)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(3)-(ra(3)+rab(3)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(3)+rp(3))
       ENDDO
    ENDDO
    ENDDO

    !
    !   compute v_{lxa,lya,lza,lxb,lyb,lzb} given v_{lxp,lyp,lzp} and alpha(ls,lxa,lxb,1)
    !   use a three step procedure
    !

    lxyz=0
    DO lzp=0,6
       coef_tyz=0.0_dp
       DO lyp=0,6-lzp
          coef_ttz=0.0_dp
          DO lxp=0,6-lzp-lyp
             lxyz=lxyz+1
             DO lxb=0,4
             DO lxa=0,2
                coef_ttz(lxa,lxb)=coef_ttz(lxa,lxb)+coef_xyz(lxyz)*alpha(lxp,lxa,lxb,1)
             ENDDO
             ENDDO

          ENDDO ! lxp

          DO lyb=0,4
          DO lya=0,2
             DO lxb=0,4-lyb
             DO lxa=0,2-lya
                coef_tyz(lxa,lxb,lya,lyb)=coef_tyz(lxa,lxb,lya,lyb)+coef_ttz(lxa,lxb)*alpha(lyp,lya,lyb,2)
             ENDDO
             ENDDO
          ENDDO
          ENDDO

       ENDDO !lyp

       DO lzb=0,4
       DO lza=0,2
          DO lyb=0,4-lzb
          DO lya=0,2-lza
             DO lxb=MAX(lb_min_local-lzb-lyb,0),4-lzb-lyb
             jco=coset(lxb,lyb,lzb)
             DO lxa=MAX(la_min_local-lza-lya,0),2-lza-lya
                ico=coset(lxa,lya,lza)
                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,lya,lyb)*alpha(lzp,lza,lzb,3)
             ENDDO
             ENDDO
          ENDDO
          ENDDO
       ENDDO
       ENDDO

    ENDDO

    END SUBROUTINE xyz_to_vab_2_4! *****************************************************************************
    SUBROUTINE xyz_to_vab_3_0

    !INTEGER :: my_rank,ierror
    !REAL(kind=dp) :: alpha(0:(3+0),0:3,0:0,3)

    INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp
    REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb

    REAL(kind=dp) :: alpha_1(0:(3+0),0:3,0:0)
    REAL(kind=dp) :: alpha_2(0:(3+0),0:3,0:0)
    REAL(kind=dp) :: alpha_3(0:(3+0),0:3,0:0)
    REAL(kind=dp) :: coef_ttz(0:3,0:0)
    REAL(kind=dp) :: coef_tyz(0:3,0:0,0:3,0:0)
    coef_xyz=coef_xyz*prefactor
!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max_local + lb_max_local
!   (current implementation is l**7)
!
!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***
  !RR    DO iaxis=1,3
  alpha_1(:,:,:)=0.0_dp
  alpha_2(:,:,:)=0.0_dp
  alpha_3(:,:,:)=0.0_dp

         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,0,0)=alpha_1(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,1,0)=alpha_1(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,1,0)=alpha_1(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,2,0)=alpha_1(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,2,0)=alpha_1(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,2,0)=alpha_1(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(3,3,0)=alpha_1(3,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(3,dp)
             a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,3,0)=alpha_1(2,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,3,0)=alpha_1(1,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(3,dp)
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,3,0)=alpha_1(0,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,0,0)=alpha_2(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,1,0)=alpha_2(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,1,0)=alpha_2(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,2,0)=alpha_2(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,2,0)=alpha_2(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,2,0)=alpha_2(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(3,3,0)=alpha_2(3,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(3,dp)
             a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,3,0)=alpha_2(2,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,3,0)=alpha_2(1,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(3,dp)
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,3,0)=alpha_2(0,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,0,0)=alpha_3(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,1,0)=alpha_3(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,1,0)=alpha_3(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,2,0)=alpha_3(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,2,0)=alpha_3(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,2,0)=alpha_3(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(3,3,0)=alpha_3(3,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(3,dp)
             a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,3,0)=alpha_3(2,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,3,0)=alpha_3(1,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(3,dp)
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,3,0)=alpha_3(0,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
  !RR    ENDDO
    !
    !   compute v_{lxa,lya,lza,lxb,lyb,lzb} given v_{lxp,lyp,lzp} and alpha(ls,lxa,lxb,1)
    !   use a three step procedure
    !

    lxyz=0
    DO lzp=0,3+0
       coef_tyz=0.0_dp
       DO lyp=0,3+0-lzp
          coef_ttz=0.0_dp
          DO lxp=0,3+0-lzp-lyp
             lxyz=lxyz+1
                coef_ttz(0,0)=coef_ttz(0,0)+coef_xyz(lxyz)*alpha_1(lxp,0,0)
                coef_ttz(1,0)=coef_ttz(1,0)+coef_xyz(lxyz)*alpha_1(lxp,1,0)
                coef_ttz(2,0)=coef_ttz(2,0)+coef_xyz(lxyz)*alpha_1(lxp,2,0)
                coef_ttz(3,0)=coef_ttz(3,0)+coef_xyz(lxyz)*alpha_1(lxp,3,0)
          ENDDO

     ! Improve locality
     !temp(:,:) = alpha(lyp,0:3,0:0,2)
     ! RR do lyb = 0, lb_max_local
     ! RR do lya = 0, la_max_local
             !DO lxb=0,0-lyb
             !DO lxa=0,3-lya
                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
                coef_tyz(1,0,0,0)=coef_tyz(1,0,0,0)+coef_ttz(1,0)*alpha_2(lyp,0,0)
                coef_tyz(2,0,0,0)=coef_tyz(2,0,0,0)+coef_ttz(2,0)*alpha_2(lyp,0,0)
                coef_tyz(3,0,0,0)=coef_tyz(3,0,0,0)+coef_ttz(3,0)*alpha_2(lyp,0,0)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,0-lyb
             !DO lxa=0,3-lya
                coef_tyz(0,0,1,0)=coef_tyz(0,0,1,0)+coef_ttz(0,0)*alpha_2(lyp,1,0)
                coef_tyz(1,0,1,0)=coef_tyz(1,0,1,0)+coef_ttz(1,0)*alpha_2(lyp,1,0)
                coef_tyz(2,0,1,0)=coef_tyz(2,0,1,0)+coef_ttz(2,0)*alpha_2(lyp,1,0)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,0-lyb
             !DO lxa=0,3-lya
                coef_tyz(0,0,2,0)=coef_tyz(0,0,2,0)+coef_ttz(0,0)*alpha_2(lyp,2,0)
                coef_tyz(1,0,2,0)=coef_tyz(1,0,2,0)+coef_ttz(1,0)*alpha_2(lyp,2,0)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,0-lyb
             !DO lxa=0,3-lya
                coef_tyz(0,0,3,0)=coef_tyz(0,0,3,0)+coef_ttz(0,0)*alpha_2(lyp,3,0)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
       ENDDO


          !DO lyb=0,0-lzb
          !DO lya=0,3-lza
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=la_min_local,3

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-0-1,0),2

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-0-2,0),1

                ico=coset(lxa,2,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,0,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,3,0)
                ico=17

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,3,0)*alpha_3(lzp,0,0)



          !ENDDO
          !ENDDO
          !DO lyb=0,0-lzb
          !DO lya=0,3-lza
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-1-0,0),2

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-1-1,0),1

                ico=coset(lxa,1,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,1,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,2,1)
                ico=18

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,1,0)



          !ENDDO
          !ENDDO
          !DO lyb=0,0-lzb
          !DO lya=0,3-lza
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-2-0,0),1

                ico=coset(lxa,0,2)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,2,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,1,2)
                ico=19

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,2,0)



          !ENDDO
          !ENDDO
          !DO lyb=0,0-lzb
          !DO lya=0,3-lza
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,0,3)
                ico=20

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,3,0)



          !ENDDO
          !ENDDO

   ENDDO

    END SUBROUTINE xyz_to_vab_3_0

! *****************************************************************************
    SUBROUTINE xyz_to_vab_3_1

    !INTEGER :: my_rank,ierror
    !REAL(kind=dp) :: alpha(0:(3+1),0:3,0:1,3)

    INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp
    REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb

    REAL(kind=dp) :: alpha_1(0:(3+1),0:3,0:1)
    REAL(kind=dp) :: alpha_2(0:(3+1),0:3,0:1)
    REAL(kind=dp) :: alpha_3(0:(3+1),0:3,0:1)
    REAL(kind=dp) :: coef_ttz(0:3,0:1)
    REAL(kind=dp) :: coef_tyz(0:3,0:1,0:3,0:1)
    coef_xyz=coef_xyz*prefactor
!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max_local + lb_max_local
!   (current implementation is l**7)
!
!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***
  !RR    DO iaxis=1,3
  alpha_1(:,:,:)=0.0_dp
  alpha_2(:,:,:)=0.0_dp
  alpha_3(:,:,:)=0.0_dp

         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,0,0)=alpha_1(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,0,1)=alpha_1(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,0,1)=alpha_1(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,1,0)=alpha_1(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,1,0)=alpha_1(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,1,1)=alpha_1(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,1,1)=alpha_1(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,1,1)=alpha_1(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,1,1)=alpha_1(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,2,0)=alpha_1(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,2,0)=alpha_1(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,2,0)=alpha_1(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(3,2,1)=alpha_1(3,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(2,2,1)=alpha_1(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,2,1)=alpha_1(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,2,1)=alpha_1(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,2,1)=alpha_1(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,2,1)=alpha_1(0,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(3,3,0)=alpha_1(3,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(3,dp)
             a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,3,0)=alpha_1(2,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,3,0)=alpha_1(1,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(3,dp)
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,3,0)=alpha_1(0,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(4,3,1)=alpha_1(4,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(3,3,1)=alpha_1(3,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(3,dp)
             a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(3,3,1)=alpha_1(3,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(2,3,1)=alpha_1(2,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,3,1)=alpha_1(2,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(1,3,1)=alpha_1(1,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(3,dp)
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,3,1)=alpha_1(1,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(1)-(ra(1)+rab(1)))
          !ENDDO
          ! DO l=0,lxb
             alpha_1(0,3,1)=alpha_1(0,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,0,0)=alpha_2(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,0,1)=alpha_2(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,0,1)=alpha_2(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,1,0)=alpha_2(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,1,0)=alpha_2(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,1,1)=alpha_2(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,1,1)=alpha_2(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,1,1)=alpha_2(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,1,1)=alpha_2(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,2,0)=alpha_2(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,2,0)=alpha_2(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,2,0)=alpha_2(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(3,2,1)=alpha_2(3,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(2,2,1)=alpha_2(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,2,1)=alpha_2(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,2,1)=alpha_2(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,2,1)=alpha_2(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,2,1)=alpha_2(0,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(3,3,0)=alpha_2(3,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(3,dp)
             a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,3,0)=alpha_2(2,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,3,0)=alpha_2(1,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(3,dp)
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,3,0)=alpha_2(0,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(4,3,1)=alpha_2(4,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(3,3,1)=alpha_2(3,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(3,dp)
             a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(3,3,1)=alpha_2(3,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(2,3,1)=alpha_2(2,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,3,1)=alpha_2(2,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(1,3,1)=alpha_2(1,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(3,dp)
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,3,1)=alpha_2(1,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(2)-(ra(2)+rab(2)))
          !ENDDO
          ! DO l=0,lxb
             alpha_2(0,3,1)=alpha_2(0,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,0,0)=alpha_3(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,0,1)=alpha_3(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,0,1)=alpha_3(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,1,0)=alpha_3(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,1,0)=alpha_3(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,1,1)=alpha_3(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,1,1)=alpha_3(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,1,1)=alpha_3(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,1,1)=alpha_3(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,2,0)=alpha_3(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,2,0)=alpha_3(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,2,0)=alpha_3(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(3,2,1)=alpha_3(3,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(2,2,1)=alpha_3(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,2,1)=alpha_3(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,2,1)=alpha_3(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,2,1)=alpha_3(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,2,1)=alpha_3(0,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(3,3,0)=alpha_3(3,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(3,dp)
             a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,3,0)=alpha_3(2,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,3,0)=alpha_3(1,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(3,dp)
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,3,0)=alpha_3(0,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(4,3,1)=alpha_3(4,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(3,3,1)=alpha_3(3,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(3,dp)
             a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(3,3,1)=alpha_3(3,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(2,3,1)=alpha_3(2,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,3,1)=alpha_3(2,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(1,3,1)=alpha_3(1,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(3,dp)
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,3,1)=alpha_3(1,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*(rp(3)-(ra(3)+rab(3)))
          !ENDDO
          ! DO l=0,lxb
             alpha_3(0,3,1)=alpha_3(0,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
  !RR    ENDDO
    !
    !   compute v_{lxa,lya,lza,lxb,lyb,lzb} given v_{lxp,lyp,lzp} and alpha(ls,lxa,lxb,1)
    !   use a three step procedure
    !

    lxyz=0
    DO lzp=0,3+1
       coef_tyz=0.0_dp
       DO lyp=0,3+1-lzp
          coef_ttz=0.0_dp
          DO lxp=0,3+1-lzp-lyp
             lxyz=lxyz+1
                coef_ttz(0,0)=coef_ttz(0,0)+coef_xyz(lxyz)*alpha_1(lxp,0,0)
                coef_ttz(1,0)=coef_ttz(1,0)+coef_xyz(lxyz)*alpha_1(lxp,1,0)
                coef_ttz(2,0)=coef_ttz(2,0)+coef_xyz(lxyz)*alpha_1(lxp,2,0)
                coef_ttz(3,0)=coef_ttz(3,0)+coef_xyz(lxyz)*alpha_1(lxp,3,0)
                coef_ttz(0,1)=coef_ttz(0,1)+coef_xyz(lxyz)*alpha_1(lxp,0,1)
                coef_ttz(1,1)=coef_ttz(1,1)+coef_xyz(lxyz)*alpha_1(lxp,1,1)
                coef_ttz(2,1)=coef_ttz(2,1)+coef_xyz(lxyz)*alpha_1(lxp,2,1)
                coef_ttz(3,1)=coef_ttz(3,1)+coef_xyz(lxyz)*alpha_1(lxp,3,1)
          ENDDO

     ! Improve locality
     !temp(:,:) = alpha(lyp,0:3,0:1,2)
     ! RR do lyb = 0, lb_max_local
     ! RR do lya = 0, la_max_local
             !DO lxb=0,1-lyb
             !DO lxa=0,3-lya
                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
                coef_tyz(1,0,0,0)=coef_tyz(1,0,0,0)+coef_ttz(1,0)*alpha_2(lyp,0,0)
                coef_tyz(2,0,0,0)=coef_tyz(2,0,0,0)+coef_ttz(2,0)*alpha_2(lyp,0,0)
                coef_tyz(3,0,0,0)=coef_tyz(3,0,0,0)+coef_ttz(3,0)*alpha_2(lyp,0,0)
             !ENDDO
             !DO lxa=0,3-lya
                coef_tyz(0,1,0,0)=coef_tyz(0,1,0,0)+coef_ttz(0,1)*alpha_2(lyp,0,0)
                coef_tyz(1,1,0,0)=coef_tyz(1,1,0,0)+coef_ttz(1,1)*alpha_2(lyp,0,0)
                coef_tyz(2,1,0,0)=coef_tyz(2,1,0,0)+coef_ttz(2,1)*alpha_2(lyp,0,0)
                coef_tyz(3,1,0,0)=coef_tyz(3,1,0,0)+coef_ttz(3,1)*alpha_2(lyp,0,0)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,1-lyb
             !DO lxa=0,3-lya
                coef_tyz(0,0,1,0)=coef_tyz(0,0,1,0)+coef_ttz(0,0)*alpha_2(lyp,1,0)
                coef_tyz(1,0,1,0)=coef_tyz(1,0,1,0)+coef_ttz(1,0)*alpha_2(lyp,1,0)
                coef_tyz(2,0,1,0)=coef_tyz(2,0,1,0)+coef_ttz(2,0)*alpha_2(lyp,1,0)
             !ENDDO
             !DO lxa=0,3-lya
                coef_tyz(0,1,1,0)=coef_tyz(0,1,1,0)+coef_ttz(0,1)*alpha_2(lyp,1,0)
                coef_tyz(1,1,1,0)=coef_tyz(1,1,1,0)+coef_ttz(1,1)*alpha_2(lyp,1,0)
                coef_tyz(2,1,1,0)=coef_tyz(2,1,1,0)+coef_ttz(2,1)*alpha_2(lyp,1,0)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,1-lyb
             !DO lxa=0,3-lya
                coef_tyz(0,0,2,0)=coef_tyz(0,0,2,0)+coef_ttz(0,0)*alpha_2(lyp,2,0)
                coef_tyz(1,0,2,0)=coef_tyz(1,0,2,0)+coef_ttz(1,0)*alpha_2(lyp,2,0)
             !ENDDO
             !DO lxa=0,3-lya
                coef_tyz(0,1,2,0)=coef_tyz(0,1,2,0)+coef_ttz(0,1)*alpha_2(lyp,2,0)
                coef_tyz(1,1,2,0)=coef_tyz(1,1,2,0)+coef_ttz(1,1)*alpha_2(lyp,2,0)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,1-lyb
             !DO lxa=0,3-lya
                coef_tyz(0,0,3,0)=coef_tyz(0,0,3,0)+coef_ttz(0,0)*alpha_2(lyp,3,0)
             !ENDDO
             !DO lxa=0,3-lya
                coef_tyz(0,1,3,0)=coef_tyz(0,1,3,0)+coef_ttz(0,1)*alpha_2(lyp,3,0)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
     ! RR do lya = 0, la_max_local
             !DO lxb=0,1-lyb
             !DO lxa=0,3-lya
                coef_tyz(0,0,0,1)=coef_tyz(0,0,0,1)+coef_ttz(0,0)*alpha_2(lyp,0,1)
                coef_tyz(1,0,0,1)=coef_tyz(1,0,0,1)+coef_ttz(1,0)*alpha_2(lyp,0,1)
                coef_tyz(2,0,0,1)=coef_tyz(2,0,0,1)+coef_ttz(2,0)*alpha_2(lyp,0,1)
                coef_tyz(3,0,0,1)=coef_tyz(3,0,0,1)+coef_ttz(3,0)*alpha_2(lyp,0,1)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,1-lyb
             !DO lxa=0,3-lya
                coef_tyz(0,0,1,1)=coef_tyz(0,0,1,1)+coef_ttz(0,0)*alpha_2(lyp,1,1)
                coef_tyz(1,0,1,1)=coef_tyz(1,0,1,1)+coef_ttz(1,0)*alpha_2(lyp,1,1)
                coef_tyz(2,0,1,1)=coef_tyz(2,0,1,1)+coef_ttz(2,0)*alpha_2(lyp,1,1)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,1-lyb
             !DO lxa=0,3-lya
                coef_tyz(0,0,2,1)=coef_tyz(0,0,2,1)+coef_ttz(0,0)*alpha_2(lyp,2,1)
                coef_tyz(1,0,2,1)=coef_tyz(1,0,2,1)+coef_ttz(1,0)*alpha_2(lyp,2,1)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,1-lyb
             !DO lxa=0,3-lya
                coef_tyz(0,0,3,1)=coef_tyz(0,0,3,1)+coef_ttz(0,0)*alpha_2(lyp,3,1)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
       ENDDO


          !DO lyb=0,1-lzb
          !DO lya=0,3-lza
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

              DO lxa=la_min_local,3

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO

          !ENDDO
          !ENDDO
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

              DO lxa=MAX(la_min_local-0-1,0),2

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO

          !ENDDO
          !ENDDO
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

              DO lxa=MAX(la_min_local-0-2,0),1

                ico=coset(lxa,2,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO

          !ENDDO
          !ENDDO
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,3,0)
                ico=17

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,3,0)*alpha_3(lzp,0,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,3-lza
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

              DO lxa=la_min_local,3

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

              DO lxa=MAX(la_min_local-0-1,0),2

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,0,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

              DO lxa=MAX(la_min_local-0-2,0),1

                ico=coset(lxa,2,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,1)*alpha_3(lzp,0,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

                lxa = 0
                !ico=coset(0,3,0)
                ico=17

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,3,1)*alpha_3(lzp,0,0)



          !ENDDO
          !ENDDO
          !DO lyb=0,1-lzb
          !DO lya=0,3-lza
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

              DO lxa=MAX(la_min_local-1-0,0),2

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,0)

             ENDDO

             ENDDO

          !ENDDO
          !ENDDO
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

              DO lxa=MAX(la_min_local-1-1,0),1

                ico=coset(lxa,1,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,1,0)

             ENDDO

             ENDDO

          !ENDDO
          !ENDDO
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,2,1)
                ico=18

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,1,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,3-lza
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

              DO lxa=MAX(la_min_local-1-0,0),2

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,1,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

              DO lxa=MAX(la_min_local-1-1,0),1

                ico=coset(lxa,1,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,1,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

                lxa = 0
                !ico=coset(0,2,1)
                ico=18

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,1)*alpha_3(lzp,1,0)



          !ENDDO
          !ENDDO
          !DO lyb=0,1-lzb
          !DO lya=0,3-lza
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

              DO lxa=MAX(la_min_local-2-0,0),1

                ico=coset(lxa,0,2)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,2,0)

             ENDDO

             ENDDO

          !ENDDO
          !ENDDO
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,1,2)
                ico=19

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,2,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,3-lza
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

              DO lxa=MAX(la_min_local-2-0,0),1

                ico=coset(lxa,0,2)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,2,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

                lxa = 0
                !ico=coset(0,1,2)
                ico=19

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,2,0)



          !ENDDO
          !ENDDO
          !DO lyb=0,1-lzb
          !DO lya=0,3-lza
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,0,3)
                ico=20

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,3,0)


             ENDDO

          !ENDDO
          !ENDDO
          !DO lya=0,3-lza
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

                lxa = 0
                !ico=coset(0,0,3)
                ico=20

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,3,0)



          !ENDDO
          !ENDDO
          !DO lyb=0,1-lzb
          !DO lya=0,3-lza
             lxb = 0
             !jco=coset(0,0,1)
             jco=4

              DO lxa=la_min_local,3

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,1)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,1)
             jco=4

              DO lxa=MAX(la_min_local-0-1,0),2

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,1)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,1)
             jco=4

              DO lxa=MAX(la_min_local-0-2,0),1

                ico=coset(lxa,2,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,0,1)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,1)
             jco=4

                lxa = 0
                !ico=coset(0,3,0)
                ico=17

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,3,0)*alpha_3(lzp,0,1)



          !ENDDO
          !ENDDO
          !DO lyb=0,1-lzb
          !DO lya=0,3-lza
             lxb = 0
             !jco=coset(0,0,1)
             jco=4

              DO lxa=MAX(la_min_local-1-0,0),2

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,1)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,1)
             jco=4

              DO lxa=MAX(la_min_local-1-1,0),1

                ico=coset(lxa,1,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,1,1)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,1)
             jco=4

                lxa = 0
                !ico=coset(0,2,1)
                ico=18

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,1,1)



          !ENDDO
          !ENDDO
          !DO lyb=0,1-lzb
          !DO lya=0,3-lza
             lxb = 0
             !jco=coset(0,0,1)
             jco=4

              DO lxa=MAX(la_min_local-2-0,0),1

                ico=coset(lxa,0,2)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,2,1)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,1)
             jco=4

                lxa = 0
                !ico=coset(0,1,2)
                ico=19

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,2,1)



          !ENDDO
          !ENDDO
          !DO lyb=0,1-lzb
          !DO lya=0,3-lza
             lxb = 0
             !jco=coset(0,0,1)
             jco=4

                lxa = 0
                !ico=coset(0,0,3)
                ico=20

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,3,1)



          !ENDDO
          !ENDDO

   ENDDO

    END SUBROUTINE xyz_to_vab_3_1

! *****************************************************************************
    SUBROUTINE xyz_to_vab_3_2

    INTEGER :: my_rank,ierror

    REAL(kind=dp) :: alpha(0:(3+2),0:3,0:2,3)
    REAL(kind=dp) :: coef_ttz(0:3,0:2)
    REAL(kind=dp) :: coef_tyz(0:3,0:2,0:3,0:2)

    coef_xyz=coef_xyz*prefactor

!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max + lb_max
!   (current implementation is l**7)
!
!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***

!**** RR : Values of la_max and lb_max
!    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)
!    WRITE(100+my_rank,*) "la=",la_max_local,"lb=",lb_max_local
!****
    alpha(:,:,:,:)=0.0_dp

    DO lxa=0,3
    DO lxb=0,2
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,1)=alpha(lxa-l+lxb-k,lxa,lxb,1)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(1)-(ra(1)+rab(1)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(1)+rp(1))
       ENDDO
    ENDDO
    ENDDO
    DO lxa=0,3
    DO lxb=0,2
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,2)=alpha(lxa-l+lxb-k,lxa,lxb,2)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(2)-(ra(2)+rab(2)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(2)+rp(2))
       ENDDO
    ENDDO
    ENDDO
    DO lxa=0,3
    DO lxb=0,2
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,3)=alpha(lxa-l+lxb-k,lxa,lxb,3)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(3)-(ra(3)+rab(3)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(3)+rp(3))
       ENDDO
    ENDDO
    ENDDO

    !
    !   compute v_{lxa,lya,lza,lxb,lyb,lzb} given v_{lxp,lyp,lzp} and alpha(ls,lxa,lxb,1)
    !   use a three step procedure
    !

    lxyz=0
    DO lzp=0,5
       coef_tyz=0.0_dp
       DO lyp=0,5-lzp
          coef_ttz=0.0_dp
          DO lxp=0,5-lzp-lyp
             lxyz=lxyz+1
             DO lxb=0,2
             DO lxa=0,3
                coef_ttz(lxa,lxb)=coef_ttz(lxa,lxb)+coef_xyz(lxyz)*alpha(lxp,lxa,lxb,1)
             ENDDO
             ENDDO

          ENDDO ! lxp

          DO lyb=0,2
          DO lya=0,3
             DO lxb=0,2-lyb
             DO lxa=0,3-lya
                coef_tyz(lxa,lxb,lya,lyb)=coef_tyz(lxa,lxb,lya,lyb)+coef_ttz(lxa,lxb)*alpha(lyp,lya,lyb,2)
             ENDDO
             ENDDO
          ENDDO
          ENDDO

       ENDDO !lyp

       DO lzb=0,2
       DO lza=0,3
          DO lyb=0,2-lzb
          DO lya=0,3-lza
             DO lxb=MAX(lb_min_local-lzb-lyb,0),2-lzb-lyb
             jco=coset(lxb,lyb,lzb)
             DO lxa=MAX(la_min_local-lza-lya,0),3-lza-lya
                ico=coset(lxa,lya,lza)
                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,lya,lyb)*alpha(lzp,lza,lzb,3)
             ENDDO
             ENDDO
          ENDDO
          ENDDO
       ENDDO
       ENDDO

    ENDDO

    END SUBROUTINE xyz_to_vab_3_2! *****************************************************************************
    SUBROUTINE xyz_to_vab_3_3

    INTEGER :: my_rank,ierror

    REAL(kind=dp) :: alpha(0:(3+3),0:3,0:3,3)
    REAL(kind=dp) :: coef_ttz(0:3,0:3)
    REAL(kind=dp) :: coef_tyz(0:3,0:3,0:3,0:3)

    coef_xyz=coef_xyz*prefactor

!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max + lb_max
!   (current implementation is l**7)
!
!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***

!**** RR : Values of la_max and lb_max
!    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)
!    WRITE(100+my_rank,*) "la=",la_max_local,"lb=",lb_max_local
!****
    alpha(:,:,:,:)=0.0_dp

    DO lxa=0,3
    DO lxb=0,3
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,1)=alpha(lxa-l+lxb-k,lxa,lxb,1)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(1)-(ra(1)+rab(1)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(1)+rp(1))
       ENDDO
    ENDDO
    ENDDO
    DO lxa=0,3
    DO lxb=0,3
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,2)=alpha(lxa-l+lxb-k,lxa,lxb,2)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(2)-(ra(2)+rab(2)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(2)+rp(2))
       ENDDO
    ENDDO
    ENDDO
    DO lxa=0,3
    DO lxb=0,3
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,3)=alpha(lxa-l+lxb-k,lxa,lxb,3)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(3)-(ra(3)+rab(3)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(3)+rp(3))
       ENDDO
    ENDDO
    ENDDO

    !
    !   compute v_{lxa,lya,lza,lxb,lyb,lzb} given v_{lxp,lyp,lzp} and alpha(ls,lxa,lxb,1)
    !   use a three step procedure
    !

    lxyz=0
    DO lzp=0,6
       coef_tyz=0.0_dp
       DO lyp=0,6-lzp
          coef_ttz=0.0_dp
          DO lxp=0,6-lzp-lyp
             lxyz=lxyz+1
             DO lxb=0,3
             DO lxa=0,3
                coef_ttz(lxa,lxb)=coef_ttz(lxa,lxb)+coef_xyz(lxyz)*alpha(lxp,lxa,lxb,1)
             ENDDO
             ENDDO

          ENDDO ! lxp

          DO lyb=0,3
          DO lya=0,3
             DO lxb=0,3-lyb
             DO lxa=0,3-lya
                coef_tyz(lxa,lxb,lya,lyb)=coef_tyz(lxa,lxb,lya,lyb)+coef_ttz(lxa,lxb)*alpha(lyp,lya,lyb,2)
             ENDDO
             ENDDO
          ENDDO
          ENDDO

       ENDDO !lyp

       DO lzb=0,3
       DO lza=0,3
          DO lyb=0,3-lzb
          DO lya=0,3-lza
             DO lxb=MAX(lb_min_local-lzb-lyb,0),3-lzb-lyb
             jco=coset(lxb,lyb,lzb)
             DO lxa=MAX(la_min_local-lza-lya,0),3-lza-lya
                ico=coset(lxa,lya,lza)
                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,lya,lyb)*alpha(lzp,lza,lzb,3)
             ENDDO
             ENDDO
          ENDDO
          ENDDO
       ENDDO
       ENDDO

    ENDDO

    END SUBROUTINE xyz_to_vab_3_3! *****************************************************************************
    SUBROUTINE xyz_to_vab_3_4

    INTEGER :: my_rank,ierror

    REAL(kind=dp) :: alpha(0:(3+4),0:3,0:4,3)
    REAL(kind=dp) :: coef_ttz(0:3,0:4)
    REAL(kind=dp) :: coef_tyz(0:3,0:4,0:3,0:4)

    coef_xyz=coef_xyz*prefactor

!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max + lb_max
!   (current implementation is l**7)
!
!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***

!**** RR : Values of la_max and lb_max
!    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)
!    WRITE(100+my_rank,*) "la=",la_max_local,"lb=",lb_max_local
!****
    alpha(:,:,:,:)=0.0_dp

    DO lxa=0,3
    DO lxb=0,4
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,1)=alpha(lxa-l+lxb-k,lxa,lxb,1)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(1)-(ra(1)+rab(1)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(1)+rp(1))
       ENDDO
    ENDDO
    ENDDO
    DO lxa=0,3
    DO lxb=0,4
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,2)=alpha(lxa-l+lxb-k,lxa,lxb,2)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(2)-(ra(2)+rab(2)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(2)+rp(2))
       ENDDO
    ENDDO
    ENDDO
    DO lxa=0,3
    DO lxb=0,4
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,3)=alpha(lxa-l+lxb-k,lxa,lxb,3)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(3)-(ra(3)+rab(3)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(3)+rp(3))
       ENDDO
    ENDDO
    ENDDO

    !
    !   compute v_{lxa,lya,lza,lxb,lyb,lzb} given v_{lxp,lyp,lzp} and alpha(ls,lxa,lxb,1)
    !   use a three step procedure
    !

    lxyz=0
    DO lzp=0,7
       coef_tyz=0.0_dp
       DO lyp=0,7-lzp
          coef_ttz=0.0_dp
          DO lxp=0,7-lzp-lyp
             lxyz=lxyz+1
             DO lxb=0,4
             DO lxa=0,3
                coef_ttz(lxa,lxb)=coef_ttz(lxa,lxb)+coef_xyz(lxyz)*alpha(lxp,lxa,lxb,1)
             ENDDO
             ENDDO

          ENDDO ! lxp

          DO lyb=0,4
          DO lya=0,3
             DO lxb=0,4-lyb
             DO lxa=0,3-lya
                coef_tyz(lxa,lxb,lya,lyb)=coef_tyz(lxa,lxb,lya,lyb)+coef_ttz(lxa,lxb)*alpha(lyp,lya,lyb,2)
             ENDDO
             ENDDO
          ENDDO
          ENDDO

       ENDDO !lyp

       DO lzb=0,4
       DO lza=0,3
          DO lyb=0,4-lzb
          DO lya=0,3-lza
             DO lxb=MAX(lb_min_local-lzb-lyb,0),4-lzb-lyb
             jco=coset(lxb,lyb,lzb)
             DO lxa=MAX(la_min_local-lza-lya,0),3-lza-lya
                ico=coset(lxa,lya,lza)
                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,lya,lyb)*alpha(lzp,lza,lzb,3)
             ENDDO
             ENDDO
          ENDDO
          ENDDO
       ENDDO
       ENDDO

    ENDDO

    END SUBROUTINE xyz_to_vab_3_4! *****************************************************************************
    SUBROUTINE xyz_to_vab_4_0

    !INTEGER :: my_rank,ierror
    !REAL(kind=dp) :: alpha(0:(4+0),0:4,0:0,3)

    INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp
    REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb

    REAL(kind=dp) :: alpha_1(0:(4+0),0:4,0:0)
    REAL(kind=dp) :: alpha_2(0:(4+0),0:4,0:0)
    REAL(kind=dp) :: alpha_3(0:(4+0),0:4,0:0)
    REAL(kind=dp) :: coef_ttz(0:4,0:0)
    REAL(kind=dp) :: coef_tyz(0:4,0:0,0:4,0:0)
    coef_xyz=coef_xyz*prefactor
!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max_local + lb_max_local
!   (current implementation is l**7)
!
!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***
  !RR    DO iaxis=1,3
  alpha_1(:,:,:)=0.0_dp
  alpha_2(:,:,:)=0.0_dp
  alpha_3(:,:,:)=0.0_dp

         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,0,0)=alpha_1(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,1,0)=alpha_1(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,1,0)=alpha_1(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,2,0)=alpha_1(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,2,0)=alpha_1(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,2,0)=alpha_1(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(3,3,0)=alpha_1(3,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(3,dp)
             a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,3,0)=alpha_1(2,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,3,0)=alpha_1(1,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(3,dp)
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,3,0)=alpha_1(0,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(4,4,0)=alpha_1(4,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(4,dp)
             a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(3,4,0)=alpha_1(3,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(3,dp)/REAL(2,dp)
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(2,4,0)=alpha_1(2,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(2,dp)/REAL(3,dp)
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(1,4,0)=alpha_1(1,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(4,dp)
          a=a*(-ra(1)+rp(1))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_1(0,4,0)=alpha_1(0,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,0,0)=alpha_2(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,1,0)=alpha_2(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,1,0)=alpha_2(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,2,0)=alpha_2(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,2,0)=alpha_2(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,2,0)=alpha_2(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(3,3,0)=alpha_2(3,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(3,dp)
             a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,3,0)=alpha_2(2,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,3,0)=alpha_2(1,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(3,dp)
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,3,0)=alpha_2(0,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(4,4,0)=alpha_2(4,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(4,dp)
             a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(3,4,0)=alpha_2(3,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(3,dp)/REAL(2,dp)
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(2,4,0)=alpha_2(2,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(2,dp)/REAL(3,dp)
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(1,4,0)=alpha_2(1,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(4,dp)
          a=a*(-ra(2)+rp(2))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_2(0,4,0)=alpha_2(0,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,0,0)=alpha_3(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,1,0)=alpha_3(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,1,0)=alpha_3(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,2,0)=alpha_3(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,2,0)=alpha_3(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,2,0)=alpha_3(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(3,3,0)=alpha_3(3,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(3,dp)
             a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,3,0)=alpha_3(2,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,3,0)=alpha_3(1,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(3,dp)
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,3,0)=alpha_3(0,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
         binomial_k_lxa=1.0_dp
         a=1.0_dp
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(4,4,0)=alpha_3(4,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             binomial_k_lxa=binomial_k_lxa*REAL(4,dp)
             a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(3,4,0)=alpha_3(3,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(3,dp)/REAL(2,dp)
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(2,4,0)=alpha_3(2,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(2,dp)/REAL(3,dp)
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(1,4,0)=alpha_3(1,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
          binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(4,dp)
          a=a*(-ra(3)+rp(3))
         ! ENDDO
         ! DO k=0,lxa
          binomial_l_lxb=1.0_dp
          b=1.0_dp
          ! DO l=0,lxb
             alpha_3(0,4,0)=alpha_3(0,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
          !ENDDO
             ! Removed  here as lxa-k == 0
         ! ENDDO
  !RR    ENDDO
    !
    !   compute v_{lxa,lya,lza,lxb,lyb,lzb} given v_{lxp,lyp,lzp} and alpha(ls,lxa,lxb,1)
    !   use a three step procedure
    !

    lxyz=0
    DO lzp=0,4+0
       coef_tyz=0.0_dp
       DO lyp=0,4+0-lzp
          coef_ttz=0.0_dp
          DO lxp=0,4+0-lzp-lyp
             lxyz=lxyz+1
                coef_ttz(0,0)=coef_ttz(0,0)+coef_xyz(lxyz)*alpha_1(lxp,0,0)
                coef_ttz(1,0)=coef_ttz(1,0)+coef_xyz(lxyz)*alpha_1(lxp,1,0)
                coef_ttz(2,0)=coef_ttz(2,0)+coef_xyz(lxyz)*alpha_1(lxp,2,0)
                coef_ttz(3,0)=coef_ttz(3,0)+coef_xyz(lxyz)*alpha_1(lxp,3,0)
                coef_ttz(4,0)=coef_ttz(4,0)+coef_xyz(lxyz)*alpha_1(lxp,4,0)
          ENDDO

     ! Improve locality
     !temp(:,:) = alpha(lyp,0:4,0:0,2)
     ! RR do lyb = 0, lb_max_local
     ! RR do lya = 0, la_max_local
             !DO lxb=0,0-lyb
             !DO lxa=0,4-lya
                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
                coef_tyz(1,0,0,0)=coef_tyz(1,0,0,0)+coef_ttz(1,0)*alpha_2(lyp,0,0)
                coef_tyz(2,0,0,0)=coef_tyz(2,0,0,0)+coef_ttz(2,0)*alpha_2(lyp,0,0)
                coef_tyz(3,0,0,0)=coef_tyz(3,0,0,0)+coef_ttz(3,0)*alpha_2(lyp,0,0)
                coef_tyz(4,0,0,0)=coef_tyz(4,0,0,0)+coef_ttz(4,0)*alpha_2(lyp,0,0)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,0-lyb
             !DO lxa=0,4-lya
                coef_tyz(0,0,1,0)=coef_tyz(0,0,1,0)+coef_ttz(0,0)*alpha_2(lyp,1,0)
                coef_tyz(1,0,1,0)=coef_tyz(1,0,1,0)+coef_ttz(1,0)*alpha_2(lyp,1,0)
                coef_tyz(2,0,1,0)=coef_tyz(2,0,1,0)+coef_ttz(2,0)*alpha_2(lyp,1,0)
                coef_tyz(3,0,1,0)=coef_tyz(3,0,1,0)+coef_ttz(3,0)*alpha_2(lyp,1,0)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,0-lyb
             !DO lxa=0,4-lya
                coef_tyz(0,0,2,0)=coef_tyz(0,0,2,0)+coef_ttz(0,0)*alpha_2(lyp,2,0)
                coef_tyz(1,0,2,0)=coef_tyz(1,0,2,0)+coef_ttz(1,0)*alpha_2(lyp,2,0)
                coef_tyz(2,0,2,0)=coef_tyz(2,0,2,0)+coef_ttz(2,0)*alpha_2(lyp,2,0)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,0-lyb
             !DO lxa=0,4-lya
                coef_tyz(0,0,3,0)=coef_tyz(0,0,3,0)+coef_ttz(0,0)*alpha_2(lyp,3,0)
                coef_tyz(1,0,3,0)=coef_tyz(1,0,3,0)+coef_ttz(1,0)*alpha_2(lyp,3,0)
             !ENDDO
             !ENDDO
  ! RR enddo
             !DO lxb=0,0-lyb
             !DO lxa=0,4-lya
                coef_tyz(0,0,4,0)=coef_tyz(0,0,4,0)+coef_ttz(0,0)*alpha_2(lyp,4,0)
             !ENDDO
             !ENDDO
  ! RR enddo

  ! RR enddo
       ENDDO


          !DO lyb=0,0-lzb
          !DO lya=0,4-lza
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=la_min_local,4

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-0-1,0),3

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-0-2,0),2

                ico=coset(lxa,2,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,0,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-0-3,0),1

                ico=coset(lxa,3,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,3,0)*alpha_3(lzp,0,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,4,0)
                ico=31

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,4,0)*alpha_3(lzp,0,0)



          !ENDDO
          !ENDDO
          !DO lyb=0,0-lzb
          !DO lya=0,4-lza
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-1-0,0),3

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-1-1,0),2

                ico=coset(lxa,1,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,1,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-1-2,0),1

                ico=coset(lxa,2,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,1,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,3,1)
                ico=32

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,3,0)*alpha_3(lzp,1,0)



          !ENDDO
          !ENDDO
          !DO lyb=0,0-lzb
          !DO lya=0,4-lza
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-2-0,0),2

                ico=coset(lxa,0,2)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,2,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-2-1,0),1

                ico=coset(lxa,1,2)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,2,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,2,2)
                ico=33

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,2,0)



          !ENDDO
          !ENDDO
          !DO lyb=0,0-lzb
          !DO lya=0,4-lza
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-3-0,0),1

                ico=coset(lxa,0,3)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,3,0)

             ENDDO


          !ENDDO
          !ENDDO
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,1,3)
                ico=34

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,3,0)



          !ENDDO
          !ENDDO
          !DO lyb=0,0-lzb
          !DO lya=0,4-lza
             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,0,4)
                ico=35

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,4,0)



          !ENDDO
          !ENDDO

   ENDDO

    END SUBROUTINE xyz_to_vab_4_0

! *****************************************************************************
    SUBROUTINE xyz_to_vab_4_1

    INTEGER :: my_rank,ierror

    REAL(kind=dp) :: alpha(0:(4+1),0:4,0:1,3)
    REAL(kind=dp) :: coef_ttz(0:4,0:1)
    REAL(kind=dp) :: coef_tyz(0:4,0:1,0:4,0:1)

    coef_xyz=coef_xyz*prefactor

!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max + lb_max
!   (current implementation is l**7)
!
!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***

!**** RR : Values of la_max and lb_max
!    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)
!    WRITE(100+my_rank,*) "la=",la_max_local,"lb=",lb_max_local
!****
    alpha(:,:,:,:)=0.0_dp

    DO lxa=0,4
    DO lxb=0,1
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,1)=alpha(lxa-l+lxb-k,lxa,lxb,1)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(1)-(ra(1)+rab(1)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(1)+rp(1))
       ENDDO
    ENDDO
    ENDDO
    DO lxa=0,4
    DO lxb=0,1
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,2)=alpha(lxa-l+lxb-k,lxa,lxb,2)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(2)-(ra(2)+rab(2)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(2)+rp(2))
       ENDDO
    ENDDO
    ENDDO
    DO lxa=0,4
    DO lxb=0,1
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,3)=alpha(lxa-l+lxb-k,lxa,lxb,3)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(3)-(ra(3)+rab(3)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(3)+rp(3))
       ENDDO
    ENDDO
    ENDDO

    !
    !   compute v_{lxa,lya,lza,lxb,lyb,lzb} given v_{lxp,lyp,lzp} and alpha(ls,lxa,lxb,1)
    !   use a three step procedure
    !

    lxyz=0
    DO lzp=0,5
       coef_tyz=0.0_dp
       DO lyp=0,5-lzp
          coef_ttz=0.0_dp
          DO lxp=0,5-lzp-lyp
             lxyz=lxyz+1
             DO lxb=0,1
             DO lxa=0,4
                coef_ttz(lxa,lxb)=coef_ttz(lxa,lxb)+coef_xyz(lxyz)*alpha(lxp,lxa,lxb,1)
             ENDDO
             ENDDO

          ENDDO ! lxp

          DO lyb=0,1
          DO lya=0,4
             DO lxb=0,1-lyb
             DO lxa=0,4-lya
                coef_tyz(lxa,lxb,lya,lyb)=coef_tyz(lxa,lxb,lya,lyb)+coef_ttz(lxa,lxb)*alpha(lyp,lya,lyb,2)
             ENDDO
             ENDDO
          ENDDO
          ENDDO

       ENDDO !lyp

       DO lzb=0,1
       DO lza=0,4
          DO lyb=0,1-lzb
          DO lya=0,4-lza
             DO lxb=MAX(lb_min_local-lzb-lyb,0),1-lzb-lyb
             jco=coset(lxb,lyb,lzb)
             DO lxa=MAX(la_min_local-lza-lya,0),4-lza-lya
                ico=coset(lxa,lya,lza)
                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,lya,lyb)*alpha(lzp,lza,lzb,3)
             ENDDO
             ENDDO
          ENDDO
          ENDDO
       ENDDO
       ENDDO

    ENDDO

    END SUBROUTINE xyz_to_vab_4_1! *****************************************************************************
    SUBROUTINE xyz_to_vab_4_2

    INTEGER :: my_rank,ierror

    REAL(kind=dp) :: alpha(0:(4+2),0:4,0:2,3)
    REAL(kind=dp) :: coef_ttz(0:4,0:2)
    REAL(kind=dp) :: coef_tyz(0:4,0:2,0:4,0:2)

    coef_xyz=coef_xyz*prefactor

!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max + lb_max
!   (current implementation is l**7)
!
!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***

!**** RR : Values of la_max and lb_max
!    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)
!    WRITE(100+my_rank,*) "la=",la_max_local,"lb=",lb_max_local
!****
    alpha(:,:,:,:)=0.0_dp

    DO lxa=0,4
    DO lxb=0,2
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,1)=alpha(lxa-l+lxb-k,lxa,lxb,1)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(1)-(ra(1)+rab(1)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(1)+rp(1))
       ENDDO
    ENDDO
    ENDDO
    DO lxa=0,4
    DO lxb=0,2
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,2)=alpha(lxa-l+lxb-k,lxa,lxb,2)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(2)-(ra(2)+rab(2)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(2)+rp(2))
       ENDDO
    ENDDO
    ENDDO
    DO lxa=0,4
    DO lxb=0,2
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,3)=alpha(lxa-l+lxb-k,lxa,lxb,3)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(3)-(ra(3)+rab(3)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(3)+rp(3))
       ENDDO
    ENDDO
    ENDDO

    !
    !   compute v_{lxa,lya,lza,lxb,lyb,lzb} given v_{lxp,lyp,lzp} and alpha(ls,lxa,lxb,1)
    !   use a three step procedure
    !

    lxyz=0
    DO lzp=0,6
       coef_tyz=0.0_dp
       DO lyp=0,6-lzp
          coef_ttz=0.0_dp
          DO lxp=0,6-lzp-lyp
             lxyz=lxyz+1
             DO lxb=0,2
             DO lxa=0,4
                coef_ttz(lxa,lxb)=coef_ttz(lxa,lxb)+coef_xyz(lxyz)*alpha(lxp,lxa,lxb,1)
             ENDDO
             ENDDO

          ENDDO ! lxp

          DO lyb=0,2
          DO lya=0,4
             DO lxb=0,2-lyb
             DO lxa=0,4-lya
                coef_tyz(lxa,lxb,lya,lyb)=coef_tyz(lxa,lxb,lya,lyb)+coef_ttz(lxa,lxb)*alpha(lyp,lya,lyb,2)
             ENDDO
             ENDDO
          ENDDO
          ENDDO

       ENDDO !lyp

       DO lzb=0,2
       DO lza=0,4
          DO lyb=0,2-lzb
          DO lya=0,4-lza
             DO lxb=MAX(lb_min_local-lzb-lyb,0),2-lzb-lyb
             jco=coset(lxb,lyb,lzb)
             DO lxa=MAX(la_min_local-lza-lya,0),4-lza-lya
                ico=coset(lxa,lya,lza)
                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,lya,lyb)*alpha(lzp,lza,lzb,3)
             ENDDO
             ENDDO
          ENDDO
          ENDDO
       ENDDO
       ENDDO

    ENDDO

    END SUBROUTINE xyz_to_vab_4_2! *****************************************************************************
    SUBROUTINE xyz_to_vab_4_3

    INTEGER :: my_rank,ierror

    REAL(kind=dp) :: alpha(0:(4+3),0:4,0:3,3)
    REAL(kind=dp) :: coef_ttz(0:4,0:3)
    REAL(kind=dp) :: coef_tyz(0:4,0:3,0:4,0:3)

    coef_xyz=coef_xyz*prefactor

!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max + lb_max
!   (current implementation is l**7)
!
!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***

!**** RR : Values of la_max and lb_max
!    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)
!    WRITE(100+my_rank,*) "la=",la_max_local,"lb=",lb_max_local
!****
    alpha(:,:,:,:)=0.0_dp

    DO lxa=0,4
    DO lxb=0,3
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,1)=alpha(lxa-l+lxb-k,lxa,lxb,1)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(1)-(ra(1)+rab(1)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(1)+rp(1))
       ENDDO
    ENDDO
    ENDDO
    DO lxa=0,4
    DO lxb=0,3
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,2)=alpha(lxa-l+lxb-k,lxa,lxb,2)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(2)-(ra(2)+rab(2)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(2)+rp(2))
       ENDDO
    ENDDO
    ENDDO
    DO lxa=0,4
    DO lxb=0,3
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,3)=alpha(lxa-l+lxb-k,lxa,lxb,3)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(3)-(ra(3)+rab(3)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(3)+rp(3))
       ENDDO
    ENDDO
    ENDDO

    !
    !   compute v_{lxa,lya,lza,lxb,lyb,lzb} given v_{lxp,lyp,lzp} and alpha(ls,lxa,lxb,1)
    !   use a three step procedure
    !

    lxyz=0
    DO lzp=0,7
       coef_tyz=0.0_dp
       DO lyp=0,7-lzp
          coef_ttz=0.0_dp
          DO lxp=0,7-lzp-lyp
             lxyz=lxyz+1
             DO lxb=0,3
             DO lxa=0,4
                coef_ttz(lxa,lxb)=coef_ttz(lxa,lxb)+coef_xyz(lxyz)*alpha(lxp,lxa,lxb,1)
             ENDDO
             ENDDO

          ENDDO ! lxp

          DO lyb=0,3
          DO lya=0,4
             DO lxb=0,3-lyb
             DO lxa=0,4-lya
                coef_tyz(lxa,lxb,lya,lyb)=coef_tyz(lxa,lxb,lya,lyb)+coef_ttz(lxa,lxb)*alpha(lyp,lya,lyb,2)
             ENDDO
             ENDDO
          ENDDO
          ENDDO

       ENDDO !lyp

       DO lzb=0,3
       DO lza=0,4
          DO lyb=0,3-lzb
          DO lya=0,4-lza
             DO lxb=MAX(lb_min_local-lzb-lyb,0),3-lzb-lyb
             jco=coset(lxb,lyb,lzb)
             DO lxa=MAX(la_min_local-lza-lya,0),4-lza-lya
                ico=coset(lxa,lya,lza)
                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,lya,lyb)*alpha(lzp,lza,lzb,3)
             ENDDO
             ENDDO
          ENDDO
          ENDDO
       ENDDO
       ENDDO

    ENDDO

    END SUBROUTINE xyz_to_vab_4_3! *****************************************************************************
    SUBROUTINE xyz_to_vab_4_4

    INTEGER :: my_rank,ierror

    REAL(kind=dp) :: alpha(0:(4+4),0:4,0:4,3)
    REAL(kind=dp) :: coef_ttz(0:4,0:4)
    REAL(kind=dp) :: coef_tyz(0:4,0:4,0:4,0:4)

    coef_xyz=coef_xyz*prefactor

!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max + lb_max
!   (current implementation is l**7)
!
!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***

!**** RR : Values of la_max and lb_max
!    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)
!    WRITE(100+my_rank,*) "la=",la_max_local,"lb=",lb_max_local
!****
    alpha(:,:,:,:)=0.0_dp

    DO lxa=0,4
    DO lxb=0,4
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,1)=alpha(lxa-l+lxb-k,lxa,lxb,1)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(1)-(ra(1)+rab(1)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(1)+rp(1))
       ENDDO
    ENDDO
    ENDDO
    DO lxa=0,4
    DO lxb=0,4
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,2)=alpha(lxa-l+lxb-k,lxa,lxb,2)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(2)-(ra(2)+rab(2)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(2)+rp(2))
       ENDDO
    ENDDO
    ENDDO
    DO lxa=0,4
    DO lxb=0,4
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,3)=alpha(lxa-l+lxb-k,lxa,lxb,3)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(3)-(ra(3)+rab(3)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(3)+rp(3))
       ENDDO
    ENDDO
    ENDDO

    !
    !   compute v_{lxa,lya,lza,lxb,lyb,lzb} given v_{lxp,lyp,lzp} and alpha(ls,lxa,lxb,1)
    !   use a three step procedure
    !

    lxyz=0
    DO lzp=0,8
       coef_tyz=0.0_dp
       DO lyp=0,8-lzp
          coef_ttz=0.0_dp
          DO lxp=0,8-lzp-lyp
             lxyz=lxyz+1
             DO lxb=0,4
             DO lxa=0,4
                coef_ttz(lxa,lxb)=coef_ttz(lxa,lxb)+coef_xyz(lxyz)*alpha(lxp,lxa,lxb,1)
             ENDDO
             ENDDO

          ENDDO ! lxp

          DO lyb=0,4
          DO lya=0,4
             DO lxb=0,4-lyb
             DO lxa=0,4-lya
                coef_tyz(lxa,lxb,lya,lyb)=coef_tyz(lxa,lxb,lya,lyb)+coef_ttz(lxa,lxb)*alpha(lyp,lya,lyb,2)
             ENDDO
             ENDDO
          ENDDO
          ENDDO

       ENDDO !lyp

       DO lzb=0,4
       DO lza=0,4
          DO lyb=0,4-lzb
          DO lya=0,4-lza
             DO lxb=MAX(lb_min_local-lzb-lyb,0),4-lzb-lyb
             jco=coset(lxb,lyb,lzb)
             DO lxa=MAX(la_min_local-lza-lya,0),4-lza-lya
                ico=coset(lxa,lya,lza)
                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,lya,lyb)*alpha(lzp,lza,lzb,3)
             ENDDO
             ENDDO
          ENDDO
          ENDDO
       ENDDO
       ENDDO

    ENDDO

    END SUBROUTINE xyz_to_vab_4_4
	    SUBROUTINE call_to_xyz_to_vab


            SELECT CASE(la_max_local)
                 CASE (0)
                      SELECT CASE(lb_max_local)
                               CASE(0) 
                                    CALL xyz_to_vab_0_0
                               CASE(1) 
                                    CALL xyz_to_vab_0_1
                               CASE(2) 
                                    CALL xyz_to_vab_0_2
                               CASE(3) 
                                    CALL xyz_to_vab_0_3
                               CASE(4) 
                                    CALL xyz_to_vab_0_4
                               CASE DEFAULT
                                    CALL xyz_to_vab 
                       END SELECT
                 CASE (1)
                      SELECT CASE(lb_max_local)
                               CASE(0) 
                                    CALL xyz_to_vab_1_0
                               CASE(1) 
                                    CALL xyz_to_vab_1_1
                               CASE(2) 
                                    CALL xyz_to_vab_1_2
                               CASE(3) 
                                    CALL xyz_to_vab_1_3
                               CASE(4) 
                                    CALL xyz_to_vab_1_4
                               CASE DEFAULT
                                    CALL xyz_to_vab 
                       END SELECT
                 CASE (2)
                      SELECT CASE(lb_max_local)
                               CASE(0) 
                                    CALL xyz_to_vab_2_0
                               CASE(1) 
                                    CALL xyz_to_vab_2_1
                               CASE(2) 
                                    CALL xyz_to_vab_2_2
                               CASE(3) 
                                    CALL xyz_to_vab_2_3
                               CASE(4) 
                                    CALL xyz_to_vab_2_4
                               CASE DEFAULT
                                    CALL xyz_to_vab 
                       END SELECT
                 CASE (3)
                      SELECT CASE(lb_max_local)
                               CASE(0) 
                                    CALL xyz_to_vab_3_0
                               CASE(1) 
                                    CALL xyz_to_vab_3_1
                               CASE(2) 
                                    CALL xyz_to_vab_3_2
                               CASE(3) 
                                    CALL xyz_to_vab_3_3
                               CASE(4) 
                                    CALL xyz_to_vab_3_4
                               CASE DEFAULT
                                    CALL xyz_to_vab 
                       END SELECT
                 CASE (4)
                      SELECT CASE(lb_max_local)
                               CASE(0) 
                                    CALL xyz_to_vab_4_0
                               CASE(1) 
                                    CALL xyz_to_vab_4_1
                               CASE(2) 
                                    CALL xyz_to_vab_4_2
                               CASE(3) 
                                    CALL xyz_to_vab_4_3
                               CASE(4) 
                                    CALL xyz_to_vab_4_4
                               CASE DEFAULT
                                    CALL xyz_to_vab 
                       END SELECT
                 CASE DEFAULT
                    CALL xyz_to_vab 
             END SELECT
	   

	    END SUBROUTINE call_to_xyz_to_vab
	
