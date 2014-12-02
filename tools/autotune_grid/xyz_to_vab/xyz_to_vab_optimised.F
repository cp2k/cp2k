! *****************************************************************************
!> \\brief provides specific implementations for common cases of xyz_to_vab routine.
!> 
!> \\note
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
!>      The script to regenerate this file can be found at autotune_grid/xyz_to_vab/
!> \par History
!>      07.2012 Fixed a bug with the dimension of coef_xyz
!>      06.2012 Template cleanup 
!>      05.2012 created
!>
!> \\author Ruyman Reyes
! *****************************************************************************

! Caller routine
      SUBROUTINE call_to_xyz_to_vab(prefactor, coef_xyz, lp, &
                                    la_max_local, lb_max_local, &
                                    rp, ra, rab, vab, coset, & 
                                    la_min_local, lb_min_local, &
                                    maxl, lvab, hvab)
      
              USE kinds,                        ONLY: dp
      
             REAL(KIND=dp), INTENT(IN) :: prefactor
             REAL(KIND=dp), DIMENSION(((lp+1)*(lp+2)*(lp+3))/6), &
                                  INTENT(INOUT) :: coef_xyz
             INTEGER, INTENT(IN) :: lp, la_max_local, lb_max_local,&
                                    la_min_local, lb_min_local
             INTEGER, INTENT(IN) :: maxl, lvab, hvab
             REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: rp, ra, rab
             REAL(KIND=dp), DIMENSION(lvab,hvab), INTENT(INOUT)  :: vab
             INTEGER,  DIMENSION(-1:maxl,-1:maxl,-1:maxl), &
                       INTENT(IN)     :: coset

      
          SELECT CASE(la_max_local)
               CASE (0)
                    SELECT CASE(lb_max_local)
                             CASE(0) 
                               CALL xyz_to_vab_0_0(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, & 
                                                  lb_min_local, maxl, lvab, hvab)
  
                             CASE(1) 
                               CALL xyz_to_vab_0_1(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, & 
                                                  lb_min_local, maxl, lvab, hvab)
  
                             CASE(2) 
                               CALL xyz_to_vab_0_2(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, & 
                                                  lb_min_local, maxl, lvab, hvab)
  
                             CASE(3) 
                               CALL xyz_to_vab_0_3(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, & 
                                                  lb_min_local, maxl, lvab, hvab)
  
                             CASE(4) 
                               CALL xyz_to_vab_0_4(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, & 
                                                  lb_min_local, maxl, lvab, hvab)
  
                             CASE DEFAULT
  			          CALL xyz_to_vab(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, &
                                                  lb_min_local, maxl, lvab, hvab)
                     END SELECT
               CASE (1)
                    SELECT CASE(lb_max_local)
                             CASE(0) 
                               CALL xyz_to_vab_1_0(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, & 
                                                  lb_min_local, maxl, lvab, hvab)
  
                             CASE(1) 
                               CALL xyz_to_vab_1_1(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, & 
                                                  lb_min_local, maxl, lvab, hvab)
  
                             CASE(2) 
                               CALL xyz_to_vab_1_2(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, & 
                                                  lb_min_local, maxl, lvab, hvab)
  
                             CASE(3) 
                               CALL xyz_to_vab_1_3(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, & 
                                                  lb_min_local, maxl, lvab, hvab)
  
                             CASE(4) 
                               CALL xyz_to_vab_1_4(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, & 
                                                  lb_min_local, maxl, lvab, hvab)
  
                             CASE DEFAULT
  			          CALL xyz_to_vab(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, &
                                                  lb_min_local, maxl, lvab, hvab)
                     END SELECT
               CASE (2)
                    SELECT CASE(lb_max_local)
                             CASE(0) 
                               CALL xyz_to_vab_2_0(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, & 
                                                  lb_min_local, maxl, lvab, hvab)
  
                             CASE(1) 
                               CALL xyz_to_vab_2_1(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, & 
                                                  lb_min_local, maxl, lvab, hvab)
  
                             CASE(2) 
                               CALL xyz_to_vab_2_2(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, & 
                                                  lb_min_local, maxl, lvab, hvab)
  
                             CASE(3) 
                               CALL xyz_to_vab_2_3(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, & 
                                                  lb_min_local, maxl, lvab, hvab)
  
                             CASE(4) 
                               CALL xyz_to_vab_2_4(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, & 
                                                  lb_min_local, maxl, lvab, hvab)
  
                             CASE DEFAULT
  			          CALL xyz_to_vab(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, &
                                                  lb_min_local, maxl, lvab, hvab)
                     END SELECT
               CASE (3)
                    SELECT CASE(lb_max_local)
                             CASE(0) 
                               CALL xyz_to_vab_3_0(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, & 
                                                  lb_min_local, maxl, lvab, hvab)
  
                             CASE(1) 
                               CALL xyz_to_vab_3_1(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, & 
                                                  lb_min_local, maxl, lvab, hvab)
  
                             CASE(2) 
                               CALL xyz_to_vab_3_2(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, & 
                                                  lb_min_local, maxl, lvab, hvab)
  
                             CASE(3) 
                               CALL xyz_to_vab_3_3(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, & 
                                                  lb_min_local, maxl, lvab, hvab)
  
                             CASE(4) 
                               CALL xyz_to_vab_3_4(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, & 
                                                  lb_min_local, maxl, lvab, hvab)
  
                             CASE DEFAULT
  			          CALL xyz_to_vab(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, &
                                                  lb_min_local, maxl, lvab, hvab)
                     END SELECT
               CASE (4)
                    SELECT CASE(lb_max_local)
                             CASE(0) 
                               CALL xyz_to_vab_4_0(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, & 
                                                  lb_min_local, maxl, lvab, hvab)
  
                             CASE(1) 
                               CALL xyz_to_vab_4_1(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, & 
                                                  lb_min_local, maxl, lvab, hvab)
  
                             CASE(2) 
                               CALL xyz_to_vab_4_2(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, & 
                                                  lb_min_local, maxl, lvab, hvab)
  
                             CASE(3) 
                               CALL xyz_to_vab_4_3(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, & 
                                                  lb_min_local, maxl, lvab, hvab)
  
                             CASE(4) 
                               CALL xyz_to_vab_4_4(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, & 
                                                  lb_min_local, maxl, lvab, hvab)
  
                             CASE DEFAULT
  			          CALL xyz_to_vab(prefactor, coef_xyz, lp, &
                                                  la_max_local, lb_max_local, &
                                                  rp, ra, rab, vab, coset, la_min_local, &
                                                  lb_min_local, maxl, lvab, hvab)
                     END SELECT
               CASE DEFAULT
  			CALL xyz_to_vab(prefactor, coef_xyz, lp, &
  					 la_max_local, lb_max_local, &
  					  rp, ra, rab, vab, coset, la_min_local, lb_min_local, maxl, lvab, hvab)
  
           END SELECT
  	   
      
      	    END SUBROUTINE call_to_xyz_to_vab
      
!!! Cloned routines
          ! Fully unrolled version
    ! *****************************************************************************
    SUBROUTINE xyz_to_vab_0_0 (prefactor, coef_xyz, & 
                           lp, la_max_local, lb_max_local, rp, ra, rab, vab, coset, &
                           la_min_local, lb_min_local, maxl, lvab, hvab)
       

       USE kinds,                        ONLY: dp

       ! PARAMETERS
       REAL(KIND=dp), INTENT(IN) :: prefactor
       INTEGER, INTENT(IN) :: lp, la_max_local, lb_max_local, &
                              lb_min_local, la_min_local, maxl, &
                              lvab,hvab

       INTEGER,  DIMENSION(-1:maxl,-1:maxl,-1:maxl), &
                        INTENT(IN)     :: coset

       REAL(KIND=dp),  DIMENSION(((lp+1)*(lp+2)*(lp+3))/6), &
                        INTENT(INOUT) :: coef_xyz
       REAL(KIND=dp),  DIMENSION(lvab,hvab), INTENT(INOUT) :: vab
       REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: rp, ra, rab
   
       ! LOCAL
       INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp, iaxis, l, k, lya,lyb,lza,lzb
       REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb,const_rp,const_ra
   

    REAL(KIND=dp) :: alpha_1(0:(0+0),0:0,0:0)
    REAL(KIND=dp) :: alpha_2(0:(0+0),0:0,0:0)
    REAL(KIND=dp) :: alpha_3(0:(0+0),0:0,0:0)
    coef_xyz=coef_xyz*prefactor


             jco=1 ! coset(0,0,0)
                ico=1 ! coset(0,0,0)
                vab(ico,jco)=vab(ico,jco)+coef_xyz(1)

    END SUBROUTINE xyz_to_vab_0_0



          ! Fully unrolled version
    ! *****************************************************************************
    SUBROUTINE xyz_to_vab_0_1 (prefactor, coef_xyz, & 
                           lp, la_max_local, lb_max_local, rp, ra, rab, vab, coset, &
                           la_min_local, lb_min_local, maxl, lvab, hvab)
       

       USE kinds,                        ONLY: dp

       ! PARAMETERS
       REAL(KIND=dp), INTENT(IN) :: prefactor
       INTEGER, INTENT(IN) :: lp, la_max_local, lb_max_local, &
                              lb_min_local, la_min_local, maxl, &
                              lvab,hvab

       INTEGER,  DIMENSION(-1:maxl,-1:maxl,-1:maxl), &
                        INTENT(IN)     :: coset

       REAL(KIND=dp),  DIMENSION(((lp+1)*(lp+2)*(lp+3))/6), &
                        INTENT(INOUT) :: coef_xyz
       REAL(KIND=dp),  DIMENSION(lvab,hvab), INTENT(INOUT) :: vab
       REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: rp, ra, rab
   
       ! LOCAL
       INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp, iaxis, l, k, lya,lyb,lza,lzb
       REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb,const_rp,const_ra
   

    REAL(KIND=dp) :: alpha_1(0:(0+1),0:0,0:1)
    REAL(KIND=dp) :: alpha_2(0:(0+1),0:0,0:1)
    REAL(KIND=dp) :: alpha_3(0:(0+1),0:0,0:1)
    REAL(KIND=dp) :: coef_ttz(0:0,0:1)
    REAL(KIND=dp) :: coef_tyz(0:0,0:1,0:0,0:1)
    coef_xyz=coef_xyz*prefactor

  alpha_1(:,:,:)=0.0_dp
  alpha_2(:,:,:)=0.0_dp
  alpha_3(:,:,:)=0.0_dp

         const_rp = rp(1)-(ra(1)+rab(1))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,0,0)=alpha_1(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,0,1)=alpha_1(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(0,0,1)=alpha_1(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         const_rp = rp(2)-(ra(2)+rab(2))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,0,0)=alpha_2(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,0,1)=alpha_2(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(0,0,1)=alpha_2(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         const_rp = rp(3)-(ra(3)+rab(3))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,0,0)=alpha_3(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,0,1)=alpha_3(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(0,0,1)=alpha_3(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0

    lxyz=0
       coef_tyz=0.0_dp
       DO lyp=0,0+1-0
          coef_ttz=0.0_dp
          DO lxp=0,0+1-0-lyp
             lxyz=lxyz+1
                coef_ttz(0,0)=coef_ttz(0,0)+coef_xyz(lxyz)*alpha_1(lxp,0,0)
                coef_ttz(0,1)=coef_ttz(0,1)+coef_xyz(lxyz)*alpha_1(lxp,0,1)
          ENDDO

                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
                coef_tyz(0,1,0,0)=coef_tyz(0,1,0,0)+coef_ttz(0,1)*alpha_2(lyp,0,0)

                coef_tyz(0,0,0,1)=coef_tyz(0,0,0,1)+coef_ttz(0,0)*alpha_2(lyp,0,1)

       ENDDO


              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(0,0,0)


             ENDDO
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(0,0,0)


             lxb = 0
             !jco=coset(0,0,1)
             jco=4

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(0,0,1)



       coef_tyz=0.0_dp
       DO lyp=0,0+1-1
          coef_ttz=0.0_dp
          DO lxp=0,0+1-1-lyp
             lxyz=lxyz+1
                coef_ttz(0,0)=coef_ttz(0,0)+coef_xyz(lxyz)*alpha_1(lxp,0,0)
                coef_ttz(0,1)=coef_ttz(0,1)+coef_xyz(lxyz)*alpha_1(lxp,0,1)
          ENDDO

                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
                coef_tyz(0,1,0,0)=coef_tyz(0,1,0,0)+coef_ttz(0,1)*alpha_2(lyp,0,0)

                coef_tyz(0,0,0,1)=coef_tyz(0,0,0,1)+coef_ttz(0,0)*alpha_2(lyp,0,1)

       ENDDO


              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(1,0,0)


             ENDDO
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(1,0,0)


             lxb = 0
             !jco=coset(0,0,1)
             jco=4

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(1,0,1)




    END SUBROUTINE xyz_to_vab_0_1



          ! Fully unrolled version
    ! *****************************************************************************
    SUBROUTINE xyz_to_vab_0_2 (prefactor, coef_xyz, & 
                           lp, la_max_local, lb_max_local, rp, ra, rab, vab, coset, &
                           la_min_local, lb_min_local, maxl, lvab, hvab)
       

       USE kinds,                        ONLY: dp

       ! PARAMETERS
       REAL(KIND=dp), INTENT(IN) :: prefactor
       INTEGER, INTENT(IN) :: lp, la_max_local, lb_max_local, &
                              lb_min_local, la_min_local, maxl, &
                              lvab,hvab

       INTEGER,  DIMENSION(-1:maxl,-1:maxl,-1:maxl), &
                        INTENT(IN)     :: coset

       REAL(KIND=dp),  DIMENSION(((lp+1)*(lp+2)*(lp+3))/6), &
                        INTENT(INOUT) :: coef_xyz
       REAL(KIND=dp),  DIMENSION(lvab,hvab), INTENT(INOUT) :: vab
       REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: rp, ra, rab
   
       ! LOCAL
       INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp, iaxis, l, k, lya,lyb,lza,lzb
       REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb,const_rp,const_ra
   

    REAL(KIND=dp) :: alpha_1(0:(0+2),0:0,0:2)
    REAL(KIND=dp) :: alpha_2(0:(0+2),0:0,0:2)
    REAL(KIND=dp) :: alpha_3(0:(0+2),0:0,0:2)
    REAL(KIND=dp) :: coef_ttz(0:0,0:2)
    REAL(KIND=dp) :: coef_tyz(0:0,0:2,0:0,0:2)
    coef_xyz=coef_xyz*prefactor

  alpha_1(:,:,:)=0.0_dp
  alpha_2(:,:,:)=0.0_dp
  alpha_3(:,:,:)=0.0_dp

         const_rp = rp(1)-(ra(1)+rab(1))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,0,0)=alpha_1(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,0,1)=alpha_1(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(0,0,1)=alpha_1(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,0,2)=alpha_1(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_1(1,0,2)=alpha_1(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_1(0,0,2)=alpha_1(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         const_rp = rp(2)-(ra(2)+rab(2))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,0,0)=alpha_2(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,0,1)=alpha_2(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(0,0,1)=alpha_2(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,0,2)=alpha_2(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_2(1,0,2)=alpha_2(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_2(0,0,2)=alpha_2(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         const_rp = rp(3)-(ra(3)+rab(3))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,0,0)=alpha_3(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,0,1)=alpha_3(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(0,0,1)=alpha_3(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,0,2)=alpha_3(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_3(1,0,2)=alpha_3(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_3(0,0,2)=alpha_3(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0

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

                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
                coef_tyz(0,1,0,0)=coef_tyz(0,1,0,0)+coef_ttz(0,1)*alpha_2(lyp,0,0)
                coef_tyz(0,2,0,0)=coef_tyz(0,2,0,0)+coef_ttz(0,2)*alpha_2(lyp,0,0)
                coef_tyz(0,0,0,1)=coef_tyz(0,0,0,1)+coef_ttz(0,0)*alpha_2(lyp,0,1)
                coef_tyz(0,1,0,1)=coef_tyz(0,1,0,1)+coef_ttz(0,1)*alpha_2(lyp,0,1)
                coef_tyz(0,0,0,2)=coef_tyz(0,0,0,2)+coef_ttz(0,0)*alpha_2(lyp,0,2)
       ENDDO


              DO lxb=lb_min_local,2
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,0)


             ENDDO
              DO lxb=MAX(lb_min_local-0-1,0),1
             jco=coset(lxb,1,0)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,0)


             ENDDO
             lxb = 0
             !jco=coset(0,2,0)
             jco=8

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,0,0)


              DO lxb=MAX(lb_min_local-1-0,0),1
             jco=coset(lxb,0,1)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,1)


             ENDDO
             lxb = 0
             !jco=coset(0,1,1)
             jco=9

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,1)


             lxb = 0
             !jco=coset(0,0,2)
             jco=10

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,2)



   ENDDO

    END SUBROUTINE xyz_to_vab_0_2



          ! Fully unrolled version
    ! *****************************************************************************
    SUBROUTINE xyz_to_vab_0_3 (prefactor, coef_xyz, & 
                           lp, la_max_local, lb_max_local, rp, ra, rab, vab, coset, &
                           la_min_local, lb_min_local, maxl, lvab, hvab)
       

       USE kinds,                        ONLY: dp

       ! PARAMETERS
       REAL(KIND=dp), INTENT(IN) :: prefactor
       INTEGER, INTENT(IN) :: lp, la_max_local, lb_max_local, &
                              lb_min_local, la_min_local, maxl, &
                              lvab,hvab

       INTEGER,  DIMENSION(-1:maxl,-1:maxl,-1:maxl), &
                        INTENT(IN)     :: coset

       REAL(KIND=dp),  DIMENSION(((lp+1)*(lp+2)*(lp+3))/6), &
                        INTENT(INOUT) :: coef_xyz
       REAL(KIND=dp),  DIMENSION(lvab,hvab), INTENT(INOUT) :: vab
       REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: rp, ra, rab
   
       ! LOCAL
       INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp, iaxis, l, k, lya,lyb,lza,lzb
       REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb,const_rp,const_ra
   

    REAL(KIND=dp) :: alpha_1(0:(0+3),0:0,0:3)
    REAL(KIND=dp) :: alpha_2(0:(0+3),0:0,0:3)
    REAL(KIND=dp) :: alpha_3(0:(0+3),0:0,0:3)
    REAL(KIND=dp) :: coef_ttz(0:0,0:3)
    REAL(KIND=dp) :: coef_tyz(0:0,0:3,0:0,0:3)
    coef_xyz=coef_xyz*prefactor

  alpha_1(:,:,:)=0.0_dp
  alpha_2(:,:,:)=0.0_dp
  alpha_3(:,:,:)=0.0_dp

         const_rp = rp(1)-(ra(1)+rab(1))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,0,0)=alpha_1(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,0,1)=alpha_1(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(0,0,1)=alpha_1(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,0,2)=alpha_1(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_1(1,0,2)=alpha_1(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_1(0,0,2)=alpha_1(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(3,0,3)=alpha_1(3,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*const_rp
             alpha_1(2,0,3)=alpha_1(2,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(1,0,3)=alpha_1(1,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*const_rp
             alpha_1(0,0,3)=alpha_1(0,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         const_rp = rp(2)-(ra(2)+rab(2))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,0,0)=alpha_2(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,0,1)=alpha_2(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(0,0,1)=alpha_2(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,0,2)=alpha_2(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_2(1,0,2)=alpha_2(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_2(0,0,2)=alpha_2(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(3,0,3)=alpha_2(3,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*const_rp
             alpha_2(2,0,3)=alpha_2(2,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(1,0,3)=alpha_2(1,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*const_rp
             alpha_2(0,0,3)=alpha_2(0,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         const_rp = rp(3)-(ra(3)+rab(3))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,0,0)=alpha_3(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,0,1)=alpha_3(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(0,0,1)=alpha_3(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,0,2)=alpha_3(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_3(1,0,2)=alpha_3(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_3(0,0,2)=alpha_3(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(3,0,3)=alpha_3(3,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*const_rp
             alpha_3(2,0,3)=alpha_3(2,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(1,0,3)=alpha_3(1,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*const_rp
             alpha_3(0,0,3)=alpha_3(0,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0

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

                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
                coef_tyz(0,1,0,0)=coef_tyz(0,1,0,0)+coef_ttz(0,1)*alpha_2(lyp,0,0)
                coef_tyz(0,2,0,0)=coef_tyz(0,2,0,0)+coef_ttz(0,2)*alpha_2(lyp,0,0)
                coef_tyz(0,3,0,0)=coef_tyz(0,3,0,0)+coef_ttz(0,3)*alpha_2(lyp,0,0)
                coef_tyz(0,0,0,1)=coef_tyz(0,0,0,1)+coef_ttz(0,0)*alpha_2(lyp,0,1)
                coef_tyz(0,1,0,1)=coef_tyz(0,1,0,1)+coef_ttz(0,1)*alpha_2(lyp,0,1)
                coef_tyz(0,2,0,1)=coef_tyz(0,2,0,1)+coef_ttz(0,2)*alpha_2(lyp,0,1)
                coef_tyz(0,0,0,2)=coef_tyz(0,0,0,2)+coef_ttz(0,0)*alpha_2(lyp,0,2)
                coef_tyz(0,1,0,2)=coef_tyz(0,1,0,2)+coef_ttz(0,1)*alpha_2(lyp,0,2)
                coef_tyz(0,0,0,3)=coef_tyz(0,0,0,3)+coef_ttz(0,0)*alpha_2(lyp,0,3)
       ENDDO


              DO lxb=lb_min_local,3
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,0)


             ENDDO
              DO lxb=MAX(lb_min_local-0-1,0),2
             jco=coset(lxb,1,0)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,0)


             ENDDO
              DO lxb=MAX(lb_min_local-0-2,0),1
             jco=coset(lxb,2,0)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,0,0)


             ENDDO
             lxb = 0
             !jco=coset(0,3,0)
             jco=17

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,3)*alpha_3(lzp,0,0)


              DO lxb=MAX(lb_min_local-1-0,0),2
             jco=coset(lxb,0,1)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,1)


             ENDDO
              DO lxb=MAX(lb_min_local-1-1,0),1
             jco=coset(lxb,1,1)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,1)


             ENDDO
             lxb = 0
             !jco=coset(0,2,1)
             jco=18

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,0,1)


              DO lxb=MAX(lb_min_local-2-0,0),1
             jco=coset(lxb,0,2)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,2)


             ENDDO
             lxb = 0
             !jco=coset(0,1,2)
             jco=19

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,2)


             lxb = 0
             !jco=coset(0,0,3)
             jco=20

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,3)



   ENDDO

    END SUBROUTINE xyz_to_vab_0_3



          ! Fully unrolled version
    ! *****************************************************************************
    SUBROUTINE xyz_to_vab_0_4 (prefactor, coef_xyz, & 
                           lp, la_max_local, lb_max_local, rp, ra, rab, vab, coset, &
                           la_min_local, lb_min_local, maxl, lvab, hvab)
       

       USE kinds,                        ONLY: dp

       ! PARAMETERS
       REAL(KIND=dp), INTENT(IN) :: prefactor
       INTEGER, INTENT(IN) :: lp, la_max_local, lb_max_local, &
                              lb_min_local, la_min_local, maxl, &
                              lvab,hvab

       INTEGER,  DIMENSION(-1:maxl,-1:maxl,-1:maxl), &
                        INTENT(IN)     :: coset

       REAL(KIND=dp),  DIMENSION(((lp+1)*(lp+2)*(lp+3))/6), &
                        INTENT(INOUT) :: coef_xyz
       REAL(KIND=dp),  DIMENSION(lvab,hvab), INTENT(INOUT) :: vab
       REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: rp, ra, rab
   
       ! LOCAL
       INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp, iaxis, l, k, lya,lyb,lza,lzb
       REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb,const_rp,const_ra
   

    REAL(KIND=dp) :: alpha_1(0:(0+4),0:0,0:4)
    REAL(KIND=dp) :: alpha_2(0:(0+4),0:0,0:4)
    REAL(KIND=dp) :: alpha_3(0:(0+4),0:0,0:4)
    REAL(KIND=dp) :: coef_ttz(0:0,0:4)
    REAL(KIND=dp) :: coef_tyz(0:0,0:4,0:0,0:4)
    coef_xyz=coef_xyz*prefactor

  alpha_1(:,:,:)=0.0_dp
  alpha_2(:,:,:)=0.0_dp
  alpha_3(:,:,:)=0.0_dp

         const_rp = rp(1)-(ra(1)+rab(1))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,0,0)=alpha_1(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,0,1)=alpha_1(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(0,0,1)=alpha_1(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,0,2)=alpha_1(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_1(1,0,2)=alpha_1(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_1(0,0,2)=alpha_1(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(3,0,3)=alpha_1(3,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*const_rp
             alpha_1(2,0,3)=alpha_1(2,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(1,0,3)=alpha_1(1,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*const_rp
             alpha_1(0,0,3)=alpha_1(0,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(4,0,4)=alpha_1(4,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(4,dp)
                 b=b*const_rp
             alpha_1(3,0,4)=alpha_1(3,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_1(2,0,4)=alpha_1(2,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)/REAL(3,dp)
                 b=b*const_rp
             alpha_1(1,0,4)=alpha_1(1,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(4,dp)
                 b=b*const_rp
             alpha_1(0,0,4)=alpha_1(0,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         const_rp = rp(2)-(ra(2)+rab(2))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,0,0)=alpha_2(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,0,1)=alpha_2(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(0,0,1)=alpha_2(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,0,2)=alpha_2(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_2(1,0,2)=alpha_2(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_2(0,0,2)=alpha_2(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(3,0,3)=alpha_2(3,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*const_rp
             alpha_2(2,0,3)=alpha_2(2,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(1,0,3)=alpha_2(1,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*const_rp
             alpha_2(0,0,3)=alpha_2(0,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(4,0,4)=alpha_2(4,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(4,dp)
                 b=b*const_rp
             alpha_2(3,0,4)=alpha_2(3,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_2(2,0,4)=alpha_2(2,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)/REAL(3,dp)
                 b=b*const_rp
             alpha_2(1,0,4)=alpha_2(1,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(4,dp)
                 b=b*const_rp
             alpha_2(0,0,4)=alpha_2(0,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         const_rp = rp(3)-(ra(3)+rab(3))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,0,0)=alpha_3(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,0,1)=alpha_3(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(0,0,1)=alpha_3(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,0,2)=alpha_3(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_3(1,0,2)=alpha_3(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_3(0,0,2)=alpha_3(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(3,0,3)=alpha_3(3,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*const_rp
             alpha_3(2,0,3)=alpha_3(2,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(1,0,3)=alpha_3(1,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*const_rp
             alpha_3(0,0,3)=alpha_3(0,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(4,0,4)=alpha_3(4,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(4,dp)
                 b=b*const_rp
             alpha_3(3,0,4)=alpha_3(3,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_3(2,0,4)=alpha_3(2,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)/REAL(3,dp)
                 b=b*const_rp
             alpha_3(1,0,4)=alpha_3(1,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(4,dp)
                 b=b*const_rp
             alpha_3(0,0,4)=alpha_3(0,0,4)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0

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

                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
                coef_tyz(0,1,0,0)=coef_tyz(0,1,0,0)+coef_ttz(0,1)*alpha_2(lyp,0,0)
                coef_tyz(0,2,0,0)=coef_tyz(0,2,0,0)+coef_ttz(0,2)*alpha_2(lyp,0,0)
                coef_tyz(0,3,0,0)=coef_tyz(0,3,0,0)+coef_ttz(0,3)*alpha_2(lyp,0,0)
                coef_tyz(0,4,0,0)=coef_tyz(0,4,0,0)+coef_ttz(0,4)*alpha_2(lyp,0,0)
                coef_tyz(0,0,0,1)=coef_tyz(0,0,0,1)+coef_ttz(0,0)*alpha_2(lyp,0,1)
                coef_tyz(0,1,0,1)=coef_tyz(0,1,0,1)+coef_ttz(0,1)*alpha_2(lyp,0,1)
                coef_tyz(0,2,0,1)=coef_tyz(0,2,0,1)+coef_ttz(0,2)*alpha_2(lyp,0,1)
                coef_tyz(0,3,0,1)=coef_tyz(0,3,0,1)+coef_ttz(0,3)*alpha_2(lyp,0,1)
                coef_tyz(0,0,0,2)=coef_tyz(0,0,0,2)+coef_ttz(0,0)*alpha_2(lyp,0,2)
                coef_tyz(0,1,0,2)=coef_tyz(0,1,0,2)+coef_ttz(0,1)*alpha_2(lyp,0,2)
                coef_tyz(0,2,0,2)=coef_tyz(0,2,0,2)+coef_ttz(0,2)*alpha_2(lyp,0,2)
                coef_tyz(0,0,0,3)=coef_tyz(0,0,0,3)+coef_ttz(0,0)*alpha_2(lyp,0,3)
                coef_tyz(0,1,0,3)=coef_tyz(0,1,0,3)+coef_ttz(0,1)*alpha_2(lyp,0,3)
                coef_tyz(0,0,0,4)=coef_tyz(0,0,0,4)+coef_ttz(0,0)*alpha_2(lyp,0,4)
       ENDDO


              DO lxb=lb_min_local,4
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,0)


             ENDDO
              DO lxb=MAX(lb_min_local-0-1,0),3
             jco=coset(lxb,1,0)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,0)


             ENDDO
              DO lxb=MAX(lb_min_local-0-2,0),2
             jco=coset(lxb,2,0)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,0,0)


             ENDDO
              DO lxb=MAX(lb_min_local-0-3,0),1
             jco=coset(lxb,3,0)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,3)*alpha_3(lzp,0,0)


             ENDDO
             lxb = 0
             !jco=coset(0,4,0)
             jco=31

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,4)*alpha_3(lzp,0,0)


              DO lxb=MAX(lb_min_local-1-0,0),3
             jco=coset(lxb,0,1)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,1)


             ENDDO
              DO lxb=MAX(lb_min_local-1-1,0),2
             jco=coset(lxb,1,1)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,1)


             ENDDO
              DO lxb=MAX(lb_min_local-1-2,0),1
             jco=coset(lxb,2,1)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,0,1)


             ENDDO
             lxb = 0
             !jco=coset(0,3,1)
             jco=32

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,3)*alpha_3(lzp,0,1)


              DO lxb=MAX(lb_min_local-2-0,0),2
             jco=coset(lxb,0,2)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,2)


             ENDDO
              DO lxb=MAX(lb_min_local-2-1,0),1
             jco=coset(lxb,1,2)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,2)


             ENDDO
             lxb = 0
             !jco=coset(0,2,2)
             jco=33

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,0,2)


              DO lxb=MAX(lb_min_local-3-0,0),1
             jco=coset(lxb,0,3)

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,3)


             ENDDO
             lxb = 0
             !jco=coset(0,1,3)
             jco=34

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,3)


             lxb = 0
             !jco=coset(0,0,4)
             jco=35

                lxa = 0
                !ico=coset(0,0,0)
                ico=1

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,4)



   ENDDO

    END SUBROUTINE xyz_to_vab_0_4



          ! Fully unrolled version
    ! *****************************************************************************
    SUBROUTINE xyz_to_vab_1_0 (prefactor, coef_xyz, & 
                           lp, la_max_local, lb_max_local, rp, ra, rab, vab, coset, &
                           la_min_local, lb_min_local, maxl, lvab, hvab)
       

       USE kinds,                        ONLY: dp

       ! PARAMETERS
       REAL(KIND=dp), INTENT(IN) :: prefactor
       INTEGER, INTENT(IN) :: lp, la_max_local, lb_max_local, &
                              lb_min_local, la_min_local, maxl, &
                              lvab,hvab

       INTEGER,  DIMENSION(-1:maxl,-1:maxl,-1:maxl), &
                        INTENT(IN)     :: coset

       REAL(KIND=dp),  DIMENSION(((lp+1)*(lp+2)*(lp+3))/6), &
                        INTENT(INOUT) :: coef_xyz
       REAL(KIND=dp),  DIMENSION(lvab,hvab), INTENT(INOUT) :: vab
       REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: rp, ra, rab
   
       ! LOCAL
       INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp, iaxis, l, k, lya,lyb,lza,lzb
       REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb,const_rp,const_ra
   

    REAL(KIND=dp) :: alpha_1(0:(1+0),0:1,0:0)
    REAL(KIND=dp) :: alpha_2(0:(1+0),0:1,0:0)
    REAL(KIND=dp) :: alpha_3(0:(1+0),0:1,0:0)
    REAL(KIND=dp) :: coef_ttz(0:1,0:0)
    REAL(KIND=dp) :: coef_tyz(0:1,0:0,0:1,0:0)
    coef_xyz=coef_xyz*prefactor

  alpha_1(:,:,:)=0.0_dp
  alpha_2(:,:,:)=0.0_dp
  alpha_3(:,:,:)=0.0_dp

         const_ra = (-ra(1)+rp(1))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,0,0)=alpha_1(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,1,0)=alpha_1(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,1,0)=alpha_1(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         const_ra = (-ra(2)+rp(2))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,0,0)=alpha_2(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,1,0)=alpha_2(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,1,0)=alpha_2(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         const_ra = (-ra(3)+rp(3))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,0,0)=alpha_3(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,1,0)=alpha_3(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,1,0)=alpha_3(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0

    lxyz=0
       coef_tyz=0.0_dp
       DO lyp=0,1+0-0
          coef_ttz=0.0_dp
          DO lxp=0,1+0-0-lyp
             lxyz=lxyz+1
                coef_ttz(0,0)=coef_ttz(0,0)+coef_xyz(lxyz)*alpha_1(lxp,0,0)
                coef_ttz(1,0)=coef_ttz(1,0)+coef_xyz(lxyz)*alpha_1(lxp,1,0)
          ENDDO

                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
                coef_tyz(1,0,0,0)=coef_tyz(1,0,0,0)+coef_ttz(1,0)*alpha_2(lyp,0,0)
                coef_tyz(0,0,1,0)=coef_tyz(0,0,1,0)+coef_ttz(0,0)*alpha_2(lyp,1,0)

       ENDDO


             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=la_min_local,1
                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(0,0,0)

             ENDDO

             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(0,0,0)


             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(0,1,0)



       coef_tyz=0.0_dp
       DO lyp=0,1+0-1
          coef_ttz=0.0_dp
          DO lxp=0,1+0-1-lyp
             lxyz=lxyz+1
                coef_ttz(0,0)=coef_ttz(0,0)+coef_xyz(lxyz)*alpha_1(lxp,0,0)
                coef_ttz(1,0)=coef_ttz(1,0)+coef_xyz(lxyz)*alpha_1(lxp,1,0)
          ENDDO

                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
                coef_tyz(1,0,0,0)=coef_tyz(1,0,0,0)+coef_ttz(1,0)*alpha_2(lyp,0,0)
                coef_tyz(0,0,1,0)=coef_tyz(0,0,1,0)+coef_ttz(0,0)*alpha_2(lyp,1,0)

       ENDDO


             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=la_min_local,1
                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(1,0,0)

             ENDDO

             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(1,0,0)


             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(1,1,0)




    END SUBROUTINE xyz_to_vab_1_0



          ! Fully unrolled version
    ! *****************************************************************************
    SUBROUTINE xyz_to_vab_1_1 (prefactor, coef_xyz, & 
                           lp, la_max_local, lb_max_local, rp, ra, rab, vab, coset, &
                           la_min_local, lb_min_local, maxl, lvab, hvab)
       

       USE kinds,                        ONLY: dp

       ! PARAMETERS
       REAL(KIND=dp), INTENT(IN) :: prefactor
       INTEGER, INTENT(IN) :: lp, la_max_local, lb_max_local, &
                              lb_min_local, la_min_local, maxl, &
                              lvab,hvab

       INTEGER,  DIMENSION(-1:maxl,-1:maxl,-1:maxl), &
                        INTENT(IN)     :: coset

       REAL(KIND=dp),  DIMENSION(((lp+1)*(lp+2)*(lp+3))/6), &
                        INTENT(INOUT) :: coef_xyz
       REAL(KIND=dp),  DIMENSION(lvab,hvab), INTENT(INOUT) :: vab
       REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: rp, ra, rab
   
       ! LOCAL
       INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp, iaxis, l, k, lya,lyb,lza,lzb
       REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb,const_rp,const_ra
   

    REAL(KIND=dp) :: alpha_1(0:(1+1),0:1,0:1)
    REAL(KIND=dp) :: alpha_2(0:(1+1),0:1,0:1)
    REAL(KIND=dp) :: alpha_3(0:(1+1),0:1,0:1)
    REAL(KIND=dp) :: coef_ttz(0:1,0:1)
    REAL(KIND=dp) :: coef_tyz(0:1,0:1,0:1,0:1)
    coef_xyz=coef_xyz*prefactor

  alpha_1(:,:,:)=0.0_dp
  alpha_2(:,:,:)=0.0_dp
  alpha_3(:,:,:)=0.0_dp

         const_rp = rp(1)-(ra(1)+rab(1))
         const_ra = (-ra(1)+rp(1))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,0,0)=alpha_1(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,0,1)=alpha_1(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(0,0,1)=alpha_1(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,1,0)=alpha_1(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,1,0)=alpha_1(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,1,1)=alpha_1(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(1,1,1)=alpha_1(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,1,1)=alpha_1(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(0,1,1)=alpha_1(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         const_rp = rp(2)-(ra(2)+rab(2))
         const_ra = (-ra(2)+rp(2))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,0,0)=alpha_2(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,0,1)=alpha_2(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(0,0,1)=alpha_2(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,1,0)=alpha_2(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,1,0)=alpha_2(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,1,1)=alpha_2(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(1,1,1)=alpha_2(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,1,1)=alpha_2(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(0,1,1)=alpha_2(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         const_rp = rp(3)-(ra(3)+rab(3))
         const_ra = (-ra(3)+rp(3))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,0,0)=alpha_3(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,0,1)=alpha_3(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(0,0,1)=alpha_3(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,1,0)=alpha_3(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,1,0)=alpha_3(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,1,1)=alpha_3(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(1,1,1)=alpha_3(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,1,1)=alpha_3(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(0,1,1)=alpha_3(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0

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

                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
                coef_tyz(1,0,0,0)=coef_tyz(1,0,0,0)+coef_ttz(1,0)*alpha_2(lyp,0,0)
                coef_tyz(0,1,0,0)=coef_tyz(0,1,0,0)+coef_ttz(0,1)*alpha_2(lyp,0,0)
                coef_tyz(1,1,0,0)=coef_tyz(1,1,0,0)+coef_ttz(1,1)*alpha_2(lyp,0,0)
                coef_tyz(0,0,1,0)=coef_tyz(0,0,1,0)+coef_ttz(0,0)*alpha_2(lyp,1,0)
                coef_tyz(0,1,1,0)=coef_tyz(0,1,1,0)+coef_ttz(0,1)*alpha_2(lyp,1,0)
                coef_tyz(0,0,0,1)=coef_tyz(0,0,0,1)+coef_ttz(0,0)*alpha_2(lyp,0,1)
                coef_tyz(1,0,0,1)=coef_tyz(1,0,0,1)+coef_ttz(1,0)*alpha_2(lyp,0,1)
                coef_tyz(0,0,1,1)=coef_tyz(0,0,1,1)+coef_ttz(0,0)*alpha_2(lyp,1,1)
       ENDDO


              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,0)


             ENDDO
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,0)

             ENDDO

             lxb = 0
             !jco=coset(0,1,0)
             jco=3

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,0,0)


              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,0)


             ENDDO
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,1,0)


             lxb = 0
             !jco=coset(0,0,1)
             jco=4

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,1)

             ENDDO

             lxb = 0
             !jco=coset(0,0,1)
             jco=4

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,1)


             lxb = 0
             !jco=coset(0,0,1)
             jco=4

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,1)



   ENDDO

    END SUBROUTINE xyz_to_vab_1_1



          ! Fully unrolled version
    ! *****************************************************************************
    SUBROUTINE xyz_to_vab_1_2 (prefactor, coef_xyz, & 
                           lp, la_max_local, lb_max_local, rp, ra, rab, vab, coset, &
                           la_min_local, lb_min_local, maxl, lvab, hvab)
       

       USE kinds,                        ONLY: dp

       ! PARAMETERS
       REAL(KIND=dp), INTENT(IN) :: prefactor
       INTEGER, INTENT(IN) :: lp, la_max_local, lb_max_local, &
                              lb_min_local, la_min_local, maxl, &
                              lvab,hvab

       INTEGER,  DIMENSION(-1:maxl,-1:maxl,-1:maxl), &
                        INTENT(IN)     :: coset

       REAL(KIND=dp),  DIMENSION(((lp+1)*(lp+2)*(lp+3))/6), &
                        INTENT(INOUT) :: coef_xyz
       REAL(KIND=dp),  DIMENSION(lvab,hvab), INTENT(INOUT) :: vab
       REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: rp, ra, rab
   
       ! LOCAL
       INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp, iaxis, l, k, lya,lyb,lza,lzb
       REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb,const_rp,const_ra
   

    REAL(KIND=dp) :: alpha_1(0:(1+2),0:1,0:2)
    REAL(KIND=dp) :: alpha_2(0:(1+2),0:1,0:2)
    REAL(KIND=dp) :: alpha_3(0:(1+2),0:1,0:2)
    REAL(KIND=dp) :: coef_ttz(0:1,0:2)
    REAL(KIND=dp) :: coef_tyz(0:1,0:2,0:1,0:2)
    coef_xyz=coef_xyz*prefactor

  alpha_1(:,:,:)=0.0_dp
  alpha_2(:,:,:)=0.0_dp
  alpha_3(:,:,:)=0.0_dp

         const_rp = rp(1)-(ra(1)+rab(1))
         const_ra = (-ra(1)+rp(1))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,0,0)=alpha_1(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,0,1)=alpha_1(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(0,0,1)=alpha_1(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,0,2)=alpha_1(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_1(1,0,2)=alpha_1(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_1(0,0,2)=alpha_1(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,1,0)=alpha_1(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,1,0)=alpha_1(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,1,1)=alpha_1(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(1,1,1)=alpha_1(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,1,1)=alpha_1(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(0,1,1)=alpha_1(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(3,1,2)=alpha_1(3,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_1(2,1,2)=alpha_1(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_1(1,1,2)=alpha_1(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,1,2)=alpha_1(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_1(1,1,2)=alpha_1(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_1(0,1,2)=alpha_1(0,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         const_rp = rp(2)-(ra(2)+rab(2))
         const_ra = (-ra(2)+rp(2))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,0,0)=alpha_2(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,0,1)=alpha_2(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(0,0,1)=alpha_2(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,0,2)=alpha_2(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_2(1,0,2)=alpha_2(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_2(0,0,2)=alpha_2(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,1,0)=alpha_2(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,1,0)=alpha_2(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,1,1)=alpha_2(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(1,1,1)=alpha_2(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,1,1)=alpha_2(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(0,1,1)=alpha_2(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(3,1,2)=alpha_2(3,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_2(2,1,2)=alpha_2(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_2(1,1,2)=alpha_2(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,1,2)=alpha_2(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_2(1,1,2)=alpha_2(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_2(0,1,2)=alpha_2(0,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         const_rp = rp(3)-(ra(3)+rab(3))
         const_ra = (-ra(3)+rp(3))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,0,0)=alpha_3(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,0,1)=alpha_3(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(0,0,1)=alpha_3(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,0,2)=alpha_3(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_3(1,0,2)=alpha_3(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_3(0,0,2)=alpha_3(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,1,0)=alpha_3(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,1,0)=alpha_3(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,1,1)=alpha_3(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(1,1,1)=alpha_3(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,1,1)=alpha_3(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(0,1,1)=alpha_3(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(3,1,2)=alpha_3(3,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_3(2,1,2)=alpha_3(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_3(1,1,2)=alpha_3(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,1,2)=alpha_3(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_3(1,1,2)=alpha_3(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_3(0,1,2)=alpha_3(0,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0

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

                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
                coef_tyz(1,0,0,0)=coef_tyz(1,0,0,0)+coef_ttz(1,0)*alpha_2(lyp,0,0)
                coef_tyz(0,1,0,0)=coef_tyz(0,1,0,0)+coef_ttz(0,1)*alpha_2(lyp,0,0)
                coef_tyz(1,1,0,0)=coef_tyz(1,1,0,0)+coef_ttz(1,1)*alpha_2(lyp,0,0)
                coef_tyz(0,2,0,0)=coef_tyz(0,2,0,0)+coef_ttz(0,2)*alpha_2(lyp,0,0)
                coef_tyz(1,2,0,0)=coef_tyz(1,2,0,0)+coef_ttz(1,2)*alpha_2(lyp,0,0)
                coef_tyz(0,0,1,0)=coef_tyz(0,0,1,0)+coef_ttz(0,0)*alpha_2(lyp,1,0)
                coef_tyz(0,1,1,0)=coef_tyz(0,1,1,0)+coef_ttz(0,1)*alpha_2(lyp,1,0)
                coef_tyz(0,2,1,0)=coef_tyz(0,2,1,0)+coef_ttz(0,2)*alpha_2(lyp,1,0)
                coef_tyz(0,0,0,1)=coef_tyz(0,0,0,1)+coef_ttz(0,0)*alpha_2(lyp,0,1)
                coef_tyz(1,0,0,1)=coef_tyz(1,0,0,1)+coef_ttz(1,0)*alpha_2(lyp,0,1)
                coef_tyz(0,1,0,1)=coef_tyz(0,1,0,1)+coef_ttz(0,1)*alpha_2(lyp,0,1)
                coef_tyz(1,1,0,1)=coef_tyz(1,1,0,1)+coef_ttz(1,1)*alpha_2(lyp,0,1)
                coef_tyz(0,0,1,1)=coef_tyz(0,0,1,1)+coef_ttz(0,0)*alpha_2(lyp,1,1)
                coef_tyz(0,1,1,1)=coef_tyz(0,1,1,1)+coef_ttz(0,1)*alpha_2(lyp,1,1)
                coef_tyz(0,0,0,2)=coef_tyz(0,0,0,2)+coef_ttz(0,0)*alpha_2(lyp,0,2)
                coef_tyz(1,0,0,2)=coef_tyz(1,0,0,2)+coef_ttz(1,0)*alpha_2(lyp,0,2)
                coef_tyz(0,0,1,2)=coef_tyz(0,0,1,2)+coef_ttz(0,0)*alpha_2(lyp,1,2)
       ENDDO


              DO lxb=lb_min_local,2
             jco=coset(lxb,0,0)

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO
              DO lxb=lb_min_local,2
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,0)


             ENDDO
              DO lxb=MAX(lb_min_local-0-1,0),1
             jco=coset(lxb,1,0)

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO
              DO lxb=MAX(lb_min_local-0-1,0),1
             jco=coset(lxb,1,0)

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,0,0)


             ENDDO
             lxb = 0
             !jco=coset(0,2,0)
             jco=8

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,0,0)

             ENDDO

             lxb = 0
             !jco=coset(0,2,0)
             jco=8

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,2)*alpha_3(lzp,0,0)


              DO lxb=lb_min_local,2
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,0)


             ENDDO
              DO lxb=MAX(lb_min_local-0-1,0),1
             jco=coset(lxb,1,0)

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,1,0)


             ENDDO
             lxb = 0
             !jco=coset(0,2,0)
             jco=8

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,1,0)


              DO lxb=MAX(lb_min_local-1-0,0),1
             jco=coset(lxb,0,1)

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,1)

             ENDDO

             ENDDO
              DO lxb=MAX(lb_min_local-1-0,0),1
             jco=coset(lxb,0,1)

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,1)


             ENDDO
             lxb = 0
             !jco=coset(0,1,1)
             jco=9

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,1)

             ENDDO

             lxb = 0
             !jco=coset(0,1,1)
             jco=9

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,0,1)


              DO lxb=MAX(lb_min_local-1-0,0),1
             jco=coset(lxb,0,1)

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,1)


             ENDDO
             lxb = 0
             !jco=coset(0,1,1)
             jco=9

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,1,1)


             lxb = 0
             !jco=coset(0,0,2)
             jco=10

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,2)

             ENDDO

             lxb = 0
             !jco=coset(0,0,2)
             jco=10

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,2)


             lxb = 0
             !jco=coset(0,0,2)
             jco=10

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,2)



   ENDDO

    END SUBROUTINE xyz_to_vab_1_2



          ! Fully unrolled version
    ! *****************************************************************************
    SUBROUTINE xyz_to_vab_1_3 (prefactor, coef_xyz, & 
                           lp, la_max_local, lb_max_local, rp, ra, rab, vab, coset, &
                           la_min_local, lb_min_local, maxl, lvab, hvab)
       

       USE kinds,                        ONLY: dp

       ! PARAMETERS
       REAL(KIND=dp), INTENT(IN) :: prefactor
       INTEGER, INTENT(IN) :: lp, la_max_local, lb_max_local, &
                              lb_min_local, la_min_local, maxl, &
                              lvab,hvab

       INTEGER,  DIMENSION(-1:maxl,-1:maxl,-1:maxl), &
                        INTENT(IN)     :: coset

       REAL(KIND=dp),  DIMENSION(((lp+1)*(lp+2)*(lp+3))/6), &
                        INTENT(INOUT) :: coef_xyz
       REAL(KIND=dp),  DIMENSION(lvab,hvab), INTENT(INOUT) :: vab
       REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: rp, ra, rab
   
       ! LOCAL
       INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp, iaxis, l, k, lya,lyb,lza,lzb
       REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb,const_rp,const_ra
   

    REAL(KIND=dp) :: alpha_1(0:(1+3),0:1,0:3)
    REAL(KIND=dp) :: alpha_2(0:(1+3),0:1,0:3)
    REAL(KIND=dp) :: alpha_3(0:(1+3),0:1,0:3)
    REAL(KIND=dp) :: coef_ttz(0:1,0:3)
    REAL(KIND=dp) :: coef_tyz(0:1,0:3,0:1,0:3)
    coef_xyz=coef_xyz*prefactor

  alpha_1(:,:,:)=0.0_dp
  alpha_2(:,:,:)=0.0_dp
  alpha_3(:,:,:)=0.0_dp

         const_rp = rp(1)-(ra(1)+rab(1))
         const_ra = (-ra(1)+rp(1))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,0,0)=alpha_1(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,0,1)=alpha_1(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(0,0,1)=alpha_1(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,0,2)=alpha_1(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_1(1,0,2)=alpha_1(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_1(0,0,2)=alpha_1(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(3,0,3)=alpha_1(3,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*const_rp
             alpha_1(2,0,3)=alpha_1(2,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(1,0,3)=alpha_1(1,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*const_rp
             alpha_1(0,0,3)=alpha_1(0,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,1,0)=alpha_1(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,1,0)=alpha_1(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,1,1)=alpha_1(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(1,1,1)=alpha_1(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,1,1)=alpha_1(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(0,1,1)=alpha_1(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(3,1,2)=alpha_1(3,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_1(2,1,2)=alpha_1(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_1(1,1,2)=alpha_1(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,1,2)=alpha_1(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_1(1,1,2)=alpha_1(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_1(0,1,2)=alpha_1(0,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(4,1,3)=alpha_1(4,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*const_rp
             alpha_1(3,1,3)=alpha_1(3,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(2,1,3)=alpha_1(2,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*const_rp
             alpha_1(1,1,3)=alpha_1(1,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(3,1,3)=alpha_1(3,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*const_rp
             alpha_1(2,1,3)=alpha_1(2,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(1,1,3)=alpha_1(1,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*const_rp
             alpha_1(0,1,3)=alpha_1(0,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         const_rp = rp(2)-(ra(2)+rab(2))
         const_ra = (-ra(2)+rp(2))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,0,0)=alpha_2(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,0,1)=alpha_2(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(0,0,1)=alpha_2(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,0,2)=alpha_2(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_2(1,0,2)=alpha_2(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_2(0,0,2)=alpha_2(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(3,0,3)=alpha_2(3,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*const_rp
             alpha_2(2,0,3)=alpha_2(2,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(1,0,3)=alpha_2(1,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*const_rp
             alpha_2(0,0,3)=alpha_2(0,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,1,0)=alpha_2(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,1,0)=alpha_2(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,1,1)=alpha_2(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(1,1,1)=alpha_2(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,1,1)=alpha_2(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(0,1,1)=alpha_2(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(3,1,2)=alpha_2(3,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_2(2,1,2)=alpha_2(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_2(1,1,2)=alpha_2(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,1,2)=alpha_2(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_2(1,1,2)=alpha_2(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_2(0,1,2)=alpha_2(0,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(4,1,3)=alpha_2(4,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*const_rp
             alpha_2(3,1,3)=alpha_2(3,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(2,1,3)=alpha_2(2,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*const_rp
             alpha_2(1,1,3)=alpha_2(1,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(3,1,3)=alpha_2(3,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*const_rp
             alpha_2(2,1,3)=alpha_2(2,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(1,1,3)=alpha_2(1,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*const_rp
             alpha_2(0,1,3)=alpha_2(0,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         const_rp = rp(3)-(ra(3)+rab(3))
         const_ra = (-ra(3)+rp(3))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,0,0)=alpha_3(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,0,1)=alpha_3(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(0,0,1)=alpha_3(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,0,2)=alpha_3(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_3(1,0,2)=alpha_3(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_3(0,0,2)=alpha_3(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(3,0,3)=alpha_3(3,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*const_rp
             alpha_3(2,0,3)=alpha_3(2,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(1,0,3)=alpha_3(1,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*const_rp
             alpha_3(0,0,3)=alpha_3(0,0,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,1,0)=alpha_3(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,1,0)=alpha_3(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,1,1)=alpha_3(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(1,1,1)=alpha_3(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,1,1)=alpha_3(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(0,1,1)=alpha_3(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(3,1,2)=alpha_3(3,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_3(2,1,2)=alpha_3(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_3(1,1,2)=alpha_3(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,1,2)=alpha_3(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_3(1,1,2)=alpha_3(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_3(0,1,2)=alpha_3(0,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(4,1,3)=alpha_3(4,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*const_rp
             alpha_3(3,1,3)=alpha_3(3,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(2,1,3)=alpha_3(2,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*const_rp
             alpha_3(1,1,3)=alpha_3(1,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(3,1,3)=alpha_3(3,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(3,dp)
                 b=b*const_rp
             alpha_3(2,1,3)=alpha_3(2,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(1,1,3)=alpha_3(1,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(3,dp)
                 b=b*const_rp
             alpha_3(0,1,3)=alpha_3(0,1,3)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0

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

                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
                coef_tyz(1,0,0,0)=coef_tyz(1,0,0,0)+coef_ttz(1,0)*alpha_2(lyp,0,0)
                coef_tyz(0,1,0,0)=coef_tyz(0,1,0,0)+coef_ttz(0,1)*alpha_2(lyp,0,0)
                coef_tyz(1,1,0,0)=coef_tyz(1,1,0,0)+coef_ttz(1,1)*alpha_2(lyp,0,0)
                coef_tyz(0,2,0,0)=coef_tyz(0,2,0,0)+coef_ttz(0,2)*alpha_2(lyp,0,0)
                coef_tyz(1,2,0,0)=coef_tyz(1,2,0,0)+coef_ttz(1,2)*alpha_2(lyp,0,0)
                coef_tyz(0,3,0,0)=coef_tyz(0,3,0,0)+coef_ttz(0,3)*alpha_2(lyp,0,0)
                coef_tyz(1,3,0,0)=coef_tyz(1,3,0,0)+coef_ttz(1,3)*alpha_2(lyp,0,0)
                coef_tyz(0,0,1,0)=coef_tyz(0,0,1,0)+coef_ttz(0,0)*alpha_2(lyp,1,0)
                coef_tyz(0,1,1,0)=coef_tyz(0,1,1,0)+coef_ttz(0,1)*alpha_2(lyp,1,0)
                coef_tyz(0,2,1,0)=coef_tyz(0,2,1,0)+coef_ttz(0,2)*alpha_2(lyp,1,0)
                coef_tyz(0,3,1,0)=coef_tyz(0,3,1,0)+coef_ttz(0,3)*alpha_2(lyp,1,0)
                coef_tyz(0,0,0,1)=coef_tyz(0,0,0,1)+coef_ttz(0,0)*alpha_2(lyp,0,1)
                coef_tyz(1,0,0,1)=coef_tyz(1,0,0,1)+coef_ttz(1,0)*alpha_2(lyp,0,1)
                coef_tyz(0,1,0,1)=coef_tyz(0,1,0,1)+coef_ttz(0,1)*alpha_2(lyp,0,1)
                coef_tyz(1,1,0,1)=coef_tyz(1,1,0,1)+coef_ttz(1,1)*alpha_2(lyp,0,1)
                coef_tyz(0,2,0,1)=coef_tyz(0,2,0,1)+coef_ttz(0,2)*alpha_2(lyp,0,1)
                coef_tyz(1,2,0,1)=coef_tyz(1,2,0,1)+coef_ttz(1,2)*alpha_2(lyp,0,1)
                coef_tyz(0,0,1,1)=coef_tyz(0,0,1,1)+coef_ttz(0,0)*alpha_2(lyp,1,1)
                coef_tyz(0,1,1,1)=coef_tyz(0,1,1,1)+coef_ttz(0,1)*alpha_2(lyp,1,1)
                coef_tyz(0,2,1,1)=coef_tyz(0,2,1,1)+coef_ttz(0,2)*alpha_2(lyp,1,1)
                coef_tyz(0,0,0,2)=coef_tyz(0,0,0,2)+coef_ttz(0,0)*alpha_2(lyp,0,2)
                coef_tyz(1,0,0,2)=coef_tyz(1,0,0,2)+coef_ttz(1,0)*alpha_2(lyp,0,2)
                coef_tyz(0,1,0,2)=coef_tyz(0,1,0,2)+coef_ttz(0,1)*alpha_2(lyp,0,2)
                coef_tyz(1,1,0,2)=coef_tyz(1,1,0,2)+coef_ttz(1,1)*alpha_2(lyp,0,2)
                coef_tyz(0,0,1,2)=coef_tyz(0,0,1,2)+coef_ttz(0,0)*alpha_2(lyp,1,2)
                coef_tyz(0,1,1,2)=coef_tyz(0,1,1,2)+coef_ttz(0,1)*alpha_2(lyp,1,2)
                coef_tyz(0,0,0,3)=coef_tyz(0,0,0,3)+coef_ttz(0,0)*alpha_2(lyp,0,3)
                coef_tyz(1,0,0,3)=coef_tyz(1,0,0,3)+coef_ttz(1,0)*alpha_2(lyp,0,3)
                coef_tyz(0,0,1,3)=coef_tyz(0,0,1,3)+coef_ttz(0,0)*alpha_2(lyp,1,3)
       ENDDO


              DO lxb=lb_min_local,3
             jco=coset(lxb,0,0)

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO
              DO lxb=lb_min_local,3
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,0)


             ENDDO
              DO lxb=MAX(lb_min_local-0-1,0),2
             jco=coset(lxb,1,0)

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO
              DO lxb=MAX(lb_min_local-0-1,0),2
             jco=coset(lxb,1,0)

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,0,0)


             ENDDO
              DO lxb=MAX(lb_min_local-0-2,0),1
             jco=coset(lxb,2,0)

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO
              DO lxb=MAX(lb_min_local-0-2,0),1
             jco=coset(lxb,2,0)

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,2)*alpha_3(lzp,0,0)


             ENDDO
             lxb = 0
             !jco=coset(0,3,0)
             jco=17

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,3)*alpha_3(lzp,0,0)

             ENDDO

             lxb = 0
             !jco=coset(0,3,0)
             jco=17

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,3)*alpha_3(lzp,0,0)


              DO lxb=lb_min_local,3
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,0)


             ENDDO
              DO lxb=MAX(lb_min_local-0-1,0),2
             jco=coset(lxb,1,0)

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,1,0)


             ENDDO
              DO lxb=MAX(lb_min_local-0-2,0),1
             jco=coset(lxb,2,0)

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,1,0)


             ENDDO
             lxb = 0
             !jco=coset(0,3,0)
             jco=17

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,3)*alpha_3(lzp,1,0)


              DO lxb=MAX(lb_min_local-1-0,0),2
             jco=coset(lxb,0,1)

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,1)

             ENDDO

             ENDDO
              DO lxb=MAX(lb_min_local-1-0,0),2
             jco=coset(lxb,0,1)

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,1)


             ENDDO
              DO lxb=MAX(lb_min_local-1-1,0),1
             jco=coset(lxb,1,1)

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,1)

             ENDDO

             ENDDO
              DO lxb=MAX(lb_min_local-1-1,0),1
             jco=coset(lxb,1,1)

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,0,1)


             ENDDO
             lxb = 0
             !jco=coset(0,2,1)
             jco=18

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,0,1)

             ENDDO

             lxb = 0
             !jco=coset(0,2,1)
             jco=18

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,2)*alpha_3(lzp,0,1)


              DO lxb=MAX(lb_min_local-1-0,0),2
             jco=coset(lxb,0,1)

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,1)


             ENDDO
              DO lxb=MAX(lb_min_local-1-1,0),1
             jco=coset(lxb,1,1)

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,1,1)


             ENDDO
             lxb = 0
             !jco=coset(0,2,1)
             jco=18

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,1,1)


              DO lxb=MAX(lb_min_local-2-0,0),1
             jco=coset(lxb,0,2)

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,2)

             ENDDO

             ENDDO
              DO lxb=MAX(lb_min_local-2-0,0),1
             jco=coset(lxb,0,2)

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,2)


             ENDDO
             lxb = 0
             !jco=coset(0,1,2)
             jco=19

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,2)

             ENDDO

             lxb = 0
             !jco=coset(0,1,2)
             jco=19

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,0,2)


              DO lxb=MAX(lb_min_local-2-0,0),1
             jco=coset(lxb,0,2)

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,2)


             ENDDO
             lxb = 0
             !jco=coset(0,1,2)
             jco=19

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,1,2)


             lxb = 0
             !jco=coset(0,0,3)
             jco=20

              DO lxa=la_min_local,1

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,3)

             ENDDO

             lxb = 0
             !jco=coset(0,0,3)
             jco=20

                lxa = 0
                !ico=coset(0,1,0)
                ico=3

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,3)


             lxb = 0
             !jco=coset(0,0,3)
             jco=20

                lxa = 0
                !ico=coset(0,0,1)
                ico=4

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,3)



   ENDDO

    END SUBROUTINE xyz_to_vab_1_3



          ! Partial unrolled version
    ! *****************************************************************************
    SUBROUTINE xyz_to_vab_1_4 (prefactor, coef_xyz, & 
                           lp, la_max_local, lb_max_local, rp, ra, rab, vab, coset,&
                           la_min_local, lb_min_local, maxl, lvab, hvab)
 

      USE kinds,                        ONLY: dp
       ! PARAMETERS
       INTEGER, INTENT(IN) :: lp, la_max_local, lb_max_local, & 
                               maxl, lb_min_local, la_min_local, &
                               lvab,hvab

       INTEGER,  DIMENSION(-1:maxl,-1:maxl,-1:maxl), &
                       INTENT(IN)     :: coset

       REAL(KIND=dp), INTENT(IN) :: prefactor
       
       REAL(KIND=dp),  DIMENSION(((lp+1)*(lp+2)*(lp+3))/6),&
                         INTENT(INOUT) :: coef_xyz
       REAL(KIND=dp), DIMENSION(lvab,hvab), INTENT(INOUT) :: vab
       REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: rp, ra, rab
   
       ! LOCAL
       INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp, iaxis, l, k, lya,lyb,lza,lzb
       REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb,const_rp,const_ra
 

    REAL(KIND=dp) :: alpha(0:(1+4),0:1,0:4,3)
    REAL(KIND=dp) :: coef_ttz(0:1,0:4)
    REAL(KIND=dp) :: coef_tyz(0:1,0:4,0:1,0:4)
    
    coef_xyz=coef_xyz*prefactor

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

    END SUBROUTINE xyz_to_vab_1_4

          ! Fully unrolled version
    ! *****************************************************************************
    SUBROUTINE xyz_to_vab_2_0 (prefactor, coef_xyz, & 
                           lp, la_max_local, lb_max_local, rp, ra, rab, vab, coset, &
                           la_min_local, lb_min_local, maxl, lvab, hvab)
       

       USE kinds,                        ONLY: dp

       ! PARAMETERS
       REAL(KIND=dp), INTENT(IN) :: prefactor
       INTEGER, INTENT(IN) :: lp, la_max_local, lb_max_local, &
                              lb_min_local, la_min_local, maxl, &
                              lvab,hvab

       INTEGER,  DIMENSION(-1:maxl,-1:maxl,-1:maxl), &
                        INTENT(IN)     :: coset

       REAL(KIND=dp),  DIMENSION(((lp+1)*(lp+2)*(lp+3))/6), &
                        INTENT(INOUT) :: coef_xyz
       REAL(KIND=dp),  DIMENSION(lvab,hvab), INTENT(INOUT) :: vab
       REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: rp, ra, rab
   
       ! LOCAL
       INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp, iaxis, l, k, lya,lyb,lza,lzb
       REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb,const_rp,const_ra
   

    REAL(KIND=dp) :: alpha_1(0:(2+0),0:2,0:0)
    REAL(KIND=dp) :: alpha_2(0:(2+0),0:2,0:0)
    REAL(KIND=dp) :: alpha_3(0:(2+0),0:2,0:0)
    REAL(KIND=dp) :: coef_ttz(0:2,0:0)
    REAL(KIND=dp) :: coef_tyz(0:2,0:0,0:2,0:0)
    coef_xyz=coef_xyz*prefactor

  alpha_1(:,:,:)=0.0_dp
  alpha_2(:,:,:)=0.0_dp
  alpha_3(:,:,:)=0.0_dp

         const_ra = (-ra(1)+rp(1))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,0,0)=alpha_1(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,1,0)=alpha_1(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,1,0)=alpha_1(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,2,0)=alpha_1(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,2,0)=alpha_1(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,2,0)=alpha_1(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         const_ra = (-ra(2)+rp(2))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,0,0)=alpha_2(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,1,0)=alpha_2(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,1,0)=alpha_2(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,2,0)=alpha_2(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,2,0)=alpha_2(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,2,0)=alpha_2(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         const_ra = (-ra(3)+rp(3))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,0,0)=alpha_3(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,1,0)=alpha_3(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,1,0)=alpha_3(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,2,0)=alpha_3(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,2,0)=alpha_3(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,2,0)=alpha_3(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0

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

                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
                coef_tyz(1,0,0,0)=coef_tyz(1,0,0,0)+coef_ttz(1,0)*alpha_2(lyp,0,0)
                coef_tyz(2,0,0,0)=coef_tyz(2,0,0,0)+coef_ttz(2,0)*alpha_2(lyp,0,0)
                coef_tyz(0,0,1,0)=coef_tyz(0,0,1,0)+coef_ttz(0,0)*alpha_2(lyp,1,0)
                coef_tyz(1,0,1,0)=coef_tyz(1,0,1,0)+coef_ttz(1,0)*alpha_2(lyp,1,0)
                coef_tyz(0,0,2,0)=coef_tyz(0,0,2,0)+coef_ttz(0,0)*alpha_2(lyp,2,0)
       ENDDO


             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=la_min_local,2

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,0)

             ENDDO

             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-0-1,0),1

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,0)

             ENDDO

             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,2,0)
                ico=8

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,0,0)


             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-1-0,0),1

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,0)

             ENDDO

             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,1,1)
                ico=9

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,1,0)


             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,0,2)
                ico=10

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,2,0)



   ENDDO

    END SUBROUTINE xyz_to_vab_2_0



          ! Fully unrolled version
    ! *****************************************************************************
    SUBROUTINE xyz_to_vab_2_1 (prefactor, coef_xyz, & 
                           lp, la_max_local, lb_max_local, rp, ra, rab, vab, coset, &
                           la_min_local, lb_min_local, maxl, lvab, hvab)
       

       USE kinds,                        ONLY: dp

       ! PARAMETERS
       REAL(KIND=dp), INTENT(IN) :: prefactor
       INTEGER, INTENT(IN) :: lp, la_max_local, lb_max_local, &
                              lb_min_local, la_min_local, maxl, &
                              lvab,hvab

       INTEGER,  DIMENSION(-1:maxl,-1:maxl,-1:maxl), &
                        INTENT(IN)     :: coset

       REAL(KIND=dp),  DIMENSION(((lp+1)*(lp+2)*(lp+3))/6), &
                        INTENT(INOUT) :: coef_xyz
       REAL(KIND=dp),  DIMENSION(lvab,hvab), INTENT(INOUT) :: vab
       REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: rp, ra, rab
   
       ! LOCAL
       INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp, iaxis, l, k, lya,lyb,lza,lzb
       REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb,const_rp,const_ra
   

    REAL(KIND=dp) :: alpha_1(0:(2+1),0:2,0:1)
    REAL(KIND=dp) :: alpha_2(0:(2+1),0:2,0:1)
    REAL(KIND=dp) :: alpha_3(0:(2+1),0:2,0:1)
    REAL(KIND=dp) :: coef_ttz(0:2,0:1)
    REAL(KIND=dp) :: coef_tyz(0:2,0:1,0:2,0:1)
    coef_xyz=coef_xyz*prefactor

  alpha_1(:,:,:)=0.0_dp
  alpha_2(:,:,:)=0.0_dp
  alpha_3(:,:,:)=0.0_dp

         const_rp = rp(1)-(ra(1)+rab(1))
         const_ra = (-ra(1)+rp(1))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,0,0)=alpha_1(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,0,1)=alpha_1(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(0,0,1)=alpha_1(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,1,0)=alpha_1(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,1,0)=alpha_1(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,1,1)=alpha_1(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(1,1,1)=alpha_1(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,1,1)=alpha_1(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(0,1,1)=alpha_1(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,2,0)=alpha_1(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,2,0)=alpha_1(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,2,0)=alpha_1(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(3,2,1)=alpha_1(3,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(2,2,1)=alpha_1(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,2,1)=alpha_1(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(1,2,1)=alpha_1(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,2,1)=alpha_1(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(0,2,1)=alpha_1(0,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         const_rp = rp(2)-(ra(2)+rab(2))
         const_ra = (-ra(2)+rp(2))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,0,0)=alpha_2(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,0,1)=alpha_2(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(0,0,1)=alpha_2(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,1,0)=alpha_2(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,1,0)=alpha_2(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,1,1)=alpha_2(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(1,1,1)=alpha_2(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,1,1)=alpha_2(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(0,1,1)=alpha_2(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,2,0)=alpha_2(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,2,0)=alpha_2(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,2,0)=alpha_2(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(3,2,1)=alpha_2(3,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(2,2,1)=alpha_2(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,2,1)=alpha_2(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(1,2,1)=alpha_2(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,2,1)=alpha_2(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(0,2,1)=alpha_2(0,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         const_rp = rp(3)-(ra(3)+rab(3))
         const_ra = (-ra(3)+rp(3))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,0,0)=alpha_3(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,0,1)=alpha_3(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(0,0,1)=alpha_3(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,1,0)=alpha_3(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,1,0)=alpha_3(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,1,1)=alpha_3(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(1,1,1)=alpha_3(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,1,1)=alpha_3(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(0,1,1)=alpha_3(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,2,0)=alpha_3(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,2,0)=alpha_3(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,2,0)=alpha_3(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(3,2,1)=alpha_3(3,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(2,2,1)=alpha_3(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,2,1)=alpha_3(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(1,2,1)=alpha_3(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,2,1)=alpha_3(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(0,2,1)=alpha_3(0,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0

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

                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
                coef_tyz(1,0,0,0)=coef_tyz(1,0,0,0)+coef_ttz(1,0)*alpha_2(lyp,0,0)
                coef_tyz(2,0,0,0)=coef_tyz(2,0,0,0)+coef_ttz(2,0)*alpha_2(lyp,0,0)
                coef_tyz(0,1,0,0)=coef_tyz(0,1,0,0)+coef_ttz(0,1)*alpha_2(lyp,0,0)
                coef_tyz(1,1,0,0)=coef_tyz(1,1,0,0)+coef_ttz(1,1)*alpha_2(lyp,0,0)
                coef_tyz(2,1,0,0)=coef_tyz(2,1,0,0)+coef_ttz(2,1)*alpha_2(lyp,0,0)
                coef_tyz(0,0,1,0)=coef_tyz(0,0,1,0)+coef_ttz(0,0)*alpha_2(lyp,1,0)
                coef_tyz(1,0,1,0)=coef_tyz(1,0,1,0)+coef_ttz(1,0)*alpha_2(lyp,1,0)
                coef_tyz(0,1,1,0)=coef_tyz(0,1,1,0)+coef_ttz(0,1)*alpha_2(lyp,1,0)
                coef_tyz(1,1,1,0)=coef_tyz(1,1,1,0)+coef_ttz(1,1)*alpha_2(lyp,1,0)
                coef_tyz(0,0,2,0)=coef_tyz(0,0,2,0)+coef_ttz(0,0)*alpha_2(lyp,2,0)
                coef_tyz(0,1,2,0)=coef_tyz(0,1,2,0)+coef_ttz(0,1)*alpha_2(lyp,2,0)
                coef_tyz(0,0,0,1)=coef_tyz(0,0,0,1)+coef_ttz(0,0)*alpha_2(lyp,0,1)
                coef_tyz(1,0,0,1)=coef_tyz(1,0,0,1)+coef_ttz(1,0)*alpha_2(lyp,0,1)
                coef_tyz(2,0,0,1)=coef_tyz(2,0,0,1)+coef_ttz(2,0)*alpha_2(lyp,0,1)
                coef_tyz(0,0,1,1)=coef_tyz(0,0,1,1)+coef_ttz(0,0)*alpha_2(lyp,1,1)
                coef_tyz(1,0,1,1)=coef_tyz(1,0,1,1)+coef_ttz(1,0)*alpha_2(lyp,1,1)
                coef_tyz(0,0,2,1)=coef_tyz(0,0,2,1)+coef_ttz(0,0)*alpha_2(lyp,2,1)
       ENDDO


              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

              DO lxa=la_min_local,2

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

              DO lxa=MAX(la_min_local-0-1,0),1

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,2,0)
                ico=8

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,0,0)


             ENDDO
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

              DO lxa=la_min_local,2

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,0)

             ENDDO

             lxb = 0
             !jco=coset(0,1,0)
             jco=3

              DO lxa=MAX(la_min_local-0-1,0),1

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,0,0)

             ENDDO

             lxb = 0
             !jco=coset(0,1,0)
             jco=3

                lxa = 0
                !ico=coset(0,2,0)
                ico=8

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,1)*alpha_3(lzp,0,0)


              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

              DO lxa=MAX(la_min_local-1-0,0),1

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,0)

             ENDDO

             ENDDO
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,1,1)
                ico=9

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,1,0)


             ENDDO
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

              DO lxa=MAX(la_min_local-1-0,0),1

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,1,0)

             ENDDO

             lxb = 0
             !jco=coset(0,1,0)
             jco=3

                lxa = 0
                !ico=coset(0,1,1)
                ico=9

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,1,0)


              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,0,2)
                ico=10

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,2,0)


             ENDDO
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

                lxa = 0
                !ico=coset(0,0,2)
                ico=10

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,2,0)


             lxb = 0
             !jco=coset(0,0,1)
             jco=4

              DO lxa=la_min_local,2

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,1)

             ENDDO

             lxb = 0
             !jco=coset(0,0,1)
             jco=4

              DO lxa=MAX(la_min_local-0-1,0),1

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,1)

             ENDDO

             lxb = 0
             !jco=coset(0,0,1)
             jco=4

                lxa = 0
                !ico=coset(0,2,0)
                ico=8

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,0,1)


             lxb = 0
             !jco=coset(0,0,1)
             jco=4

              DO lxa=MAX(la_min_local-1-0,0),1

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,1)

             ENDDO

             lxb = 0
             !jco=coset(0,0,1)
             jco=4

                lxa = 0
                !ico=coset(0,1,1)
                ico=9

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,1,1)


             lxb = 0
             !jco=coset(0,0,1)
             jco=4

                lxa = 0
                !ico=coset(0,0,2)
                ico=10

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,2,1)



   ENDDO

    END SUBROUTINE xyz_to_vab_2_1



          ! Fully unrolled version
    ! *****************************************************************************
    SUBROUTINE xyz_to_vab_2_2 (prefactor, coef_xyz, & 
                           lp, la_max_local, lb_max_local, rp, ra, rab, vab, coset, &
                           la_min_local, lb_min_local, maxl, lvab, hvab)
       

       USE kinds,                        ONLY: dp

       ! PARAMETERS
       REAL(KIND=dp), INTENT(IN) :: prefactor
       INTEGER, INTENT(IN) :: lp, la_max_local, lb_max_local, &
                              lb_min_local, la_min_local, maxl, &
                              lvab,hvab

       INTEGER,  DIMENSION(-1:maxl,-1:maxl,-1:maxl), &
                        INTENT(IN)     :: coset

       REAL(KIND=dp),  DIMENSION(((lp+1)*(lp+2)*(lp+3))/6), &
                        INTENT(INOUT) :: coef_xyz
       REAL(KIND=dp),  DIMENSION(lvab,hvab), INTENT(INOUT) :: vab
       REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: rp, ra, rab
   
       ! LOCAL
       INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp, iaxis, l, k, lya,lyb,lza,lzb
       REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb,const_rp,const_ra
   

    REAL(KIND=dp) :: alpha_1(0:(2+2),0:2,0:2)
    REAL(KIND=dp) :: alpha_2(0:(2+2),0:2,0:2)
    REAL(KIND=dp) :: alpha_3(0:(2+2),0:2,0:2)
    REAL(KIND=dp) :: coef_ttz(0:2,0:2)
    REAL(KIND=dp) :: coef_tyz(0:2,0:2,0:2,0:2)
    coef_xyz=coef_xyz*prefactor

  alpha_1(:,:,:)=0.0_dp
  alpha_2(:,:,:)=0.0_dp
  alpha_3(:,:,:)=0.0_dp

         const_rp = rp(1)-(ra(1)+rab(1))
         const_ra = (-ra(1)+rp(1))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,0,0)=alpha_1(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,0,1)=alpha_1(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(0,0,1)=alpha_1(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,0,2)=alpha_1(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_1(1,0,2)=alpha_1(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_1(0,0,2)=alpha_1(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,1,0)=alpha_1(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,1,0)=alpha_1(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,1,1)=alpha_1(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(1,1,1)=alpha_1(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,1,1)=alpha_1(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(0,1,1)=alpha_1(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(3,1,2)=alpha_1(3,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_1(2,1,2)=alpha_1(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_1(1,1,2)=alpha_1(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,1,2)=alpha_1(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_1(1,1,2)=alpha_1(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_1(0,1,2)=alpha_1(0,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,2,0)=alpha_1(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,2,0)=alpha_1(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,2,0)=alpha_1(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(3,2,1)=alpha_1(3,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(2,2,1)=alpha_1(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,2,1)=alpha_1(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(1,2,1)=alpha_1(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,2,1)=alpha_1(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(0,2,1)=alpha_1(0,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(4,2,2)=alpha_1(4,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_1(3,2,2)=alpha_1(3,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_1(2,2,2)=alpha_1(2,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(3,2,2)=alpha_1(3,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_1(2,2,2)=alpha_1(2,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_1(1,2,2)=alpha_1(1,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,2,2)=alpha_1(2,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_1(1,2,2)=alpha_1(1,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_1(0,2,2)=alpha_1(0,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         const_rp = rp(2)-(ra(2)+rab(2))
         const_ra = (-ra(2)+rp(2))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,0,0)=alpha_2(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,0,1)=alpha_2(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(0,0,1)=alpha_2(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,0,2)=alpha_2(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_2(1,0,2)=alpha_2(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_2(0,0,2)=alpha_2(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,1,0)=alpha_2(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,1,0)=alpha_2(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,1,1)=alpha_2(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(1,1,1)=alpha_2(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,1,1)=alpha_2(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(0,1,1)=alpha_2(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(3,1,2)=alpha_2(3,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_2(2,1,2)=alpha_2(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_2(1,1,2)=alpha_2(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,1,2)=alpha_2(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_2(1,1,2)=alpha_2(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_2(0,1,2)=alpha_2(0,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,2,0)=alpha_2(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,2,0)=alpha_2(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,2,0)=alpha_2(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(3,2,1)=alpha_2(3,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(2,2,1)=alpha_2(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,2,1)=alpha_2(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(1,2,1)=alpha_2(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,2,1)=alpha_2(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(0,2,1)=alpha_2(0,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(4,2,2)=alpha_2(4,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_2(3,2,2)=alpha_2(3,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_2(2,2,2)=alpha_2(2,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(3,2,2)=alpha_2(3,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_2(2,2,2)=alpha_2(2,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_2(1,2,2)=alpha_2(1,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,2,2)=alpha_2(2,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_2(1,2,2)=alpha_2(1,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_2(0,2,2)=alpha_2(0,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         const_rp = rp(3)-(ra(3)+rab(3))
         const_ra = (-ra(3)+rp(3))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,0,0)=alpha_3(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,0,1)=alpha_3(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(0,0,1)=alpha_3(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,0,2)=alpha_3(2,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_3(1,0,2)=alpha_3(1,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_3(0,0,2)=alpha_3(0,0,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,1,0)=alpha_3(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,1,0)=alpha_3(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,1,1)=alpha_3(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(1,1,1)=alpha_3(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,1,1)=alpha_3(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(0,1,1)=alpha_3(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(3,1,2)=alpha_3(3,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_3(2,1,2)=alpha_3(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_3(1,1,2)=alpha_3(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,1,2)=alpha_3(2,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_3(1,1,2)=alpha_3(1,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_3(0,1,2)=alpha_3(0,1,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,2,0)=alpha_3(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,2,0)=alpha_3(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,2,0)=alpha_3(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(3,2,1)=alpha_3(3,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(2,2,1)=alpha_3(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,2,1)=alpha_3(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(1,2,1)=alpha_3(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,2,1)=alpha_3(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(0,2,1)=alpha_3(0,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(4,2,2)=alpha_3(4,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_3(3,2,2)=alpha_3(3,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_3(2,2,2)=alpha_3(2,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(3,2,2)=alpha_3(3,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_3(2,2,2)=alpha_3(2,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_3(1,2,2)=alpha_3(1,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,2,2)=alpha_3(2,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(2,dp)
                 b=b*const_rp
             alpha_3(1,2,2)=alpha_3(1,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(1,dp)/REAL(2,dp)
                 b=b*const_rp
             alpha_3(0,2,2)=alpha_3(0,2,2)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0

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

                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
                coef_tyz(1,0,0,0)=coef_tyz(1,0,0,0)+coef_ttz(1,0)*alpha_2(lyp,0,0)
                coef_tyz(2,0,0,0)=coef_tyz(2,0,0,0)+coef_ttz(2,0)*alpha_2(lyp,0,0)
                coef_tyz(0,1,0,0)=coef_tyz(0,1,0,0)+coef_ttz(0,1)*alpha_2(lyp,0,0)
                coef_tyz(1,1,0,0)=coef_tyz(1,1,0,0)+coef_ttz(1,1)*alpha_2(lyp,0,0)
                coef_tyz(2,1,0,0)=coef_tyz(2,1,0,0)+coef_ttz(2,1)*alpha_2(lyp,0,0)
                coef_tyz(0,2,0,0)=coef_tyz(0,2,0,0)+coef_ttz(0,2)*alpha_2(lyp,0,0)
                coef_tyz(1,2,0,0)=coef_tyz(1,2,0,0)+coef_ttz(1,2)*alpha_2(lyp,0,0)
                coef_tyz(2,2,0,0)=coef_tyz(2,2,0,0)+coef_ttz(2,2)*alpha_2(lyp,0,0)
                coef_tyz(0,0,1,0)=coef_tyz(0,0,1,0)+coef_ttz(0,0)*alpha_2(lyp,1,0)
                coef_tyz(1,0,1,0)=coef_tyz(1,0,1,0)+coef_ttz(1,0)*alpha_2(lyp,1,0)
                coef_tyz(0,1,1,0)=coef_tyz(0,1,1,0)+coef_ttz(0,1)*alpha_2(lyp,1,0)
                coef_tyz(1,1,1,0)=coef_tyz(1,1,1,0)+coef_ttz(1,1)*alpha_2(lyp,1,0)
                coef_tyz(0,2,1,0)=coef_tyz(0,2,1,0)+coef_ttz(0,2)*alpha_2(lyp,1,0)
                coef_tyz(1,2,1,0)=coef_tyz(1,2,1,0)+coef_ttz(1,2)*alpha_2(lyp,1,0)
                coef_tyz(0,0,2,0)=coef_tyz(0,0,2,0)+coef_ttz(0,0)*alpha_2(lyp,2,0)
                coef_tyz(0,1,2,0)=coef_tyz(0,1,2,0)+coef_ttz(0,1)*alpha_2(lyp,2,0)
                coef_tyz(0,2,2,0)=coef_tyz(0,2,2,0)+coef_ttz(0,2)*alpha_2(lyp,2,0)
                coef_tyz(0,0,0,1)=coef_tyz(0,0,0,1)+coef_ttz(0,0)*alpha_2(lyp,0,1)
                coef_tyz(1,0,0,1)=coef_tyz(1,0,0,1)+coef_ttz(1,0)*alpha_2(lyp,0,1)
                coef_tyz(2,0,0,1)=coef_tyz(2,0,0,1)+coef_ttz(2,0)*alpha_2(lyp,0,1)
                coef_tyz(0,1,0,1)=coef_tyz(0,1,0,1)+coef_ttz(0,1)*alpha_2(lyp,0,1)
                coef_tyz(1,1,0,1)=coef_tyz(1,1,0,1)+coef_ttz(1,1)*alpha_2(lyp,0,1)
                coef_tyz(2,1,0,1)=coef_tyz(2,1,0,1)+coef_ttz(2,1)*alpha_2(lyp,0,1)
                coef_tyz(0,0,1,1)=coef_tyz(0,0,1,1)+coef_ttz(0,0)*alpha_2(lyp,1,1)
                coef_tyz(1,0,1,1)=coef_tyz(1,0,1,1)+coef_ttz(1,0)*alpha_2(lyp,1,1)
                coef_tyz(0,1,1,1)=coef_tyz(0,1,1,1)+coef_ttz(0,1)*alpha_2(lyp,1,1)
                coef_tyz(1,1,1,1)=coef_tyz(1,1,1,1)+coef_ttz(1,1)*alpha_2(lyp,1,1)
                coef_tyz(0,0,2,1)=coef_tyz(0,0,2,1)+coef_ttz(0,0)*alpha_2(lyp,2,1)
                coef_tyz(0,1,2,1)=coef_tyz(0,1,2,1)+coef_ttz(0,1)*alpha_2(lyp,2,1)
                coef_tyz(0,0,0,2)=coef_tyz(0,0,0,2)+coef_ttz(0,0)*alpha_2(lyp,0,2)
                coef_tyz(1,0,0,2)=coef_tyz(1,0,0,2)+coef_ttz(1,0)*alpha_2(lyp,0,2)
                coef_tyz(2,0,0,2)=coef_tyz(2,0,0,2)+coef_ttz(2,0)*alpha_2(lyp,0,2)
                coef_tyz(0,0,1,2)=coef_tyz(0,0,1,2)+coef_ttz(0,0)*alpha_2(lyp,1,2)
                coef_tyz(1,0,1,2)=coef_tyz(1,0,1,2)+coef_ttz(1,0)*alpha_2(lyp,1,2)
                coef_tyz(0,0,2,2)=coef_tyz(0,0,2,2)+coef_ttz(0,0)*alpha_2(lyp,2,2)
       ENDDO


              DO lxb=lb_min_local,2
             jco=coset(lxb,0,0)

              DO lxa=la_min_local,2

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO
              DO lxb=lb_min_local,2
             jco=coset(lxb,0,0)

              DO lxa=MAX(la_min_local-0-1,0),1

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO
              DO lxb=lb_min_local,2
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,2,0)
                ico=8

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,0,0)


             ENDDO
              DO lxb=MAX(lb_min_local-0-1,0),1
             jco=coset(lxb,1,0)

              DO lxa=la_min_local,2

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO
              DO lxb=MAX(lb_min_local-0-1,0),1
             jco=coset(lxb,1,0)

              DO lxa=MAX(la_min_local-0-1,0),1

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO
              DO lxb=MAX(lb_min_local-0-1,0),1
             jco=coset(lxb,1,0)

                lxa = 0
                !ico=coset(0,2,0)
                ico=8

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,1)*alpha_3(lzp,0,0)


             ENDDO
             lxb = 0
             !jco=coset(0,2,0)
             jco=8

              DO lxa=la_min_local,2

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,0,0)

             ENDDO

             lxb = 0
             !jco=coset(0,2,0)
             jco=8

              DO lxa=MAX(la_min_local-0-1,0),1

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,2)*alpha_3(lzp,0,0)

             ENDDO

             lxb = 0
             !jco=coset(0,2,0)
             jco=8

                lxa = 0
                !ico=coset(0,2,0)
                ico=8

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,2)*alpha_3(lzp,0,0)


              DO lxb=lb_min_local,2
             jco=coset(lxb,0,0)

              DO lxa=MAX(la_min_local-1-0,0),1

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,0)

             ENDDO

             ENDDO
              DO lxb=lb_min_local,2
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,1,1)
                ico=9

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,1,0)


             ENDDO
              DO lxb=MAX(lb_min_local-0-1,0),1
             jco=coset(lxb,1,0)

              DO lxa=MAX(la_min_local-1-0,0),1

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,1,0)

             ENDDO

             ENDDO
              DO lxb=MAX(lb_min_local-0-1,0),1
             jco=coset(lxb,1,0)

                lxa = 0
                !ico=coset(0,1,1)
                ico=9

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,1,0)


             ENDDO
             lxb = 0
             !jco=coset(0,2,0)
             jco=8

              DO lxa=MAX(la_min_local-1-0,0),1

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,1,0)

             ENDDO

             lxb = 0
             !jco=coset(0,2,0)
             jco=8

                lxa = 0
                !ico=coset(0,1,1)
                ico=9

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,2)*alpha_3(lzp,1,0)


              DO lxb=lb_min_local,2
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,0,2)
                ico=10

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,2,0)


             ENDDO
              DO lxb=MAX(lb_min_local-0-1,0),1
             jco=coset(lxb,1,0)

                lxa = 0
                !ico=coset(0,0,2)
                ico=10

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,2,0)


             ENDDO
             lxb = 0
             !jco=coset(0,2,0)
             jco=8

                lxa = 0
                !ico=coset(0,0,2)
                ico=10

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,2)*alpha_3(lzp,2,0)


              DO lxb=MAX(lb_min_local-1-0,0),1
             jco=coset(lxb,0,1)

              DO lxa=la_min_local,2

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,1)

             ENDDO

             ENDDO
              DO lxb=MAX(lb_min_local-1-0,0),1
             jco=coset(lxb,0,1)

              DO lxa=MAX(la_min_local-0-1,0),1

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,1)

             ENDDO

             ENDDO
              DO lxb=MAX(lb_min_local-1-0,0),1
             jco=coset(lxb,0,1)

                lxa = 0
                !ico=coset(0,2,0)
                ico=8

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,0,1)


             ENDDO
             lxb = 0
             !jco=coset(0,1,1)
             jco=9

              DO lxa=la_min_local,2

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,1)

             ENDDO

             lxb = 0
             !jco=coset(0,1,1)
             jco=9

              DO lxa=MAX(la_min_local-0-1,0),1

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,0,1)

             ENDDO

             lxb = 0
             !jco=coset(0,1,1)
             jco=9

                lxa = 0
                !ico=coset(0,2,0)
                ico=8

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,1)*alpha_3(lzp,0,1)


              DO lxb=MAX(lb_min_local-1-0,0),1
             jco=coset(lxb,0,1)

              DO lxa=MAX(la_min_local-1-0,0),1

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,1)

             ENDDO

             ENDDO
              DO lxb=MAX(lb_min_local-1-0,0),1
             jco=coset(lxb,0,1)

                lxa = 0
                !ico=coset(0,1,1)
                ico=9

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,1,1)


             ENDDO
             lxb = 0
             !jco=coset(0,1,1)
             jco=9

              DO lxa=MAX(la_min_local-1-0,0),1

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,1,1)

             ENDDO

             lxb = 0
             !jco=coset(0,1,1)
             jco=9

                lxa = 0
                !ico=coset(0,1,1)
                ico=9

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,1,1)


              DO lxb=MAX(lb_min_local-1-0,0),1
             jco=coset(lxb,0,1)

                lxa = 0
                !ico=coset(0,0,2)
                ico=10

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,2,1)


             ENDDO
             lxb = 0
             !jco=coset(0,1,1)
             jco=9

                lxa = 0
                !ico=coset(0,0,2)
                ico=10

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,2,1)


             lxb = 0
             !jco=coset(0,0,2)
             jco=10

              DO lxa=la_min_local,2

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,2)

             ENDDO

             lxb = 0
             !jco=coset(0,0,2)
             jco=10

              DO lxa=MAX(la_min_local-0-1,0),1

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,2)

             ENDDO

             lxb = 0
             !jco=coset(0,0,2)
             jco=10

                lxa = 0
                !ico=coset(0,2,0)
                ico=8

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,0,2)


             lxb = 0
             !jco=coset(0,0,2)
             jco=10

              DO lxa=MAX(la_min_local-1-0,0),1

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,2)

             ENDDO

             lxb = 0
             !jco=coset(0,0,2)
             jco=10

                lxa = 0
                !ico=coset(0,1,1)
                ico=9

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,1,2)


             lxb = 0
             !jco=coset(0,0,2)
             jco=10

                lxa = 0
                !ico=coset(0,0,2)
                ico=10

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,2,2)



   ENDDO

    END SUBROUTINE xyz_to_vab_2_2



          ! Partial unrolled version
    ! *****************************************************************************
    SUBROUTINE xyz_to_vab_2_3 (prefactor, coef_xyz, & 
                           lp, la_max_local, lb_max_local, rp, ra, rab, vab, coset,&
                           la_min_local, lb_min_local, maxl, lvab, hvab)
 

      USE kinds,                        ONLY: dp
       ! PARAMETERS
       INTEGER, INTENT(IN) :: lp, la_max_local, lb_max_local, & 
                               maxl, lb_min_local, la_min_local, &
                               lvab,hvab

       INTEGER,  DIMENSION(-1:maxl,-1:maxl,-1:maxl), &
                       INTENT(IN)     :: coset

       REAL(KIND=dp), INTENT(IN) :: prefactor
       
       REAL(KIND=dp),  DIMENSION(((lp+1)*(lp+2)*(lp+3))/6),&
                         INTENT(INOUT) :: coef_xyz
       REAL(KIND=dp), DIMENSION(lvab,hvab), INTENT(INOUT) :: vab
       REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: rp, ra, rab
   
       ! LOCAL
       INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp, iaxis, l, k, lya,lyb,lza,lzb
       REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb,const_rp,const_ra
 

    REAL(KIND=dp) :: alpha(0:(2+3),0:2,0:3,3)
    REAL(KIND=dp) :: coef_ttz(0:2,0:3)
    REAL(KIND=dp) :: coef_tyz(0:2,0:3,0:2,0:3)
    
    coef_xyz=coef_xyz*prefactor

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

    END SUBROUTINE xyz_to_vab_2_3

          ! Partial unrolled version
    ! *****************************************************************************
    SUBROUTINE xyz_to_vab_2_4 (prefactor, coef_xyz, & 
                           lp, la_max_local, lb_max_local, rp, ra, rab, vab, coset,&
                           la_min_local, lb_min_local, maxl, lvab, hvab)
 

      USE kinds,                        ONLY: dp
       ! PARAMETERS
       INTEGER, INTENT(IN) :: lp, la_max_local, lb_max_local, & 
                               maxl, lb_min_local, la_min_local, &
                               lvab,hvab

       INTEGER,  DIMENSION(-1:maxl,-1:maxl,-1:maxl), &
                       INTENT(IN)     :: coset

       REAL(KIND=dp), INTENT(IN) :: prefactor
       
       REAL(KIND=dp),  DIMENSION(((lp+1)*(lp+2)*(lp+3))/6),&
                         INTENT(INOUT) :: coef_xyz
       REAL(KIND=dp), DIMENSION(lvab,hvab), INTENT(INOUT) :: vab
       REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: rp, ra, rab
   
       ! LOCAL
       INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp, iaxis, l, k, lya,lyb,lza,lzb
       REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb,const_rp,const_ra
 

    REAL(KIND=dp) :: alpha(0:(2+4),0:2,0:4,3)
    REAL(KIND=dp) :: coef_ttz(0:2,0:4)
    REAL(KIND=dp) :: coef_tyz(0:2,0:4,0:2,0:4)
    
    coef_xyz=coef_xyz*prefactor

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

    END SUBROUTINE xyz_to_vab_2_4

          ! Fully unrolled version
    ! *****************************************************************************
    SUBROUTINE xyz_to_vab_3_0 (prefactor, coef_xyz, & 
                           lp, la_max_local, lb_max_local, rp, ra, rab, vab, coset, &
                           la_min_local, lb_min_local, maxl, lvab, hvab)
       

       USE kinds,                        ONLY: dp

       ! PARAMETERS
       REAL(KIND=dp), INTENT(IN) :: prefactor
       INTEGER, INTENT(IN) :: lp, la_max_local, lb_max_local, &
                              lb_min_local, la_min_local, maxl, &
                              lvab,hvab

       INTEGER,  DIMENSION(-1:maxl,-1:maxl,-1:maxl), &
                        INTENT(IN)     :: coset

       REAL(KIND=dp),  DIMENSION(((lp+1)*(lp+2)*(lp+3))/6), &
                        INTENT(INOUT) :: coef_xyz
       REAL(KIND=dp),  DIMENSION(lvab,hvab), INTENT(INOUT) :: vab
       REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: rp, ra, rab
   
       ! LOCAL
       INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp, iaxis, l, k, lya,lyb,lza,lzb
       REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb,const_rp,const_ra
   

    REAL(KIND=dp) :: alpha_1(0:(3+0),0:3,0:0)
    REAL(KIND=dp) :: alpha_2(0:(3+0),0:3,0:0)
    REAL(KIND=dp) :: alpha_3(0:(3+0),0:3,0:0)
    REAL(KIND=dp) :: coef_ttz(0:3,0:0)
    REAL(KIND=dp) :: coef_tyz(0:3,0:0,0:3,0:0)
    coef_xyz=coef_xyz*prefactor

  alpha_1(:,:,:)=0.0_dp
  alpha_2(:,:,:)=0.0_dp
  alpha_3(:,:,:)=0.0_dp

         const_ra = (-ra(1)+rp(1))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,0,0)=alpha_1(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,1,0)=alpha_1(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,1,0)=alpha_1(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,2,0)=alpha_1(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,2,0)=alpha_1(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,2,0)=alpha_1(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(3,3,0)=alpha_1(3,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(3,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,3,0)=alpha_1(2,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,3,0)=alpha_1(1,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(3,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,3,0)=alpha_1(0,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         const_ra = (-ra(2)+rp(2))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,0,0)=alpha_2(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,1,0)=alpha_2(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,1,0)=alpha_2(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,2,0)=alpha_2(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,2,0)=alpha_2(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,2,0)=alpha_2(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(3,3,0)=alpha_2(3,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(3,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,3,0)=alpha_2(2,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,3,0)=alpha_2(1,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(3,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,3,0)=alpha_2(0,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         const_ra = (-ra(3)+rp(3))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,0,0)=alpha_3(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,1,0)=alpha_3(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,1,0)=alpha_3(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,2,0)=alpha_3(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,2,0)=alpha_3(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,2,0)=alpha_3(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(3,3,0)=alpha_3(3,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(3,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,3,0)=alpha_3(2,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,3,0)=alpha_3(1,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(3,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,3,0)=alpha_3(0,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0

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

                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
                coef_tyz(1,0,0,0)=coef_tyz(1,0,0,0)+coef_ttz(1,0)*alpha_2(lyp,0,0)
                coef_tyz(2,0,0,0)=coef_tyz(2,0,0,0)+coef_ttz(2,0)*alpha_2(lyp,0,0)
                coef_tyz(3,0,0,0)=coef_tyz(3,0,0,0)+coef_ttz(3,0)*alpha_2(lyp,0,0)
                coef_tyz(0,0,1,0)=coef_tyz(0,0,1,0)+coef_ttz(0,0)*alpha_2(lyp,1,0)
                coef_tyz(1,0,1,0)=coef_tyz(1,0,1,0)+coef_ttz(1,0)*alpha_2(lyp,1,0)
                coef_tyz(2,0,1,0)=coef_tyz(2,0,1,0)+coef_ttz(2,0)*alpha_2(lyp,1,0)
                coef_tyz(0,0,2,0)=coef_tyz(0,0,2,0)+coef_ttz(0,0)*alpha_2(lyp,2,0)
                coef_tyz(1,0,2,0)=coef_tyz(1,0,2,0)+coef_ttz(1,0)*alpha_2(lyp,2,0)
                coef_tyz(0,0,3,0)=coef_tyz(0,0,3,0)+coef_ttz(0,0)*alpha_2(lyp,3,0)
       ENDDO


             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=la_min_local,3

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,0)

             ENDDO

             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-0-1,0),2

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,0)

             ENDDO

             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-0-2,0),1

                ico=coset(lxa,2,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,0,0)

             ENDDO

             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,3,0)
                ico=17

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,3,0)*alpha_3(lzp,0,0)


             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-1-0,0),2

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,0)

             ENDDO

             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-1-1,0),1

                ico=coset(lxa,1,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,1,0)

             ENDDO

             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,2,1)
                ico=18

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,1,0)


             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-2-0,0),1

                ico=coset(lxa,0,2)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,2,0)

             ENDDO

             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,1,2)
                ico=19

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,2,0)


             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,0,3)
                ico=20

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,3,0)



   ENDDO

    END SUBROUTINE xyz_to_vab_3_0



          ! Fully unrolled version
    ! *****************************************************************************
    SUBROUTINE xyz_to_vab_3_1 (prefactor, coef_xyz, & 
                           lp, la_max_local, lb_max_local, rp, ra, rab, vab, coset, &
                           la_min_local, lb_min_local, maxl, lvab, hvab)
       

       USE kinds,                        ONLY: dp

       ! PARAMETERS
       REAL(KIND=dp), INTENT(IN) :: prefactor
       INTEGER, INTENT(IN) :: lp, la_max_local, lb_max_local, &
                              lb_min_local, la_min_local, maxl, &
                              lvab,hvab

       INTEGER,  DIMENSION(-1:maxl,-1:maxl,-1:maxl), &
                        INTENT(IN)     :: coset

       REAL(KIND=dp),  DIMENSION(((lp+1)*(lp+2)*(lp+3))/6), &
                        INTENT(INOUT) :: coef_xyz
       REAL(KIND=dp),  DIMENSION(lvab,hvab), INTENT(INOUT) :: vab
       REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: rp, ra, rab
   
       ! LOCAL
       INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp, iaxis, l, k, lya,lyb,lza,lzb
       REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb,const_rp,const_ra
   

    REAL(KIND=dp) :: alpha_1(0:(3+1),0:3,0:1)
    REAL(KIND=dp) :: alpha_2(0:(3+1),0:3,0:1)
    REAL(KIND=dp) :: alpha_3(0:(3+1),0:3,0:1)
    REAL(KIND=dp) :: coef_ttz(0:3,0:1)
    REAL(KIND=dp) :: coef_tyz(0:3,0:1,0:3,0:1)
    coef_xyz=coef_xyz*prefactor

  alpha_1(:,:,:)=0.0_dp
  alpha_2(:,:,:)=0.0_dp
  alpha_3(:,:,:)=0.0_dp

         const_rp = rp(1)-(ra(1)+rab(1))
         const_ra = (-ra(1)+rp(1))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,0,0)=alpha_1(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,0,1)=alpha_1(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(0,0,1)=alpha_1(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,1,0)=alpha_1(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,1,0)=alpha_1(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,1,1)=alpha_1(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(1,1,1)=alpha_1(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,1,1)=alpha_1(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(0,1,1)=alpha_1(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,2,0)=alpha_1(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,2,0)=alpha_1(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,2,0)=alpha_1(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(3,2,1)=alpha_1(3,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(2,2,1)=alpha_1(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,2,1)=alpha_1(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(1,2,1)=alpha_1(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,2,1)=alpha_1(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(0,2,1)=alpha_1(0,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(3,3,0)=alpha_1(3,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(3,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,3,0)=alpha_1(2,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,3,0)=alpha_1(1,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(3,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,3,0)=alpha_1(0,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(4,3,1)=alpha_1(4,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(3,3,1)=alpha_1(3,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(3,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(3,3,1)=alpha_1(3,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(2,3,1)=alpha_1(2,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,3,1)=alpha_1(2,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(1,3,1)=alpha_1(1,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(3,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,3,1)=alpha_1(1,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_1(0,3,1)=alpha_1(0,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         const_rp = rp(2)-(ra(2)+rab(2))
         const_ra = (-ra(2)+rp(2))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,0,0)=alpha_2(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,0,1)=alpha_2(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(0,0,1)=alpha_2(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,1,0)=alpha_2(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,1,0)=alpha_2(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,1,1)=alpha_2(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(1,1,1)=alpha_2(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,1,1)=alpha_2(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(0,1,1)=alpha_2(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,2,0)=alpha_2(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,2,0)=alpha_2(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,2,0)=alpha_2(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(3,2,1)=alpha_2(3,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(2,2,1)=alpha_2(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,2,1)=alpha_2(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(1,2,1)=alpha_2(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,2,1)=alpha_2(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(0,2,1)=alpha_2(0,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(3,3,0)=alpha_2(3,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(3,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,3,0)=alpha_2(2,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,3,0)=alpha_2(1,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(3,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,3,0)=alpha_2(0,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(4,3,1)=alpha_2(4,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(3,3,1)=alpha_2(3,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(3,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(3,3,1)=alpha_2(3,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(2,3,1)=alpha_2(2,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,3,1)=alpha_2(2,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(1,3,1)=alpha_2(1,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(3,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,3,1)=alpha_2(1,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_2(0,3,1)=alpha_2(0,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         const_rp = rp(3)-(ra(3)+rab(3))
         const_ra = (-ra(3)+rp(3))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,0,0)=alpha_3(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,0,1)=alpha_3(1,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(0,0,1)=alpha_3(0,0,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,1,0)=alpha_3(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,1,0)=alpha_3(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,1,1)=alpha_3(2,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(1,1,1)=alpha_3(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,1,1)=alpha_3(1,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(0,1,1)=alpha_3(0,1,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,2,0)=alpha_3(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,2,0)=alpha_3(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,2,0)=alpha_3(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(3,2,1)=alpha_3(3,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(2,2,1)=alpha_3(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,2,1)=alpha_3(2,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(1,2,1)=alpha_3(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,2,1)=alpha_3(1,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(0,2,1)=alpha_3(0,2,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(3,3,0)=alpha_3(3,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(3,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,3,0)=alpha_3(2,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,3,0)=alpha_3(1,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(3,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,3,0)=alpha_3(0,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(4,3,1)=alpha_3(4,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(3,3,1)=alpha_3(3,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(3,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(3,3,1)=alpha_3(3,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(2,3,1)=alpha_3(2,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,3,1)=alpha_3(2,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(1,3,1)=alpha_3(1,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(3,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,3,1)=alpha_3(1,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! binomial_l_lxb=binomial_l_lxb*1
                 b=b*const_rp
             alpha_3(0,3,1)=alpha_3(0,3,1)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0

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

                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
                coef_tyz(1,0,0,0)=coef_tyz(1,0,0,0)+coef_ttz(1,0)*alpha_2(lyp,0,0)
                coef_tyz(2,0,0,0)=coef_tyz(2,0,0,0)+coef_ttz(2,0)*alpha_2(lyp,0,0)
                coef_tyz(3,0,0,0)=coef_tyz(3,0,0,0)+coef_ttz(3,0)*alpha_2(lyp,0,0)
                coef_tyz(0,1,0,0)=coef_tyz(0,1,0,0)+coef_ttz(0,1)*alpha_2(lyp,0,0)
                coef_tyz(1,1,0,0)=coef_tyz(1,1,0,0)+coef_ttz(1,1)*alpha_2(lyp,0,0)
                coef_tyz(2,1,0,0)=coef_tyz(2,1,0,0)+coef_ttz(2,1)*alpha_2(lyp,0,0)
                coef_tyz(3,1,0,0)=coef_tyz(3,1,0,0)+coef_ttz(3,1)*alpha_2(lyp,0,0)
                coef_tyz(0,0,1,0)=coef_tyz(0,0,1,0)+coef_ttz(0,0)*alpha_2(lyp,1,0)
                coef_tyz(1,0,1,0)=coef_tyz(1,0,1,0)+coef_ttz(1,0)*alpha_2(lyp,1,0)
                coef_tyz(2,0,1,0)=coef_tyz(2,0,1,0)+coef_ttz(2,0)*alpha_2(lyp,1,0)
                coef_tyz(0,1,1,0)=coef_tyz(0,1,1,0)+coef_ttz(0,1)*alpha_2(lyp,1,0)
                coef_tyz(1,1,1,0)=coef_tyz(1,1,1,0)+coef_ttz(1,1)*alpha_2(lyp,1,0)
                coef_tyz(2,1,1,0)=coef_tyz(2,1,1,0)+coef_ttz(2,1)*alpha_2(lyp,1,0)
                coef_tyz(0,0,2,0)=coef_tyz(0,0,2,0)+coef_ttz(0,0)*alpha_2(lyp,2,0)
                coef_tyz(1,0,2,0)=coef_tyz(1,0,2,0)+coef_ttz(1,0)*alpha_2(lyp,2,0)
                coef_tyz(0,1,2,0)=coef_tyz(0,1,2,0)+coef_ttz(0,1)*alpha_2(lyp,2,0)
                coef_tyz(1,1,2,0)=coef_tyz(1,1,2,0)+coef_ttz(1,1)*alpha_2(lyp,2,0)
                coef_tyz(0,0,3,0)=coef_tyz(0,0,3,0)+coef_ttz(0,0)*alpha_2(lyp,3,0)
                coef_tyz(0,1,3,0)=coef_tyz(0,1,3,0)+coef_ttz(0,1)*alpha_2(lyp,3,0)
                coef_tyz(0,0,0,1)=coef_tyz(0,0,0,1)+coef_ttz(0,0)*alpha_2(lyp,0,1)
                coef_tyz(1,0,0,1)=coef_tyz(1,0,0,1)+coef_ttz(1,0)*alpha_2(lyp,0,1)
                coef_tyz(2,0,0,1)=coef_tyz(2,0,0,1)+coef_ttz(2,0)*alpha_2(lyp,0,1)
                coef_tyz(3,0,0,1)=coef_tyz(3,0,0,1)+coef_ttz(3,0)*alpha_2(lyp,0,1)
                coef_tyz(0,0,1,1)=coef_tyz(0,0,1,1)+coef_ttz(0,0)*alpha_2(lyp,1,1)
                coef_tyz(1,0,1,1)=coef_tyz(1,0,1,1)+coef_ttz(1,0)*alpha_2(lyp,1,1)
                coef_tyz(2,0,1,1)=coef_tyz(2,0,1,1)+coef_ttz(2,0)*alpha_2(lyp,1,1)
                coef_tyz(0,0,2,1)=coef_tyz(0,0,2,1)+coef_ttz(0,0)*alpha_2(lyp,2,1)
                coef_tyz(1,0,2,1)=coef_tyz(1,0,2,1)+coef_ttz(1,0)*alpha_2(lyp,2,1)
                coef_tyz(0,0,3,1)=coef_tyz(0,0,3,1)+coef_ttz(0,0)*alpha_2(lyp,3,1)
       ENDDO


              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

              DO lxa=la_min_local,3

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

              DO lxa=MAX(la_min_local-0-1,0),2

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

              DO lxa=MAX(la_min_local-0-2,0),1

                ico=coset(lxa,2,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,0,0)

             ENDDO

             ENDDO
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,3,0)
                ico=17

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,3,0)*alpha_3(lzp,0,0)


             ENDDO
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

              DO lxa=la_min_local,3

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,0,0)

             ENDDO

             lxb = 0
             !jco=coset(0,1,0)
             jco=3

              DO lxa=MAX(la_min_local-0-1,0),2

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,0,0)

             ENDDO

             lxb = 0
             !jco=coset(0,1,0)
             jco=3

              DO lxa=MAX(la_min_local-0-2,0),1

                ico=coset(lxa,2,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,1)*alpha_3(lzp,0,0)

             ENDDO

             lxb = 0
             !jco=coset(0,1,0)
             jco=3

                lxa = 0
                !ico=coset(0,3,0)
                ico=17

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,3,1)*alpha_3(lzp,0,0)


              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

              DO lxa=MAX(la_min_local-1-0,0),2

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,0)

             ENDDO

             ENDDO
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

              DO lxa=MAX(la_min_local-1-1,0),1

                ico=coset(lxa,1,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,1,0)

             ENDDO

             ENDDO
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,2,1)
                ico=18

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,1,0)


             ENDDO
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

              DO lxa=MAX(la_min_local-1-0,0),2

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,1,0)

             ENDDO

             lxb = 0
             !jco=coset(0,1,0)
             jco=3

              DO lxa=MAX(la_min_local-1-1,0),1

                ico=coset(lxa,1,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,1,0)

             ENDDO

             lxb = 0
             !jco=coset(0,1,0)
             jco=3

                lxa = 0
                !ico=coset(0,2,1)
                ico=18

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,1)*alpha_3(lzp,1,0)


              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

              DO lxa=MAX(la_min_local-2-0,0),1

                ico=coset(lxa,0,2)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,2,0)

             ENDDO

             ENDDO
              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,1,2)
                ico=19

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,2,0)


             ENDDO
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

              DO lxa=MAX(la_min_local-2-0,0),1

                ico=coset(lxa,0,2)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,2,0)

             ENDDO

             lxb = 0
             !jco=coset(0,1,0)
             jco=3

                lxa = 0
                !ico=coset(0,1,2)
                ico=19

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,1)*alpha_3(lzp,2,0)


              DO lxb=lb_min_local,1
             jco=coset(lxb,0,0)

                lxa = 0
                !ico=coset(0,0,3)
                ico=20

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,3,0)


             ENDDO
             lxb = 0
             !jco=coset(0,1,0)
             jco=3

                lxa = 0
                !ico=coset(0,0,3)
                ico=20

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,1)*alpha_3(lzp,3,0)


             lxb = 0
             !jco=coset(0,0,1)
             jco=4

              DO lxa=la_min_local,3

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,1)

             ENDDO

             lxb = 0
             !jco=coset(0,0,1)
             jco=4

              DO lxa=MAX(la_min_local-0-1,0),2

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,1)

             ENDDO

             lxb = 0
             !jco=coset(0,0,1)
             jco=4

              DO lxa=MAX(la_min_local-0-2,0),1

                ico=coset(lxa,2,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,0,1)

             ENDDO

             lxb = 0
             !jco=coset(0,0,1)
             jco=4

                lxa = 0
                !ico=coset(0,3,0)
                ico=17

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,3,0)*alpha_3(lzp,0,1)


             lxb = 0
             !jco=coset(0,0,1)
             jco=4

              DO lxa=MAX(la_min_local-1-0,0),2

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,1)

             ENDDO

             lxb = 0
             !jco=coset(0,0,1)
             jco=4

              DO lxa=MAX(la_min_local-1-1,0),1

                ico=coset(lxa,1,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,1,1)

             ENDDO

             lxb = 0
             !jco=coset(0,0,1)
             jco=4

                lxa = 0
                !ico=coset(0,2,1)
                ico=18

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,1,1)


             lxb = 0
             !jco=coset(0,0,1)
             jco=4

              DO lxa=MAX(la_min_local-2-0,0),1

                ico=coset(lxa,0,2)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,2,1)

             ENDDO

             lxb = 0
             !jco=coset(0,0,1)
             jco=4

                lxa = 0
                !ico=coset(0,1,2)
                ico=19

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,2,1)


             lxb = 0
             !jco=coset(0,0,1)
             jco=4

                lxa = 0
                !ico=coset(0,0,3)
                ico=20

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,3,1)



   ENDDO

    END SUBROUTINE xyz_to_vab_3_1



          ! Partial unrolled version
    ! *****************************************************************************
    SUBROUTINE xyz_to_vab_3_2 (prefactor, coef_xyz, & 
                           lp, la_max_local, lb_max_local, rp, ra, rab, vab, coset,&
                           la_min_local, lb_min_local, maxl, lvab, hvab)
 

      USE kinds,                        ONLY: dp
       ! PARAMETERS
       INTEGER, INTENT(IN) :: lp, la_max_local, lb_max_local, & 
                               maxl, lb_min_local, la_min_local, &
                               lvab,hvab

       INTEGER,  DIMENSION(-1:maxl,-1:maxl,-1:maxl), &
                       INTENT(IN)     :: coset

       REAL(KIND=dp), INTENT(IN) :: prefactor
       
       REAL(KIND=dp),  DIMENSION(((lp+1)*(lp+2)*(lp+3))/6),&
                         INTENT(INOUT) :: coef_xyz
       REAL(KIND=dp), DIMENSION(lvab,hvab), INTENT(INOUT) :: vab
       REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: rp, ra, rab
   
       ! LOCAL
       INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp, iaxis, l, k, lya,lyb,lza,lzb
       REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb,const_rp,const_ra
 

    REAL(KIND=dp) :: alpha(0:(3+2),0:3,0:2,3)
    REAL(KIND=dp) :: coef_ttz(0:3,0:2)
    REAL(KIND=dp) :: coef_tyz(0:3,0:2,0:3,0:2)
    
    coef_xyz=coef_xyz*prefactor

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

    END SUBROUTINE xyz_to_vab_3_2

          ! Partial unrolled version
    ! *****************************************************************************
    SUBROUTINE xyz_to_vab_3_3 (prefactor, coef_xyz, & 
                           lp, la_max_local, lb_max_local, rp, ra, rab, vab, coset,&
                           la_min_local, lb_min_local, maxl, lvab, hvab)
 

      USE kinds,                        ONLY: dp
       ! PARAMETERS
       INTEGER, INTENT(IN) :: lp, la_max_local, lb_max_local, & 
                               maxl, lb_min_local, la_min_local, &
                               lvab,hvab

       INTEGER,  DIMENSION(-1:maxl,-1:maxl,-1:maxl), &
                       INTENT(IN)     :: coset

       REAL(KIND=dp), INTENT(IN) :: prefactor
       
       REAL(KIND=dp),  DIMENSION(((lp+1)*(lp+2)*(lp+3))/6),&
                         INTENT(INOUT) :: coef_xyz
       REAL(KIND=dp), DIMENSION(lvab,hvab), INTENT(INOUT) :: vab
       REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: rp, ra, rab
   
       ! LOCAL
       INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp, iaxis, l, k, lya,lyb,lza,lzb
       REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb,const_rp,const_ra
 

    REAL(KIND=dp) :: alpha(0:(3+3),0:3,0:3,3)
    REAL(KIND=dp) :: coef_ttz(0:3,0:3)
    REAL(KIND=dp) :: coef_tyz(0:3,0:3,0:3,0:3)
    
    coef_xyz=coef_xyz*prefactor

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

    END SUBROUTINE xyz_to_vab_3_3

          ! Partial unrolled version
    ! *****************************************************************************
    SUBROUTINE xyz_to_vab_3_4 (prefactor, coef_xyz, & 
                           lp, la_max_local, lb_max_local, rp, ra, rab, vab, coset,&
                           la_min_local, lb_min_local, maxl, lvab, hvab)
 

      USE kinds,                        ONLY: dp
       ! PARAMETERS
       INTEGER, INTENT(IN) :: lp, la_max_local, lb_max_local, & 
                               maxl, lb_min_local, la_min_local, &
                               lvab,hvab

       INTEGER,  DIMENSION(-1:maxl,-1:maxl,-1:maxl), &
                       INTENT(IN)     :: coset

       REAL(KIND=dp), INTENT(IN) :: prefactor
       
       REAL(KIND=dp),  DIMENSION(((lp+1)*(lp+2)*(lp+3))/6),&
                         INTENT(INOUT) :: coef_xyz
       REAL(KIND=dp), DIMENSION(lvab,hvab), INTENT(INOUT) :: vab
       REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: rp, ra, rab
   
       ! LOCAL
       INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp, iaxis, l, k, lya,lyb,lza,lzb
       REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb,const_rp,const_ra
 

    REAL(KIND=dp) :: alpha(0:(3+4),0:3,0:4,3)
    REAL(KIND=dp) :: coef_ttz(0:3,0:4)
    REAL(KIND=dp) :: coef_tyz(0:3,0:4,0:3,0:4)
    
    coef_xyz=coef_xyz*prefactor

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

    END SUBROUTINE xyz_to_vab_3_4

          ! Fully unrolled version
    ! *****************************************************************************
    SUBROUTINE xyz_to_vab_4_0 (prefactor, coef_xyz, & 
                           lp, la_max_local, lb_max_local, rp, ra, rab, vab, coset, &
                           la_min_local, lb_min_local, maxl, lvab, hvab)
       

       USE kinds,                        ONLY: dp

       ! PARAMETERS
       REAL(KIND=dp), INTENT(IN) :: prefactor
       INTEGER, INTENT(IN) :: lp, la_max_local, lb_max_local, &
                              lb_min_local, la_min_local, maxl, &
                              lvab,hvab

       INTEGER,  DIMENSION(-1:maxl,-1:maxl,-1:maxl), &
                        INTENT(IN)     :: coset

       REAL(KIND=dp),  DIMENSION(((lp+1)*(lp+2)*(lp+3))/6), &
                        INTENT(INOUT) :: coef_xyz
       REAL(KIND=dp),  DIMENSION(lvab,hvab), INTENT(INOUT) :: vab
       REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: rp, ra, rab
   
       ! LOCAL
       INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp, iaxis, l, k, lya,lyb,lza,lzb
       REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb,const_rp,const_ra
   

    REAL(KIND=dp) :: alpha_1(0:(4+0),0:4,0:0)
    REAL(KIND=dp) :: alpha_2(0:(4+0),0:4,0:0)
    REAL(KIND=dp) :: alpha_3(0:(4+0),0:4,0:0)
    REAL(KIND=dp) :: coef_ttz(0:4,0:0)
    REAL(KIND=dp) :: coef_tyz(0:4,0:0,0:4,0:0)
    coef_xyz=coef_xyz*prefactor

  alpha_1(:,:,:)=0.0_dp
  alpha_2(:,:,:)=0.0_dp
  alpha_3(:,:,:)=0.0_dp

         const_ra = (-ra(1)+rp(1))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,0,0)=alpha_1(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,1,0)=alpha_1(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,1,0)=alpha_1(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,2,0)=alpha_1(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,2,0)=alpha_1(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,2,0)=alpha_1(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(3,3,0)=alpha_1(3,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(3,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,3,0)=alpha_1(2,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,3,0)=alpha_1(1,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(3,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,3,0)=alpha_1(0,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(4,4,0)=alpha_1(4,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(4,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(3,4,0)=alpha_1(3,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(3,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(2,4,0)=alpha_1(2,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)/REAL(3,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(1,4,0)=alpha_1(1,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(4,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_1(0,4,0)=alpha_1(0,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         const_ra = (-ra(2)+rp(2))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,0,0)=alpha_2(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,1,0)=alpha_2(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,1,0)=alpha_2(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,2,0)=alpha_2(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,2,0)=alpha_2(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,2,0)=alpha_2(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(3,3,0)=alpha_2(3,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(3,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,3,0)=alpha_2(2,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,3,0)=alpha_2(1,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(3,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,3,0)=alpha_2(0,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(4,4,0)=alpha_2(4,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(4,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(3,4,0)=alpha_2(3,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(3,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(2,4,0)=alpha_2(2,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)/REAL(3,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(1,4,0)=alpha_2(1,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(4,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_2(0,4,0)=alpha_2(0,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         const_ra = (-ra(3)+rp(3))
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,0,0)=alpha_3(0,0,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,1,0)=alpha_3(1,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,1,0)=alpha_3(0,1,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,2,0)=alpha_3(2,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,2,0)=alpha_3(1,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,2,0)=alpha_3(0,2,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(3,3,0)=alpha_3(3,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(3,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,3,0)=alpha_3(2,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,3,0)=alpha_3(1,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(3,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,3,0)=alpha_3(0,3,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0
         binomial_k_lxa=1.0_dp
         a=1.0_dp
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(4,4,0)=alpha_3(4,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(4,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(3,4,0)=alpha_3(3,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(3,dp)/REAL(2,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(2,4,0)=alpha_3(2,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(2,dp)/REAL(3,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(1,4,0)=alpha_3(1,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             binomial_k_lxa=binomial_k_lxa*REAL(1,dp)/REAL(4,dp)
             a=a*const_ra
          binomial_l_lxb=1.0_dp
          b=1.0_dp
             alpha_3(0,4,0)=alpha_3(0,4,0)+ &
                               binomial_k_lxa*binomial_l_lxb*a*b
                 ! Removed here as lxb - l = 0
             ! Removed  here as lxa-k == 0

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

                coef_tyz(0,0,0,0)=coef_tyz(0,0,0,0)+coef_ttz(0,0)*alpha_2(lyp,0,0)
                coef_tyz(1,0,0,0)=coef_tyz(1,0,0,0)+coef_ttz(1,0)*alpha_2(lyp,0,0)
                coef_tyz(2,0,0,0)=coef_tyz(2,0,0,0)+coef_ttz(2,0)*alpha_2(lyp,0,0)
                coef_tyz(3,0,0,0)=coef_tyz(3,0,0,0)+coef_ttz(3,0)*alpha_2(lyp,0,0)
                coef_tyz(4,0,0,0)=coef_tyz(4,0,0,0)+coef_ttz(4,0)*alpha_2(lyp,0,0)
                coef_tyz(0,0,1,0)=coef_tyz(0,0,1,0)+coef_ttz(0,0)*alpha_2(lyp,1,0)
                coef_tyz(1,0,1,0)=coef_tyz(1,0,1,0)+coef_ttz(1,0)*alpha_2(lyp,1,0)
                coef_tyz(2,0,1,0)=coef_tyz(2,0,1,0)+coef_ttz(2,0)*alpha_2(lyp,1,0)
                coef_tyz(3,0,1,0)=coef_tyz(3,0,1,0)+coef_ttz(3,0)*alpha_2(lyp,1,0)
                coef_tyz(0,0,2,0)=coef_tyz(0,0,2,0)+coef_ttz(0,0)*alpha_2(lyp,2,0)
                coef_tyz(1,0,2,0)=coef_tyz(1,0,2,0)+coef_ttz(1,0)*alpha_2(lyp,2,0)
                coef_tyz(2,0,2,0)=coef_tyz(2,0,2,0)+coef_ttz(2,0)*alpha_2(lyp,2,0)
                coef_tyz(0,0,3,0)=coef_tyz(0,0,3,0)+coef_ttz(0,0)*alpha_2(lyp,3,0)
                coef_tyz(1,0,3,0)=coef_tyz(1,0,3,0)+coef_ttz(1,0)*alpha_2(lyp,3,0)
                coef_tyz(0,0,4,0)=coef_tyz(0,0,4,0)+coef_ttz(0,0)*alpha_2(lyp,4,0)
       ENDDO


             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=la_min_local,4

                ico=coset(lxa,0,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,0,0)

             ENDDO

             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-0-1,0),3

                ico=coset(lxa,1,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,0,0)

             ENDDO

             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-0-2,0),2

                ico=coset(lxa,2,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,0,0)

             ENDDO

             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-0-3,0),1

                ico=coset(lxa,3,0)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,3,0)*alpha_3(lzp,0,0)

             ENDDO

             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,4,0)
                ico=31

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,4,0)*alpha_3(lzp,0,0)


             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-1-0,0),3

                ico=coset(lxa,0,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,1,0)

             ENDDO

             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-1-1,0),2

                ico=coset(lxa,1,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,1,0)

             ENDDO

             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-1-2,0),1

                ico=coset(lxa,2,1)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,1,0)

             ENDDO

             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,3,1)
                ico=32

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,3,0)*alpha_3(lzp,1,0)


             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-2-0,0),2

                ico=coset(lxa,0,2)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,2,0)

             ENDDO

             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-2-1,0),1

                ico=coset(lxa,1,2)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,2,0)

             ENDDO

             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,2,2)
                ico=33

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,2,0)*alpha_3(lzp,2,0)


             lxb = 0
             !jco=coset(0,0,0)
             jco=1

              DO lxa=MAX(la_min_local-3-0,0),1

                ico=coset(lxa,0,3)

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,3,0)

             ENDDO

             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,1,3)
                ico=34

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,1,0)*alpha_3(lzp,3,0)


             lxb = 0
             !jco=coset(0,0,0)
             jco=1

                lxa = 0
                !ico=coset(0,0,4)
                ico=35

                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,0,0)*alpha_3(lzp,4,0)



   ENDDO

    END SUBROUTINE xyz_to_vab_4_0



          ! Partial unrolled version
    ! *****************************************************************************
    SUBROUTINE xyz_to_vab_4_1 (prefactor, coef_xyz, & 
                           lp, la_max_local, lb_max_local, rp, ra, rab, vab, coset,&
                           la_min_local, lb_min_local, maxl, lvab, hvab)
 

      USE kinds,                        ONLY: dp
       ! PARAMETERS
       INTEGER, INTENT(IN) :: lp, la_max_local, lb_max_local, & 
                               maxl, lb_min_local, la_min_local, &
                               lvab,hvab

       INTEGER,  DIMENSION(-1:maxl,-1:maxl,-1:maxl), &
                       INTENT(IN)     :: coset

       REAL(KIND=dp), INTENT(IN) :: prefactor
       
       REAL(KIND=dp),  DIMENSION(((lp+1)*(lp+2)*(lp+3))/6),&
                         INTENT(INOUT) :: coef_xyz
       REAL(KIND=dp), DIMENSION(lvab,hvab), INTENT(INOUT) :: vab
       REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: rp, ra, rab
   
       ! LOCAL
       INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp, iaxis, l, k, lya,lyb,lza,lzb
       REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb,const_rp,const_ra
 

    REAL(KIND=dp) :: alpha(0:(4+1),0:4,0:1,3)
    REAL(KIND=dp) :: coef_ttz(0:4,0:1)
    REAL(KIND=dp) :: coef_tyz(0:4,0:1,0:4,0:1)
    
    coef_xyz=coef_xyz*prefactor

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

    END SUBROUTINE xyz_to_vab_4_1

          ! Partial unrolled version
    ! *****************************************************************************
    SUBROUTINE xyz_to_vab_4_2 (prefactor, coef_xyz, & 
                           lp, la_max_local, lb_max_local, rp, ra, rab, vab, coset,&
                           la_min_local, lb_min_local, maxl, lvab, hvab)
 

      USE kinds,                        ONLY: dp
       ! PARAMETERS
       INTEGER, INTENT(IN) :: lp, la_max_local, lb_max_local, & 
                               maxl, lb_min_local, la_min_local, &
                               lvab,hvab

       INTEGER,  DIMENSION(-1:maxl,-1:maxl,-1:maxl), &
                       INTENT(IN)     :: coset

       REAL(KIND=dp), INTENT(IN) :: prefactor
       
       REAL(KIND=dp),  DIMENSION(((lp+1)*(lp+2)*(lp+3))/6),&
                         INTENT(INOUT) :: coef_xyz
       REAL(KIND=dp), DIMENSION(lvab,hvab), INTENT(INOUT) :: vab
       REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: rp, ra, rab
   
       ! LOCAL
       INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp, iaxis, l, k, lya,lyb,lza,lzb
       REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb,const_rp,const_ra
 

    REAL(KIND=dp) :: alpha(0:(4+2),0:4,0:2,3)
    REAL(KIND=dp) :: coef_ttz(0:4,0:2)
    REAL(KIND=dp) :: coef_tyz(0:4,0:2,0:4,0:2)
    
    coef_xyz=coef_xyz*prefactor

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

    END SUBROUTINE xyz_to_vab_4_2

          ! Partial unrolled version
    ! *****************************************************************************
    SUBROUTINE xyz_to_vab_4_3 (prefactor, coef_xyz, & 
                           lp, la_max_local, lb_max_local, rp, ra, rab, vab, coset,&
                           la_min_local, lb_min_local, maxl, lvab, hvab)
 

      USE kinds,                        ONLY: dp
       ! PARAMETERS
       INTEGER, INTENT(IN) :: lp, la_max_local, lb_max_local, & 
                               maxl, lb_min_local, la_min_local, &
                               lvab,hvab

       INTEGER,  DIMENSION(-1:maxl,-1:maxl,-1:maxl), &
                       INTENT(IN)     :: coset

       REAL(KIND=dp), INTENT(IN) :: prefactor
       
       REAL(KIND=dp),  DIMENSION(((lp+1)*(lp+2)*(lp+3))/6),&
                         INTENT(INOUT) :: coef_xyz
       REAL(KIND=dp), DIMENSION(lvab,hvab), INTENT(INOUT) :: vab
       REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: rp, ra, rab
   
       ! LOCAL
       INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp, iaxis, l, k, lya,lyb,lza,lzb
       REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb,const_rp,const_ra
 

    REAL(KIND=dp) :: alpha(0:(4+3),0:4,0:3,3)
    REAL(KIND=dp) :: coef_ttz(0:4,0:3)
    REAL(KIND=dp) :: coef_tyz(0:4,0:3,0:4,0:3)
    
    coef_xyz=coef_xyz*prefactor

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

    END SUBROUTINE xyz_to_vab_4_3

          ! Partial unrolled version
    ! *****************************************************************************
    SUBROUTINE xyz_to_vab_4_4 (prefactor, coef_xyz, & 
                           lp, la_max_local, lb_max_local, rp, ra, rab, vab, coset,&
                           la_min_local, lb_min_local, maxl, lvab, hvab)
 

      USE kinds,                        ONLY: dp
       ! PARAMETERS
       INTEGER, INTENT(IN) :: lp, la_max_local, lb_max_local, & 
                               maxl, lb_min_local, la_min_local, &
                               lvab,hvab

       INTEGER,  DIMENSION(-1:maxl,-1:maxl,-1:maxl), &
                       INTENT(IN)     :: coset

       REAL(KIND=dp), INTENT(IN) :: prefactor
       
       REAL(KIND=dp),  DIMENSION(((lp+1)*(lp+2)*(lp+3))/6),&
                         INTENT(INOUT) :: coef_xyz
       REAL(KIND=dp), DIMENSION(lvab,hvab), INTENT(INOUT) :: vab
       REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: rp, ra, rab
   
       ! LOCAL
       INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp, iaxis, l, k, lya,lyb,lza,lzb
       REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb,const_rp,const_ra
 

    REAL(KIND=dp) :: alpha(0:(4+4),0:4,0:4,3)
    REAL(KIND=dp) :: coef_ttz(0:4,0:4)
    REAL(KIND=dp) :: coef_tyz(0:4,0:4,0:4,0:4)
    
    coef_xyz=coef_xyz*prefactor

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


!!! End of cloned routines

      ! Original routine 
      SUBROUTINE xyz_to_vab(prefactor, coef_xyz, lp, la_max_local, &
                             lb_max_local, rp, ra, rab, vab, coset, &
                             la_min_local, lb_min_local, maxl, &
                             lvab, hvab)
      
              USE kinds,                        ONLY: dp
              
         implicit none
      
             ! PARAMETERS
             REAL(KIND=dp), INTENT(IN) :: prefactor
             INTEGER, INTENT(IN) :: lp, la_max_local, lb_max_local, &
                                        la_min_local, lb_min_local, &
                                        maxl, lvab, hvab
             INTEGER,  DIMENSION(-1:maxl,-1:maxl,-1:maxl), &
                       INTENT(IN)     :: coset


             REAL(KIND=dp), DIMENSION(((lp+1)*(lp+2)*(lp+3))/6), &
                             INTENT(INOUT) :: coef_xyz
             REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: rp, ra, rab
             REAL(KIND=dp), DIMENSION(lvab,hvab), INTENT(INOUT) :: vab
      
             ! LOCAL
             INTEGER :: ico,jco,lxa,lxb,lxp,lxyz,lyp,lzp, iaxis, l, k, &
                        lya,lyb,lza,lzb
             REAL(KIND=dp) :: a,b,binomial_k_lxa,binomial_l_lxb
      
          REAL(KIND=dp), DIMENSION(:,:,:,:), ALLOCATABLE :: alpha
          REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: coef_ttz
          REAL(KIND=dp), DIMENSION(:,:,:,:), ALLOCATABLE :: coef_tyz
      
      
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
          ALLOCATE(alpha(0:lp,0:la_max_local,0:lb_max_local,3))
          alpha(:,:,:,:)=0.0_dp
          DO iaxis=1,3
          DO lxa=0,la_max_local
          DO lxb=0,lb_max_local
             binomial_k_lxa=1.0_dp
             a=1.0_dp
             DO k=0,lxa
              binomial_l_lxb=1.0_dp
              b=1.0_dp
              DO l=0,lxb
                 alpha(lxa-l+lxb-k,lxa,lxb,iaxis)=alpha(lxa-l+lxb-k,lxa,lxb,iaxis)+ &
                                   binomial_k_lxa*binomial_l_lxb*a*b
                 binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
                 b=b*(rp(iaxis)-(ra(iaxis)+rab(iaxis)))
              ENDDO
              binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
              a=a*(-ra(iaxis)+rp(iaxis))
             ENDDO
          ENDDO
          ENDDO
          ENDDO
      
          !
          !   compute v_{lxa,lya,lza,lxb,lyb,lzb} given v_{lxp,lyp,lzp} and alpha(ls,lxa,lxb,1)
          !   use a three step procedure
          !
      
          ALLOCATE(coef_ttz(0:la_max_local,0:lb_max_local))
          ALLOCATE(coef_tyz(0:la_max_local,0:lb_max_local,0:la_max_local,0:lb_max_local))
          lxyz=0
          DO lzp=0,lp
             coef_tyz=0.0_dp
             DO lyp=0,lp-lzp
                coef_ttz=0.0_dp
                DO lxp=0,lp-lzp-lyp
                   lxyz=lxyz+1
                   DO lxb=0,lb_max_local
                   DO lxa=0,la_max_local
                      coef_ttz(lxa,lxb)=coef_ttz(lxa,lxb)+coef_xyz(lxyz)*alpha(lxp,lxa,lxb,1)
                   ENDDO
                   ENDDO
      
                ENDDO
      
                DO lyb=0,lb_max_local
                DO lya=0,la_max_local
                   DO lxb=0,lb_max_local-lyb
                   DO lxa=0,la_max_local-lya
                      coef_tyz(lxa,lxb,lya,lyb)=coef_tyz(lxa,lxb,lya,lyb)+coef_ttz(lxa,lxb)*alpha(lyp,lya,lyb,2)
                   ENDDO
                   ENDDO
                ENDDO
                ENDDO
      
             ENDDO
      
             DO lzb=0,lb_max_local
             DO lza=0,la_max_local
                DO lyb=0,lb_max_local-lzb
                DO lya=0,la_max_local-lza
                   DO lxb=MAX(lb_min_local-lzb-lyb,0),lb_max_local-lzb-lyb
                   jco=coset(lxb,lyb,lzb)
                   DO lxa=MAX(la_min_local-lza-lya,0),la_max_local-lza-lya
                      ico=coset(lxa,lya,lza)
                      vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,lya,lyb)*alpha(lzp,lza,lzb,3)
                   ENDDO
                   ENDDO
                ENDDO
                ENDDO
             ENDDO
             ENDDO
      
          ENDDO
      
          DEALLOCATE(coef_tyz)
          DEALLOCATE(coef_ttz)
          DEALLOCATE(alpha)
      
      
      
            end subroutine xyz_to_vab
      
      
      

