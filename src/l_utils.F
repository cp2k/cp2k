! ***************************************************************
MODULE l_utils
  USE kinds,                           ONLY: dp
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

   IMPLICIT NONE

   public :: l_info_type
   public :: init_l_info
   public :: return_l_info
   public :: destroy_l_info
   
   TYPE l_set_info_type
        integer                        :: lx_max,lxy_max,lxyz_max
        integer, pointer, dimension(:) :: ly_max,lz_max

        integer, pointer, dimension(:,:)  :: map
        real(KIND = dp), pointer, dimension(:,:) :: polx,poly,polz,alpha, &
                                             dpy,dpz
        real(KIND = dp), pointer, dimension(:)   :: pzyx
        integer,  pointer, dimension(:,:,:,:,:,:)   :: ipzyx
   END TYPE

   TYPE l_info_type
        integer :: lmax
        integer :: cmax
        integer :: nthread
        TYPE(l_set_info_type), pointer, dimension(:,:,:,:,:) :: data
   END TYPE

CONTAINS

   SUBROUTINE return_l_info(l_info,la_min,la_max,lb_min,lb_max,ithread, &
                            lx_max,lxy_max,lxyz_max,ly_max,lz_max,&
                            map,polx,poly,polz,dpy,dpz,alpha,pzyx,ipzyx,cmax)
    type(l_info_type)                        :: l_info
    integer                                  :: la_min, la_max, lb_min, &
                                                lb_max, ithread, lx_max, &
                                                lxy_max, lxyz_max
    integer, dimension(:), pointer           :: ly_max, lz_max
    integer, dimension(:, :), pointer        :: map
    real(KIND=dp), dimension(:, :), pointer  :: polx, poly, polz, dpy, dpz, &
                                                alpha
    real(KIND=dp), dimension(:), pointer     :: pzyx
    integer, dimension(:, :, :, :, :, :), &
      pointer                                :: ipzyx
    integer                                  :: cmax

      if (la_max.gt.l_info%lmax .or. lb_max.gt.l_info%lmax) then
        CALL stop_program("return_l_info","called with l too large")
      endif
      if (ithread.ge.l_info%nthread) then
        CALL stop_program("return_l_info","called with ithread too large")
      endif
      lx_max=l_info%data(la_min,la_max,lb_min,lb_max,ithread)%lx_max
      lxy_max=l_info%data(la_min,la_max,lb_min,lb_max,ithread)%lxy_max
      lxyz_max=l_info%data(la_min,la_max,lb_min,lb_max,ithread)%lxyz_max
      ly_max=>l_info%data(la_min,la_max,lb_min,lb_max,ithread)%ly_max
      lz_max=>l_info%data(la_min,la_max,lb_min,lb_max,ithread)%lz_max
      map=>l_info%data(la_min,la_max,lb_min,lb_max,ithread)%map
      polx=>l_info%data(la_min,la_max,lb_min,lb_max,ithread)%polx
      poly=>l_info%data(la_min,la_max,lb_min,lb_max,ithread)%poly
      polz=>l_info%data(la_min,la_max,lb_min,lb_max,ithread)%polz
      dpy=>l_info%data(la_min,la_max,lb_min,lb_max,ithread)%dpy
      dpz=>l_info%data(la_min,la_max,lb_min,lb_max,ithread)%dpz
      alpha=>l_info%data(la_min,la_max,lb_min,lb_max,ithread)%alpha
      pzyx=>l_info%data(la_min,la_max,lb_min,lb_max,ithread)%pzyx
      ipzyx=>l_info%data(la_min,la_max,lb_min,lb_max,ithread)%ipzyx
      cmax=l_info%cmax
   END SUBROUTINE
! needs a nthread agrument because it provides pre-allocated datastreactures for the collocate routines
! hence, in some way, these data structures need the 'thread-private' attribute
! really not so elegant...
! the nthread argument is optional, the right number of threads should be determined automatically
   SUBROUTINE init_l_info(l_info,lmax,cmax,nthread)
    type(l_info_type)                        :: l_info
    integer                                  :: lmax, cmax
    integer, optional                        :: nthread

    integer                                  :: handle, ithread, la_max, &
                                                la_min, lb_max, lb_min, &
                                                nthread_l
!$  integer                                  :: omp_get_max_threads

      CALL timeset("init_l_info","I","",handle)
      IF (PRESENT(nthread)) THEN
         nthread_l=nthread
      ELSE
         nthread_l=1
!$       nthread_l=omp_get_max_threads()
      ENDIF

      l_info%lmax=lmax
      l_info%cmax=cmax
      l_info%nthread=nthread_l
      allocate(l_info%data(0:lmax,0:lmax,0:lmax,0:lmax,0:nthread_l-1))
      DO ithread=0,nthread_l-1
      DO lb_max=0,lmax
      DO lb_min=0,lb_max
      DO la_max=0,lmax
      DO la_min=0,la_max
         call init_l_set_info_type(la_min,la_max,lb_min,lb_max,cmax, &
                                   l_info%data(la_min,la_max,lb_min,lb_max,ithread))       
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      CALL timestop(0.0_dp,handle)
   END SUBROUTINE

   SUBROUTINE destroy_l_info(l_info)
    type(l_info_type)                        :: l_info

    integer                                  :: istat, ithread, la_max, &
                                                la_min, lb_max, lb_min

      DO ithread=0,l_info%nthread-1
      DO lb_max=0,l_info%lmax
      DO lb_min=0,lb_max
      DO la_max=0,l_info%lmax
      DO la_min=0,la_max
         call destroy_l_set_info_type(l_info%data(la_min,la_max,lb_min,lb_max,ithread))
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      deallocate(l_info%data,STAT=istat)
      IF (istat /= 0) CALL stop_memory("destroy_l_info","l_info%data")
   END SUBROUTINE

   SUBROUTINE init_l_set_info_type(la_min,la_max,lb_min,lb_max,cmax,l_set_info)
    integer                                  :: la_min, la_max, lb_min, &
                                                lb_max, cmax
    TYPE(l_set_info_type)                    :: l_set_info

    integer                                  :: coef_max, lx, lxa, lxb, lxy, &
                                                lxyz, lya, lyb, lza, lzb

      coef_max=la_max+lb_max+1

      l_set_info%lxyz_max=0
      l_set_info%lxy_max=0
      l_set_info%lx_max=0
      DO lxa=0,la_max
      DO lxb=0,lb_max
         l_set_info%lx_max=l_set_info%lx_max+1
         DO lya=0,la_max-lxa
         DO lyb=0,lb_max-lxb
            l_set_info%lxy_max=l_set_info%lxy_max+1
            DO lza=max(la_min-lxa-lya,0),la_max-lxa-lya
            DO lzb=max(lb_min-lxb-lyb,0),lb_max-lxb-lyb
               l_set_info%lxyz_max=l_set_info%lxyz_max+1
            ENDDO
            ENDDO
         ENDDO
         ENDDO
      ENDDO
      ENDDO
      allocate(l_set_info%ly_max(l_set_info%lx_max),l_set_info%lz_max(l_set_info%lxy_max))
      allocate(l_set_info%map(-cmax:cmax,1:3),&
      l_set_info%pzyx(l_set_info%lxyz_max),&
      l_set_info%ipzyx(0:la_max,0:la_max,0:la_max,&
                       0:lb_max,0:lb_max,0:lb_max),&
      l_set_info%polz(l_set_info%lxyz_max,-cmax:cmax))
      allocate(l_set_info%poly(l_set_info%lxy_max,-cmax:cmax),&
      l_set_info%polx(coef_max,-cmax:cmax),&
      l_set_info%alpha(coef_max,l_set_info%lx_max))
      allocate(l_set_info%dpy(0:lb_max,0:la_max),&
      l_set_info%dpz(0:lb_max,0:la_max))
      l_set_info%ly_max(:)=0
      l_set_info%lz_max(:)=0
      lx=0
      lxy=0
      lxyz=0
      DO lxa=0,la_max
      DO lxb=0,lb_max
         lx=lx+1
         DO lya=0,la_max-lxa
         DO lyb=0,lb_max-lxb
            lxy=lxy+1
            l_set_info%ly_max(lx)=l_set_info%ly_max(lx)+1
            DO lza=max(la_min-lxa-lya,0),la_max-lxa-lya
            DO lzb=max(lb_min-lxb-lyb,0),lb_max-lxb-lyb
               lxyz=lxyz+1
               l_set_info%lz_max(lxy)=l_set_info%lz_max(lxy)+1
               l_set_info%ipzyx(lxa,lya,lza,lxb,lyb,lzb)=lxyz 
            ENDDO
            ENDDO
         ENDDO
         ENDDO
      ENDDO
      ENDDO

   END SUBROUTINE

   SUBROUTINE destroy_l_set_info_type(l_set_info)
    TYPE(l_set_info_type)                    :: l_set_info

    INTEGER                                  :: istat

      deallocate(l_set_info%ly_max,l_set_info%lz_max,l_set_info%map,&
                 l_set_info%polz,l_set_info%poly,l_set_info%polx,&
                 l_set_info%dpy,l_set_info%dpz,l_set_info%pzyx,&
                 l_set_info%ipzyx,l_set_info%alpha,STAT=istat)
      IF (istat /= 0) CALL stop_memory("destroy_l_set_info_type","l_set_info%")
   END SUBROUTINE
    
END MODULE
