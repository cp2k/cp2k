! ***************************************************************
MODULE L_UTILS
   USE termination,           ONLY: stop_program
   IMPLICIT NONE

   public :: l_info_type
   public :: init_l_info
   public :: return_l_info
   public :: destroy_l_info
   
   TYPE l_set_info_type
        integer                        :: lx_max,lxy_max,lxyz_max
        integer, pointer, dimension(:) :: ly_max,lz_max
   END TYPE

   TYPE l_info_type
        integer :: lmax
        TYPE(l_set_info_type), pointer, dimension(:,:,:,:) :: data
   END TYPE

CONTAINS

   SUBROUTINE return_l_info(l_info,la_min,la_max,lb_min,lb_max, &
                            lx_max,lxy_max,lxyz_max,ly_max,lz_max)
      IMPLICIT NONE
      type(l_info_type) :: l_info
      integer           :: la_min,la_max,lb_min,lb_max
      integer           :: lx_max,lxy_max,lxyz_max
      integer, pointer, dimension(:) :: ly_max,lz_max
      if (la_max.gt.l_info%lmax .or. lb_max.gt.l_info%lmax) then
        CALL stop_program("return_l_info","called with l too large")
      endif
      lx_max=l_info%data(la_min,la_max,lb_min,lb_max)%lx_max
      lxy_max=l_info%data(la_min,la_max,lb_min,lb_max)%lxy_max
      lxyz_max=l_info%data(la_min,la_max,lb_min,lb_max)%lxyz_max
      ly_max=>l_info%data(la_min,la_max,lb_min,lb_max)%ly_max
      lz_max=>l_info%data(la_min,la_max,lb_min,lb_max)%lz_max
   END SUBROUTINE

   SUBROUTINE init_l_info(l_info,lmax)
      IMPLICIT NONE
      integer :: lmax
      type(l_info_type) :: l_info
      integer :: la_min,la_max,lb_min,lb_max
      l_info%lmax=lmax
      allocate(l_info%data(0:lmax,0:lmax,0:lmax,0:lmax))
      DO lb_max=0,lmax
      DO lb_min=0,lb_max
      DO la_max=0,lmax
      DO la_min=0,la_max
         call init_l_set_info_type(la_min,la_max,lb_min,lb_max, &
                                   l_info%data(la_min,la_max,lb_min,lb_max))       
      ENDDO
      ENDDO
      ENDDO
      ENDDO
   END SUBROUTINE

   SUBROUTINE destroy_l_info(l_info)
      IMPLICIT NONE
      type(l_info_type) :: l_info
      integer :: la_min,la_max,lb_min,lb_max

      DO lb_max=0,l_info%lmax
      DO lb_min=0,lb_max
      DO la_max=0,l_info%lmax
      DO la_min=0,la_max
         call destroy_l_set_info_type(l_info%data(la_min,la_max,lb_min,lb_max))
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      deallocate(l_info%data)
   END SUBROUTINE

   SUBROUTINE init_l_set_info_type(la_min,la_max,lb_min,lb_max,l_set_info)
      IMPLICIT NONE
      integer               :: la_min,lb_min,la_max,lb_max
      TYPE(l_set_info_type) :: l_set_info
      integer lxa,lxb,lya,lyb,lza,lzb,lx,lxy

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
      l_set_info%ly_max(:)=0
      l_set_info%lz_max(:)=0
      lx=0
      lxy=0
      DO lxa=0,la_max
      DO lxb=0,lb_max
         lx=lx+1
         DO lya=0,la_max-lxa
         DO lyb=0,lb_max-lxb
            lxy=lxy+1
            l_set_info%ly_max(lx)=l_set_info%ly_max(lx)+1
            DO lza=max(la_min-lxa-lya,0),la_max-lxa-lya
            DO lzb=max(lb_min-lxb-lyb,0),lb_max-lxb-lyb
               l_set_info%lz_max(lxy)=l_set_info%lz_max(lxy)+1
            ENDDO
            ENDDO
         ENDDO
         ENDDO
      ENDDO
      ENDDO

   END SUBROUTINE

   SUBROUTINE destroy_l_set_info_type(l_set_info)
      TYPE(l_set_info_type) :: l_set_info
      deallocate(l_set_info%ly_max,l_set_info%lz_max)
   END SUBROUTINE
    
END MODULE
