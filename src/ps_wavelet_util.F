


!!****** cp2k/ ps_wavelet_util [1.0] *
!!
!!   NAME
!!     ps_wavelet_util
!!
!!   FUNCTION
!!     Performs a wavelet based solution of the Poisson equation.
!!
!!   AUTHOR
!!     Florian Schiffmann (09.2007,fschiff)
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE ps_wavelet_util


  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fourpi
  USE ps_wavelet_fft3d,                ONLY: fourier_dim,&
                                             ctrig,&
                                             fftstp

  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_types,                        ONLY: pw_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE message_passing,                 ONLY: mp_cart_rank,&
                                             mp_alltoall,&
                                             mp_bcast
  USE ps_wavelet_base
#include "cp_common_uses.h"


  IMPLICIT NONE
  
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ps_wavelet_util'

! *** Public data types ***

  PUBLIC :: PSolver,&
            PS_dim4allocation,&
            P_FFT_dimensions,&
            S_FFT_dimensions,&
            F_FFT_dimensions

  CONTAINS



!!****h* BigDFT/PSolver
!! NAME
!!    PSolver
!!
!! FUNCTION
!!    Calculate the Poisson equation $\nabla^2 V(x,y,z)=-4 \pi \rho(x,y,z)$
!!    from a given $\rho$, for different boundary conditions an for different data distributions.
!!    Following the boundary conditions, it applies the Poisson Kernel previously calculated.
!!
!! SYNOPSIS
!!    geocode  Indicates the boundary conditions (BC) of the problem:
!!            'F' free BC, isolated systems.
!!                The program calculates the solution as if the given density is
!!                "alone" in R^3 space.
!!            'S' surface BC, isolated in y direction, periodic in xz plane                
!!                The given density is supposed to be periodic in the xz plane,
!!                so the dimensions in these direction mus be compatible with the FFT
!!                Beware of the fact that the isolated direction is y!
!!            'P' periodic BC.
!!                The density is supposed to be periodic in all the three directions,
!!                then all the dimensions must be compatible with the FFT.
!!                No need for setting up the kernel.
!!    datacode Indicates the distribution of the data of the input/output array:
!!            'G' global data. Each process has the whole array of the density 
!!                which will be overwritten with the whole array of the potential
!!            'D' distributed data. Each process has only the needed part of the density
!!                and of the potential. The data distribution is such that each processor
!!                has the xy planes needed for the calculation AND for the evaluation of the 
!!                gradient, needed for XC part, and for the White-Bird correction, which
!!                may lead up to 8 planes more on each side. Due to this fact, the information
!!                between the processors may overlap.
!!    nproc       number of processors
!!    iproc       label of the process,from 0 to nproc-1
!!    n01,n02,n03 global dimension in the three directions. They are the same no matter if the 
!!                datacode is in 'G' or in 'D' position.
!!    ixc         eXchange-Correlation code. Indicates the XC functional to be used 
!!                for calculating XC energies and potential. 
!!                ixc=0 indicates that no XC terms are computed. The XC functional codes follow
!!                the ABINIT convention.
!!    hx,hy,hz    grid spacings. For the isolated BC case for the moment they are supposed to 
!!                be equal in the three directions
!!    rhopot      main input/output array.
!!                On input, it represents the density values on the grid points
!!                On output, it is the Hartree potential, namely the solution of the Poisson 
!!                equation PLUS (when ixc/=0) the XC potential PLUS (again for ixc/=0) the 
!!                pot_ion array. The output is non overlapping, in the sense that it does not
!!                consider the points that are related to gradient and WB calculation
!!    karray      kernel of the poisson equation. It is provided in distributed case, with
!!                dimensions that are related to the output of the PS_dim4allocation routine
!!                it MUST be created by following the same geocode as the Poisson Solver.
!!    pot_ion     additional external potential that is added to the output, 
!!                when the XC parameter ixc/=0. It is always provided in the distributed form,
!!                clearly without the overlapping terms which are needed only for the XC part
!!    eh,exc,vxc  Hartree energy, XC energy and integral of $\rho V_{xc}$ respectively
!!    offset      value of the potential at the point 1,1,1 of the grid.
!!                To be used only in the periodic case, ignored for other boundary conditions.
!! WARNING
!!    The dimensions of the arrays must be compatible with geocode, datacode, nproc, 
!!    ixc and iproc. Since the arguments of these routines are indicated with the *, it
!!    is IMPERATIVE to use the PS_dim4allocation routine for calculation arrays sizes.
!!
!! AUTHOR
!!    Luigi Genovese
!! CREATION DATE
!!    February 2007
!!
!! SOURCE
!!
SUBROUTINE PSolver(geocode,datacode,iproc,nproc,n01,n02,n03,hx,hy,hz,&
     rhopot,karray,eh,exc,vxc,offset,pw_grid,error)
!  INCLUDE 'mpif.h'
  CHARACTER(len=1), INTENT(in) :: geocode
  CHARACTER(len=1), INTENT(in) :: datacode
  INTEGER, INTENT(in) :: iproc,nproc,n01,n02,n03
  REAL(KIND=dp), INTENT(in) :: hx,hy,hz,offset
  REAL(KIND=dp), DIMENSION(*), INTENT(in) :: karray
  REAL(KIND=dp), INTENT(out) :: eh,exc,vxc
  REAL(KIND=dp), DIMENSION(*), INTENT(inout):: rhopot
  TYPE(pw_grid_type), POINTER         :: pw_grid
  TYPE(cp_error_type), INTENT(inout)       :: error
  !local variables
  INTEGER, PARAMETER :: nordgr=4 !the order of the finite-difference gradient (fixed)
  INTEGER :: m1,m2,m3,md1,md2,md3,n1,n2,n3,nd1,nd2,nd3,ixc
  INTEGER :: i_all,i_stat,ierr,ind,ind2,ind3
  INTEGER :: i1,i2,i3,j2,istart,iend,i3start,jend,jproc,i3xcsh
  INTEGER :: nxc,nwbl,nwbr,nxt,nwb,nxcl,nxcr,nlim
  REAL(KIND=dp) :: ehartreeLOC,eexcuLOC,vexcuLOC,red_fact
  REAL(KIND=dp) :: hgrid,scal,newoffset,correction,pot,factor
  REAL(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: zf,zfionxc,pot_ion
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: gather_arr
  REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: energies_mpi

  
ixc =0
    ALLOCATE(pot_ion(n01,n02,n03))
    pot_ion=0.0_dp
  !calculate the dimensions wrt the geocode
  IF (geocode == 'P') THEN
     CALL P_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
  ELSE IF (geocode == 'S') THEN
     CALL S_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
  ELSE IF (geocode == 'F') THEN
     CALL F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
  ELSE
     STOP 'PSolver: geometry code not admitted'
  END IF
  !array allocations
  i_all=0
  ALLOCATE(zf(md1,md3,md2/nproc),stat=i_stat)
  i_all=i_all+i_stat
  ALLOCATE(zfionxc(md1,md3,md2/nproc),stat=i_stat)
  i_all=i_all+i_stat
  IF (nproc > 1) THEN
     IF (datacode == 'G') THEN
        ALLOCATE(gather_arr(0:nproc-1,2),stat=i_stat)
        i_all=i_all+i_stat
     END IF
     ALLOCATE(energies_mpi(6),stat=i_stat)
     i_all=i_all+i_stat
  END IF
  IF (i_all /= 0) THEN
     PRINT *,"PSolver: Problem of memory allocation"
     STOP
  END IF
!  CALL timing(iproc,'Exchangecorr  ','ON')
  !dimension for exchange-correlation (different in the global or distributed case)
  !let us calculate the dimension of the portion of the rhopot array to be passed 
  !to the xc routine
  !this portion will depend on the need of calculating the gradient or not, 
  !and whether the White-Bird correction must be inserted or not 
  !(absent only in the LB ixc=13 case)
  
  !nxc is the effective part of the third dimension that is being processed
  !nxt is the dimension of the part of rhopot that must be passed to the gradient routine
  !nwb is the dimension of the part of rhopot in the wb-postprocessing routine
  !note: nxc <= nwb <= nxt
  !the dimension are related by the values of nwbl and nwbr
  !      nxc+nxcl+nxcr-2 = nwb
  !      nwb+nwbl+nwbr = nxt
  istart=iproc*(md2/nproc)
  iend=min((iproc+1)*md2/nproc,m2)
  nxc=iend-istart
  nwbl=0
  nwbr=0
  nxcl=1
  nxcr=1

  nwb=nxcl+nxc+nxcr-2
  nxt=nwbr+nwb+nwbl

  if (datacode=='G') then
     !starting address of rhopot in the case of global i/o
     i3start=istart+2-nxcl-nwbl
  else if (datacode == 'D') then
     !distributed i/o
     i3start=1
  else
     stop 'PSolver: datacode not admitted'
  end if

  !calculate the actual limit of the array for the zero padded FFT
  if (geocode == 'P') then
     nlim=n2
  else if (geocode == 'S') then
     nlim=n2
  else if (geocode == 'F') then
     nlim=n2/2
  end if

 !!$  print *,'density must go from',min(istart+1,m2),'to',iend,'with n2/2=',n2/2
 !!$  print *,'        it goes from',i3start+nwbl+nxcl-1,'to',i3start+nxc-1
 
   IF (istart+1 <= m2) THEN 
      red_fact=1._dp
      call xc_energy(geocode,m1,m2,m3,md1,md2,md3,nxc,nwb,nxt,nwbl,nwbr,nxcl,nxcr,&
           ixc,hx,hy,hz,rhopot,pot_ion,zf,zfionxc,&
           eexcuLOC,vexcuLOC,iproc,nproc,red_fact)
   ELSE IF (istart+1 <= nlim) THEN !this condition assures that we have perform good zero padding
      DO i2=istart+1,MIN(nlim,istart+md2/nproc)
         j2=i2-istart
         DO i3=1,md3
            DO i1=1,md1
               zf(i1,i3,j2)=0.d0
               zfionxc(i1,i3,j2)=0.d0 !this is not needed, only if pot is updated in Solver
            END DO
         END DO
      END DO
      eexcuLOC=0.d0
      vexcuLOC=0.d0
   ELSE
      eexcuLOC=0.d0
      vexcuLOC=0.d0
   END IF
      
  !this routine builds the values for each process of the potential (zf), multiplying by scal 
  if(geocode == 'P') then
     !no powers of hgrid because they are incorporated in the plane wave treatment
     scal=1.d0/real(n1*n2*n3,KIND=dp)
     call P_PoissonSolver(n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc,zf,&
          scal,hx,hy,hz,pw_grid%para%rs_group,error)

     !offset correction for the periodic treatment
     if (iproc == 0) newoffset=0._dp
     call mp_bcast(newoffset,0,pw_grid%para%rs_group)
     correction=offset-newoffset
     factor=0.5d0*hx*hy*hz

  else if (geocode == 'S') then
     !only one power of hgrid 
     scal=hy/real(n1*n2*n3,KIND=dp)
     call S_PoissonSolver(n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc,karray,zf,&
          scal,hx,hy,hz,pw_grid%para%rs_group,error)!,ehartreeLOC)
     correction=0.d0
     factor=0.5d0*hx*hy*hz
  else if (geocode == 'F') then
     hgrid=max(hx,hy,hz)
     scal=hgrid**3/real(n1*n2*n3,KIND=dp)
     call F_PoissonSolver(n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc,karray,zf,&
          scal,hgrid,pw_grid%para%rs_group,error)!,ehartreeLOC)
     correction=0.d0
     factor=0.5d0*hgrid**3

  end if
  
!  call timing(iproc,'PSolv_comput  ','ON')

  !the value of the shift depends on the distributed i/o or not
  if (datacode=='G') then
     i3xcsh=istart
  else if (datacode=='D') then
     i3xcsh=nxcl+nwbl-1
  end if

  ehartreeLOC=0.d0
  do j2=1,nxc!i2=istart+1,iend !the index i2 must be changed in the distributed i/o case
     i2=j2+i3xcsh!j2=i2-istart
     ind3=(i2-1)*n01*n02
     do i3=1,m3
        ind2=(i3-1)*n01+ind3
        do i1=1,m1
           ind=i1+ind2
           pot=zf(i1,i3,j2)+correction
           ehartreeLOC=ehartreeLOC+rhopot(ind)*pot
           rhopot(ind)=pot+zfionxc(i1,i3,j2)
        end do
     end do
  end do

  IF(geocode=='F')THEN
     red_fact=1._dp
  ELSE
     red_fact=-fourpi
  END IF
  CALL xc_energy(geocode,m1,m2,m3,md1,md2,md3,nxc,nwb,nxt,nwbl,nwbr,nxcl,nxcr,&
       ixc,hx,hy,hz,zf,pot_ion,rhopot,zfionxc,&
       eexcuLOC,vexcuLOC,iproc,nproc,red_fact)
  
  ehartreeLOC=ehartreeLOC*2.d0*factor


  deallocate(zf,stat=i_all)
  deallocate(zfionxc,stat=i_stat)
  if (i_all+i_stat /= 0) then
     write(*,*)' psolver: problem of memory deallocation'
     stop
  end if

!  call timing(iproc,'PSolv_comput  ','OF')

  !gathering the data to obtain the distribution array
  !evaluating the total ehartree,eexcu,vexcu
!!!  if (nproc.gt.1) then
!!!     call timing(iproc,'PSolv_commun  ','ON')
!!!     energies_mpi(1)=ehartreeLOC
!!!     energies_mpi(2)=eexcuLOC
!!!     energies_mpi(3)=vexcuLOC
!!!     call mpi_allreduce(energies_mpi(1),energies_mpi(4),3,MPI_double_precision,  &
!!!          MPI_SUM,pw_grid%para%rs_group,ierr)
!!!     eh=energies_mpi(4)
!!!     exc=energies_mpi(5)
!!!     vxc=energies_mpi(6)
!!!     call timing(iproc,'PSolv_commun  ','OF')
!!!
!!!     if (datacode == 'G') then
!!!        !building the array of the data to be sent from each process
!!!        !and the array of the displacement
!!!        call timing(iproc,'PSolv_comput  ','ON')
!!!        do jproc=0,nproc-1
!!!           istart=min(jproc*(md2/nproc),m2-1)
!!!           jend=max(min(md2/nproc,m2-md2/nproc*jproc),0)
!!!           gather_arr(jproc,1)=m1*m3*jend
!!!           gather_arr(jproc,2)=m1*m3*istart
!!!        end do
!!!        call timing(iproc,'PSolv_comput  ','OF')
!!!
!!!        !gather all the results in the same rhopot array
!!!        istart=min(iproc*(md2/nproc),m2-1)
!!!        call timing(iproc,'PSolv_commun  ','ON')
!!!        call mp_allgather(rhopot(1+n01*n02*istart),gather_arr(iproc,1)&
!!!             rhopot,gather_arr(:,1),gather_arr(:,2),MPI_double_precision,MPI_COMM_WORLD,ierr)
!!!        call timing(iproc,'PSolv_commun  ','OF')
!!!     end if
!!!
!!!  else
!!!     eh=ehartreeLOC
!!!     exc=eexcuLOC
!!!     vxc=vexcuLOC
!!!  end if

  IF (ALLOCATED(gather_arr)) DEALLOCATE(gather_arr)
  IF (ALLOCATED(energies_mpi)) DEALLOCATE(energies_mpi)
  

END SUBROUTINE PSolver




!!****h* BigDFT/PS_dim4allocation
!! NAME
!!    PS_dim4allocation
!!
!! FUNCTION
!!    Calculate the dimensions needed for the allocation of the arrays 
!!    related to the Poisson Solver
!!
!! SYNOPSIS
!!    geocode  Indicates the boundary conditions (BC) of the problem:
!!            'F' free BC, isolated systems.
!!                The program calculates the solution as if the given density is
!!                "alone" in R^3 space.
!!            'S' surface BC, isolated in y direction, periodic in xz plane                
!!                The given density is supposed to be periodic in the xz plane,
!!                so the dimensions in these direction mus be compatible with the FFT
!!                Beware of the fact that the isolated direction is y!
!!            'P' periodic BC.
!!                The density is supposed to be periodic in all the three directions,
!!                then all the dimensions must be compatible with the FFT.
!!                No need for setting up the kernel.
!!    datacode Indicates the distribution of the data of the input/output array:
!!            'G' global data. Each process has the whole array of the density 
!!                which will be overwritten with the whole array of the potential
!!            'D' distributed data. Each process has only the needed part of the density
!!                and of the potential. The data distribution is such that each processor
!!                has the xy planes needed for the calculation AND for the evaluation of the 
!!                gradient, needed for XC part, and for the White-Bird correction, which
!!                may lead up to 8 planes more on each side. Due to this fact, the information
!!                between the processors may overlap.
!!    iproc,nproc number of process, number of processes
!!    n01,n02,n03 dimensions of the real space grid to be hit with the Poisson Solver
!!    ixc         eXchange-Correlation code. Indicates the XC functional to be used 
!!                for calculating XC energies and potential. 
!!                ixc=0 indicates that no XC terms are computed. The XC functional codes follow
!!                the ABINIT convention.
!!    n3d         third dimension of the density. For distributed data, it takes into account 
!!                the enlarging needed for calculating the XC functionals.
!!                For global data it it simply equal to n03. 
!!                When there are too many processes and there is no room for the density n3d=0
!!    n3p         third dimension for the potential. The same as n3d, but without 
!!                taking into account the enlargment for the XC part. For non-GGA XC, n3p=n3d.
!!    n3pi        Dimension of the pot_ion array, always with distributed data. 
!!                For distributed data n3pi=n3p
!!    i3xcsh      Shift of the density that must be performed to enter in the 
!!                non-overlapping region. Useful for recovering the values of the potential
!!                when using GGA XC functionals. If the density starts from rhopot(1,1,1),
!!                the potential starts from rhopot(1,1,i3xcsh+1). 
!!                For non-GGA XCs and for global distribution data i3xcsh=0
!!    i3s         Starting point of the density effectively treated by each processor 
!!                in the third direction.
!!                It takes into account also the XC enlarging. The array rhopot will correspond
!!                To the planes of third coordinate from i3s to i3s+n3d-1. 
!!                The potential to the planes from i3s+i3xcsh to i3s+i3xcsh+n3p-1
!!                The array pot_ion to the planes from i3s+i3xcsh to i3s+i3xcsh+n3pi-1
!!                For global disposition i3s is equal to distributed case with i3xcsh=0.
!!
!!
!! WARNING
!!    The XC enlarging due to GGA part is not present for surfaces and 
!!    periodic boundary condition. This is related to the fact that the calculation of the
!!    gradient and the White-Bird correction are not yet implemented for non-isolated systems
!!
!! AUTHOR
!!    Luigi Genovese
!! CREATION DATE
!!    February 2007
!!
!! SOURCE
!!
subroutine PS_dim4allocation(geocode,datacode,iproc,nproc,n01,n02,n03,ixc,&
     n3d,n3p,n3pi,i3xcsh,i3s)
  implicit none
  character(len=1), intent(in) :: geocode
  character(len=1), intent(in) :: datacode
  integer, intent(in) :: iproc,nproc,n01,n02,n03,ixc
  integer, intent(out) :: n3d,n3p,n3pi,i3xcsh,i3s
  !local variables
  integer, parameter :: nordgr=4
  integer :: m1,m2,m3,md1,md2,md3,n1,n2,n3,nd1,nd2,nd3
  integer :: istart,iend,nxc,nwb,nxt,nxcl,nxcr,nwbl,nwbr


  !calculate the dimensions wrt the geocode
  if (geocode == 'P') then
     call P_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
  else if (geocode == 'S') then
     call S_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
  else if (geocode == 'F') then
     call F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
  else
     stop 'PS_dim4allocation: geometry code not admitted'
  end if

  !formal start and end of the slice
  istart=iproc*(md2/nproc)
  iend=min((iproc+1)*md2/nproc,m2)

  IF (datacode == 'D') THEN
     IF (istart <= m2-1) THEN
        nxc=iend-istart
        IF (ixc >= 11 .AND. ixc <= 16 .AND. geocode == 'F') THEN
           IF (ixc==13) THEN
              !now the dimension of the part required for the gradient
              nwbl=MIN(istart,nordgr)
              nwbr=MIN(m2-iend,nordgr)
              nxcl=1
              nxcr=1
           ELSE
              !now the dimension of the part required for the gradient
              IF(istart<=nordgr) THEN
                 nxcl=istart+1
                 nwbl=0
              ELSE
                 nxcl=nordgr+1
                 nwbl=MIN(nordgr,istart-nordgr)
              END IF
              IF(iend>=m2-nordgr+1) THEN
                 nxcr=m2-iend+1
                 nwbr=0
              ELSE
                 nxcr=nordgr+1
                 nwbr=min(nordgr,m2-nordgr-iend)
              END IF
           END IF
        ELSE !(for the moment GGA is not implemented in the non free BC)
           nwbl=0
           nwbr=0
           nxcl=1
           nxcr=1
        END IF
        nwb=nxcl+nxc+nxcr-2
        nxt=nwbr+nwb+nwbl

        i3xcsh=nxcl+nwbl-1
        i3s=istart+1-i3xcsh
     ELSE
        nxc=0
        nxt=0
        i3xcsh=0
        i3s=m2
     end if
     n3p=nxc
     n3d=nxt
     n3pi=n3p
  ELSE IF (datacode == 'G') THEN
     n3d=n03
     n3p=n03
     i3xcsh=0
     i3s=min(istart,m2-1)+1
     n3pi=max(iend-istart,0)
  ELSE
     print *,datacode
     stop 'PS_dim4allocation: data code not admitted'
  end if

!!$  print *,'P4,iproc',iproc,'nxc,ncxl,ncxr,nwbl,nwbr',nxc,nxcl,nxcr,nwbl,nwbr,&
!!$       'ixc,n3d,n3p,i3xcsh,i3s',ixc,n3d,n3p,i3xcsh,i3s

end subroutine PS_dim4allocation


!!***
!!****h* BigDFT/P_FFT_dimensions
!! NAME
!!   P_FFT_dimensions
!!
!! FUNCTION
!!    Calculate four sets of dimension needed for the calculation of the
!!    convolution for the periodic system
!!
!! SYNOPSIS
!!    n01,n02,n03 original real dimensions (input)
!!
!!    m1,m2,m3 original real dimension, with m2 and m3 exchanged
!!
!!    n1,n2,n3 the first FFT dimensions, for the moment supposed to be even
!!
!!    md1,md2,md3 the n1,n2,n3 dimensions. They contain the real unpadded space,
!!                properly enlarged to be compatible with the FFT dimensions n_i.
!!                md2 is further enlarged to be a multiple of nproc
!!
!!    nd1,nd2,nd3 fourier dimensions for which the kernel is injective,
!!                formally 1/8 of the fourier grid. Here the dimension nd3 is
!!                enlarged to be a multiple of nproc
!!
!! WARNING
!!    This four sets of dimensions are actually redundant (mi=n0i), 
!!    due to the backward-compatibility
!!    with the other geometries of the Poisson Solver.
!!    The dimensions 2 and 3 are exchanged.
!!
!! AUTHOR
!!    Luigi Genovese
!! CREATION DATE
!!    October 2006
!!
!! SOURCE
!!
subroutine P_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
 implicit none
 integer, intent(in) :: n01,n02,n03,nproc
 integer, intent(out) :: m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3
 integer :: l1,l2,l3,l1A,l3A

 !dimensions of the density in the real space
 m1=n01
 m2=n03
 m3=n02

 ! real space grid dimension (suitable for number of processors)
 l1=m1
 l2=m2
 l3=m3 !beware of the half dimension
    call fourier_dim(l1,n1)
    if (n1 == m1) then
    else
       print *,'the FFT in the x direction is not allowed'
       print *,'n01 dimension',n01
       stop
    end if
    l1=l1+1
    call fourier_dim(l2,n2)
    if (n2 == m2) then
    else
       print *,'the FFT in the z direction is not allowed'
       print *,'n03 dimension',n03
       stop
    end if
    call fourier_dim(l3,n3)
    if (n3 == m3) then
    else
       print *,'the FFT in the y direction is not allowed'
       print *,'n02 dimension',n02
       stop
    end if

 !dimensions that contain the unpadded real space,
 ! compatible with the number of processes
 md1=n1
 md2=n2
 md3=n3
151 if (nproc*(md2/nproc).lt.n2) then
    md2=md2+1
    goto 151
 endif


 !dimensions of the kernel, 1/8 of the total volume,
 !compatible with nproc
 nd1=n1/2+1 
 nd2=n2/2+1
 nd3=n3/2+1
250 if (modulo(nd3,nproc) .ne. 0) then
    nd3=nd3+1
    goto 250
 endif

end subroutine P_FFT_dimensions


!!****h* BigDFT/S_FFT_dimensions
!! NAME
!!   S_FFT_dimensions
!!
!! FUNCTION
!!    Calculate four sets of dimension needed for the calculation of the
!!    convolution for the surface system
!!
!! SYNOPSIS
!!    n01,n02,n03 original real dimensions (input)
!!
!!    m1,m2,m3 original real dimension, with 2 and 3 exchanged
!!
!!    n1,n2 the first FFT dimensions, for the moment supposed to be even
!!    n3    the double of the first FFT even dimension greater than m3
!!          (improved for the HalFFT procedure)
!!
!!    md1,md2     the n1,n2 dimensions. 
!!    md3         half of n3 dimension. They contain the real unpadded space,
!!                properly enlarged to be compatible with the FFT dimensions n_i.
!!                md2 is further enlarged to be a multiple of nproc
!!
!!    nd1,nd2,nd3 fourier dimensions for which the kernel FFT is injective,
!!                formally 1/8 of the fourier grid. Here the dimension nd3 is
!!                enlarged to be a multiple of nproc
!!
!! WARNING
!!    This four sets of dimensions are actually redundant (mi=n0i), 
!!    due to the backward-compatibility
!!    with the Poisson Solver with other geometries.
!!    Dimensions n02 and n03 were exchanged
!!
!! AUTHOR
!!    Luigi Genovese
!! CREATION DATE
!!    October 2006
!!
!! SOURCE
!!
subroutine S_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
 implicit none
 integer, intent(in) :: n01,n02,n03,nproc
 integer, intent(out) :: m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3
 INTEGER :: l1,l2,l3,l1A,l3A,handle
 CHARACTER(len=*), PARAMETER :: routineN = 'S_FFT_dimensions'
 !dimensions of the density in the real space

 CALL timeset(routineN,"I","",handle)
 m1=n01
 m2=n03
 m3=n02

 ! real space grid dimension (suitable for number of processors)
 l1=m1
 l2=m2
 l3=m3 !beware of the half dimension
    call fourier_dim(l1,n1)
    if (n1 == m1) then
    else
       print *,'the FFT in the x direction is not allowed'
       print *,'n01 dimension',n01
       stop
    end if
    l1=l1+1
    call fourier_dim(l2,n2)
    if (n2 == m2) then
    else
       print *,'the FFT in the z direction is not allowed'
       print *,'n03 dimension',n03
       stop
    end if
 do
    call fourier_dim(l3,n3)
    if (modulo(n3,2) == 0) then
       exit
    end if
    l3=l3+1
 end do
 n3=2*n3

 !dimensions that contain the unpadded real space,
 ! compatible with the number of processes
 md1=n1
 md2=n2
 md3=n3/2
151 if (nproc*(md2/nproc).lt.n2) then
    md2=md2+1
    goto 151
 endif


 !dimensions of the kernel, 1/8 of the total volume,
 !compatible with nproc

 !these two dimensions are like that since they are even
 nd1=n1/2+1
 nd2=n2/2+1

 nd3=n3/2+1
250 if (modulo(nd3,nproc) .ne. 0) then
    nd3=nd3+1
    goto 250
 endif 
    CALL timestop(0.0_dp,handle)

end subroutine S_FFT_dimensions
!!***

!!****h* BigDFT/F_FFT_dimensions
!! NAME
!!   F_FFT_pardimensions
!!
!! FUNCTION
!!    Calculate four sets of dimension needed for the calculation of the
!!    zero-padded convolution
!!
!! SYNOPSIS
!!    n01,n02,n03 original real dimensions (input)
!!
!!    m1,m2,m3 original real dimension with the dimension 2 and 3 exchanged
!!
!!    n1,n2 the first FFT even dimensions greater that 2*m1, 2*m2
!!    n3    the double of the first FFT even dimension greater than m3
!!          (improved for the HalFFT procedure)
!!
!!    md1,md2,md3 half of n1,n2,n3 dimension. They contain the real unpadded space,
!!                properly enlarged to be compatible with the FFT dimensions n_i.
!!                md2 is further enlarged to be a multiple of nproc
!!
!!    nd1,nd2,nd3 fourier dimensions for which the kernel FFT is injective,
!!                formally 1/8 of the fourier grid. Here the dimension nd3 is
!!                enlarged to be a multiple of nproc
!!
!! WARNING
!!    The dimension m2 and m3 correspond to n03 and n02 respectively
!!    this is needed since the convolution routine manage arrays of dimension
!!    (md1,md3,md2/nproc)
!!
!! AUTHOR
!!    Luigi Genovese
!! CREATION DATE
!!    February 2006
!!
!! SOURCE
!!
SUBROUTINE F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
 implicit none
 integer, intent(in) :: n01,n02,n03,nproc
 integer, intent(out) :: m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3
 integer :: l1,l2,l3,l1A,l3A

 !dimensions of the density in the real space, inverted for convenience
 m1=n01
 m2=n03
 m3=n02
 ! real space grid dimension (suitable for number of processors)
 l1=2*m1
 l2=2*m2
 l3=m3 !beware of the half dimension
 DO
    CALL fourier_dim(l1,n1)
    IF (MODULO(n1,2) == 0) THEN
       EXIT
    end if
    l1=l1+1
 END DO
 DO
    CALL fourier_dim(l2,n2)
    IF (MODULO(n2,2) == 0) THEN
       EXIT
    END IF
    l2=l2+1
 END DO
 DO
    CALL fourier_dim(l3,n3)
    IF (MODULO(n3,2) == 0) THEN
       EXIT
    END IF
    l3=l3+1
 END DO
 n3=2*n3
 
 !dimensions that contain the unpadded real space,
 ! compatible with the number of processes
 md1=n1/2
 md2=n2/2
 md3=n3/2
151 IF (nproc*(md2/nproc).LT.n2/2) THEN
    md2=md2+1
    GOTO 151
 ENDIF
 
 !dimensions of the kernel, 1/8 of the total volume,
 !compatible with nproc
 nd1=n1/2+1
 nd2=n2/2+1
 nd3=n3/2+1
 
250 IF (MODULO(nd3,nproc) .NE. 0) THEN
    nd3=nd3+1
    GOTO 250
 ENDIF
 
END SUBROUTINE F_FFT_dimensions




SUBROUTINE xc_energy(geocode,m1,m2,m3,md1,md2,md3,nxc,nwb,nxt,nwbl,nwbr,&
     nxcl,nxcr,ixc,hx,hy,hz,rhopot,pot_ion,zf,zfionxc,exc,vxc,iproc,nproc,factor)

  IMPLICIT NONE

  !Arguments----------------------
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: m1,m2,m3,nxc,nwb,nxcl,nxcr,nxt,md1,md2,md3,ixc,iproc,nproc
  integer, intent(in) :: nwbl,nwbr
  REAL(KIND=dp), INTENT(in) :: hx,hy,hz,factor
  REAL(KIND=dp), DIMENSION(md1,md3,md2/nproc), INTENT(inout) :: rhopot
  real(KIND=dp), dimension(m1,m3,nxc), intent(in) :: pot_ion
  real(KIND=dp), dimension(md1,md3,md2/nproc), intent(out) :: zf,zfionxc
  real(KIND=dp), intent(out) :: exc,vxc

  !Local variables----------------
  REAL(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: exci,d2vxci
  real(KIND=dp), dimension(:,:,:,:), allocatable :: vxci,dvxci,dvxcdgr
  real(KIND=dp), dimension(:,:,:,:,:), allocatable :: gradient
  REAL(KIND=dp) :: elocal,vlocal,rho,pot,potion,hgrid,facpotion
  integer :: npts,i_all,nspden,order,offset,i_stat
  integer :: i1,i2,i3,j1,j2,j3,jp2,jpp2,jppp2
  integer :: ndvxc,nvxcdgr,ngr2

  !Body

  zf = 0.0_dp
  zfionxc = 0.0_dp
vlocal=0.0_dp
potion=0.0
DO jp2=1,nxc
   j2=offset+jp2+nxcl-2
   jpp2=jp2     
   jppp2=jp2+nxcl-1
   DO j3=1,m3
      DO j1=1,m1
         zf(j1,j3,jp2)=factor*rhopot(j1,j3,jp2)
         zfionxc(j1,j3,jp2)=potion+vlocal
      END DO
      DO j1=m1+1,md1
         zf(j1,j3,jp2)=0.d0
         zfionxc(j1,j3,jp2)=0.d0
      END DO
   END DO
   DO j3=m3+1,md3
      DO j1=1,md1
         zf(j1,j3,jp2)=0.d0
         zfionxc(j1,j3,jp2)=0.d0
      END DO
   END DO
END DO
DO jp2=nxc+1,md2/nproc
   DO j3=1,md3
      DO j1=1,md1
         zf(j1,j3,jp2)=0.d0
         zfionxc(j1,j3,jp2)=0.d0
      END DO
   END DO
END DO

      !the two factor is due to the 
  !need of using the density of states in abinit routines


  !De-allocations


end subroutine xc_energy
END MODULE ps_wavelet_util
