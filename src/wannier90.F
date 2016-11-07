!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Outtakes from Wannier90 code
!> \par History
!>      06.2016 created [JGH]
!> \author JGH
! **************************************************************************************************
!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
!                                                            !
!                       WANNIER90                            !
!                                                            !
!          The Maximally-Localised Generalised               !
!                 Wannier Functions Code                     !
!                                                            !
! Wannier90 v2.0 authors:                                    !
!           Arash A. Mostofi   (Imperial College London)     !
!           Jonathan R. Yates  (University of Oxford)        !
!           Giovanni Pizzi     (EPFL, Switzerland)           !
!           Ivo Souza          (Universidad del Pais Vasco)  !
!                                                            !
! Contributors:                                              !
!          Young-Su Lee        (KIST, S. Korea)              !
!          Matthew Shelley     (Imperial College London)     !
!          Nicolas Poilvert    (Penn State University)       !
!          Raffaello Bianco    (Paris 6 and CNRS)            !
!          Gabriele Sclauzero  (ETH Zurich)                  !
!                                                            !
!  Please cite                                               !
!                                                            !
!  [ref] A. A. Mostofi, J. R. Yates, Y.-S. Lee, I. Souza,    !
!        D. Vanderbilt and N. Marzari, "Wannier90: A Tool    !
!        for Obtaining Maximally Localised Wannier           !
!        Functions", Computer Physics Communications,        !
!        178, 685 (2008)                                     !
!                                                            !
!  in any publications arising from the use of this code.    !
!                                                            !
!                                                            !
!  Wannier90 is based on Wannier77, written by N. Marzari,   !
!  I. Souza and D. Vanderbilt. For the method please cite    !
!                                                            !
!  [ref] N. Marzari and D. Vanderbilt,                       !
!        Phys. Rev. B 56 12847 (1997)                        !
!                                                            !
!  [ref] I. Souza, N. Marzari and D. Vanderbilt,             !
!        Phys. Rev. B 65 035109 (2001)                       !
!                                                            !
!                                                            !
! Copyright (C) 2007-13 Jonathan Yates, Arash Mostofi,       !
!                Giovanni Pizzi, Young-Su Lee,               !
!                Nicola Marzari, Ivo Souza, David Vanderbilt !
!                                                            !
! This file is distributed under the terms of the GNU        !
! General Public License. See the file `LICENSE' in          !
! the root directory of the present distribution, or         !
! http://www.gnu.org/copyleft/gpl.txt .                      !
!                                                            !
!------------------------------------------------------------!

MODULE wannier90
   USE kinds,                           ONLY: dp
   USE physcon,                         ONLY: bohr
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'wannier90'

  !Input
  INTEGER            :: iprint
  CHARACTER(len=20)  :: length_unit
  INTEGER            :: stdout

  !parameters dervied from input
  INTEGER            :: num_kpts
  REAL(kind=dp)      :: real_lattice(3,3)
  REAL(kind=dp)      :: recip_lattice(3,3)
  REAL(kind=dp)      :: cell_volume
  REAL(kind=dp), ALLOCATABLE     ::kpt_cart(:,:) !kpoints in cartesians
  REAL(kind=dp)      :: lenconfac
  INTEGER            :: num_exclude_bands

  REAL(kind=dp)      :: kmesh_tol
  INTEGER            :: num_shells
  INTEGER            :: mp_grid(3)
  INTEGER            :: search_shells
  INTEGER, ALLOCATABLE :: shell_list(:)
  REAL(kind=dp), ALLOCATABLE     :: kpt_latt(:,:) !kpoints in lattice vecs


  ! kmesh parameters (set in kmesh)
  INTEGER                     :: nnh           ! the number of b-directions (bka)
  INTEGER                     :: nntot         ! total number of neighbours for each k-point
  INTEGER       , ALLOCATABLE :: nnlist(:,:)   ! list of neighbours for each k-point
  INTEGER       , ALLOCATABLE :: neigh(:,:)
  INTEGER       , ALLOCATABLE :: nncell(:,:,:) ! gives BZ of each neighbour of each k-point
  REAL(kind=dp)               :: wbtot
  REAL(kind=dp) , ALLOCATABLE :: wb(:)         ! weights associated with neighbours of each k-point
  REAL(kind=dp) , ALLOCATABLE :: bk(:,:,:)     ! the b-vectors that go from each k-point to its neighbours
  REAL(kind=dp) , ALLOCATABLE :: bka(:,:)      ! the b-directions from 1st k-point to its neighbours

  ! The maximum number of shells we need to satisfy B1 condition in kmesh
  INTEGER, PARAMETER :: max_shells=6
  INTEGER, PARAMETER :: num_nnmax=12

  INTEGER, PARAMETER :: nsupcell=5
  INTEGER            :: lmn(3,(2*nsupcell+1)**3)

  REAL(kind=dp), PARAMETER    :: eps5  = 1.0e-5_dp
  REAL(kind=dp), PARAMETER    :: eps6  = 1.0e-6_dp
  REAL(kind=dp), PARAMETER    :: eps8  = 1.0e-8_dp

   PUBLIC :: w90_write_header, wannier_setup

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param mp_grid_loc ...
!> \param num_kpts_loc ...
!> \param real_lattice_loc ...
!> \param recip_lattice_loc ...
!> \param kpt_latt_loc ...
!> \param nntot_loc ...
!> \param nnlist_loc ...
!> \param nncell_loc ...
!> \param iounit ...
! **************************************************************************************************
SUBROUTINE wannier_setup(mp_grid_loc,num_kpts_loc,&
     real_lattice_loc,recip_lattice_loc,kpt_latt_loc, &
     nntot_loc,nnlist_loc,nncell_loc,iounit)

      INTEGER, DIMENSION(3), INTENT(in)                  :: mp_grid_loc
      INTEGER, INTENT(in)                                :: num_kpts_loc
      REAL(kind=dp), DIMENSION(3, 3), INTENT(in)         :: real_lattice_loc, recip_lattice_loc
      REAL(kind=dp), DIMENSION(3, num_kpts_loc), &
         INTENT(in)                                      :: kpt_latt_loc
      INTEGER, INTENT(out)                               :: nntot_loc
      INTEGER, DIMENSION(num_kpts_loc, num_nnmax), &
         INTENT(out)                                     :: nnlist_loc
      INTEGER, DIMENSION(3, num_kpts_loc, num_nnmax), &
         INTENT(out)                                     :: nncell_loc
      INTEGER, INTENT(in)                                :: iounit

      INTEGER                                            :: nkp

  ! interface uses atomic units
  length_unit='bohr'
  lenconfac=1.0_dp/bohr
  stdout = iounit

  CALL w90_write_header(stdout)

  WRITE(stdout,'(a/)') ' Setting up k-point neighbours...'

  ! copy local data into module variables
  mp_grid=mp_grid_loc
  num_kpts=num_kpts_loc
  real_lattice=real_lattice_loc
  recip_lattice=recip_lattice_loc
  ALLOCATE ( kpt_latt(3,num_kpts) )
  ALLOCATE ( kpt_cart(3,num_kpts) )
  kpt_latt(1:3,1:nkp)=kpt_latt_loc(1:3,1:nkp)
  DO nkp=1,num_kpts
     kpt_cart(:,nkp)=MATMUL(kpt_latt(:,nkp),recip_lattice(:,:))
  END DO

  num_shells = 0
  ALLOCATE( shell_list(max_shells) )

  cell_volume = real_lattice(1,1)*(real_lattice(2,2)*real_lattice(3,3)-real_lattice(3,2)*real_lattice(2,3)) +&
                real_lattice(1,2)*(real_lattice(2,3)*real_lattice(3,1)-real_lattice(3,3)*real_lattice(2,1)) +&
                real_lattice(1,3)*(real_lattice(2,1)*real_lattice(3,2)-real_lattice(3,1)*real_lattice(2,2))
  iprint = 1
  search_shells = 12
  kmesh_tol = 0.000001_dp
  num_exclude_bands=0

  CALL w90_param_write(stdout)

  CALL w90_kmesh_get()

  nntot_loc = nntot
  nnlist_loc = 0
  nnlist_loc(:,1:nntot) =  nnlist(:,1:nntot)
  nncell_loc = 0
  nncell_loc(:,:,1:nntot) =  nncell(:,:,1:nntot)

  DEALLOCATE(bk,bka,wb)
  DEALLOCATE(nncell,neigh,nnlist)
  DEALLOCATE(kpt_latt,kpt_cart,shell_list)

  WRITE(stdout,'(/a/)') ' Finished setting up k-point neighbours.'

END SUBROUTINE wannier_setup
! **************************************************************************************************
!> \brief ...
!> \param stdout ...
! **************************************************************************************************
   SUBROUTINE w90_write_header(stdout)
      INTEGER, INTENT(IN)                                :: stdout

      WRITE (stdout, *)
      WRITE (stdout, *) '            +---------------------------------------------------+'
      WRITE (stdout, *) '            |                                                   |'
      WRITE (stdout, *) '            |                   WANNIER90                       |'
      WRITE (stdout, *) '            |                                                   |'
      WRITE (stdout, *) '            +---------------------------------------------------+'
      WRITE (stdout, *) '            |                                                   |'
      WRITE (stdout, *) '            |        Welcome to the Maximally-Localized         |'
      WRITE (stdout, *) '            |        Generalized Wannier Functions code         |'
      WRITE (stdout, *) '            |            http://www.wannier.org                 |'
      WRITE (stdout, *) '            |                                                   |'
      WRITE (stdout, *) '            |  Wannier90 v2.0 Authors:                          |'
      WRITE (stdout, *) '            |    Arash A. Mostofi  (Imperial College London)    |'
      WRITE (stdout, *) '            |    Giovanni Pizzi    (EPFL)                       |'
      WRITE (stdout, *) '            |    Ivo Souza         (Universidad del Pais Vasco) |'
      WRITE (stdout, *) '            |    Jonathan R. Yates (University of Oxford)       |'
      WRITE (stdout, *) '            |                                                   |'
      WRITE (stdout, *) '            |  Wannier90 Contributors:                          |'
      WRITE (stdout, *) '            |    Young-Su Lee       (KIST, S. Korea)            |'
      WRITE (stdout, *) '            |    Matthew Shelley    (Imperial College London)   |'
      WRITE (stdout, *) '            |    Nicolas Poilvert   (Penn State University)     |'
      WRITE (stdout, *) '            |    Raffaello Bianco   (Paris 6 and CNRS)          |'
      WRITE (stdout, *) '            |    Gabriele Sclauzero (ETH Zurich)                |'
      WRITE (stdout, *) '            |                                                   |'
      WRITE (stdout, *) '            |  Wannier77 Authors:                               |'
      WRITE (stdout, *) '            |    Nicola Marzari    (EPFL)                       |'
      WRITE (stdout, *) '            |    Ivo Souza         (Universidad del Pais Vasco) |'
      WRITE (stdout, *) '            |    David Vanderbilt  (Rutgers University)         |'
      WRITE (stdout, *) '            |                                                   |'
      WRITE (stdout, *) '            |  Please cite                                      |'
      WRITE (stdout, *) '            |                                                   |'
      WRITE (stdout, *) '            |  [ref] "Wannier90: A Tool for Obtaining Maximally |'
      WRITE (stdout, *) '            |         Localised Wannier Functions"              |'
      WRITE (stdout, *) '            |        A. A. Mostofi, J. R. Yates, Y.-S. Lee,     |'
      WRITE (stdout, *) '            |        I. Souza, D. Vanderbilt and N. Marzari     |'
      WRITE (stdout, *) '            |        Comput. Phys. Commun. 178, 685 (2008)      |'
      WRITE (stdout, *) '            |                                                   |'
      WRITE (stdout, *) '            |  in any publications arising from the use of      |'
      WRITE (stdout, *) '            |  this code. For the method please cite            |'
      WRITE (stdout, *) '            |                                                   |'
      WRITE (stdout, *) '            |  [ref] "Maximally Localized Generalised Wannier   |'
      WRITE (stdout, *) '            |         Functions for Composite Energy Bands"     |'
      WRITE (stdout, *) '            |         N. Marzari and D. Vanderbilt              |'
      WRITE (stdout, *) '            |         Phys. Rev. B 56 12847 (1997)              |'
      WRITE (stdout, *) '            |                                                   |'
      WRITE (stdout, *) '            |  [ref] "Maximally Localized Wannier Functions     |'
      WRITE (stdout, *) '            |         for Entangled Energy Bands"               |'
      WRITE (stdout, *) '            |         I. Souza, N. Marzari and D. Vanderbilt    |'
      WRITE (stdout, *) '            |         Phys. Rev. B 65 035109 (2001)             |'
      WRITE (stdout, *) '            |                                                   |'
      WRITE (stdout, *) '            |                                                   |'
      WRITE (stdout, *) '            | Copyright (c) 1996-2015                           |'
      WRITE (stdout, *) '            |        Arash A. Mostofi, Jonathan R. Yates,       |'
      WRITE (stdout, *) '            |        Young-Su Lee, Giovanni Pizzi, Ivo Souza,   |'
      WRITE (stdout, *) '            |        David Vanderbilt and Nicola Marzari        |'
      WRITE (stdout, *) '            |                                                   |'
      WRITE (stdout, *) '            |        Release: 2.0.1   2nd April 2015            |'
      WRITE (stdout, *) '            |                                                   |'
      WRITE (stdout, *) '            | This program is free software; you can            |'
      WRITE (stdout, *) '            | redistribute it and/or modify it under the terms  |'
      WRITE (stdout, *) '            | of the GNU General Public License as published by |'
      WRITE (stdout, *) '            | the Free Software Foundation; either version 2 of |'
      WRITE (stdout, *) '            | the License, or (at your option) any later version|'
      WRITE (stdout, *) '            |                                                   |'
      WRITE (stdout, *) '            | This program is distributed in the hope that it   |'
      WRITE (stdout, *) '            | will be useful, but WITHOUT ANY WARRANTY; without |'
      WRITE (stdout, *) '            | even the implied warranty of MERCHANTABILITY or   |'
      WRITE (stdout, *) '            | FITNESS FOR A PARTICULAR PURPOSE. See the GNU     |'
      WRITE (stdout, *) '            | General Public License for more details.          |'
      WRITE (stdout, *) '            |                                                   |'
      WRITE (stdout, *) '            | You should have received a copy of the GNU General|'
      WRITE (stdout, *) '            | Public License along with this program; if not,   |'
      WRITE (stdout, *) '            | write to the Free Software Foundation, Inc.,      |'
      WRITE (stdout, *) '            | 675 Mass Ave, Cambridge, MA 02139, USA.           |'
      WRITE (stdout, *) '            |                                                   |'
      WRITE (stdout, *) '            +---------------------------------------------------+'
      WRITE (stdout, *) ''

   END SUBROUTINE w90_write_header

! **************************************************************************************************
!> \brief ...
!> \param stdout ...
! **************************************************************************************************
  SUBROUTINE w90_param_write(stdout)
      INTEGER, INTENT(IN)                                :: stdout

      INTEGER                                            :: i, nkp

    ! System
    WRITE(stdout,'(36x,a6)') '------'
    WRITE(stdout,'(36x,a6)') 'SYSTEM'
    WRITE(stdout,'(36x,a6)') '------'
    WRITE(stdout,*)
    WRITE(stdout,'(28x,a22)') 'Lattice Vectors (Bohr)'
    WRITE(stdout,'(20x,a3,2x,3F11.6)') 'a_1',(real_lattice(1,I)*lenconfac, i=1,3)
    WRITE(stdout,'(20x,a3,2x,3F11.6)') 'a_2',(real_lattice(2,I)*lenconfac, i=1,3)
    WRITE(stdout,'(20x,a3,2x,3F11.6)') 'a_3',(real_lattice(3,I)*lenconfac, i=1,3)
    WRITE(stdout,*)
    WRITE(stdout,'(19x,a17,3x,f11.5)',advance='no') &
         'Unit Cell Volume:',cell_volume*lenconfac**3
    WRITE(stdout,'(2x,a8)') '(Bohr^3)'
    WRITE(stdout,*)
    WRITE(stdout,'(22x,a34)') 'Reciprocal-Space Vectors (Bohr^-1)'
    WRITE(stdout,'(20x,a3,2x,3F11.6)') 'b_1',(recip_lattice(1,I)/lenconfac, i=1,3)
    WRITE(stdout,'(20x,a3,2x,3F11.6)') 'b_2',(recip_lattice(2,I)/lenconfac, i=1,3)
    WRITE(stdout,'(20x,a3,2x,3F11.6)') 'b_3',(recip_lattice(3,I)/lenconfac, i=1,3)
    WRITE(stdout,*)   ' '
    WRITE(stdout,*) ' '
    ! K-points
    WRITE(stdout,'(32x,a)') '------------'
    WRITE(stdout,'(32x,a)') 'K-POINT GRID'
    WRITE(stdout,'(32x,a)') '------------'
    WRITE(stdout,*) ' '
    WRITE(stdout,'(13x,a,i3,1x,a1,i3,1x,a1,i3,6x,a,i5)') 'Grid size =',mp_grid(1),'x',mp_grid(2),'x',mp_grid(3),&
         'Total points =',num_kpts
    WRITE(stdout,*) ' '
    WRITE(stdout,'(1x,a)') '*----------------------------------------------------------------------------*'
    WRITE(stdout,'(1x,a)') '| k-point      Fractional Coordinate        Cartesian Coordinate (Bohr^-1)   |'
    WRITE(stdout,'(1x,a)') '+----------------------------------------------------------------------------+'
    DO nkp=1,num_kpts
       WRITE(stdout,'(1x,a1,i6,1x,3F10.5,3x,a1,1x,3F10.5,4x,a1)') '|',nkp,kpt_latt(:,nkp),'|',kpt_cart(:,nkp)/lenconfac,'|'
    END DO
    WRITE(stdout,'(1x,a)') '*----------------------------------------------------------------------------*'
    WRITE(stdout,*) ' '

  END SUBROUTINE w90_param_write

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
  SUBROUTINE w90_kmesh_get()

    ! Variables that are private

      REAL(kind=dp), PARAMETER                           :: eta = 99999999.0_dp    

      INTEGER :: counter, i, ifneg, ifound, ifpos, j, l, loop, loop_b, loop_s, m, &
         multi(search_shells), n, na, nap, ndnn, ndnntot, ndnnx, nkp, nkp2, nlist, nn, nnsh, &
         nnshell(num_kpts,search_shells), nnx, shell
      REAL(kind=dp) :: bb1, bbn, bk_local(3,num_nnmax,num_kpts), bweight(max_shells), ddelta, &
         dist, dnn(search_shells), dnn0, dnn1, vkpp(3), vkpp2(3), wb_local(num_nnmax)
      REAL(kind=dp), ALLOCATABLE                         :: bvec_tmp(:,:)

     WRITE(stdout,'(/1x,a)') &
       '*---------------------------------- K-MESH ----------------------------------*'
 
     ! Sort the cell neighbours so we loop in order of distance from the home shell
     CALL w90_kmesh_supercell_sort
 
     ! find the distance between k-point 1 and its nearest-neighbour shells
     ! if we have only one k-point, the n-neighbours are its periodic images
 
     dnn0 = 0.0_dp
     dnn1 = eta
     ndnntot = 0
     DO nlist = 1, search_shells
        DO nkp = 1, num_kpts
           DO loop=1,(2*nsupcell+1)**3
              l=lmn(1,loop);m=lmn(2,loop);n=lmn(3,loop)
              !
              vkpp=kpt_cart(:,nkp)+MATMUL(lmn(:,loop),recip_lattice)
              dist= SQRT( (kpt_cart(1,1)-vkpp(1))**2 &
                   + (kpt_cart(2,1)-vkpp(2))**2 + (kpt_cart(3,1)-vkpp(3))**2 )
              !
              IF ( (dist.gt.kmesh_tol) .AND. (dist.gt.dnn0 + kmesh_tol) ) THEN
                 IF(dist.lt.dnn1-kmesh_tol) THEN
                    dnn1=dist  ! found a closer shell
                    counter=0
                 END IF
                 IF(dist.gt.(dnn1-kmesh_tol) .AND. dist.lt.(dnn1+kmesh_tol)) THEN
                    counter=counter+1 ! count the multiplicity of the shell
                 END IF
              END IF
           ENDDO
        ENDDO
        IF (dnn1.lt.eta - kmesh_tol) ndnntot = ndnntot + 1
        dnn(nlist) = dnn1
        multi(nlist)=counter
        dnn0 = dnn1
        dnn1 = eta
     ENDDO
 
     WRITE(stdout,'(1x,a)') '+----------------------------------------------------------------------------+'
     WRITE(stdout,'(1x,a)') '|                    Distance to Nearest-Neighbour Shells                    |'
     WRITE(stdout,'(1x,a)') '|                    ------------------------------------                    |'
     WRITE(stdout,'(1x,a)') '|          Shell             Distance (Bohr^-1)         Multiplicity         |'
     WRITE(stdout,'(1x,a)') '|          -----             ------------------         ------------         |'
     DO ndnn = 1, ndnntot
        WRITE(stdout,'(1x,a,11x,i3,17x,f10.6,19x,i4,12x,a)') '|',ndnn,dnn(ndnn)/lenconfac,multi(ndnn),'|'
     ENDDO
     WRITE(stdout,'(1x,a)') '+----------------------------------------------------------------------------+'
 
     ! Get the shell weights to satisfy the B1 condition
     CALL kmesh_shell_automatic(multi,dnn,bweight)

     WRITE(stdout,'(1x,a)',advance='no') '| The following shells are used: '
     DO ndnn=1,num_shells
        IF (ndnn.eq.num_shells) THEN
           WRITE(stdout,'(i3,1x)',advance='no') shell_list(ndnn)
        ELSE
           WRITE(stdout,'(i3,",")',advance='no') shell_list(ndnn)
        ENDIF
     ENDDO
     DO l=1,11-num_shells
        WRITE(stdout,'(4x)',advance='no')
     ENDDO
     WRITE(stdout,'("|")')

     nntot=0
     DO loop_s=1,num_shells
        nntot=nntot+multi(shell_list(loop_s))
     END DO
     IF(nntot>num_nnmax) THEN
     WRITE(stdout,'(a,i2,a)') ' **WARNING: kmesh has found >',num_nnmax,' nearest neighbours**'
     WRITE(stdout,'(a)') ' '
     WRITE(stdout,'(a)') ' This is probably caused by an error in your unit cell specification'
     WRITE(stdout,'(a)') ' '
 
     ALLOCATE( bvec_tmp(3,MAXVAL(multi)))
     bvec_tmp=0.0_dp
     counter=0
     DO shell=1,search_shells
        CALL kmesh_get_bvectors(multi(shell),1,dnn(shell),bvec_tmp(:,1:multi(shell)))
        DO loop=1,multi(shell)
           counter=counter+1
           WRITE(stdout,'(a,I4,1x,a,2x,3f12.6,2x,a,2x,f12.6,a)') ' | b-vector  ',counter,': (', &
                bvec_tmp(:,loop)/lenconfac,')',dnn(shell)/lenconfac,'  |'
        END DO
     END DO
     WRITE(stdout,'(a)') ' '
     DEALLOCATE( bvec_tmp)
     CPABORT('kmesh_get: something wrong, found too many nearest neighbours')
     END IF
 
     ALLOCATE(nnlist(num_kpts,nntot))
     ALLOCATE(neigh(num_kpts,nntot/2))
     ALLOCATE(nncell(3,num_kpts,nntot))
 
     ALLOCATE(wb(nntot))
     ALLOCATE(bka(3,nntot/2))
     ALLOCATE(bk(3,nntot,num_kpts))
 
     nnx=0
     DO loop_s=1,num_shells
        DO loop_b=1,multi(shell_list(loop_s))
           nnx=nnx+1
           wb_local(nnx)=bweight(loop_s)
        END DO
     END DO
 
     ! Now build up the list of nearest-neighbour shells for each k-point.
     ! nnlist(nkp,1...nnx) points to the nnx neighbours (ordered along increa
     ! shells) of the k-point nkp. nncell(i,nkp,nnth) tells us in which BZ is
     ! nnth nearest-neighbour of the k-point nkp. Construct the nnx b-vectors
     ! go from k-point nkp to each neighbour bk(1:3,nkp,1...nnx). 
     ! Comment: Now we have bk(3,nntot,num_kps) 09/04/2006
 
     WRITE(stdout,'(1x,a)') '+----------------------------------------------------------------------------+'
     WRITE(stdout,'(1x,a)') '|                        Shell   # Nearest-Neighbours                        |'
     WRITE(stdout,'(1x,a)') '|                        -----   --------------------                        |'
     !
     ! Standard routine
     !
      nnshell=0
    DO nkp = 1, num_kpts
        nnx = 0
        ok: DO ndnnx = 1, num_shells
           ndnn = shell_list(ndnnx)
           DO loop=1,(2*nsupcell+1)**3
              l=lmn(1,loop);m=lmn(2,loop);n=lmn(3,loop)
              vkpp2=MATMUL(lmn(:,loop),recip_lattice)
              DO nkp2 = 1, num_kpts
                 vkpp=vkpp2+kpt_cart(:,nkp2)
                 dist= SQRT( (kpt_cart(1,nkp)-vkpp(1))**2 &
                      + (kpt_cart(2,nkp)-vkpp(2))**2 + (kpt_cart(3,nkp)-vkpp(3)) **2)
                 IF ( (dist.ge.dnn(ndnn)*(1-kmesh_tol)) .AND. (dist.le.dnn(ndnn)*(1+kmesh_tol)) ) THEN
                    nnx = nnx + 1
                    nnshell(nkp,ndnn) = nnshell(nkp,ndnn) + 1
                    nnlist(nkp,nnx) = nkp2
                    nncell(1,nkp,nnx) = l
                    nncell(2,nkp,nnx) = m
                    nncell(3,nkp,nnx) = n
                    bk_local(:,nnx,nkp) = vkpp(:) - kpt_cart(:,nkp)
                 ENDIF
                 !if we have the right number of neighbours we can exit
                 IF(nnshell(nkp,ndnn)==multi(ndnn)) CYCLE ok
              ENDDO
           ENDDO
           ! check to see if too few neighbours here
        END DO ok
 
     END DO
 
     DO ndnnx=1, num_shells
        ndnn = shell_list(ndnnx)
        WRITE(stdout,'(1x,a,24x,i3,13x,i3,33x,a)') '|',ndnn,nnshell(1,ndnn),'|'
     END DO
     WRITE(stdout,'(1x,"+",76("-"),"+")')
 
 
     DO nkp = 1, num_kpts
        nnx = 0
        DO ndnnx = 1, num_shells
           ndnn = shell_list(ndnnx)
           DO nnsh = 1, nnshell(nkp,ndnn)
              bb1 = 0.0_dp
              bbn = 0.0_dp
              nnx = nnx + 1
              DO i = 1, 3
                 bb1 = bb1 + bk_local(i,nnx,1) * bk_local(i,nnx,1)
                 bbn = bbn + bk_local(i,nnx,nkp) * bk_local(i,nnx,nkp)
              ENDDO
              IF (ABS(SQRT(bb1)-SQRT(bbn)).gt.kmesh_tol) THEN
                 WRITE(stdout,'(1x,2f10.6)') bb1,bbn
                 CPABORT('Non-symmetric k-point neighbours!')
              ENDIF
           ENDDO
        ENDDO
     ENDDO
 
 
     ! now check that the completeness relation is satisfied for every kpoint
     ! We know it is true for kpt=1; but we check the rest to be safe.
     ! Eq. B1 in Appendix  B PRB 56 12847 (1997)
 
        DO nkp = 1, num_kpts
           DO i = 1, 3
              DO j = 1, 3
                 ddelta = 0.0_dp
                 nnx = 0
                 DO ndnnx = 1, num_shells
                    ndnn = shell_list(ndnnx)
                    DO nnsh = 1, nnshell(1,ndnn)
                       nnx = nnx + 1
                       ddelta = ddelta + wb_local(nnx) * bk_local(i,nnx,nkp) * bk_local(j,nnx,nkp)
                    ENDDO
                 ENDDO
                 IF ( (i.eq.j) .AND. (ABS(ddelta-1.0_dp).gt.kmesh_tol) ) THEN
                    WRITE(stdout,'(1x,2i3,f12.8)') i,j,ddelta
                    CPABORT('Eq. (B1) not satisfied in kmesh_get (1)')
                 ENDIF
                 IF ( (i.ne.j) .AND. (ABS(ddelta).gt.kmesh_tol) ) THEN
                    WRITE(stdout,'(1x,2i3,f12.8)') i,j,ddelta
                    CPABORT('Eq. (B1) not satisfied in kmesh_get (2)')
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
 
     WRITE(stdout,'(1x,a)') '| Completeness relation is fully satisfied [Eq. (B1), PRB 56, 12847 (1997)]  |'
     WRITE(stdout,'(1x,"+",76("-"),"+")')
 
     !
     wbtot = 0.0_dp
     nnx = 0
     DO ndnnx = 1, num_shells
        ndnn = shell_list(ndnnx)
        DO nnsh = 1, nnshell(1,ndnn)
           nnx = nnx + 1
           wbtot = wbtot + wb_local(nnx)
        ENDDO
     ENDDO
 
     nnh = nntot / 2
     ! make list of bka vectors from neighbours of first k-point
     ! delete any inverse vectors as you collect them
     na = 0
     DO nn = 1, nntot
        ifound = 0
        IF (na.ne.0) THEN
           DO nap = 1, na
              CALL utility_compar(bka(1,nap),bk_local(1,nn,1),ifpos,ifneg)
              IF (ifneg.eq.1) ifound = 1
           ENDDO
        ENDIF
        IF (ifound.eq.0) THEN
           !         found new vector to add to set
           na = na + 1
           bka(1,na) = bk_local(1,nn,1)
           bka(2,na) = bk_local(2,nn,1)
           bka(3,na) = bk_local(3,nn,1)
        ENDIF
     ENDDO
     IF (na.ne.nnh) CPABORT('Did not find right number of bk directions')
 
     WRITE(stdout,'(1x,a)') '|                 b_k Vectors (Bohr^-1) and Weights (Bohr^2)                 |'
     WRITE(stdout,'(1x,a)') '|                 ------------------------------------------                 |'
     WRITE(stdout,'(1x,a)') '|            No.         b_k(x)      b_k(y)      b_k(z)        w_b           |'
     WRITE(stdout,'(1x,a)') '|            ---        --------------------------------     --------        |'
     DO i = 1, nntot
        WRITE (stdout,'(1x,"|",11x,i3,5x,3f12.6,3x,f10.6,8x,"|")') &
             i,(bk_local(j,i,1)/lenconfac,j=1,3),wb_local(i)*lenconfac**2
     ENDDO
     WRITE(stdout,'(1x,"+",76("-"),"+")')
     WRITE(stdout,'(1x,a)') '|                           b_k Directions (Bohr^-1)                         |'
     WRITE(stdout,'(1x,a)') '|                           ------------------------                         |'
     WRITE(stdout,'(1x,a)') '|            No.           x           y           z                         |'
     WRITE(stdout,'(1x,a)') '|            ---        --------------------------------                     |'
     DO i = 1, nnh
        WRITE(stdout,'(1x,"|",11x,i3,5x,3f12.6,21x,"|")') i,(bka(j,i)/lenconfac,j=1,3)
     ENDDO
     WRITE(stdout,'(1x,"+",76("-"),"+")')
     WRITE(stdout,*) ' '


    ! find index array
    DO nkp = 1, num_kpts
       DO na = 1, nnh
          ! first, zero the index array so we can check it gets filled
          neigh(nkp,na) = 0
          ! now search through list of neighbours of this k-point
          DO nn = 1, nntot
             CALL utility_compar(bka(1,na),bk_local(1,nn,nkp),ifpos,ifneg)
             IF (ifpos.eq.1) neigh(nkp,na)=nn
          ENDDO
          ! check found
          IF (neigh(nkp,na).eq.0) THEN
             WRITE(stdout,*) ' nkp,na=',nkp,na
             CPABORT('kmesh_get: failed to find neighbours for this kpoint')
          ENDIF
       ENDDO
    ENDDO


    !fill in the global arrays from the local ones
    DO loop=1,nntot
       wb(loop)=wb_local(loop)
    END DO

    DO loop_s=1,num_kpts
       DO loop=1,nntot
          bk(:,loop,loop_s)=bk_local(:,loop,loop_s)
       END DO
    END DO

  END SUBROUTINE w90_kmesh_get

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
  SUBROUTINE w90_kmesh_supercell_sort
    !==================================================================!
    !                                                                  !
    ! We look for kpoint neighbours in a large supercell of reciprocal !
    ! unit cells. Done sequentially this is very slow.                 !
    ! Here we order the cells by the distance from the origin          !
    ! Doing the search in this order gives a dramatic speed up         !
    !                                                                  !
    !==================================================================!  
      INTEGER                                            :: counter, indx(1), l, &
                                                            lmn_cp( 3,(2*nsupcell+1)**3), loop, m, &
                                                            n
      REAL(kind=dp)                                      :: dist((2*nsupcell+1)**3), &
                                                            dist_cp((2*nsupcell+1)**3), pos(3)

    counter=1
    lmn(:,counter)=0
    dist(counter)=0.0_dp
    DO l = -nsupcell,  nsupcell
       DO m = -nsupcell,   nsupcell
          DO n = -nsupcell,  nsupcell
             IF(l==0 .AND. m==0 .AND. n==0) CYCLE
             counter=counter+1
             lmn(1,counter)=l;lmn(2,counter)=m;lmn(3,counter)=n
             pos=MATMUL(lmn(:,counter),recip_lattice)
             dist(counter)=SQRT(DOT_PRODUCT(pos,pos))
          END DO
       END DO
    END DO

    DO loop=(2*nsupcell+1)**3,1,-1
       indx=internal_maxloc(dist)
       dist_cp(loop)=dist(indx(1))
       lmn_cp(:,loop)=lmn(:,indx(1))
       dist(indx(1))=-1.0_dp
    END DO

    lmn=lmn_cp
    dist=dist_cp

  END SUBROUTINE w90_kmesh_supercell_sort

! **************************************************************************************************
!> \brief ...
!> \param dist ...
!> \retval internal_maxloc ...
! **************************************************************************************************
     FUNCTION internal_maxloc(dist)
    !=========================================================================!
    !                                                                         !
    !  A predictable maxloc.                                                  !
    !                                                                         !
    !=========================================================================!


      REAL(kind=dp), INTENT(in)                          :: dist((2*nsupcell+1)**3)
      INTEGER                                            :: internal_maxloc

      INTEGER                                            :: counter, guess(1), &
                                                            list((2*nsupcell+1)**3), loop

    list=0
    counter=1

    guess=MAXLOC(dist)
    list(1)=guess(1)
    ! look for any degenerate values
    DO loop=1,(2*nsupcell+1)**3
       IF (loop==guess(1)) CYCLE
       IF ( ABS(dist(loop)-dist(guess(1))) < eps8 ) THEN
          counter=counter+1
          list(counter)=loop
       ENDIF
    END DO
    ! and always return the lowest index
    internal_maxloc=MINVAL(list(1:counter))

    END FUNCTION internal_maxloc

! **************************************************************************************************
!> \brief ...
!> \param multi ...
!> \param dnn ...
!> \param bweight ...
! **************************************************************************************************
    SUBROUTINE kmesh_shell_automatic(multi,dnn,bweight)
    !==========================================================================!
    !                                                                          !
    ! Find the correct set of shells to satisfy B1                             !
    !  The stratagy is:                                                        !
    !        Take the bvectors from the next shell                             !
    !        Reject them if they are parallel to exisiting b vectors           !
    !        Test to see if we satisfy B1, if not add another shell and repeat !
    !                                                                          !
    !==========================================================================!  

      INTEGER, INTENT(in)                                :: multi(search_shells)
      REAL(kind=dp), INTENT(in)                          :: dnn(search_shells)
      REAL(kind=dp), INTENT(out)                         :: bweight(max_shells)

      INTEGER, PARAMETER                                 :: lwork = max_shells*10
      REAL(kind=dp), PARAMETER :: TARGET(6) = (/1.0_dp,1.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp/)

      INTEGER                                            :: cur_shell, info, loop_b, loop_bn, &
                                                            loop_i, loop_j, loop_s, shell
      LOGICAL                                            :: b1sat, lpar
      REAL(kind=dp)                                      :: delta, work(lwork)
      REAL(kind=dp), ALLOCATABLE                         :: bvector(:,:,:)
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: singv
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: amat, smat, umat, vmat

    ALLOCATE( bvector(3,MAXVAL(multi),max_shells))
    bvector=0.0_dp;bweight=0.0_dp

    WRITE(stdout,'(1x,a)') '| The b-vectors are chosen automatically                                     |'

    b1sat=.FALSE.
    DO shell=1,search_shells
       cur_shell=num_shells+1

       ! get the b vectors for the new shell
       CALL kmesh_get_bvectors(multi(shell),1,dnn(shell),bvector(:,1:multi(shell),cur_shell))

       ! We check that the new shell is not parrallel to an existing shell (cosine=1)
       lpar=.FALSE.
       IF(num_shells>0) THEN
          DO loop_bn=1,multi(shell)
             DO loop_s=1,num_shells
                DO loop_b=1,multi(shell_list(loop_s))
                   delta=DOT_PRODUCT(bvector(:,loop_bn,cur_shell),bvector(:,loop_b,loop_s))/ &
                        SQRT(DOT_PRODUCT(bvector(:,loop_bn,cur_shell),bvector(:,loop_bn,cur_shell))* &
                        DOT_PRODUCT(bvector(:,loop_b,loop_s),bvector(:,loop_b,loop_s)))
                   IF(ABS(ABS(delta)-1.0_dp)<eps6) lpar=.TRUE.
                END DO
             END DO
          END DO
       END IF

       IF(lpar) THEN
          IF(iprint>=3) THEN
             WRITE(stdout,'(1x,a)') '| This shell is linearly dependent on existing shells: Trying next shell     |'
          END IF
          CYCLE
       END IF

       num_shells=num_shells+1
       shell_list(num_shells)=shell

       ALLOCATE(amat(max_shells,num_shells))
       ALLOCATE(umat(max_shells,max_shells))
       ALLOCATE(vmat(num_shells,num_shells))
       ALLOCATE(smat(num_shells,max_shells))
       ALLOCATE(singv(num_shells))
       amat=0.0_dp;umat=0.0_dp;vmat=0.0_dp;smat=0.0_dp;singv=0.0_dp

       amat=0.0_dp
       DO loop_s=1,num_shells
          DO loop_b=1,multi(shell_list(loop_s))
             amat(1,loop_s)=amat(1,loop_s)+bvector(1,loop_b,loop_s)*bvector(1,loop_b,loop_s)
             amat(2,loop_s)=amat(2,loop_s)+bvector(2,loop_b,loop_s)*bvector(2,loop_b,loop_s)
             amat(3,loop_s)=amat(3,loop_s)+bvector(3,loop_b,loop_s)*bvector(3,loop_b,loop_s)
             amat(4,loop_s)=amat(4,loop_s)+bvector(1,loop_b,loop_s)*bvector(2,loop_b,loop_s)
             amat(5,loop_s)=amat(5,loop_s)+bvector(2,loop_b,loop_s)*bvector(3,loop_b,loop_s)
             amat(6,loop_s)=amat(6,loop_s)+bvector(3,loop_b,loop_s)*bvector(1,loop_b,loop_s)
          END DO
       END DO

       info=0
       CALL dgesvd('A','A',max_shells,num_shells,amat,max_shells,singv,umat,max_shells,vmat,num_shells,work,lwork,info)
       IF(info<0) THEN
          WRITE(stdout,'(1x,a,1x,I1,1x,a)') 'kmesh_shell_automatic: Argument',ABS(info),'of dgesvd is incorrect'
          CPABORT('kmesh_shell_automatic: Problem with Singular Value Decomposition')
       ELSE IF (info>0) THEN
          CPABORT('kmesh_shell_automatic: Singular Value Decomposition did not converge')
       END IF

       IF(ANY(ABS(singv)<eps5)) THEN
         IF(num_shells==1)  THEN
             CALL cp_abort(__LOCATION__,"kmesh_shell_automatic: "//&
                 "Singular Value Decomposition has found a very small singular value.")
         ELSE
            WRITE(stdout,'(1x,a)') '| SVD found small singular value, Rejecting this shell and trying the next   |'
            b1sat=.FALSE.
            num_shells=num_shells-1
            GOTO 200
         END IF
       END IF


       smat=0.0_dp
       DO loop_s=1,num_shells
          smat(loop_s,loop_s)=1/singv(loop_s)
       END DO

       bweight(1:num_shells)=MATMUL(TRANSPOSE(vmat),MATMUL(smat,MATMUL(TRANSPOSE(umat),TARGET)))
       IF(iprint>=2) THEN
          DO loop_s=1,num_shells
             WRITE(stdout,'(1x,a,I2,a,f12.7,5x,a8,36x,a)') '| Shell: ',loop_s,&
                  ' w_b ', bweight(loop_s)*lenconfac**2,'('//TRIM(length_unit)//'^2)','|'
          END DO
       END IF

       !check b1
       b1sat=.TRUE.
       DO loop_i=1,3
          DO loop_j=loop_i,3
             delta=0.0_dp
             DO loop_s=1,num_shells
                DO loop_b=1,multi(shell_list(loop_s))
                   delta=delta+bweight(loop_s)*bvector(loop_i,loop_b,loop_s)*bvector(loop_j,loop_b,loop_s)
                END DO
             END DO
             IF(loop_i==loop_j) THEN
                IF(ABS(delta-1.0_dp)>kmesh_tol) b1sat=.FALSE.
             END IF
             IF(loop_i/=loop_j) THEN
                IF(ABS(delta)>kmesh_tol) b1sat=.FALSE.
             END IF
          END DO
       END DO

       IF(.not.b1sat) THEN
          IF(shell<search_shells .AND. iprint>=3) THEN
             WRITE(stdout,'(1x,a,24x,a1)') '| B1 condition is not satisfied: Adding another shell','|'
          ELSEIF(shell==search_shells) THEN
             WRITE(stdout,*) ' '
       WRITE(stdout,'(1x,a,i3,a)') 'Unable to satisfy B1 with any of the first ',search_shells,' shells'
       WRITE(stdout,'(1x,a)') 'Your cell might be very long, or you may have an irregular MP grid'
       WRITE(stdout,'(1x,a)') 'Try increasing the parameter search_shells in the win file (default=12)'
       WRITE(stdout,*) ' '
             CPABORT('kmesh_get_automatic')
          END IF
       END IF

200 CONTINUE

       DEALLOCATE(amat)
       DEALLOCATE(umat)
       DEALLOCATE(vmat)
       DEALLOCATE(smat)
       DEALLOCATE(singv)

       IF(b1sat) EXIT

    END DO

    IF(.NOT. b1sat)  THEN
       WRITE(stdout,*) ' '
       WRITE(stdout,'(1x,a,i3,a)') 'Unable to satisfy B1 with any of the first ',search_shells,' shells'
       WRITE(stdout,'(1x,a)') 'Your cell might be very long, or you may have an irregular MP grid'
       WRITE(stdout,'(1x,a)') 'Try increasing the parameter search_shells in the win file (default=12)'
       WRITE(stdout,*) ' '
       CPABORT('kmesh_get_automatic')
    END IF

  END SUBROUTINE kmesh_shell_automatic

! **************************************************************************************************
!> \brief ...
!> \param multi ...
!> \param kpt ...
!> \param shell_dist ...
!> \param bvector ...
! **************************************************************************************************
    SUBROUTINE kmesh_get_bvectors(multi,kpt,shell_dist,bvector)
    !==================================================================!
    !                                                                  !
    ! Returns the bvectors for a given shell and kpoint                ! 
    !                                                                  !
    !===================================================================  

      INTEGER, INTENT(in)                                :: multi, kpt
      REAL(kind=dp), INTENT(in)                          :: shell_dist
      REAL(kind=dp), INTENT(out)                         :: bvector(3,multi)

      INTEGER                                            :: loop, nkp2, num_bvec
      REAL(kind=dp)                                      :: dist, vkpp(3), vkpp2(3)

    bvector=0.0_dp

    num_bvec=0
    ok:    DO loop=1,(2*nsupcell+1)**3
       vkpp2=MATMUL(lmn(:,loop),recip_lattice)
       DO nkp2 = 1, num_kpts
          vkpp=vkpp2+kpt_cart(:,nkp2)
          dist= SQRT( (kpt_cart(1,kpt)-vkpp(1))**2 &
               + (kpt_cart(2,kpt)-vkpp(2))**2 + (kpt_cart(3,kpt)-vkpp(3)) **2)
          IF ( (dist.ge.shell_dist*(1.0_dp-kmesh_tol)) .AND. dist.le.shell_dist*(1.0_dp+kmesh_tol) ) THEN
             num_bvec = num_bvec + 1
             bvector(:,num_bvec) = vkpp(:) - kpt_cart(:,kpt)
          ENDIF
          !if we have the right number of neighbours we can exit
          IF(num_bvec==multi) CYCLE ok
       ENDDO
    ENDDO ok

    IF(num_bvec<multi) CPABORT('kmesh_get_bvector: Not enough bvectors found')

  END SUBROUTINE kmesh_get_bvectors

! **************************************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \param ifpos ...
!> \param ifneg ...
! **************************************************************************************************
  SUBROUTINE utility_compar(a,b,ifpos,ifneg)
      REAL(kind=dp), INTENT(in)                          :: a(3), b(3)
      INTEGER, INTENT(out)                               :: ifpos, ifneg

      REAL(kind=dp)                                      :: rrm, rrp

    rrp=(a(1)-b(1))**2+(a(2)-b(2))**2+(a(3)-b(3))**2
    rrm=(a(1)+b(1))**2+(a(2)+b(2))**2+(a(3)+b(3))**2
    ifpos=0
    IF (ABS(rrp).lt.eps8) ifpos=1
    ifneg=0
    IF (ABS(rrm).lt.eps8) ifneg=1

  END SUBROUTINE utility_compar

! **************************************************************************************************

END MODULE wannier90
