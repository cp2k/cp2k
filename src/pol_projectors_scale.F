!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!***** cp2k/pol_projectors_scale [1.0] *
!!
!!   NAME
!!     pol_projectors_scale
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     gloria
!!
!!   MODIFICATION HISTORY
!!
!!   NOTES
!!
!!
!!   SOURCE
!******************************************************************************

MODULE pol_projectors_scale

  USE ai_overlap,                      ONLY: overlap
  USE ao_types,                        ONLY: ao_type
  USE atomic_kinds,                    ONLY: kind_info_type
  USE kinds,                           ONLY: dbl
  USE orbital_pointers,                ONLY: nco,&
                                             ncoset
  USE particle_types,                  ONLY: particle_type
  USE termination,                     ONLY: stop_memory

  IMPLICIT NONE

  PRIVATE

  PUBLIC ::  get_proj_scale_factor

  CONTAINS

!******************************************************************************

 SUBROUTINE get_proj_scale_factor (ao, ki, part)

    TYPE(ao_type), intent(inout)             :: ao
    TYPE(kind_info_type), DIMENSION(:), &
      POINTER                                :: ki
    TYPE(particle_type), DIMENSION(:), &
      intent(in)                             :: part

    INTEGER :: first_cgf, icgf, icoef, ii, ikind, ios, ipart, iset, ishell, &
      l_max, l_min, last_cgf, maxcgf, maxco, natoms, ncgf, nco, nkind, npgf, &
      nset, nshell
    REAL(dbl), ALLOCATABLE, DIMENSION(:, :)  :: s_block, sab, work
    REAL(dbl), DIMENSION(:), POINTER         :: rpgf, zet

 nkind = size ( ki )

 maxco = 0
 maxcgf = 0

 DO ikind=1,nkind
   DO iset=1,ki(ikind)%orb_basis_set%nset
     maxcgf = MAX(maxcgf,ki(ikind)%orb_basis_set%ncgf)
     npgf = ki(ikind)%orb_basis_set%npgf(iset)
     l_max = ki(ikind)%orb_basis_set%lmax(iset)
     maxco = MAX(maxco,npgf*ncoset(l_max))
   END DO
 END DO

 ALLOCATE (sab(maxco,maxco),STAT=ios)
 IF (ios /= 0) THEN
   CALL stop_memory("get_proj_scale_factor","sab",maxco*maxco)
 END IF
 sab(:,:) = 0.0_dbl

 ALLOCATE (work(maxco,maxcgf),STAT=ios)
 IF (ios /= 0) THEN
    CALL stop_memory("get_proj_scale_factor","work",maxco*maxcgf)
 END IF
 work(:,:) = 0.0_dbl

 ALLOCATE (s_block(maxcgf,maxcgf),STAT=ios)
 IF (ios /= 0) THEN
    CALL stop_memory("get_proj_scale_factor","s_blok",maxcgf*maxcgf)
 END IF
 s_block(:,:) = 0.0_dbl

 DO ikind = 1 , size (ki)

    natoms = ki (ikind ) % natom

    DO ii = 1, natoms

       ipart = ki(ikind) % atom_list(ii)
       nset = ki(ikind) % orb_basis_set % nset

       IF (ii==1) THEN
! calculate the norm (we need to calculate it once for each atom kind) 
          DO iset = 1, nset
             l_max = ki(ikind)%orb_basis_set%lmax(iset)
             l_min = ki(ikind)%orb_basis_set%lmin(iset)
             npgf = ki(ikind)%orb_basis_set%npgf(iset)
             nshell = ki(ikind)%orb_basis_set%nshell(iset)
             rpgf => ki(ikind)%orb_basis_set%pgf_radius(1:npgf,iset)
             zet => ki(ikind)%orb_basis_set%zet(1:npgf,iset)
             first_cgf = ki(ikind) % orb_basis_set %first_cgf(1,iset)
             last_cgf = ki(ikind) % orb_basis_set %last_cgf(nshell,iset)
             ncgf = last_cgf - first_cgf + 1
             nco = npgf*ncoset(l_max)
             CALL get_norm (l_max,zet,rpgf,l_min,nco,ncgf,sab, &
                           work,s_block(first_cgf:last_cgf,first_cgf:last_cgf), &
                           ki(ikind)%orb_basis_set%cphi(1:nco,first_cgf:last_cgf))
          
          END DO
       END IF

       DO iset = 1, nset
          DO ishell=1,ki(ikind) % orb_basis_set % nshell(iset)
             first_cgf = ki(ikind) % orb_basis_set %first_cgf(ishell,iset)
             last_cgf = ki(ikind) % orb_basis_set %last_cgf(ishell,iset)
             DO icgf = first_cgf, last_cgf
                icoef  = part ( ipart ) % coef_list( icgf )
                ao % norm ( icoef ) = SQRT(1.0_dbl/s_block(icgf,icgf))
             END DO
          END DO
       END DO

    END DO
    s_block (:,:) = 0.0_dbl
 END DO

 IF (ALLOCATED(sab)) THEN
   DEALLOCATE (sab,STAT=ios)
   IF (ios /= 0) THEN
     CALL stop_memory("get_proj_scale_factor","sab")
   END IF
 END IF
 IF (ALLOCATED(s_block)) THEN
   DEALLOCATE (s_block,STAT=ios)
   IF (ios /= 0) THEN
     CALL stop_memory("get_proj_scale_factor","s_block")
   END IF
 END IF
 IF (ALLOCATED(work)) THEN
   DEALLOCATE (work,STAT=ios)
   IF (ios /= 0) THEN
     CALL stop_memory("get_proj_scale_factor","work")
   END IF
 END IF

 END SUBROUTINE get_proj_scale_factor

!---------------------------------------------------------------------------------
  SUBROUTINE get_norm(l_max,zet,rpgf,l_min,nco,ncgf,sab, &
                              work,s_block,cphi)
!---------------------------------------------------------------------------------


    INTEGER, INTENT(IN)                      :: l_max
    REAL(dbl), DIMENSION(:), POINTER         :: zet, rpgf
    INTEGER, INTENT(IN)                      :: l_min, nco, ncgf
    REAL(dbl), DIMENSION(:, :), &
      intent(inout)                          :: sab, work, s_block
    REAL(dbl), DIMENSION(:, :), intent(in)   :: cphi

    INTEGER                                  :: lda, ldb, ldc
    REAL(dbl)                                :: rijsq
    REAL(dbl), DIMENSION(3)                  :: rij

!---------------------------------------------------------------------------

   sab = 0._dbl
   rij = 0._dbl
   rijsq = 0.0_dbl 

!  calculate the primitive overlap integral

   CALL overlap(l_max,zet,rpgf,l_min,&
                        l_max,zet,rpgf,l_min,&
                         rij,rijsq,sab)

!  Contraction step (overlap matrix)

   lda = size(sab,1)
   ldb = size(cphi,1)
   ldc = size(work,1)

   CALL dgemm("N","N",nco,ncgf,nco,1.0_dbl,sab(1,1),lda,&
                       cphi(1,1), &
                       ldb,0.0_dbl,&
                       work(1,1),ldc)
   lda = size(cphi,1)
   ldb = size(work,1)
   ldc = size(s_block,1)

   CALL dgemm("T","N",ncgf,ncgf,nco,1.0_dbl,&
                       cphi,&
                       lda,&
                       work,ldb,0.0_dbl,&
                       s_block,ldc)

!   atom_block(:,:) = atom_block(:,:) + s_block(:,:)

  END SUBROUTINE get_norm
!---------------------------------------------------------------------------------
!*****************************************************************************
 END MODULE pol_projectors_scale
!*****************************************************************************
