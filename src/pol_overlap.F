!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/pol_overlap [1.0] *
!!
!!   NAME
!!     pol_overlap
!!
!!   FUNCTION
!!     Calculation of the second order contribution of the Hohenbeg-Kohn
!!     functional that depends on the overlaps between basis functions
!!
!!     E (overl) = 1/2 sum_ij  c_i*c_j (f_i + f_j) <phi_i|phi_j>
!!
!!   AUTHOR
!!     gloria (06 Jun 2001)     
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE pol_overlap

! *****************************************************************************

  USE ai_derivatives,                  ONLY: dabdr
  USE ai_overlap,                      ONLY: overlap
  USE ao_types,                        ONLY: ao_type
  USE atomic_kinds,                    ONLY: kind_info_type
  USE empirical_parameters,            ONLY: empirical_parameter_type
  USE fist_nonbond_force,              ONLY: find_image
  USE kinds,                           ONLY: dbl
  USE method_specifications,           ONLY: maxder
  USE molecule_types,                  ONLY: linklist_images,&
                                             linklist_neighbor,&
                                             particle_node_type
  USE orbital_pointers,                ONLY: ncoset
  USE particle_types,                  ONLY: particle_type
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell_param
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             allocate_matrix,&
                                             deallocate_matrix,&
                                             get_block_node,&
                                             real_matrix_p_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE
 
  PUBLIC :: force_overlap

  TYPE atom_basis_info_type
    INTEGER :: ipart, ikind
    INTEGER :: first_cgf, last_cgf
  END TYPE atom_basis_info_type

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE force_overlap ( ao, kind_info, part, pnode, box, empparm, &
                           energy, fo_coef, fo_part )

!   ***************************************************************************

    TYPE(ao_type), INTENT(INOUT)             :: ao
    TYPE(kind_info_type), DIMENSION(:), &
      INTENT(IN)                             :: kind_info
    TYPE(particle_type), dimension(:), &
      INTENT(IN)                             :: part
    TYPE(particle_node_type), dimension(:), &
      INTENT(INOUT)                          :: pnode
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(empirical_parameter_type), &
      DIMENSION(:), INTENT(IN)               :: empparm
    REAL(dbl), INTENT(OUT)                   :: energy
    REAL(dbl), DIMENSION(:), INTENT(INOUT)   :: fo_coef
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL                :: fo_part

    INTEGER :: first_cgfa, first_cgfb, i, iatom, id, ikind, im, ipart, iset, &
      istat, j, jatom, jkind, jpart, jset, la_max, la_min, last_cgfa, &
      last_cgfb, lb_max, lb_min, maxcgf, maxco, maxdim, nat, natom, ncgf, &
      ncgfa, ncgfb, ncoa, ncob, nkind, npart, npgfa, npgfb, nseta, nsetb, &
      nshella, nshellb
    INTEGER, ALLOCATABLE, &
      DIMENSION(:, :, :), SAVE               :: n_images
    LOGICAL                                  :: deriv, first_time, forces
    REAL(dbl)                                :: r2max, radius_a, radius_b, &
                                                rcut, rijsq
    REAL(dbl), ALLOCATABLE, DIMENSION(:, :)  :: sab, sabdx, sabdy, sabdz, work
    REAL(dbl), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: atom_block
    REAL(dbl), DIMENSION(3)                  :: cell_lengths, perd, quotient, &
                                                ri, rij, s, vec
    REAL(dbl), DIMENSION(3, 3)               :: h_inv, hmat
    REAL(dbl), DIMENSION(:), POINTER         :: rpgfa, rpgfb, zeta, zetb
    REAL(dbl), DIMENSION(:, :), POINTER      :: s_block
    TYPE(atom_basis_info_type), &
      DIMENSION(:), POINTER                  :: ainfo
    TYPE(linklist_images), POINTER           :: current_image
    TYPE(linklist_neighbor), POINTER         :: current_neighbor
    TYPE(real_matrix_p_type)                 :: smat

!   ---------------------------------------------------------------------------
!  initialize atomic info

  ncgf = 0
  natom = 0 
  nat = 0
  npart = size(pnode)
  nkind = size(kind_info)
  hmat = box % hmat 
  h_inv = box % h_inv
  perd = box % perd      
  energy = 0._dbl
  deriv = .FALSE.

! initialize energy
  energy=0._dbl

  DO ikind = 1, nkind
     nat = nat + size (kind_info(ikind) % atom_list)
  END DO
 
  NULLIFY (ainfo)

  IF (.NOT. ASSOCIATED (ainfo)) then
    ALLOCATE (ainfo(nat), stat=istat)
    IF (istat /= 0) THEN
      CALL stop_memory("pol_overlap","ainfo",nat)
    END IF
  END IF
 
  DO i = 1,npart
    ikind = pnode (i) % p % kind
    if (ikind == 0) CYCLE 
    natom = natom + 1
    pnode ( i ) % p % abase = natom
    ainfo ( natom ) % ipart = i
    ainfo ( natom ) % ikind = ikind
    pnode ( i ) % p % first_cgf = ncgf + 1
    ainfo ( natom ) % first_cgf =   pnode ( i ) % p % first_cgf
    ncgf = ncgf + kind_info(ikind)%orb_basis_set%ncgf
    pnode ( i ) % p % last_cgf = ncgf
    ainfo ( natom ) % last_cgf =   pnode ( i ) % p % last_cgf
  END DO                                             

  if (nat /= natom) call stop_program("pol_overlap","inconsistent number of atoms with basis")

! Allocate work storage 

  maxco = 0
  maxcgf = 0

  IF (maxder==1) THEN
    maxdim = 4
    forces = .TRUE.
  ELSE IF (maxder==0) THEN 
    maxdim = 1
    forces = .FALSE. 
  END If

  DO ikind=1,nkind
    DO iset=1,kind_info(ikind)%orb_basis_set%nset
      maxcgf = MAX(maxcgf,kind_info(ikind)%orb_basis_set%ncgf)
      npgfa = kind_info(ikind)%orb_basis_set%npgf(iset)
      la_max = kind_info(ikind)%orb_basis_set%lmax(iset) + maxder
      maxco = MAX(maxco,npgfa*ncoset(la_max))
    END DO
  END DO

  ALLOCATE (sab(maxco,maxco),STAT=istat)
  IF (istat /= 0) THEN
    CALL stop_memory("pol_overlap","sab",maxco*maxco)
  END IF
  sab(:,:) = 0.0_dbl

  ALLOCATE (sabdx(maxco,maxco),STAT=istat)
  IF (istat /= 0) THEN
    CALL stop_memory("pol_overlap","sabdx",maxco*maxco)
  END IF
  sabdx(:,:) = 0.0_dbl

  ALLOCATE (sabdy(maxco,maxco),STAT=istat)
  IF (istat /= 0) THEN
    CALL stop_memory("pol_overlap","sabdy",maxco*maxco)
  END IF
  sabdy(:,:) = 0.0_dbl

  ALLOCATE (sabdz(maxco,maxco),STAT=istat)
  IF (istat /= 0) THEN
    CALL stop_memory("pol_overlap","sabdz",maxco*maxco)
  END IF
  sabdz(:,:) = 0.0_dbl

  ALLOCATE (work(maxco,maxcgf),STAT=istat)
  IF (istat /= 0) THEN
  CALL stop_memory("pol_overlap","work",maxco*maxcgf)
  END IF
  work(:,:) = 0.0_dbl

  ALLOCATE (atom_block(maxco,maxcgf,maxdim),STAT=istat)
  IF (istat /= 0) THEN
  CALL stop_memory("pol_overlap","atom_blok",maxco*maxcgf*maxdim)
  END IF
  atom_block(:,:,:) = 0.0_dbl

  NULLIFY (smat%matrix)

  CALL build_block_overlap_matrix(pnode,part,ainfo,smat,natom,ncgf)

! initialize forces

  IF ( PRESENT ( fo_part ) ) fo_part = 0._dbl
  fo_coef = 0._dbl

! intialize images

  first_time = .NOT.ALLOCATED ( n_images )
  IF ( first_time ) THEN
     IF ( .NOT. ALLOCATED ( n_images ) ) &
          ALLOCATE ( n_images (  nkind, nkind, 3 ), STAT = istat )
     IF ( istat /= 0 ) CALL stop_memory ( 'force_env_types', &
          'n_images', nkind ** 2 * 3 )
 
     CALL get_cell_param ( box, cell_lengths )
 
     DO ikind = 1, nkind         
        DO jkind = 1,  nkind         
           rcut = kind_info( ikind ) % orb_basis_set % kind_radius + &
                  kind_info( jkind) % orb_basis_set % kind_radius 
           quotient ( : ) = rcut / cell_lengths ( : )
           DO id=1,3
              IF ( quotient ( id ) <= 0.5_dbl ) THEN
                 n_images ( ikind, jkind, id ) = 0
              ELSE
                 n_images ( ikind, jkind, id ) &
                      = CEILING ( rcut / cell_lengths ( id ) )
              END IF
           END DO
        END DO
     END DO
  END IF
                                                                                                      
!  Loop over all atoms 

  DO i = 1, size(pnode)
    
! cycle if the atoms has no basis function!

    ipart = pnode(i) %p%iatom
    ikind = part (ipart) % kind
    IF (ikind <= 0) CYCLE

    nseta = kind_info(ikind)%orb_basis_set%nset
      
    iatom = pnode (i) % p % abase
    ri = pnode (i) % p % r
      

! Loop over all neighbor atoms of the current atom "iatom" 

    current_neighbor => pnode ( i ) % sl       

    DO j = 1, pnode (i) % nsneighbor

      jpart = current_neighbor % index 
      jatom = part (jpart) % abase
      rij = current_neighbor % p % r - ri 
      jkind = current_neighbor % p % kind                    
      nsetb = kind_info(jkind)%orb_basis_set%nset
      IF (forces.AND.(iatom/=jatom)) deriv = .TRUE.

      DO iset=1,nseta

        la_max = kind_info(ikind)%orb_basis_set%lmax(iset)
        la_min = kind_info(ikind)%orb_basis_set%lmin(iset)

        npgfa = kind_info(ikind)%orb_basis_set%npgf(iset)
        nshella = kind_info(ikind)%orb_basis_set%nshell(iset)
        rpgfa => kind_info(ikind)%orb_basis_set%pgf_radius(1:npgfa,iset)
        zeta => kind_info(ikind)%orb_basis_set%zet(1:npgfa,iset)

        first_cgfa = kind_info(ikind)%orb_basis_set%first_cgf(1,iset)
        last_cgfa = kind_info(ikind)%orb_basis_set%last_cgf(nshella,iset)
        ncgfa = last_cgfa - first_cgfa + 1
        ncoa = npgfa*ncoset(la_max)

        radius_a = kind_info(ikind)%orb_basis_set%set_radius(iset)

        DO jset=1,nsetb

          radius_b = kind_info(jkind)%orb_basis_set%set_radius(jset)

          s(1) = h_inv(1,1)*rij(1) + h_inv(1,2)*rij(2) + h_inv(1,3)*rij(3)
          s(2) = h_inv(2,1)*rij(1) + h_inv(2,2)*rij(2) + h_inv(2,3)*rij(3)
          s(3) = h_inv(3,1)*rij(1) + h_inv(3,2)*rij(2) + h_inv(3,3)*rij(3)
          IF ( ABS ( s ( 1 ) ) > 0.5_dbl ) THEN
              s(1) = s(1) - perd(1) * INT(s(1)+SIGN(0.5_dbl,s(1)))
           END IF
          IF ( ABS ( s ( 2 ) ) > 0.5_dbl ) THEN
             s(2) = s(2) - perd(2) * INT(s(2)+SIGN(0.5_dbl,s(2)))
          END IF
          IF ( ABS ( s ( 3 ) ) > 0.5_dbl ) THEN
             s(3) = s(3) - perd(3) * INT(s(3)+SIGN(0.5_dbl,s(3)))
          END IF
          rij(1) = hmat(1,1)*s(1) + hmat(1,2)*s(2) + hmat(1,3)*s(3)
          rij(2) = hmat(2,1)*s(1) + hmat(2,2)*s(2) + hmat(2,3)*s(3)
          rij(3) = hmat(3,1)*s(1) + hmat(3,2)*s(2) + hmat(3,3)*s(3)
 
          rijsq = rij ( 1 ) ** 2 + rij ( 2 ) ** 2 + rij ( 3 ) ** 2  
          r2max = (radius_a + radius_b)**2

          IF (rijsq > r2max ) CYCLE

          CALL get_block_node(matrix=smat%matrix,&
                              block_row=iatom,&
                              block_col=jatom,&
                              block=s_block)

          IF (ASSOCIATED(s_block)) THEN

            atom_block(:,:,:) = 0.0_dbl

            lb_max = kind_info(jkind)%orb_basis_set%lmax(jset)
            lb_min = kind_info(jkind)%orb_basis_set%lmin(jset)

            npgfb = kind_info(jkind)%orb_basis_set%npgf(jset)
            nshellb = kind_info(jkind)%orb_basis_set%nshell(jset)
            rpgfb => kind_info(jkind)%orb_basis_set%pgf_radius(1:npgfb,jset)
            zetb => kind_info(jkind)%orb_basis_set%zet(1:npgfb,jset)

            first_cgfb = kind_info(jkind)%orb_basis_set%first_cgf(1,jset)
            last_cgfb = kind_info(jkind)%orb_basis_set%last_cgf(nshellb,jset)
            ncgfb = last_cgfb - first_cgfb + 1
            ncob = npgfb*ncoset(lb_max)

!  Calculate the overlap integrals  and the forces on the central image
            IF ( current_neighbor % nimages == 0 ) THEN

               CALL get_overlap (la_max,zeta,rpgfa,la_min,lb_max,zetb,rpgfb,lb_min, &
                        rij,rijsq,ncoa,ncob,ncgfa,ncgfb,sab,sabdx,sabdy, &
                        sabdz,work,s_block(first_cgfa:last_cgfa,first_cgfb:last_cgfb), &
                        kind_info(ikind)%orb_basis_set%cphi(1:ncoa,first_cgfa:last_cgfa), &
                        kind_info(jkind)%orb_basis_set%cphi(1:ncob,first_cgfb:last_cgfb), &
                        atom_block(first_cgfa:last_cgfa,first_cgfb:last_cgfb,1:maxdim),deriv)

            ELSE
! summing over lattice translations of neighbors
               rij = current_neighbor % p % r - ri
               s(1) = h_inv(1,1)*rij(1) + h_inv(1,2)*rij(2) + h_inv(1,3)*rij(3)
               s(2) = h_inv(2,1)*rij(1) + h_inv(2,2)*rij(2) + h_inv(2,3)*rij(3)
               s(3) = h_inv(3,1)*rij(1) + h_inv(3,2)*rij(2) + h_inv(3,3)*rij(3)

               current_image => current_neighbor % image


               DO im = 1, current_neighbor % nimages

                  CALL find_image ( s, perd, current_image % vec, hmat, rijsq, rij )

                  IF ( rijsq <= r2max ) THEN
 
                      CALL get_overlap (la_max,zeta,rpgfa,la_min,lb_max,zetb,rpgfb,lb_min, &
                           rij,rijsq,ncoa,ncob,ncgfa,ncgfb,sab,sabdx,sabdy, &
                           sabdz,work,s_block(first_cgfa:last_cgfa,first_cgfb:last_cgfb), &
                           kind_info(ikind)%orb_basis_set%cphi(1:ncoa,first_cgfa:last_cgfa), &
                           kind_info(jkind)%orb_basis_set%cphi(1:ncob,first_cgfb:last_cgfb), &
                           atom_block(first_cgfa:last_cgfa,first_cgfb:last_cgfb,1:maxdim),deriv)
                  END IF 

                  current_image => current_image % next

               END DO

            END IF

          END IF 

!         here the contribution to the forces  on atoms and coeffs given
!         from each atom_block is calculated

          CALL get_forces_on_coefs(ao,kind_info,fo_coef,energy,empparm,atom_block, &
                                      part,ipart,jpart,iset,jset)

          IF ( ( PRESENT ( fo_part ) ) .AND. forces .AND. (iatom /= jatom) )  THEN
             CALL get_forces_on_part (ao, kind_info, fo_part, empparm, &
                                       atom_block, part, ipart, jpart, iset, jset)                        
          END IF
           
        END DO

      END DO

      current_neighbor => current_neighbor % next          

    END DO
  END DO

  IF (ALLOCATED(sab)) THEN
    DEALLOCATE (sab,STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory("pol_overlap","sab")
    END IF
  END IF

  IF (ALLOCATED(sabdx)) THEN
    DEALLOCATE (sabdx,STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory("pol_overlap","sabdx")
    END IF
  END IF

  IF (ALLOCATED(sabdy)) THEN
    DEALLOCATE (sabdy,STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory("pol_overlap","sabdy")
    END IF
  END IF

  IF (ALLOCATED(sabdz)) THEN
    DEALLOCATE (sabdz,STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory("pol_overlap","sabdz")
    END IF
  END IF

  IF (ALLOCATED(work)) THEN
    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory("pol_overlap","work")
    END IF
  END IF

  IF (ALLOCATED(atom_block)) THEN
    DEALLOCATE (atom_block,STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory("pol_overlap","atom_block")
    END IF
  END IF
  
  IF (ASSOCIATED(ainfo)) THEN
    DEALLOCATE (ainfo,STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory("pol_overlap","ainfo")
    END IF
  END IF

 CALL deallocate_matrix(smat%matrix)

  END SUBROUTINE force_overlap           

!---------------------------------------------------------------------------------
  SUBROUTINE get_overlap(la_max,zeta,rpgfa,la_min,lb_max,zetb,rpgfb,lb_min, &
                              rij,rijsq,ncoa,ncob,ncgfa,ncgfb,sab,sabdx,sabdy, &
                              sabdz,work,s_block,cphia,cphib,atom_block,forces) 
!---------------------------------------------------------------------------------


    INTEGER, INTENT(IN)                      :: la_max
    REAL(dbl), DIMENSION(:), POINTER         :: zeta, rpgfa
    INTEGER, INTENT(IN)                      :: la_min, lb_max
    REAL(dbl), DIMENSION(:), POINTER         :: zetb, rpgfb
    INTEGER, INTENT(IN)                      :: lb_min
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rij
    REAL(dbl), intent(IN)                    :: rijsq
    INTEGER, INTENT(IN)                      :: ncoa, ncob, ncgfa, ncgfb
    REAL(dbl), DIMENSION(:, :), &
      intent(inout)                          :: sab, sabdx, sabdy, sabdz, &
                                                work, s_block
    REAL(dbl), DIMENSION(:, :), intent(in)   :: cphia, cphib
    REAL(dbl), DIMENSION(:, :, :), &
      intent(inout)                          :: atom_block
    LOGICAL, INTENT(IN)                      :: forces

    INTEGER                                  :: lda, ldb, ldc
    REAL(dbl)                                :: dab

!   ---------------------------------------------------------------------------

   dab = sqrt(rijsq)
   sab = 0._dbl

!  calculate the primitive overlap integral 

   CALL overlap(la_max,zeta,rpgfa,la_min,&
                        lb_max,zetb,rpgfb,lb_min,&
                         rij,rijsq,sab)

!  Contraction step (overlap matrix) 

   lda = size(sab,1)
   ldb = size(cphib,1)
   ldc = size(work,1)

   CALL dgemm("N","N",ncoa,ncgfb,ncob,1.0_dbl,sab(1,1),lda,&
                       cphib(1,1), &
                       ldb,0.0_dbl,&
                       work(1,1),ldc)
   lda = size(cphia,1)
   ldb = size(work,1)
   ldc = size(s_block,1)

   CALL dgemm("T","N",ncgfa,ncgfb,ncoa,1.0_dbl,&
                       cphia,&
                       lda,&
                       work,ldb,0.0_dbl,&
                       s_block,ldc)

! get the total overlap integral          

  atom_block(:,:,1) = atom_block(:,:,1) + s_block(:,:)

! if requested, calculate the derivatives

  if (forces) then 

   sab = 0._dbl

!  calculate the primitive overlap integral

   CALL overlap(la_max+1,zeta,rpgfa,la_min,&
                        lb_max,zetb,rpgfb,lb_min,&
                         rij,rijsq,sab)

   sabdx = 0._dbl
   sabdy = 0._dbl
   sabdz = 0._dbl

   CALL dabdr(la_max,zeta,rpgfa,la_min,&
                   lb_max,zetb,rpgfb,lb_min,&
                   dab,sab,sabdx,sabdy,sabdz)

!  contraction step for the x component

    CALL dgemm("N","N",ncoa,ncgfb,ncob,1.0_dbl,sabdx(1,1),SIZE(sabdx,1),&
                cphib, SIZE(cphib,1),0.0_dbl,work(1,1),SIZE(work,1))

    CALL dgemm("T","N",ncgfa,ncgfb,ncoa,1.0_dbl,cphia,SIZE(cphia,1),&
                work(1,1),SIZE(work,1),0.0_dbl,s_block,SIZE(s_block,1))

! get the total x derivative

    atom_block(:,:,2) = atom_block(:,:,2) + s_block(:,:)

!  contraction step for the y component

    CALL dgemm("N","N",ncoa,ncgfb,ncob,1.0_dbl,sabdy(1,1),SIZE(sabdy,1),&
                cphib, SIZE(cphib,1),0.0_dbl,work(1,1),SIZE(work,1))

    CALL dgemm("T","N",ncgfa,ncgfb,ncoa,1.0_dbl,cphia,SIZE(cphia,1),&
                work(1,1),SIZE(work,1),0.0_dbl,s_block,SIZE(s_block,1))

! get the total y derivative

    atom_block(:,:,3) = atom_block(:,:,3) + s_block(:,:)

!  contraction step for the z component

    CALL dgemm("N","N",ncoa,ncgfb,ncob,1.0_dbl,sabdz(1,1),SIZE(sabdz,1),&
                cphib, SIZE(cphib,1),0.0_dbl,work(1,1),SIZE(work,1))

    CALL dgemm("T","N",ncgfa,ncgfb,ncoa,1.0_dbl,cphia,SIZE(cphia,1),&
                work(1,1),SIZE(work,1),0.0_dbl,s_block,SIZE(s_block,1))

! get the total z derivative

    atom_block(:,:,4) = atom_block(:,:,4) + s_block(:,:)

  END IF

  END SUBROUTINE get_overlap

!--------------------------------------------------------------------------------------
  SUBROUTINE  build_block_overlap_matrix(pnode,part,ainfo,smat,natom,ncgf)
!--------------------------------------------------------------------------------------

    TYPE(particle_node_type), dimension(:), &
      intent(in)                             :: pnode
    TYPE(particle_type), dimension(:), &
      intent(in)                             :: part
    TYPE(atom_basis_info_type), &
      DIMENSION(:), INTENT(IN)               :: ainfo
    TYPE(real_matrix_p_type), intent(out)    :: smat
    INTEGER, intent(in)                      :: natom, ncgf

    INTEGER                                  :: i, iatom, ikind, ipart, j, &
                                                jatom, jpart
    TYPE(linklist_neighbor), POINTER         :: current_neighbor

  CALL allocate_matrix(matrix=smat%matrix,&
                         nblock_row=natom,&
                         nblock_col=natom,&
                         nrow=ncgf,&
                         ncol=ncgf,&
                         first_row=ainfo(:)%first_cgf,&
                         last_row=ainfo(:)%last_cgf,&
                         first_col=ainfo(:)%first_cgf,&
                         last_col=ainfo(:)%last_cgf,&
                         matrix_name="OVERLAP MATRIX",&
                         matrix_symmetry="symmetric")  

  DO i = 1, size(pnode)

    ipart = pnode(i) %p%iatom
    ikind = part (ipart) % kind

!   cycle if the particle has no polarization basis functions

    IF (ikind <= 0) CYCLE

    iatom = pnode (i) % p % abase

!   Loop over all neighbor atoms of the current atom "iatom" 
!   (the neighbor list is constructed with the black-white scheme
!   and the blocks are allocated accordingly)

    current_neighbor => pnode ( i ) % sl

    DO j = 1, pnode (i) % nsneighbor

      jpart = current_neighbor % index
      jatom = part (jpart) % abase

      CALL add_block_node ( matrix=smat%matrix, block_row=iatom, block_col=jatom ) 
        
      current_neighbor => current_neighbor % next          

    END DO  

  END DO

  END SUBROUTINE  build_block_overlap_matrix                    

!--------------------------------------------------------------------------------------
  SUBROUTINE get_forces_on_coefs(ao,kind_info,fo_coef,energy,empparm,atom_block, &
                                   part,ipart,jpart,iset,jset)
!--------------------------------------------------------------------------------------

    TYPE(ao_type), intent(inout)             :: ao
    TYPE(kind_info_type), DIMENSION(:), &
      INTENT(IN)                             :: kind_info
    REAL(dbl), DIMENSION(:), INTENT(INOUT)   :: fo_coef
    REAL(dbl), INTENT(INOUT)                 :: energy
    TYPE(empirical_parameter_type), &
      DIMENSION(:), INTENT(IN)               :: empparm
    REAL(dbl), DIMENSION(:, :, :), &
      INTENT(IN)                             :: atom_block
    TYPE(particle_type), dimension(:), &
      intent(in)                             :: part
    INTEGER, intent(in)                      :: ipart, jpart, iset, jset

    INTEGER :: first_cgfa, first_cgfb, ic, icoef, ikind, jc, jcoef, jkind, &
      last_cgfa, last_cgfb, nshella, nshellb
    REAL(dbl)                                :: hpi, hpj

  ikind = part (ipart) %  kind          
  jkind = part (jpart) %  kind             
  
  first_cgfa = kind_info (ikind) % orb_basis_set % first_cgf (1,iset)
  nshella =  kind_info (ikind) % orb_basis_set % nshell (iset)
  last_cgfa = kind_info (ikind) % orb_basis_set % last_cgf (nshella,iset)
  
  first_cgfb =  kind_info (jkind) % orb_basis_set % first_cgf (1,jset)
  nshellb =  kind_info (jkind) % orb_basis_set % nshell (jset)
  last_cgfb =  kind_info (jkind) % orb_basis_set % last_cgf (nshellb,jset)

  hpi = empparm (ikind) % hardness_param (iset) 
  hpj = empparm (jkind) % hardness_param (jset) 
  
  DO ic = first_cgfa, last_cgfa

     icoef = part (ipart) % coef_list (ic)

     DO jc = first_cgfb, last_cgfb

        jcoef = part (jpart) % coef_list (jc)

        fo_coef (icoef) = fo_coef (icoef) - 0.5_dbl * (hpi + hpj) * ao % cr (jcoef) * &
                          atom_block (ic,jc,1)

        energy = energy + 0.25_dbl * (hpi + hpj) * ao % cr (icoef) * &
                          ao % cr (jcoef) * atom_block (ic, jc, 1)

        IF (ipart/=jpart) THEN 
          fo_coef (jcoef) = fo_coef (jcoef) - 0.5_dbl * (hpi + hpj) * ao % cr (icoef) * &
                            atom_block (ic,jc,1)
          energy = energy + 0.25_dbl * (hpi + hpj) * ao % cr (icoef) * &
                            ao % cr (jcoef) * atom_block (ic, jc, 1)
        END IF
     END DO
  END DO
       
 END SUBROUTINE get_forces_on_coefs

!--------------------------------------------------------------------------------------
  SUBROUTINE get_forces_on_part(ao,kind_info,fo_part,empparm,atom_block, &
                                   part,ipart,jpart,iset,jset)
!--------------------------------------------------------------------------------------

    TYPE(ao_type), intent(inout)             :: ao
    TYPE(kind_info_type), DIMENSION(:), &
      INTENT(IN)                             :: kind_info
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: fo_part
    TYPE(empirical_parameter_type), &
      DIMENSION(:), INTENT(IN)               :: empparm
    REAL(dbl), DIMENSION(:, :, :), &
      INTENT(IN)                             :: atom_block
    TYPE(particle_type), dimension(:), &
      intent(in)                             :: part
    INTEGER, intent(in)                      :: ipart, jpart, iset, jset

    INTEGER :: first_cgfa, first_cgfb, ic, icoef, ikind, jc, jcoef, jkind, &
      last_cgfa, last_cgfb, nshella, nshellb
    REAL(dbl)                                :: hpi, hpj

  ikind = part (ipart) %  kind               
  jkind = part (jpart) %  kind                 
 
  first_cgfa =  kind_info (ikind) % orb_basis_set % first_cgf (1,iset)
  nshella =  kind_info (ikind) % orb_basis_set % nshell (iset)
  last_cgfa =  kind_info (ikind) % orb_basis_set % last_cgf (nshella,iset)
 
  first_cgfb =  kind_info (jkind) % orb_basis_set % first_cgf (1,jset)
  nshellb =  kind_info (jkind) % orb_basis_set % nshell (jset)
  last_cgfb = kind_info (jkind) % orb_basis_set % last_cgf (nshellb,jset)

  hpi = empparm (ikind) % hardness_param (iset)
  hpj = empparm (jkind) % hardness_param (jset)

  DO ic = first_cgfa, last_cgfa

     icoef = part (ipart) % coef_list (ic)

     DO jc = first_cgfb, last_cgfb

        jcoef = part (jpart) % coef_list (jc)

        fo_part (1,ipart) = fo_part (1,ipart) - 0.5_dbl *  (hpi + hpj) * ao % cr (jcoef) * &
                          ao % cr (icoef) * atom_block (ic,jc,2)
        fo_part (1,jpart) = fo_part (1,jpart) + 0.5_dbl * (hpi + hpj) * ao % cr (icoef) * &
                          ao % cr (jcoef) * atom_block (ic,jc,2)

        fo_part (2,ipart) = fo_part (2,ipart) - 0.5_dbl * (hpi + hpj) * ao % cr (jcoef) * &
                          ao % cr (icoef) * atom_block (ic,jc,3)
        fo_part (2,jpart) = fo_part (2,jpart) + 0.5_dbl * (hpi + hpj) * ao % cr (icoef) * &
                          ao % cr (jcoef) * atom_block (ic,jc,3)

        fo_part (3,ipart) = fo_part (3,ipart) - 0.5_dbl * (hpi + hpj) * ao % cr (jcoef) * &
                          ao % cr (icoef) * atom_block (ic,jc,4)
        fo_part (3,jpart) = fo_part (3,jpart) + 0.5_dbl * (hpi + hpj) * ao % cr (icoef) * &
                          ao % cr (jcoef) * atom_block (ic,jc,4)

     END DO
 
  END DO

 END SUBROUTINE get_forces_on_part

! *****************************************************************************

END MODULE pol_overlap

! *****************************************************************************
