!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/pol_overlap_pw [1.0] *
!!
!!   NAME
!!     pol_overlap_pw
!!
!!   FUNCTION
!!     Calculation of the second order contribution of the Hohenbeg-Kohn
!!     functional that depends on the overlaps between basis functions
!!
!!     E (overl) = 1/2 sum_ij  c_i*c_j (f_i + f_j) <phi_i|phi_j>
!!
!!   AUTHOR
!!     CJM (MAR 2002)     
!!
!!   MODIFICATION HISTORY
!!     NONE
!!
!!   SOURCE
!******************************************************************************
MODULE pol_overlap_pw
!******************************************************************************

  USE ai_derivatives,                  ONLY: dabdr
  USE ai_overlap,                      ONLY: overlap
  USE ao_types,                        ONLY: ao_type
  USE atomic_kinds,                    ONLY: kind_info_type
  USE coefficient_lists,               ONLY: coef_list_type
  USE coefficient_types,               ONLY: PW_COMPLEXDATA1D,&
                                             PW_COMPLEXDATA3D,&
                                             PW_REALDATA3D,&
                                             PW_REALSPACE,&
                                             PW_RECIPROCALSPACE,&
                                             coeff_allocate,&
                                             coeff_deallocate,&
                                             coeff_transform_space,&
                                             coeff_type,&
                                             coeff_zero
  USE dg_types,                        ONLY: dg_type
  USE dgs,                             ONLY: dg_get_patch,&
                                             dg_get_strucfac,&
                                             dg_sum_patch
  USE empirical_parameters,            ONLY: empirical_parameter_type
  USE fist_nonbond_force,              ONLY: find_image
  USE kinds,                           ONLY: dp
  USE method_specifications,           ONLY: maxder
  USE molecule_types,                  ONLY: linklist_images,&
                                             linklist_neighbor,&
                                             particle_node_type
  USE orbital_pointers,                ONLY: ncoset
  USE particle_types,                  ONLY: particle_type
  USE pol_projectors_pw,               ONLY: pol_projector_type
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_types,                        ONLY: pw_add,&
                                             pw_copy
  USE realspace_grid_types,            ONLY: realspace_grid_type,&
                                             rs_grid_allocate,&
                                             rs_grid_deallocate,&
                                             rs_grid_setup,&
                                             rs_pw_transfer
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell_param
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             allocate_matrix,&
                                             deallocate_matrix,&
                                             get_block_node,&
                                             real_matrix_p_type
  USE structure_factor_types,          ONLY: structure_factor_type
  USE structure_factors,               ONLY: structure_factor_allocate,&
                                             structure_factor_deallocate
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE pol_overlap,                     ONLY: maxder

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

  SUBROUTINE force_overlap ( dg, pw_grid, ao, kind_info, part, pnode, box,  &
                             empparm, grid_s, grid_b, energy, fo_coef, fo_part,  &
                             fproj_part, e_ovl_ii, e_ovl_ij )


    TYPE(dg_type), DIMENSION(:), INTENT(IN)  :: dg
    TYPE(pw_grid_type), INTENT(IN)           :: pw_grid
    TYPE(ao_type), INTENT(INOUT)             :: ao
    TYPE(kind_info_type), DIMENSION(:), &
      INTENT(IN)                             :: kind_info
    TYPE(particle_type), dimension(:), &
      INTENT(IN)                             :: part
    TYPE(particle_node_type), DIMENSION(:), &
      INTENT(INOUT)                          :: pnode
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(empirical_parameter_type), &
      DIMENSION(:), INTENT(IN)               :: empparm
    TYPE(pw_grid_type), DIMENSION(:), &
      INTENT(IN)                             :: grid_s, grid_b
    REAL(KIND=dp), INTENT(OUT)                   :: energy
    COMPLEX(KIND=dp), DIMENSION(:), &
      INTENT(INOUT), OPTIONAL                :: fo_coef
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL                :: fo_part
    TYPE(pol_projector_type), DIMENSION(:), &
      INTENT(IN), OPTIONAL                   :: fproj_part
    REAL(KIND=dp), INTENT(OUT), OPTIONAL         :: e_ovl_ii, e_ovl_ij

    INTEGER :: allocstat, first_cgfa, first_cgfb, gt, i, iatom, id, ig, igt, &
      ikind, im, ipart, iset, j, jatom, jkind, jpart, jset, la_max, la_min, &
      last_cgfa, last_cgfb, lb_max, lb_min, maxcgf, maxco, maxdim, nat, &
      natom, ncgf, ncgfa, ncgfb, ncoa, ncob, ngrids, ngt, nkind, npart, &
      npgfa, npgfb, nproj, nseta, nsetb, nshella, nshellb
    INTEGER, ALLOCATABLE, &
      DIMENSION(:, :, :), SAVE               :: n_images
    LOGICAL                                  :: deriv, first_time, forces
    REAL(KIND=dp)                                :: r2max, radius_a, radius_b, &
                                                rcut, rijsq
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)  :: sab, sabdx, sabdy, sabdz, work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: atom_block
    REAL(KIND=dp), DIMENSION(3)                  :: cell_lengths, perd, quotient, &
                                                ri, rij, s, vec
    REAL(KIND=dp), DIMENSION(3, 3)               :: h_inv, hmat
    REAL(KIND=dp), DIMENSION(:), POINTER         :: rpgfa, rpgfb, zeta, zetb
    REAL(KIND=dp), DIMENSION(:, :), POINTER      :: s_block
    TYPE(atom_basis_info_type), &
      DIMENSION(:), POINTER                  :: ainfo
    TYPE(coeff_type)                         :: coeff_dummy, rho_tot
    TYPE(linklist_images), POINTER           :: current_image
    TYPE(linklist_neighbor), POINTER         :: current_neighbor
    TYPE(real_matrix_p_type)                 :: smat
    TYPE(realspace_grid_type), ALLOCATABLE, &
      DIMENSION(:)                           :: rho_rs
    TYPE(structure_factor_type), &
      ALLOCATABLE, DIMENSION(:)              :: exp_igr

!   ---------------------------------------------------------------------------
! allocating the structure factor vector to the number of grid_types

  ngrids = SIZE ( dg )
  ALLOCATE ( exp_igr ( ngrids ), STAT=allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "force_overlap", &
   "exp_igr", ngrids )
! allocating real-space grid structures ahead of time for
! parallelization purposes
  ngt = SIZE ( grid_s )
  ALLOCATE ( rho_rs ( ngt ), STAT=allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "force_overlap", &
   "rho_rs", ngt )
  DO gt = 1, ngt
    CALL rs_grid_setup ( rho_rs ( gt ), grid_b ( gt ), -1 )
    CALL rs_grid_allocate ( rho_rs ( gt ) )
! Zero rho_rs
    rho_rs ( gt ) % r  = 0.0_dp
  ENDDO
  
  DO ig = 1, SIZE ( dg )

    nproj = dg ( ig ) % ncoefs

    gt = dg (ig ) % grid_index

    CALL structure_factor_allocate ( grid_s ( gt ) % bounds, nproj,  &
          exp_igr ( ig ), allocate_centre = .TRUE. )

! Parrallelization
    IF (  rho_rs ( gt ) % parallel .AND. rho_rs ( gt ) % group_dim ( 1 ) /= 1 )  THEN
      CALL get_center ( part, box, exp_igr ( ig ) % centre, grid_b ( gt ) % npts )
    END IF

! Fill structure factor
    CALL structure_factor_evaluate ( dg ( ig ), part, ao, box, exp_igr ( ig ), &
                                     grid_s ( gt ), grid_b ( gt ) )
  ENDDO


! allocate rho_tot
  CALL coeff_allocate ( rho_tot, pw_grid, use_data = PW_COMPLEXDATA1D ) 
  rho_tot % pw % in_space = PW_RECIPROCALSPACE
  CALL coeff_zero ( rho_tot )

! allocate coeff_dummy
  CALL coeff_allocate ( coeff_dummy, pw_grid, use_data = PW_COMPLEXDATA1D ) 
  coeff_dummy % pw % in_space = PW_RECIPROCALSPACE
  CALL coeff_zero ( coeff_dummy )


! initialize atomic info

  ncgf = 0
  natom = 0 
  nat = 0
  npart = SIZE ( pnode )
  nkind = SIZE ( kind_info )
  hmat = box % hmat 
  h_inv = box % h_inv
  perd = box % perd      
  deriv = .FALSE.

! initialize energy
  energy=0.0_dp
  e_ovl_ii=0.0_dp
  e_ovl_ij=0.0_dp

  DO ikind = 1, nkind
     nat = nat + size (kind_info(ikind) % atom_list)
  END DO
 
  NULLIFY (ainfo)

  IF (.NOT. ASSOCIATED (ainfo)) then
    ALLOCATE (ainfo(nat), stat=allocstat)
    IF (allocstat /= 0) THEN
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

  if (nat /= natom) call stop_program("force_overlap","inconsistent number of atoms with basis")

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

  ALLOCATE (sab(maxco,maxco),STAT=allocstat)
  IF (allocstat /= 0) THEN
    CALL stop_memory("force_overlap","sab",maxco*maxco)
  END IF
  sab(:,:) = 0.0_dp

  ALLOCATE (sabdx(maxco,maxco),STAT=allocstat)
  IF (allocstat /= 0) THEN
    CALL stop_memory("force_overlap","sabdx",maxco*maxco)
  END IF
  sabdx(:,:) = 0.0_dp

  ALLOCATE (sabdy(maxco,maxco),STAT=allocstat)
  IF (allocstat /= 0) THEN
    CALL stop_memory("force_overlap","sabdy",maxco*maxco)
  END IF
  sabdy(:,:) = 0.0_dp

  ALLOCATE (sabdz(maxco,maxco),STAT=allocstat)
  IF (allocstat /= 0) THEN
    CALL stop_memory("force_overlap","sabdz",maxco*maxco)
  END IF
  sabdz(:,:) = 0.0_dp

  ALLOCATE (work(maxco,maxcgf),STAT=allocstat)
  IF (allocstat /= 0) THEN
  CALL stop_memory("force_overlap","work",maxco*maxcgf)
  END IF
  work(:,:) = 0.0_dp

  ALLOCATE (atom_block(maxco,maxcgf,maxdim),STAT=allocstat)
  IF (allocstat /= 0) THEN
  CALL stop_memory("force_overlap","atom_blok",maxco*maxcgf*maxdim)
  END IF
  atom_block(:,:,:) = 0.0_dp

  NULLIFY (smat%matrix)

  CALL build_block_overlap_matrix(pnode,part,ainfo,smat,natom,ncgf)

! initialize forces

  IF ( PRESENT ( fo_part ) ) fo_part = 0.0_dp
  IF ( PRESENT ( fo_coef ) ) fo_coef = 0.0_dp

! intialize images

  first_time = .NOT.ALLOCATED ( n_images )
  IF ( first_time ) THEN
     IF ( .NOT. ALLOCATED ( n_images ) ) &
          ALLOCATE ( n_images (  nkind, nkind, 3 ), STAT = allocstat )
     IF ( allocstat /= 0 ) CALL stop_memory ( 'force_env_types', &
          'n_images', nkind ** 2 * 3 )
 
     CALL get_cell_param ( box, cell_lengths )
 
     DO ikind = 1, nkind         
        DO jkind = 1,  nkind         
           rcut = kind_info( ikind ) % orb_basis_set % kind_radius + &
                  kind_info( jkind) % orb_basis_set % kind_radius 
           quotient ( : ) = rcut / cell_lengths ( : )
           DO id=1,3
              IF ( quotient ( id ) <= 0.5_dp ) THEN
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
          IF ( ABS ( s ( 1 ) ) > 0.5_dp ) THEN
              s(1) = s(1) - perd(1) * INT(s(1)+SIGN(0.5_dp,s(1)))
           END IF
          IF ( ABS ( s ( 2 ) ) > 0.5_dp ) THEN
             s(2) = s(2) - perd(2) * INT(s(2)+SIGN(0.5_dp,s(2)))
          END IF
          IF ( ABS ( s ( 3 ) ) > 0.5_dp ) THEN
             s(3) = s(3) - perd(3) * INT(s(3)+SIGN(0.5_dp,s(3)))
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

            atom_block(:,:,:) = 0.0_dp

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

          CALL energy_evaluate ( ao, dg, rho_rs, kind_info, part, atom_block,  &
                               empparm, ipart, jpart, iset, jset, energy, grid_s, fo_coef, &
                               e_ovl_ii, e_ovl_ij ) 

          IF ( ( PRESENT ( fo_part ) ) .AND. forces )  THEN
             CALL fpart_evaluate ( ao, kind_info, fproj_part, empparm, atom_block, &
                                   part, ipart, jpart, iset, jset, fo_part )                        
          END IF
           
        END DO

      END DO

      current_neighbor => current_neighbor % next          

    END DO
  END DO

  IF (ALLOCATED(atom_block)) THEN
    DEALLOCATE (atom_block,STAT=allocstat)
    IF (allocstat /= 0) THEN
      CALL stop_memory("force_overlap","atom_block")
    END IF
  END IF

  IF (ALLOCATED(work)) THEN
    DEALLOCATE (work,STAT=allocstat)
    IF (allocstat /= 0) THEN
      CALL stop_memory("force_overlap","work")
    END IF
  END IF

  IF (ALLOCATED(sabdz)) THEN
    DEALLOCATE (sabdz,STAT=allocstat)
    IF (allocstat /= 0) THEN
      CALL stop_memory("force_overlap","sabdz")
    END IF
  END IF

  IF (ALLOCATED(sabdy)) THEN
    DEALLOCATE (sabdy,STAT=allocstat)
    IF (allocstat /= 0) THEN
      CALL stop_memory("force_overlap","sabdy")
    END IF
  END IF

  IF (ALLOCATED(sabdx)) THEN
    DEALLOCATE (sabdx,STAT=allocstat)
    IF (allocstat /= 0) THEN
      CALL stop_memory("force_overlap","sabdx")
    END IF
  END IF

  IF (ALLOCATED(sab)) THEN
    DEALLOCATE (sab,STAT=allocstat)
    IF (allocstat /= 0) THEN
      CALL stop_memory("force_overlap","sab")
    END IF
  END IF
  
  IF (ASSOCIATED(ainfo)) THEN
    DEALLOCATE (ainfo,STAT=allocstat)
    IF (allocstat /= 0) THEN
      CALL stop_memory("force_overlap","ainfo")
    END IF
  END IF

  CALL deallocate_matrix(smat%matrix)

  IF ( PRESENT ( fo_coef ) ) THEN
! Obtain real-space force on coefficients
    CALL rs_force_evaluate ( grid_b, rho_rs, coeff_dummy, rho_tot )

! Obtain g-space force on the coefficients
    rho_tot % pw % cc = CONJG ( rho_tot % pw % cc )
    fo_coef =  CMPLX ( -2.0_dp * REAL ( rho_tot % pw % cc,KIND=dp),  &
                       2.0_dp * AIMAG ( rho_tot % pw % cc ),KIND=dp)
! renormalize to 0 the force on g=0
    fo_coef (1) =  CMPLX ( 0.0_dp, 0.0_dp,KIND=dp)
  ENDIF  

  CALL coeff_deallocate ( coeff_dummy ) 

  CALL coeff_deallocate ( rho_tot ) 

  DO igt = 1, ngt
    CALL rs_grid_deallocate ( rho_rs ( igt ) )
  ENDDO

  DEALLOCATE ( rho_rs,  STAT=allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "force_overlap", &
   "rho_rs", ngt )

  DO ig = 1, SIZE ( dg )
    CALL structure_factor_deallocate ( exp_igr ( ig ) )
  ENDDO

  DEALLOCATE ( exp_igr,  STAT=allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "force_overlap", &
   "exp_igr", ngrids )

  END SUBROUTINE force_overlap           

!---------------------------------------------------------------------------------
  SUBROUTINE get_overlap(la_max,zeta,rpgfa,la_min,lb_max,zetb,rpgfb,lb_min, &
                              rij,rijsq,ncoa,ncob,ncgfa,ncgfb,sab,sabdx,sabdy, &
                              sabdz,work,s_block,cphia,cphib,atom_block,forces) 
!---------------------------------------------------------------------------------


    INTEGER, INTENT(IN)                      :: la_max
    REAL(KIND=dp), DIMENSION(:), POINTER         :: zeta, rpgfa
    INTEGER, INTENT(IN)                      :: la_min, lb_max
    REAL(KIND=dp), DIMENSION(:), POINTER         :: zetb, rpgfb
    INTEGER, INTENT(IN)                      :: lb_min
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)      :: rij
    REAL(KIND=dp), intent(IN)                    :: rijsq
    INTEGER, INTENT(IN)                      :: ncoa, ncob, ncgfa, ncgfb
    REAL(KIND=dp), DIMENSION(:, :), &
      intent(inout)                          :: sab, sabdx, sabdy, sabdz, &
                                                work, s_block
    REAL(KIND=dp), DIMENSION(:, :), intent(in)   :: cphia, cphib
    REAL(KIND=dp), DIMENSION(:, :, :), &
      intent(inout)                          :: atom_block
    LOGICAL, INTENT(IN)                      :: forces

    INTEGER                                  :: lda, ldb, ldc
    REAL(KIND=dp)                                :: dab

!   ---------------------------------------------------------------------------

   dab = sqrt(rijsq)
   sab = 0.0_dp

!  calculate the primitive overlap integral 

   CALL overlap(la_max,zeta,rpgfa,la_min,&
                        lb_max,zetb,rpgfb,lb_min,&
                         rij,rijsq,sab)

!  Contraction step (overlap matrix) 

   lda = size(sab,1)
   ldb = size(cphib,1)
   ldc = size(work,1)

   CALL dgemm("N","N",ncoa,ncgfb,ncob,1.0_dp,sab(1,1),lda,&
                       cphib(1,1), &
                       ldb,0.0_dp,&
                       work(1,1),ldc)
   lda = size(cphia,1)
   ldb = size(work,1)
   ldc = size(s_block,1)

   CALL dgemm("T","N",ncgfa,ncgfb,ncoa,1.0_dp,&
                       cphia,&
                       lda,&
                       work,ldb,0.0_dp,&
                       s_block,ldc)

! get the total overlap integral          

  atom_block(:,:,1) = atom_block(:,:,1) + s_block(:,:)

! if requested, calculate the derivatives

  if (forces) then 

   sab = 0.0_dp

!  calculate the primitive overlap integral

   CALL overlap(la_max+1,zeta,rpgfa,la_min,&
                        lb_max,zetb,rpgfb,lb_min,&
                         rij,rijsq,sab)

   sabdx = 0.0_dp
   sabdy = 0.0_dp
   sabdz = 0.0_dp

   CALL dabdr(la_max,zeta,rpgfa,la_min,&
                   lb_max,zetb,rpgfb,lb_min,&
                   dab,sab,sabdx,sabdy,sabdz)

!  contraction step for the x component

    CALL dgemm("N","N",ncoa,ncgfb,ncob,1.0_dp,sabdx(1,1),SIZE(sabdx,1),&
                cphib, SIZE(cphib,1),0.0_dp,work(1,1),SIZE(work,1))

    CALL dgemm("T","N",ncgfa,ncgfb,ncoa,1.0_dp,cphia,SIZE(cphia,1),&
                work(1,1),SIZE(work,1),0.0_dp,s_block,SIZE(s_block,1))

! get the total x derivative

    atom_block(:,:,2) = atom_block(:,:,2) + s_block(:,:)

!  contraction step for the y component

    CALL dgemm("N","N",ncoa,ncgfb,ncob,1.0_dp,sabdy(1,1),SIZE(sabdy,1),&
                cphib, SIZE(cphib,1),0.0_dp,work(1,1),SIZE(work,1))

    CALL dgemm("T","N",ncgfa,ncgfb,ncoa,1.0_dp,cphia,SIZE(cphia,1),&
                work(1,1),SIZE(work,1),0.0_dp,s_block,SIZE(s_block,1))

! get the total y derivative

    atom_block(:,:,3) = atom_block(:,:,3) + s_block(:,:)

!  contraction step for the z component

    CALL dgemm("N","N",ncoa,ncgfb,ncob,1.0_dp,sabdz(1,1),SIZE(sabdz,1),&
                cphib, SIZE(cphib,1),0.0_dp,work(1,1),SIZE(work,1))

    CALL dgemm("T","N",ncgfa,ncgfb,ncoa,1.0_dp,cphia,SIZE(cphia,1),&
                work(1,1),SIZE(work,1),0.0_dp,s_block,SIZE(s_block,1))

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
                         sparsity_id=-1, & ! XXXXxxx unknown sparsity id
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
  SUBROUTINE energy_evaluate ( ao, dg, rho_rs, kind, part, atom_block, empparm, ipart, &
                              jpart, iset, jset, energy, grid_s, fcoef, e_ovl_ii, e_ovl_ij )
!--------------------------------------------------------------------------------------

    TYPE(ao_type), INTENT(INOUT)             :: ao
    TYPE(dg_type), DIMENSION(:), INTENT(IN)  :: dg
    TYPE(realspace_grid_type), &
      DIMENSION(:), INTENT(INOUT)            :: rho_rs
    TYPE(kind_info_type), DIMENSION(:), &
      INTENT(IN)                             :: kind
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: atom_block
    TYPE(empirical_parameter_type), &
      DIMENSION(:), INTENT(IN)               :: empparm
    INTEGER, INTENT(IN)                      :: ipart, jpart, iset, jset
    REAL(KIND=dp), INTENT(INOUT)                 :: energy
    TYPE(pw_grid_type), DIMENSION(:), &
      INTENT(IN)                             :: grid_s
    COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: fcoef
    REAL(KIND=dp), INTENT(INOUT)                 :: e_ovl_ii, e_ovl_ij

    INTEGER :: first_cgfa, first_cgfb, ic, ikind, iproj, jc, jkind, jproj, &
      last_cgfa, last_cgfb, nshella, nshellb
    REAL(KIND=dp)                                :: eij, hpi, hpj

  ikind = part ( ipart ) %  kind          
  jkind = part ( jpart ) %  kind             
  
  first_cgfa = kind ( ikind ) % orb_basis_set % first_cgf ( 1,iset )
  nshella =  kind ( ikind ) % orb_basis_set % nshell ( iset )
  last_cgfa = kind ( ikind ) % orb_basis_set % last_cgf ( nshella, iset )
  
  first_cgfb =  kind ( jkind ) % orb_basis_set % first_cgf ( 1, jset )
  nshellb =  kind ( jkind ) % orb_basis_set % nshell ( jset )
  last_cgfb =  kind ( jkind ) % orb_basis_set % last_cgf ( nshellb, jset )

  hpi = empparm ( ikind ) % hardness_param ( iset ) 
  hpj = empparm ( jkind ) % hardness_param ( jset ) 
  
  DO ic = first_cgfa, last_cgfa

     iproj = part ( ipart ) % coef_list ( ic )

     DO jc = first_cgfb, last_cgfb

        eij = 0.0_dp

        jproj = part ( jpart ) % coef_list ( jc )

        eij =  0.25_dp * ( hpi + hpj ) * atom_block ( ic, jc, 1 ) * &
               ao % norm ( iproj ) * ao % norm ( jproj )

        IF ( ipart /= jpart ) THEN 

          eij = eij + 0.25_dp * ( hpi + hpj ) * atom_block ( ic, jc, 1 ) * &
                ao % norm ( iproj ) * ao % norm ( jproj )
                                      
        END IF

        energy = energy + eij * ao % cr ( iproj ) * ao % cr ( jproj )

        IF ( ipart /= jpart ) THEN 
           e_ovl_ij = e_ovl_ij + eij * ao % cr ( iproj ) * ao % cr ( jproj )
        ELSE
           e_ovl_ii = e_ovl_ii + eij * ao % cr ( iproj ) * ao % cr ( jproj )
        END IF

        IF ( PRESENT ( fcoef ) ) CALL get_patch ( eij, dg, ao, grid_s, iproj,  &
                                                 jproj, rho_rs )
     END DO
  END DO
       
 END SUBROUTINE energy_evaluate
!******************************************************************************
!!****** pol_overlap_pw/get_patch[1.0] *
!!
!!   NAME
!!     get_patch
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM 
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE get_patch ( eij, dg, ao, grid_s, iproj, jproj, rho_rs )
    REAL(KIND=dp), INTENT(IN)                    :: eij
    TYPE(dg_type), DIMENSION(:), INTENT(IN)  :: dg
    TYPE(ao_type), INTENT(INOUT)             :: ao
    TYPE(pw_grid_type), DIMENSION(:), &
      INTENT(IN)                             :: grid_s
    INTEGER, INTENT(IN)                      :: iproj, jproj
    TYPE(realspace_grid_type), &
      DIMENSION(:), INTENT(INOUT)            :: rho_rs

    INTEGER                                  :: gt, gt1, gt2, ig, ig1, ig2
    TYPE(coeff_type)                         :: rhos1, rhos2

      ig = 0
      gt = 0

      ig1 = ao % grid_type_id ( iproj )
      gt1 = ao % grid_id ( iproj )

      ig2 = ao % grid_type_id ( jproj )
      gt2 = ao % grid_id ( jproj )

      IF ( gt1 == gt2 ) gt = gt1
      IF ( ig1 == ig2 ) ig = ig1

      IF ( ig /= 0  .AND. gt /= 0 ) THEN ! both belong to same grid_type and grid_index
! Allocate small coefficient types
        CALL coeff_allocate ( rhos1,  grid_s ( gt ), use_data = PW_REALDATA3D )
        CALL coeff_allocate ( rhos2,  grid_s ( gt ), use_data = PW_REALDATA3D )
! calculate function on small boxes (we use double packing in FFT)
        CALL sum_patch_2x ( dg ( ig ), rho_rs ( gt ), ao, iproj, jproj, &
                            eij, rhos1, rhos2)
! do a little cleaning
        CALL coeff_deallocate ( rhos2 )
        CALL coeff_deallocate ( rhos1 )

      ELSE  ! different grid_type and grid_index
! Allocate small coefficient types
        CALL coeff_allocate ( rhos1,  grid_s ( gt1 ), use_data = PW_REALDATA3D )
        CALL coeff_allocate ( rhos2,  grid_s ( gt2 ), use_data = PW_REALDATA3D )

        CALL sum_patch_1x ( dg ( ig1 ), rho_rs ( gt1 ), ao, iproj, jproj, eij, rhos1 )
        CALL sum_patch_1x ( dg ( ig2 ), rho_rs ( gt2 ), ao, jproj, iproj, eij, rhos2 )
          
! do a little cleaning
        CALL coeff_deallocate ( rhos2 )
        CALL coeff_deallocate ( rhos1 )
      ENDIF
  END SUBROUTINE get_patch
!******************************************************************************
!!****** pol_overlap_pw/rs_force_evaluate [1.0] *
!!
!!   NAME
!!     rs_force_evaluate
!!
!!   FUNCTION
!!     Evaluate real space force
!!   AUTHOR
!!     CJM 
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

   SUBROUTINE rs_force_evaluate ( grid_b, rho_rs, coeff_dummy, rho_tot )
    TYPE(pw_grid_type), DIMENSION(:), &
      INTENT(IN)                             :: grid_b
    TYPE(realspace_grid_type), &
      DIMENSION(:), INTENT(INOUT)            :: rho_rs
    TYPE(coeff_type), INTENT(INOUT)          :: coeff_dummy, rho_tot

    INTEGER                                  :: igt
    TYPE(coeff_type)                         :: rhob_g, rhob_r

! Sum to reference grid 

      DO igt = 1, SIZE ( grid_b )

        IF ( ALL ( rho_rs ( igt ) % r == 0.0_dp ) ) CYCLE

! Allocate large coefficient structure based on grid_b  
        CALL coeff_allocate ( rhob_r, grid_b ( igt ), use_data = PW_COMPLEXDATA3D )
        CALL coeff_allocate ( rhob_g, grid_b ( igt ), use_data = PW_COMPLEXDATA1D )
        rhob_r % pw % in_space = PW_REALSPACE

! Transfer from real-space -> coefficient_type
        CALL rs_pw_transfer ( rho_rs ( igt ), rhob_r % pw, "FORWARD" )
 
! DO FFT from rhob_r => rhob_g 
        CALL coeff_transform_space ( rhob_r, rhob_g )

! Copy present rho_tot to coeff_dummy
        CALL pw_copy ( rho_tot % pw, coeff_dummy % pw )

! Sum to the reference grid
        CALL pw_add ( rhob_g % pw, coeff_dummy % pw, rho_tot % pw )
! do a little cleaning
        CALL coeff_deallocate ( rhob_g )
        CALL coeff_deallocate ( rhob_r )
      ENDDO
                                   
  END SUBROUTINE rs_force_evaluate
!******************************************************************************
!!****** pol_overlap_pw/fpart_evaluate [1.0] *
!!
!!   NAME
!!     fpart_evaluate
!!
!!   FUNCTION
!!     Evaluate forces on particles
!!
!!   AUTHOR
!!     CJM 
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
!--------------------------------------------------------------------------------------
SUBROUTINE fpart_evaluate ( ao, kind_info, fproj_part, empparm, atom_block, &
                              part, ipart, jpart, iset, jset, fo_part )
!--------------------------------------------------------------------------------------

    TYPE(ao_type), intent(inout)             :: ao
    TYPE(kind_info_type), DIMENSION(:), &
      INTENT(IN)                             :: kind_info
    TYPE(pol_projector_type), DIMENSION(:), &
      INTENT(IN), TARGET                     :: fproj_part
    TYPE(empirical_parameter_type), &
      DIMENSION(:), INTENT(IN)               :: empparm
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: atom_block
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    INTEGER, intent(in)                      :: ipart, jpart, iset, jset
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: fo_part

    INTEGER :: first_cgfa, first_cgfb, ic, ikind, iproj, jc, jkind, jproj, &
      last_cgfa, last_cgfb, nshella, nshellb
    REAL(KIND=dp)                                :: hpi, hpj
    REAL(KIND=dp), DIMENSION(:), POINTER         :: fpart_i, fpart_j
    REAL(KIND=dp), POINTER                       :: norm_i, norm_j, proj_i, proj_j

  ikind = part (ipart) %  kind               
  jkind = part (jpart) %  kind                 
 
  first_cgfa =  kind_info (ikind) % orb_basis_set % first_cgf (1,iset)
  nshella =  kind_info (ikind) % orb_basis_set % nshell (iset)
  last_cgfa =  kind_info (ikind) % orb_basis_set % last_cgf (nshella,iset)
 
  first_cgfb =  kind_info (jkind) % orb_basis_set % first_cgf (1,jset)
  nshellb =  kind_info (jkind) % orb_basis_set % nshell (jset)
  last_cgfb = kind_info (jkind) % orb_basis_set % last_cgf (nshellb,jset)

  hpi = empparm ( ikind ) % hardness_param ( iset )
  hpj = empparm ( jkind ) % hardness_param ( jset )

  DO ic = first_cgfa, last_cgfa

     iproj = part ( ipart ) % coef_list ( ic )

     DO jc = first_cgfb, last_cgfb

        jproj = part ( jpart ) % coef_list ( jc )

        fpart_i => fproj_part ( iproj ) % force
        fpart_j => fproj_part ( jproj ) % force

        proj_i => ao % cr ( iproj )
        proj_j => ao % cr ( jproj )

        norm_i => ao % norm ( iproj )
        norm_j => ao % norm ( jproj )

        IF ( ipart /= jpart ) THEN
          fo_part ( 1, ipart ) = fo_part ( 1, ipart ) - 0.5_dp * ( hpi + hpj ) * &
                  (norm_i * norm_j * proj_i * proj_j * atom_block ( ic, jc, 2 ) - &
                   norm_i * norm_j * fpart_i ( 1 ) * proj_j * atom_block ( ic, jc, 1 ) )

          fo_part ( 1, jpart ) = fo_part ( 1, jpart ) + 0.5_dp * ( hpi + hpj ) * &
                  (norm_i * norm_j * proj_i * proj_j * atom_block ( ic, jc, 2 ) + &
                   norm_i * norm_j * fpart_j ( 1 ) * proj_i * atom_block ( ic, jc, 1 ) )

          fo_part ( 2, ipart ) = fo_part ( 2, ipart ) - 0.5_dp * ( hpi + hpj ) * &
                  (norm_i * norm_j * proj_i * proj_j * atom_block ( ic, jc, 3 ) - &
                   norm_i * norm_j *fpart_i ( 2 ) * proj_j * atom_block ( ic, jc, 1 ) )

          fo_part ( 2, jpart ) = fo_part ( 2, jpart ) + 0.5_dp * ( hpi + hpj ) * &
                  (norm_i * norm_j * proj_i * proj_j * atom_block ( ic, jc, 3 ) + &
                   norm_i * norm_j * fpart_j ( 2 ) * proj_i * atom_block ( ic, jc, 1 ) )

          fo_part ( 3, ipart ) = fo_part ( 3, ipart ) - 0.5_dp * ( hpi + hpj ) * &
                  (norm_i * norm_j * proj_i * proj_j * atom_block ( ic, jc, 4 ) - &
                   norm_i * norm_j * fpart_i ( 3 ) * proj_j * atom_block ( ic, jc, 1 ) )

          fo_part ( 3, jpart ) = fo_part ( 3, jpart ) + 0.5_dp * ( hpi + hpj ) * &
                  (norm_i * norm_j * proj_i * proj_j * atom_block ( ic, jc, 4 ) + &
                   norm_i * norm_j * fpart_j ( 3 ) * proj_i * atom_block ( ic, jc, 1 ) )
        ELSEIF ( ipart == jpart ) THEN
          fo_part ( 1, ipart ) = fo_part ( 1, ipart ) + 0.25_dp * ( hpi + hpj ) * &
                   norm_i * norm_j * fpart_i ( 1 ) * proj_j  * atom_block ( ic, jc, 1 ) 

          fo_part ( 1, jpart ) = fo_part ( 1, jpart ) + 0.25_dp * ( hpi + hpj ) * &
                   norm_i * norm_j * fpart_j ( 1 ) * proj_i  * atom_block ( ic, jc, 1 )

          fo_part ( 2, ipart ) = fo_part ( 2, ipart ) + 0.25_dp * ( hpi + hpj ) * &
                   norm_i * norm_j * fpart_i ( 2 ) * proj_j * atom_block ( ic, jc, 1 )

          fo_part ( 2, jpart ) = fo_part ( 2, jpart ) + 0.25_dp * ( hpi + hpj ) * &
                    norm_i * norm_j *fpart_j ( 2 ) * proj_i * atom_block ( ic, jc, 1 )

          fo_part ( 3, ipart ) = fo_part ( 3, ipart ) + 0.25_dp * ( hpi + hpj ) * &
                    norm_i * norm_j * fpart_i ( 3 ) * proj_j * atom_block ( ic, jc, 1 )

          fo_part ( 3, jpart ) = fo_part ( 3, jpart ) + 0.25_dp * ( hpi + hpj ) * &
                    norm_i * norm_j * fpart_j ( 3 ) * proj_i * atom_block ( ic, jc, 1 )
        ENDIF 
     END DO
 
  END DO

END SUBROUTINE fpart_evaluate

!******************************************************************************
!!****** pol_overlap_pw/structure_factor_evaluate [1.0] *
!!
!!   NAME
!!     structure_factor_evaluate
!!
!!   FUNCTION
!!     evaluate stucture factor and map ao % eigr => exp_igr
!!
!!   AUTHOR
!!     CJM 
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE structure_factor_evaluate ( dg, part, ao, box, exp_igr, grid_s, grid_b )


    TYPE(dg_type), INTENT(IN), TARGET        :: dg
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    TYPE(ao_type), INTENT(INOUT)             :: ao
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(structure_factor_type), INTENT(OUT) :: exp_igr
    TYPE(pw_grid_type), INTENT(IN)           :: grid_s, grid_b

    INTEGER                                  :: c1, c2, handle, ip1, ip2, p1, &
                                                p2
    REAL(KIND=dp)                                :: r1( 3 ), r2( 3 )
    TYPE(coef_list_type), POINTER            :: clist

!------------------------------------------------------------------------------

  CALL timeset ( 'STRUCTURE_FACTOR_EVALUATE', 'I', ' ', handle )

! Initialize the index of p1 and p2. This is to ensure the
! correct mapping between the particle and its array location
! in exp_igr
  ip1 = 0
  ip2 = 0

  clist => dg % clist
  DO 
     CALL set_list_drho ( clist, p1, p2, ip1, ip2, c1, c2 )
     IF ( p1 == 0 .AND. p2 == 0 ) EXIT
     r1 = part ( p1 ) % r
     r2 = 0.0_dp
     IF ( p2 /= 0 ) THEN
       r2 = part ( p2 ) % r
     ENDIF
! calculate function on small boxes (we use double packing in FFT)
     CALL get_struc_fac ( r1, r2, exp_igr, box, p2, c1, c2, ip1, ip2, &
                          grid_b, grid_s, ao )

     ip1 = ip2
  END DO

END SUBROUTINE structure_factor_evaluate

!******************************************************************************
!!****** pol_overlap_pw/set_list_drho [1.0] *
!!
!!   NAME
!!     set_list_drho
!!
!!   FUNCTION
!!     Radvances in the link list of particles to the next to be calculated
!!
!!   AUTHOR
!!     CJM (18-April-2001)
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE set_list_drho ( clist, p1, p2, ip1, ip2, c1, c2 )


    TYPE(coef_list_type), POINTER            :: clist
    INTEGER, INTENT(OUT)                     :: p1, p2
    INTEGER, INTENT(INOUT)                   :: ip1, ip2
    INTEGER, INTENT(OUT)                     :: c1, c2

!------------------------------------------------------------------------------

  p1 = 0
  p2 = 0
  c1 = 0
  c2 = 0

  IF ( ASSOCIATED ( clist ) ) THEN
    ip1 = ip1 + 1
    ip2 = ip1
    p1 = clist % particle_index
    c1 = clist % coef_index
    IF ( ASSOCIATED ( clist % next ) ) THEN
      ip2 = ip2 + 1
      clist => clist % next 
      p2 = clist % particle_index
      c2 = clist % coef_index
    END IF
    clist => clist % next 
  END IF  
 

END SUBROUTINE set_list_drho

!******************************************************************************
!!****** pol_overlap_pw/get_struc_fac [1.0] *
!!
!!   NAME
!!     get_struc_fac
!!
!!   FUNCTION
!!     Calculates strucure factor and maps to ao % eigr
!!
!!   AUTHOR
!!     CJM (Mar-2002)
!!
!!   MODIFICATION HISTORY
!!     NONE
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE get_struc_fac ( r1, r2, exp_igr, box, p2, c1, c2, &
           ip1, ip2, grid_b, grid_s, ao )


    REAL(KIND=dp), DIMENSION(3), INTENT(IN)      :: r1, r2
    TYPE(structure_factor_type)              :: exp_igr
    TYPE(cell_type), INTENT(IN)              :: box
    INTEGER, INTENT(IN)                      :: p2, c1, c2, ip1, ip2
    TYPE(pw_grid_type), INTENT(IN)           :: grid_b, grid_s
    TYPE(ao_type), INTENT(INOUT)             :: ao

    COMPLEX(KIND=dp), DIMENSION(:), POINTER      :: ex1, ex2, ey1, ey2, ez1, ez2
    INTEGER, DIMENSION(:), POINTER           :: center1, center2

!------------------------------------------------------------------------------

  center1 => exp_igr % centre ( :, ip1 )
  ex1 => exp_igr % ex ( :, ip1 )
  ey1 => exp_igr % ey ( :, ip1 )
  ez1 => exp_igr % ez ( :, ip1 )
  ao % eigr ( c1 ) % centre => exp_igr % centre ( :, ip1 )
  ao % eigr ( c1 ) % ex => exp_igr % ex ( :, ip1 )
  ao % eigr ( c1 ) % ey => exp_igr % ey ( :, ip1 )
  ao % eigr ( c1 ) % ez => exp_igr % ez ( :, ip1 )
  CALL dg_get_strucfac ( box, r1, grid_s % npts, grid_b % npts, center1, &
       exp_igr % lb, ex1, ey1, ez1 )
  

  IF ( p2 /= 0 ) THEN
     center2 => exp_igr % centre ( :, ip2 )
     ex2 => exp_igr % ex ( :, ip2 )
     ey2 => exp_igr % ey ( :, ip2 )
     ez2 => exp_igr % ez ( :, ip2 )
     ao % eigr ( c2 ) % centre => exp_igr % centre ( :, ip2 )
     ao % eigr ( c2 ) % ex => exp_igr % ex ( :, ip2 )
     ao % eigr ( c2 ) % ey => exp_igr % ey ( :, ip2 )
     ao % eigr ( c2 ) % ez => exp_igr % ez ( :, ip2 )
     CALL dg_get_strucfac ( box, r2, grid_s % npts, grid_b % npts, center2, &
          exp_igr % lb, ex2, ey2, ez2 )
  END IF

END SUBROUTINE get_struc_fac

!******************************************************************************

SUBROUTINE get_center ( part, box, centre, npts )


    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    TYPE(cell_type), INTENT(IN)              :: box
    INTEGER, DIMENSION(:, :), INTENT(OUT)    :: centre
    INTEGER, DIMENSION(:), INTENT(IN)        :: npts

    INTEGER                                  :: ipart
    REAL(KIND=dp), DIMENSION(3)                  :: s

!------------------------------------------------------------------------------

  DO ipart = 1, SIZE ( part )
    ! compute the scaled coordinate of atomi
    s = MATMUL ( box % h_inv, part ( ipart ) % r )
    s = s - NINT ( s )
    ! find the closest grid point (on big grid)
    centre ( :, ipart ) = NINT ( REAL ( npts ( 1:3 ),KIND=dp) * s ( 1:3 ) )
  END DO

END SUBROUTINE get_center

!******************************************************************************
!!****** pol_overlap_pw/get_patch_2x [1.0] *
!!
!!   NAME
!!     get_patch_2x
!!
!!   FUNCTION
!!     Calculates local density in a small box
!!
!!   AUTHOR
!!     JGH (21-Mar-2001)
!!
!!   MODIFICATION HISTORY
!!     JGH (23-Mar-2001) : Switch to integer from particle list pointers
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE sum_patch_2x ( dg, rho_rs, ao, iproj, jproj, eij, rhos1, rhos2 )


    TYPE(dg_type), INTENT(IN), TARGET        :: dg
    TYPE(realspace_grid_type), INTENT(INOUT) :: rho_rs
    TYPE(ao_type), INTENT(IN)                :: ao
    INTEGER, INTENT(IN)                      :: iproj, jproj
    REAL(KIND=dp), INTENT(IN)                    :: eij
    TYPE(coeff_type), INTENT(INOUT)          :: rhos1, rhos2

    COMPLEX(KIND=dp), DIMENSION(:), POINTER      :: ex1, ex2, ey1, ey2, ez1, ez2
    INTEGER, DIMENSION(:), POINTER           :: center1, center2
    REAL(KIND=dp)                                :: q1, q2
    TYPE(coeff_type), POINTER                :: rho0

!------------------------------------------------------------------------------

  rho0 => dg % dg_rho0 % density

  center1 => ao % eigr ( iproj ) % centre ( : )
  ex1 => ao % eigr ( iproj ) % ex ( : )
  ey1 => ao % eigr ( iproj ) % ey ( : )
  ez1 => ao % eigr ( iproj ) % ez ( : )

  center2 => ao % eigr ( jproj ) % centre ( : )
  ex2 => ao % eigr ( jproj ) % ex ( : )
  ey2 => ao % eigr ( jproj ) % ey ( : )
  ez2 => ao % eigr ( jproj ) % ez ( : )

  q1 = ao % cr ( jproj ) * ao % norm ( iproj ) * eij
  q2 = ao % cr ( iproj ) * ao % norm ( jproj ) * eij

  CALL dg_get_patch ( rho0, rhos1, rhos2, q1, q2, &
                     ex1, ey1, ez1, ex2, ey2, ez2 )

  CALL dg_sum_patch ( rho_rs, rhos1, center1 )
  CALL dg_sum_patch ( rho_rs, rhos2, center2 )
 

END SUBROUTINE sum_patch_2x

!******************************************************************************
!!****** pol_overlap_pw/sum_patch_1x [1.0] *
!!
!!   NAME
!!     get_patch_1x
!!
!!   FUNCTION
!!     Calculates local density in a small box
!!
!!   AUTHOR
!!     JGH (21-Mar-2001)
!!
!!   MODIFICATION HISTORY
!!     JGH (23-Mar-2001) : Switch to integer from particle list pointers
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE sum_patch_1x ( dg, rho_rs, ao, iproj, jproj, eij, rhos1 )


    TYPE(dg_type), INTENT(IN), TARGET        :: dg
    TYPE(realspace_grid_type), INTENT(INOUT) :: rho_rs
    TYPE(ao_type), INTENT(IN)                :: ao
    INTEGER, INTENT(IN)                      :: iproj, jproj
    REAL(KIND=dp), INTENT(IN)                    :: eij
    TYPE(coeff_type), INTENT(INOUT)          :: rhos1

    COMPLEX(KIND=dp), DIMENSION(:), POINTER      :: ex1, ey1, ez1
    INTEGER, DIMENSION(:), POINTER           :: center1
    REAL(KIND=dp)                                :: q1
    TYPE(coeff_type), POINTER                :: rho0

!------------------------------------------------------------------------------

  rho0 => dg % dg_rho0 % density

  center1 => ao % eigr ( iproj ) % centre ( : )
  ex1 => ao % eigr ( iproj ) % ex ( : )
  ey1 => ao % eigr ( iproj ) % ey ( : )
  ez1 => ao % eigr ( iproj ) % ez ( : )

  q1 = ao % cr ( jproj ) * ao % norm ( iproj ) * eij

  CALL dg_get_patch ( rho0, rhos1, q1, ex1, ey1, ez1 )

  CALL dg_sum_patch ( rho_rs, rhos1, center1 )
 

END SUBROUTINE sum_patch_1x

!*****************************************************************************
END MODULE pol_overlap_pw
!*****************************************************************************
