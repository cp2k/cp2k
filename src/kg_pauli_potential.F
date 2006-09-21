!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C ) 2006 CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/kg_pauli_potential [1.0] *
!!
!!   NAME
!!     kg_pauli_potential
!!
!!   FUNCTION
!!     Calculate energy and forces for a Pauli repulsion potential
!!     approximated by overlap terms 
!!
!!   AUTHOR
!!     JGH (01.08.2006) 
!!
!!   MODIFICATION HISTORY
!!     none
!!
!***
!******************************************************************************

MODULE kg_pauli_potential
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type
  USE kg_gpw_fm_mol_types,             ONLY: fm_mol_blocks_type,&
                                             get_fm_mol_block,&
                                             get_kg_fm_mol_set,&
                                             kg_fm_mol_set_type,&
                                             kg_fm_p_type,&
                                             mol_mo_set_p_type
  USE kg_gpw_fm_mol_utils,             ONLY: copy_sparse2mol_block
  USE kg_gpw_pw_env_types,             ONLY: get_molbox_env,&
                                             kg_molbox_env_type,&
                                             kg_sub_pw_env_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_sum
  USE molecule_kind_types,             ONLY: atom_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             molecule_type
  USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                             ncoset
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_overlap,                      ONLY: atom_overlap_matrix,&
                                             block_p_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE sparse_matrix_types,             ONLY: get_block_node,&
                                             real_matrix_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kg_pauli_potential'

  PRIVATE
  PUBLIC :: kg_gpw_pauli

!******************************************************************************

CONTAINS

!******************************************************************************

SUBROUTINE kg_gpw_pauli( kg_env, do_potential, do_forces, do_overlap, error )

    TYPE(kg_environment_type), POINTER       :: kg_env
    LOGICAL, INTENT(IN)                      :: do_potential, do_forces, &
                                                do_overlap
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_gpw_pauli', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, first_atom, handle, i, ia, iao, ifa, ikind, ila, imo, &
      imol, imol_kind, imol_local, ispin, istat, j, j_pe, ja, jao, jfa, &
      jkind, jla, jmol, jmol_kind, jmol_local, last_atom, ldai, maxco, &
      maxder, maxlgto, maxsgf, nao, nao_max, naoa, naob, nat_mol, natom, &
      natom_mol, nder, nia, njb, nmol, nmol_kind, nmol_local, nseta, nsetb, &
      nspin, output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, DIMENSION(2)                    :: nmoa, nmob
    INTEGER, DIMENSION(:), POINTER :: ifirst_ao, ilast_ao, jfirst_ao, &
      jfirst_ao_bcast, jlast_ao, jlast_ao_bcast, la_max, la_min, lb_max, &
      lb_min, npgfa, npgfb, nsgfa, nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure, ionode
    REAL(KIND=dp)                            :: epsmol, oij, pauli_energy, &
                                                rij, rk
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: fora, forb
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: dvpot, pmol, sab, smol, vpot, &
                                                work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: ai_work, smat
    REAL(KIND=dp), DIMENSION(3)              :: r1, ra, rc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: r0, set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: h_block, refp, rpgfa, rpgfb, &
                                                sphi_a, sphi_b, zeta, zetb
    TYPE(atom_type), DIMENSION(:), POINTER   :: iatom_list, jatom_list
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(block_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: sint
    TYPE(cell_type), POINTER                 :: cell, cell_mol
    TYPE(cp_2d_r_p_type), DIMENSION(:), &
      POINTER                                :: mosb_bcast
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_1d_type), POINTER      :: local_molecule
    TYPE(fm_mol_blocks_type), DIMENSION(:), &
      POINTER                                :: fm_mol_blocks
    TYPE(fm_mol_blocks_type), POINTER        :: mol_block
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(kg_fm_mol_set_type), POINTER        :: fm_mol_set
    TYPE(kg_fm_p_type), POINTER              :: kg_fm
    TYPE(kg_molbox_env_type), POINTER        :: molbox
    TYPE(kg_sub_pw_env_type), POINTER        :: kg_sub_pw_env
    TYPE(mol_mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mosa, mosb
    TYPE(molecule_kind_type), POINTER        :: imolecule_kind, &
                                                jmolecule_kind, molecule_kind
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(molecule_type), POINTER             :: imolecule, jmolecule, molecule
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_p
    TYPE(section_vals_type), POINTER         :: dft_section, kg_section

!------------------------------------------------------------------------------

    IF ( .NOT. ( do_potential .OR. do_forces .OR. do_overlap ) ) RETURN

    IF ( do_potential .OR. do_forces ) THEN
      dft_section => section_vals_get_subs_vals(kg_env%input,"DFT",error=error)
      kg_section => section_vals_get_subs_vals(dft_section,"KG",error=error)
      CALL section_vals_val_get(kg_section,"PAULI_REPULSION",&
                                r_val=epsmol,error=error)
      IF ( epsmol <= 0._dp ) RETURN
    END IF

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) output_unit = cp_logger_get_default_unit_nr(logger)

    NULLIFY(local_molecule,cell,molecule_set,particle_set,para_env)
    NULLIFY(qs_env,kg_fm,atomic_kind_set,dft_control,matrix_p)

    
    CPPrecondition(ASSOCIATED(kg_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
      CALL get_kg_env( kg_env=kg_env,local_molecules=local_molecule,&
                       cell=cell,molecule_set=molecule_set,&
                       particle_set=particle_set,para_env=para_env,&
                       sub_qs_env=qs_env,kg_fm_set=kg_fm,&
                       error=error )

      NULLIFY ( matrix_h, rho, force )
      CALL get_qs_env( qs_env=qs_env,kg_sub_pw_env=kg_sub_pw_env,&
                       atomic_kind_set=atomic_kind_set,&
                       matrix_h=matrix_h,rho=rho,force=force,&
                       dft_control=dft_control,error=error )

      matrix_p => rho%rho_ao
      nspin=dft_control%nspins
      NULLIFY(mosb_bcast)
      ALLOCATE(mosb_bcast(nspin),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

      nmol = SIZE(molecule_set)
      ALLOCATE ( refp(1:3,1:nmol), STAT=istat )
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      refp = 0._dp
      ALLOCATE ( r0(1:nmol), STAT=istat )
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      r0 = 0._dp

      natom = SIZE(particle_set)
      ALLOCATE ( atom_of_kind(1:natom), STAT=istat )
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               atom_of_kind=atom_of_kind,&
                               maxco=maxco,&
                               maxlgto=maxlgto,&
                               maxsgf=maxsgf)

      IF ( do_forces ) THEN
        nder = 1
        pauli_energy = 0._dp
      ELSE
        nder = 0
      END IF
      maxder = ncoset(nder)
      ldai = ncoset(maxlgto+nder)
      CALL init_orbital_pointers(ldai)

      ALLOCATE (ai_work(ldai,ldai,ncoset(nder)),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

      ALLOCATE (sab(maxco,maxco*maxder),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      sab(:,:) = 0.0_dp

      ALLOCATE (work(maxco,maxsgf),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      work(:,:) = 0.0_dp

      ALLOCATE (sint(maxder),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DO i=1,maxder
        NULLIFY (sint(i)%block)
      END DO

      ! calculate reference point for each molecule
      nao_max = 0
      nmol_kind = SIZE(local_molecule%n_el,1)
      DO imol_kind = 1, nmol_kind

        NULLIFY(molbox,cell_mol)
        molbox => kg_sub_pw_env%molbox_env_set(imol_kind)
        CALL get_molbox_env(molbox_env=molbox,cell_mol=cell_mol)
        rk = SQRT(SUM((cell_mol%hmat(:,1) + cell_mol%hmat(:,2) + &
                       cell_mol%hmat(:,3))**2))

        NULLIFY(fm_mol_set)
        fm_mol_set => kg_fm%kg_fm_mol_set(imol_kind)
        CALL get_kg_fm_mol_set(kg_fm_mol_set=fm_mol_set,n_ao=nao)
        nao_max = MAX(nao,nao_max)

        nmol_local = local_molecule%n_el(imol_kind)
        DO imol_local = 1, nmol_local
          imol = local_molecule%list(imol_kind)%array(imol_local)

          NULLIFY(molecule)
          molecule => molecule_set ( imol )
          CALL get_molecule(molecule,molecule_kind=molecule_kind,&
                            first_atom=first_atom,last_atom=last_atom)
          natom = last_atom - first_atom + 1
          rc = 0._dp
          DO ia = first_atom, last_atom
            IF ( ia == first_atom ) r1(:)=particle_set(ia)%r
            ra(:) = particle_set(ia)%r(:) - r1(:)
            rc(:) = rc(:) + pbc(ra(:), cell)
          END DO
          refp(:,imol) = rc(:) + r1(:)
          r0(imol) = rk*0.5_dp

        END DO
      END DO
      CALL mp_sum(refp,para_env%group)

      IF ( do_overlap ) THEN
        ALLOCATE (smol(1:nmol,1:nmol),STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        smol = 0._dp
      END IF
      IF ( do_potential .OR. do_forces ) THEN
        ALLOCATE (vpot(1:nao_max,1:nao_max),STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF
      IF ( do_forces ) THEN
        ALLOCATE (pmol(1:nao_max,1:nao_max),STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE (dvpot(1:nao_max,1:nao_max),STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE (fora(1:nao_max),STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE (forb(1:nao_max),STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF


      ! loop over all molecules to find interactions
      DO jmol = 1, nmol
        jmolecule => molecule_set ( jmol )
        CALL get_molecule(jmolecule,molecule_kind=jmolecule_kind,&
                          first_atom=jfa,last_atom=jla)
        natom_mol = jla-jfa +1

        jmol_kind = jmolecule_kind%kind_number
        NULLIFY(fm_mol_set,fm_mol_blocks)
        fm_mol_set => kg_fm%kg_fm_mol_set(jmol_kind)

        CALL get_kg_fm_mol_set(kg_fm_mol_set=fm_mol_set, &
                               n_ao=naob,n_mo=nmob, &
                               fm_mol_blocks=fm_mol_blocks)

        ALLOCATE(jfirst_ao_bcast(natom_mol),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(jlast_ao_bcast(natom_mol),STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        DO ispin = 1,nspin
          ALLOCATE(mosb_bcast(ispin)%array(naob,nmob(ispin)),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END DO 

! on which processor is jmol?
        j_pe =  kg_fm%imol_pe_pos(jmol)
! to which jmol_local does it correspond on its processor?
        jmol_local = kg_fm%imol_local_name(jmol)

        IF(para_env%mepos==j_pe) THEN
           NULLIFY(mol_block,jfirst_ao,jlast_ao,mosb)
           mol_block => fm_mol_blocks(jmol_local)
           CALL get_fm_mol_block(fm_mol_block=mol_block,&
                ifirst_ao=jfirst_ao,ilast_ao=jlast_ao,mos=mosb)
! get jfirst_ao amd jlast_ao from the right processor
           jfirst_ao_bcast(1:natom_mol) = jfirst_ao(1:natom_mol)
           jlast_ao_bcast(1:natom_mol) = jlast_ao(1:natom_mol)
! get mosb for jmol (jmol_local) from the right processor
           DO ispin = 1,nspin
             DO imo = 1,nmob(ispin)
                  mosb_bcast(ispin)%array(1:naob,imo) = mosb(ispin)%mo_set%mo(1:naob,imo) 
             END DO 
           END DO 
        END IF

        ! transfer data to other pe
        CALL mp_bcast(jfirst_ao_bcast,j_pe,para_env%group)
        CALL mp_bcast(jlast_ao_bcast,j_pe,para_env%group)
        DO ispin = 1,nspin
           CALL mp_bcast(mosb_bcast(ispin)%array,j_pe,para_env%group)
        END DO 

        ! loop over all molecules (local only)
        nmol_kind = SIZE(local_molecule%n_el,1)
        DO imol_kind = 1, nmol_kind
          nmol_local = local_molecule%n_el(imol_kind)
          DO imol_local = 1, nmol_local
            imol = local_molecule%list(imol_kind)%array(imol_local)
            imolecule => molecule_set ( imol )
            CALL get_molecule(imolecule,molecule_kind=imolecule_kind,&
                            first_atom=ifa,last_atom=ila)
            CALL get_molecule_kind(imolecule_kind,atom_list=iatom_list)
            NULLIFY(fm_mol_set,fm_mol_blocks)
            fm_mol_set => kg_fm%kg_fm_mol_set(imol_kind)
            CALL get_kg_fm_mol_set(kg_fm_mol_set=fm_mol_set, &
                                   natom=nat_mol, &
                                   n_ao=naoa,n_mo=nmoa, &
                                   fm_mol_blocks=fm_mol_blocks)
     
            mol_block => fm_mol_blocks(imol_local)
            CALL get_fm_mol_block(fm_mol_block=mol_block,&
                     ifirst_ao=ifirst_ao,ilast_ao=ilast_ao,mos=mosa)
            IF ( do_potential .OR. do_forces ) THEN
              vpot = 0._dp
            END IF
            IF ( do_forces ) THEN
              pmol = 0._dp
              DO i = 1, SIZE(matrix_p)
                CALL copy_sparse2mol_block(matrix_p(i)%matrix, mol_block,&
                        dvpot, nat_mol, naoa, naoa, error)
                pmol(1:naoa,1:naoa) = pmol(1:naoa,1:naoa) + dvpot(1:naoa,1:naoa)
              END DO
            END IF

            IF ( imol == jmol ) CYCLE
            rc(:) = refp(:,imol) - refp(:,jmol)
            ra(:) = pbc(rc(:),cell)   
            rij = SQRT( SUM ( ra(:)**2 ) )
            IF ( rij < r0(imol)+r0(jmol) ) THEN
              !Molecules interact, add potential to imol 
  
!MI FIX IT !!!!!!!!!!!!!!!!!!!!
! This cannot work in parallel as it is. fm_mol_blocks are local
! to get the info in  fm_mol_blocks(jmol_local) one has to be on the 
! corresponding processor. The indexes ifirst_ao and ilast_ao should
! correspond to the local matrix.
! If a smat matrix between two molecules that are on 2 different processor
! need to be constructed, the data have to be sent 
!              mol_block => fm_mol_blocks(jmol)
!              CALL get_fm_mol_block(fm_mol_block=mol_block,&
!                   ifirst_ao=jfirst_ao,ilast_ao=jlast_ao,mos=mosb)
!MI  FIX IT !!!!!!!!!!!!!!!!!!!!
              ALLOCATE (smat(1:naoa,1:naob,maxder),STAT=istat)
              CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

              CALL get_molecule_kind(jmolecule_kind,atom_list=jatom_list)
              DO ia = ifa, ila
                iao = ia - ifa + 1
                CALL get_atomic_kind(iatom_list(iao)%atomic_kind,&
                                     orb_basis_set=orb_basis_set)
                IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE
                CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                       first_sgf=first_sgfa,&
                                       lmax=la_max,&
                                       lmin=la_min,&
                                       npgf=npgfa,&
                                       nset=nseta,&
                                       nsgf_set=nsgfa,&
                                       pgf_radius=rpgfa,&
                                       set_radius=set_radius_a,&
                                       sphi=sphi_a,&
                                       zet=zeta)
                nia=SUM(nsgfa(1:nseta))

                DO ja = jfa, jla
                  jao = ja - jfa + 1
                  CALL get_atomic_kind(jatom_list(jao)%atomic_kind,&
                                       orb_basis_set=orb_basis_set)
                  IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE
                  CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                         first_sgf=first_sgfb,&
                                         lmax=lb_max,&
                                         lmin=lb_min,&
                                         npgf=npgfb,&
                                         nset=nsetb,&
                                         nsgf_set=nsgfb,&
                                         pgf_radius=rpgfb,&
                                         set_radius=set_radius_b,&
                                         sphi=sphi_b,&
                                         zet=zetb)
                  njb=SUM(nsgfb(1:nsetb))
                  DO i=1,maxder
                    ALLOCATE (sint(i)%block(nia,njb),STAT=istat)
                    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
                    sint(i)%block(:,:) = 0._dp
                  END DO
                  rc(:) = particle_set(ia)%r(:) - particle_set(ja)%r(:)
                  ra(:) = pbc(rc(:),cell)   
                  CALL atom_overlap_matrix(ia,nseta,set_radius_a,npgfa,&
                    nsgfa,la_max,la_min,first_sgfa,rpgfa,zeta,sphi_a,&
                    ja,nsetb,set_radius_b,npgfb,nsgfb,lb_max,&
                    lb_min,first_sgfb,rpgfb,zetb,sphi_b,&
                    ra,nder,ldai,sint,sab,work,ai_work,error)

                  DO i=1,maxder
                    smat(ifirst_ao(iao):ilast_ao(iao),&
                         jfirst_ao_bcast(jao):jlast_ao_bcast(jao),i)=sint(i)%block(:,:)
                  END DO
                  DO i=1,maxder
                    DEALLOCATE (sint(i)%block,STAT=istat)
                    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
                  END DO

                END DO
              END DO
              IF ( do_overlap ) THEN
                oij = 0._dp
                DO i = 1, nspin
                  oij = oij + SUM ( &
                    MATMUL( &
                    MATMUL(TRANSPOSE(mosa(i)%mo_set%mo(1:naoa,1:nmoa(i))),&
                           smat(1:naoa,1:naob,1) ), &
                    mosb_bcast(i)%array(1:naob,1:nmob(i)) )**2 )
                END DO
                smol(imol,jmol) = SQRT(oij)
              END IF
              IF ( do_potential .OR. do_forces ) THEN
                vpot(1:naoa,1:naoa) = vpot(1:naoa,1:naoa) + &
                  MATMUL(smat(1:naoa,1:naob,1),TRANSPOSE(smat(1:naoa,1:naob,1)))
              END IF
              IF ( do_forces ) THEN
                DO i = 1, 3
                  dvpot(1:naoa,1:naoa) = &
                    MATMUL(smat(1:naoa,1:naob,1+i),TRANSPOSE(smat(1:naoa,1:naob,1))) 
                  DO j = 1, naoa
                    fora(j) = DOT_PRODUCT(dvpot(j,1:naoa),pmol(1:naoa,j))
                  END DO
                  DO ia = ifa, ila
                    iao = ia - ifa + 1
                    CALL get_atomic_kind(iatom_list(iao)%atomic_kind,&
                                     kind_number=ikind)
                    atom_a = atom_of_kind(ia)
                    force(ikind)%kg_gpw_pauli(i,atom_a) = &
                         force(ikind)%kg_gpw_pauli(i,atom_a) + &
                         2._dp*epsmol*SUM ( fora(ifirst_ao(iao):ilast_ao(iao)) )
                  END DO
                  dvpot(1:naoa,1:naob) = &
                    MATMUL(pmol(1:naoa,1:naoa),smat(1:naoa,1:naob,1)) 
                  DO j = 1, naob
                    forb(j) = -DOT_PRODUCT(dvpot(1:naoa,j),smat(1:naoa,j,1+i))
                  END DO
                  DO ja = jfa, jla
                    jao = ja - jfa + 1
                    CALL get_atomic_kind(iatom_list(jao)%atomic_kind,&
                                     kind_number=jkind)
                    atom_a = atom_of_kind(ja)
                    force(jkind)%kg_gpw_pauli(i,atom_a) = &
                         force(jkind)%kg_gpw_pauli(i,atom_a) + &
                         2._dp*epsmol*SUM ( forb(jfirst_ao_bcast(jao):jlast_ao_bcast(jao)) )
                  END DO
                END DO
              END IF

              DEALLOCATE (smat,STAT=istat)
              CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

            END IF


            IF ( do_potential ) THEN
              DO ia=ifa,ila
                iao = ia - ifa + 1
                DO ja=ia,ila
                 jao = ja - ifa + 1
                  NULLIFY(h_block)
                  CALL get_block_node(matrix=matrix_h(1)%matrix,&
                                    block_row=ia,&
                                    block_col=ja,&
                                    BLOCK=h_block)
                  h_block(:,:)=h_block(:,:) + epsmol * &
                  vpot(ifirst_ao(iao):ilast_ao(iao),ifirst_ao(jao):ilast_ao(jao))
                END DO
              END DO
            END IF
            IF ( do_forces ) THEN
              pauli_energy = pauli_energy + epsmol * &
                           SUM ( pmol(1:naoa,1:naoa)*vpot(1:naoa,1:naoa))
            END IF
          END DO  ! imol_local
        END DO  ! imol_kind
        DEALLOCATE(jfirst_ao_bcast,jlast_ao_bcast,STAT=istat)
        DO ispin = 1,nspin
           DEALLOCATE(mosb_bcast(ispin)%array,STAT=istat)
        END DO 
      END DO  ! jmol

      IF ( do_overlap ) THEN
        CALL mp_sum(smol,para_env%group)
        IF (ionode) THEN
          WRITE(output_unit,'(/,A)') " Molecular Wavefunction Overlap"
          DO ia = 1, nmol - 1
            WRITE(output_unit,'(A,i6,1000(T21,4(F9.5,"[",i4,"]")))') &
              " Molecule:",ia,(smol(ia,ja),ja,ja=ia+1,nmol)
          END DO
        END IF
      END IF
      IF ( do_forces ) THEN
        CALL mp_sum(pauli_energy,para_env%group)
        IF (ionode) THEN
          WRITE(output_unit,'(/,A,T60,F20.10)') &
                "  Pauli Potential Energy [a.u.] ",pauli_energy
        END IF
      END IF

      DEALLOCATE (mosb_bcast,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE ( refp, STAT=istat )
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE ( r0, STAT=istat )
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE ( ai_work, STAT=istat )
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE ( work, STAT=istat )
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE ( sab, STAT=istat )
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DO i=1,maxder
         NULLIFY (sint(i)%block)
       END DO
      DEALLOCATE ( sint, STAT=istat )
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      IF ( do_overlap ) THEN
        DEALLOCATE (smol,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF
      IF ( do_potential .OR. do_forces ) THEN
        DEALLOCATE (vpot,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF
      IF ( do_forces ) THEN
        DEALLOCATE (dvpot,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        DEALLOCATE (pmol,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        DEALLOCATE (fora,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        DEALLOCATE (forb,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

    END IF

    CALL timestop(0.0_dp,handle)

END SUBROUTINE kg_gpw_pauli

!******************************************************************************

END MODULE kg_pauli_potential

!******************************************************************************
