!   CP2K: A general program to perform molecular dynamics simulations         !

!!****** cp2k/qs_localization [1.0] *
!!
!!   NAME
!!     qs_localization
!!
!!   FUNCTION
!! 
!!     should be able to find properties related to the position operator in PBC
!!     e.g. location and spread of periodic orbitals
!!          Wannier functions
!!
!!   AUTHORS
!!     Joost VandeVondele with some help of Leonardo Bernasconi (01.2003)
!!     CJM (04.2003)
!!     CJM & IKUO (08.2003)
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE qs_localization
! *****************************************************************************
  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE coefficient_types,               ONLY: coeff_type,&
                                             coeff_zero
  USE cp_cfm_basic_linalg,             ONLY: cp_cfm_lu_decompose
  USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                             cp_cfm_release,&
                                             cp_cfm_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_lu_decompose
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_element,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_type,&
                                             cp_fm_write
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE global_types,                    ONLY: WANNIER_CRAZY,&
                                             WANNIER_DIRECT,&
                                             WANNIER_JACOBI,&
                                             WANNIER_NONE
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: twopi
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_max,&
                                             mp_minloc
  USE molecule_types,                  ONLY: linklist_atoms,&
                                             molecule_structure_type
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: angstrom
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_give_back_coeff,&
                                             pw_pool_init_coeff,&
                                             pw_pool_type
  USE pw_types,                        ONLY: REALDATA3D,&
                                             REALSPACE
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_integrate_potential,          ONLY: integrate_v_rspace
  USE qs_localization_methods,         ONLY: crazy_rotations,&
                                             direct_mini,&
                                             initialize_weights,&
                                             jacobi_rotations
  USE simulation_cell,                 ONLY: cell_type,&
                                             pbc
  USE sparse_matrix_types,             ONLY: deallocate_matrix,&
                                             real_matrix_p_type,&
                                             real_matrix_type,&
                                             replicate_matrix_structure,&
                                             set_matrix
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_util,                   ONLY: find_boundary,&
                                             topology_generate_bond,&
                                             topology_generate_molecule
  IMPLICIT NONE

! *** Public ***
  PUBLIC :: qs_loc_init
  PUBLIC :: qs_loc_env_type, PIPEK_MENZY, POSITION_OPERATOR
  PUBLIC :: qs_loc_compute_properties
  PUBLIC :: qs_loc_destroy

! *** Localization Data Type ***
  TYPE qs_loc_env_type
       INTEGER :: operator_type
       INTEGER :: localization_method
       LOGICAL :: wannier
       LOGICAL :: molecular_states
       INTEGER :: is_initialized
       TYPE ( real_matrix_p_type ),POINTER :: OPERATOR(:,:)
       TYPE ( cp_para_env_type ), POINTER :: para_env
       TYPE ( cell_type ), POINTER :: cell
       TYPE ( particle_type ), POINTER :: particle_set ( : )
       REAL (KIND = dp) :: weights ( 3 )
  END TYPE qs_loc_env_type
  
  INTEGER, PARAMETER :: PIPEK_MENZY=301, POSITION_OPERATOR=302

! *****************************************************************************

CONTAINS

!******************************************************************************
!
! allocates the data, and initializes the operators
!
!******************************************************************************
SUBROUTINE qs_loc_init(qs_loc_env,qs_env)
  IMPLICIT NONE

  TYPE(qs_loc_env_type), POINTER :: qs_loc_env
  TYPE(qs_environment_type) :: qs_env
! locals
  TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: s
  INTEGER :: i, j, natoms, isos, handle

  CALL timeset("qs_loc_init","I","",handle)
  NULLIFY(qs_loc_env%para_env)
  NULLIFY(qs_loc_env%cell)
  NULLIFY(qs_loc_env%particle_set)
  !
  ! start replicating the s structure, we do have
  ! 6 operators (x,y,z)*(real,imag) for POSITION_OPERATOR
  ! and we have NATOMS operators for PIPEK_MENZY
  ! notice that the operator should have the same structure as S
  ! borrow para env from S
  !
  CALL get_qs_env(qs_env=qs_env,s=s,cell=qs_loc_env%cell, &
                   particle_set = qs_loc_env % particle_set, &
                   para_env=qs_loc_env%para_env)
  CALL cp_para_env_retain(qs_loc_env%para_env) 

! dbg
  qs_loc_env % operator_type = POSITION_OPERATOR
! dbg
  IF ( qs_loc_env % operator_type == POSITION_OPERATOR )  THEN
    ALLOCATE ( qs_loc_env % operator ( 2, 3 ), stat = isos )
    IF ( isos /= 0 ) CALL stop_memory ( 'qs_localization', 'operator' )
  ELSEIF ( qs_loc_env % operator_type == PIPEK_MENZY ) THEN
    natoms = SIZE ( qs_loc_env % particle_set )
    ALLOCATE ( qs_loc_env % operator ( 1, natoms ), stat = isos )
    IF ( isos /= 0 ) CALL stop_memory ( 'qs_localization', 'operator' )
  ENDIF

  DO i=1, SIZE ( qs_loc_env % operator, 2 )
     DO j=1, SIZE ( qs_loc_env % operator, 1 )
       NULLIFY(qs_loc_env%operator(j,i)%matrix)
       CALL replicate_matrix_structure(s(1)%matrix, &
            qs_loc_env%operator(j,i)%matrix,"operator")
       CALL set_matrix(qs_loc_env%operator(j,i)%matrix,0.0_dp)
     ENDDO
  END DO


  IF ( qs_loc_env % operator_type == POSITION_OPERATOR ) THEN 
    CALL initialize_weights ( qs_loc_env % cell, qs_loc_env % weights ) 
    CALL get_boys_foster_operator ( qs_loc_env, qs_env )
  ENDIF

  IF ( qs_loc_env % operator_type == PIPEK_MENZY ) & 
  CALL get_pipek_menzy_operator ( qs_loc_env, qs_env )
  
  qs_loc_env%molecular_states = .FALSE.

! done
  qs_loc_env%is_initialized = 777
  CALL timestop(0.0_dp,handle) 
END SUBROUTINE qs_loc_init

!******************************************************************************
!
! Computes the Boys-Foster operator for periodic systems
!
!******************************************************************************
SUBROUTINE get_boys_foster_operator ( qs_loc_env, qs_env ) 

! globals
  TYPE(qs_loc_env_type), POINTER :: qs_loc_env
  TYPE(qs_environment_type) :: qs_env
! locals
  TYPE(coeff_type) :: op_rspace
  TYPE(pw_env_type), POINTER :: pw_env
  TYPE(pw_pool_type), POINTER :: auxbas_pw_pool
  INTEGER ix, iy, iz, i, j
  REAL(KIND = dp) :: x, y, z, thearg, theval

  NULLIFY ( pw_env, auxbas_pw_pool )
!
! OK, now we construct these operators the easy way, put them on the pw grid and integrate
! analytical construction should be faster.
!
! get a pw grid
  CALL get_qs_env(qs_env=qs_env,pw_env=pw_env)
  CALL pw_env_get(pw_env=pw_env,auxbas_pw_pool=auxbas_pw_pool)
  CALL pw_pool_init_coeff(auxbas_pw_pool,op_rspace,&
                          use_data=REALDATA3D,in_space=REALSPACE)

  DO i=1,3
   DO j=1,2
      ! CALL coeff_zero(op_rspace) 
      ! the definition of the position on the PW grid has to be 
      ! compatible with the one used in the integration procedure / collocation procedure
      ! integrate_v_rspace calculate_rho_elec
      DO iz=op_rspace%pw%pw_grid%bounds_local(1,3),op_rspace%pw%pw_grid%bounds_local(2,3)
      DO iy=op_rspace%pw%pw_grid%bounds_local(1,2),op_rspace%pw%pw_grid%bounds_local(2,2)
      DO ix=op_rspace%pw%pw_grid%bounds_local(1,1),op_rspace%pw%pw_grid%bounds_local(2,1)
         x=REAL(ix-op_rspace%pw%pw_grid%bounds(1,1),dp)*twopi/REAL(op_rspace%pw%pw_grid%npts(1),dp)
         y=REAL(iy-op_rspace%pw%pw_grid%bounds(1,2),dp)*twopi/REAL(op_rspace%pw%pw_grid%npts(2),dp)
         z=REAL(iz-op_rspace%pw%pw_grid%bounds(1,3),dp)*twopi/REAL(op_rspace%pw%pw_grid%npts(3),dp) 
         SELECT CASE (i)
         CASE (1)
            thearg=x
         CASE (2)
            thearg=y
         CASE (3)
            thearg=z
         END SELECT
         SELECT CASE (j)
         CASE (1)
            theval=COS(thearg)
         CASE (2)
            theval=SIN(thearg)
         END SELECT
         ! notice the volume element, slightly ugly I think
         op_rspace%pw%cr3d(ix,iy,iz)=theval*op_rspace%pw%pw_grid%dvol 
      ENDDO
      ENDDO
      ENDDO
! do the actual integration
      CALL integrate_v_rspace(op_rspace,qs_loc_env%operator(j,i),qs_env)
   ENDDO
  ENDDO

  CALL pw_pool_give_back_coeff ( auxbas_pw_pool, op_rspace )

END SUBROUTINE get_boys_foster_operator

!******************************************************************************
!
! Computes the Pipek-Menzy operator for periodic systems
!
!******************************************************************************
SUBROUTINE get_pipek_menzy_operator ( qs_loc_env, qs_env )
  IMPLICIT NONE
! globals
  TYPE(qs_loc_env_type), POINTER :: qs_loc_env
  TYPE(qs_environment_type) :: qs_env
  CALL stop_program("qs_localization","PIPEK MENZY not implemented")
END SUBROUTINE get_pipek_menzy_operator

!******************************************************************************
!
! applies the operators to obtain the necessary data
!
!******************************************************************************
SUBROUTINE qs_loc_compute_properties(qs_loc_env, molecule, vectors, center, ionode, iunit)

  TYPE(qs_loc_env_type), POINTER :: qs_loc_env
  TYPE ( molecule_structure_type ), POINTER :: molecule ( : )
  TYPE(cp_fm_type), POINTER :: vectors
  REAL ( KIND = dp ), DIMENSION(:,:), POINTER :: center 
  LOGICAL :: ionode
  INTEGER :: iunit
! locals
  TYPE(cp_fm_type), POINTER :: opvec
  TYPE(cp_fm_p_type), ALLOCATABLE :: ZIJ ( :, : )
  TYPE(cp_fm_struct_type), POINTER :: tmp_fm_struct

  INTEGER :: n,k,natoms
  INTEGER :: i,j,l,isos
  INTEGER :: handle,sweeps

  CALL timeset("qs_loc_compute_properties","I","",handle)
  NULLIFY(tmp_fm_struct,opvec)

  ! get rows and cols of the input
  CALL cp_fm_get_info(vectors,nrow_global=n,ncol_global=k)
  ! replicate the input kind of matrix
  CALL cp_fm_create(opvec,vectors%matrix_struct)
  ! create a new matrix structure k x k in size, borrowing most of info
  CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=k, &
            ncol_global=k,para_env=vectors%matrix_struct%para_env, &
            context=vectors%matrix_struct%context)

  IF (ANY(SHAPE(center).NE.(/k,6/))) CALL stop_program('centers_and_spreads','center')

! Allocating ZIJ
  IF ( qs_loc_env % operator_type == POSITION_OPERATOR ) THEN
    ALLOCATE ( ZIJ ( 2, 3 ), stat = isos )
    IF ( isos /= 0 ) CALL stop_memory ( 'qs_localization', 'ZIJ' )
  ELSEIF ( qs_loc_env % operator_type == PIPEK_MENZY ) THEN
    natoms = SIZE ( qs_loc_env % operator, 2 ) 
    ALLOCATE ( ZIJ ( 1, natoms ), stat = isos )
    IF ( isos /= 0 ) CALL stop_memory ( 'qs_localization', 'ZIJ' )
  ENDIF
! make the actual matrices
  DO i=1,SIZE ( ZIJ, 2 )
    DO j=1, SIZE ( ZIJ, 1 )
        NULLIFY ( ZIJ ( j, i ) % matrix )
        CALL cp_fm_create ( ZIJ ( j, i ) % matrix, tmp_fm_struct )
    ENDDO
  ENDDO

! release one instance of the fm_struct
  CALL cp_fm_struct_release ( tmp_fm_struct ) 

! Compute zij here
  DO i=1,SIZE ( ZIJ, 2 )
    DO j=1,SIZE ( ZIJ, 1 )
       CALL cp_sm_fm_multiply(qs_loc_env%operator(j,i)%matrix, &
                            vectors,opvec,ncol=k,para_env=qs_loc_env%para_env)
       CALL cp_fm_gemm("T","N",k,k,n,1.0_dp,vectors,opvec,0.0_dp,ZIJ(j,i)%matrix)
    ENDDO
  ENDDO

  IF ( qs_loc_env % wannier ) THEN
       SELECT CASE (qs_loc_env % localization_method) 
       CASE ( WANNIER_JACOBI )
          CALL jacobi_rotations ( qs_loc_env % weights, zij, vectors, sweeps ) 
       CASE ( WANNIER_CRAZY )
          CALL crazy_rotations  ( qs_loc_env % weights, zij, vectors, sweeps ) 
       CASE ( WANNIER_DIRECT )
          CALL direct_mini      ( qs_loc_env % weights, zij, vectors, sweeps ) 
       CASE DEFAULT
          CALL stop_program("qs_localization","unknown localization method")
       END SELECT
       IF (ionode) WRITE(6,'(A,I6,A)') " Wannier localization converged in ",sweeps," iterations"
  ELSE
       IF (ionode) WRITE(6,'(A,I6,A)') " No Wannier localization applied "
  ENDIF

  IF ( qs_loc_env % operator_type == POSITION_OPERATOR )  THEN
    CALL centers_and_spreads ( ionode, iunit, qs_loc_env, zij, center ) 
    IF (qs_loc_env % molecular_states ) &
    CALL wfc_to_molecule ( qs_loc_env, center, molecule )
  ENDIF 


! release results
  DO i=1,SIZE ( ZIJ, 2 )
    DO j=1, SIZE ( ZIJ, 1 )
       CALL cp_fm_release(ZIJ(j,i)%matrix)
    ENDDO
  ENDDO
 
  DEALLOCATE(zij, STAT = isos)
  IF(isos/=0) CALL stop_memory('qs_loc_compute_properties','zij')
  CALL cp_fm_release(opvec)

  CALL timestop(0.0_dp,handle) 

END SUBROUTINE qs_loc_compute_properties

!******************************************************************************
!
! computes spread and centers
!
!******************************************************************************

SUBROUTINE centers_and_spreads ( ionode, iunit, qs_loc_env, zij, center )
! globals
  INTEGER, INTENT ( IN ) :: iunit
  LOGICAL, INTENT ( IN ) :: ionode
  REAL ( KIND = dp ), INTENT ( OUT ) :: center ( :, : )
  TYPE ( cp_fm_p_type ), INTENT ( INOUT ) :: zij ( :, : )
  TYPE(qs_loc_env_type), POINTER :: qs_loc_env

! locals
  TYPE ( cell_type ), POINTER :: cell
  TYPE(cp_cfm_type), POINTER :: complexresult
  COMPLEX ( KIND = dp ) :: almost_determinant
  COMPLEX ( KIND = dp ) :: z
  REAL ( KIND = dp ) :: realpart, imagpart
  INTEGER :: l, i, nstates

! get nstates
  CALL cp_fm_get_info ( zij ( 1, 1 ) % matrix, nrow_global = nstates )

! initialize complexresult for boys-foster
  NULLIFY ( complexresult )
  CALL cp_cfm_create ( complexresult, zij ( 1, 1 ) % matrix % matrix_struct )
!
  
! OK print out the centers ( center ( 1:3) )  and spreads ( center ( 4:6 ) ) of these friends
  cell => qs_loc_env % cell
  DO l=1,nstates
    DO i=1,3
       CALL cp_fm_get_element(zij(1,i)%matrix,l,l,realpart)
       CALL cp_fm_get_element(zij(2,i)%matrix,l,l,imagpart)
       z=CMPLX(realpart,imagpart)
       center (l ,i )=( cell % hmat ( i, i ) / twopi ) * AIMAG ( LOG ( z ) )
       center (l, i+3 )=SQRT( -qs_loc_env % weights ( i ) * 2.0_dp * LOG ( ABS ( Z ) ) ) / twopi
    ENDDO
    IF (ionode) WRITE(iunit,'(A,6F16.8)') "He", angstrom*center(l, 1 : 6 )
  ENDDO
  
  ! now, find the spread of the orbitals
  DO i=1,3
     complexresult%local_data ( :, : ) = CMPLX (ZIJ(1,i)%matrix%local_data(:,:), &
                                         ZIJ(2,i)%matrix%local_data(:,:))
     CALL cp_cfm_lu_decompose ( complexresult, almost_determinant )
     IF (ionode) WRITE(6,'(A,I6,F16.10)') " D Eq. 11 ",i,SQRT(-2.0*nstates*LOG(ABS(almost_determinant)))
  ENDDO

  CALL cp_cfm_release(complexresult)

END SUBROUTINE centers_and_spreads

!******************************************************************************
!
! deallocates all data 
!
!******************************************************************************
SUBROUTINE qs_loc_destroy ( qs_loc_env )

  TYPE(qs_loc_env_type), POINTER :: qs_loc_env

  INTEGER :: i , j

  IF (qs_loc_env%is_initialized.ne.777) CALL stop_program("qs_localization","problems")

  DO i=1,SIZE ( qs_loc_env % operator, 2 )
   DO j=1,SIZE ( qs_loc_env % operator, 1 )
     CALL deallocate_matrix(qs_loc_env%operator(j,i)%matrix)
   ENDDO
  END DO
  CALL cp_para_env_release(qs_loc_env%para_env) 

END SUBROUTINE qs_loc_destroy

!******************************************************************************
!
! maps wfc's to molecules 
!
!******************************************************************************
SUBROUTINE wfc_to_molecule ( qs_loc_env, center, molecule ) 
  IMPLICIT NONE
! arguments
  TYPE ( qs_loc_env_type ), INTENT ( IN ) :: qs_loc_env
  TYPE ( molecule_structure_type ), POINTER :: molecule ( : )
  REAL ( KIND = dp ), INTENT ( IN ) :: center ( :, : )

! locals
  TYPE ( particle_type ), POINTER :: particle_set ( : )
  TYPE ( cell_type ), POINTER :: cell
  TYPE ( cp_para_env_type ), POINTER :: para_env
  TYPE ( linklist_atoms ), POINTER :: ll_atoms
  INTEGER :: i, j , k, isos, iatom, istate, imoltyp, imolnum
  INTEGER :: counter, imol, natom, nbond, nmol, nstate, first, last, imol_now
  INTEGER :: oldtyp,oldnum
  INTEGER, POINTER :: bond_a ( : ), bond_b ( : )
  REAL ( KIND = dp ) :: r_pbc ( 3 ), dr (3)
  REAL ( KIND = dp ), POINTER :: distance ( : ), atm_mass ( : ), atm_charge ( : )
  REAL ( KIND = dp ), POINTER :: r ( :, : )
  CHARACTER ( LEN = 20 ), POINTER :: label_atmname ( : ), atm_mol_name ( : )
  INTEGER, POINTER :: map_mol_typ ( : ), map_mol_num ( : ) ,  &
                      wfc_to_atom_map ( : ), ll_counter ( : , :),  &
                      ll_list ( : , : , : )

  INTEGER :: local_location
  INTEGER :: natom_max
  REAL ( KIND = dp ) :: mydist(2)

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  cell => qs_loc_env % cell
  particle_set => qs_loc_env % particle_set
  para_env => qs_loc_env % para_env
  nstate = SIZE ( center, 1 )

  ALLOCATE(wfc_to_atom_map(nstate),STAT=isos)
  IF(isos/=0) CALL stop_memory('wfc_to_molecule','wfc_to_atom_map')
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  nmol = SIZE ( molecule )
  natom = 0
  DO imol=1,nmol
    natom = natom + molecule(imol)%natoms_mol
  END DO

  !Determine the max number of atoms across all the procs
  natom_max = natom
  CALL mp_max(natom_max,para_env%group)

  ALLOCATE(r(3,natom_max),STAT=isos)
  IF(isos/=0) CALL stop_memory('wfc_to_molecule','r')
  ALLOCATE(label_atmname(natom_max),STAT=isos)
  IF(isos/=0) CALL stop_memory('wfc_to_molecule','label_atmname')
  ALLOCATE(atm_mol_name(natom_max),STAT=isos)
  IF(isos/=0) CALL stop_memory('wfc_to_molecule','atm_mol_name')
  ALLOCATE(distance(natom_max),STAT=isos)
  IF(isos/=0) CALL stop_memory('wfc_to_molecule','distance')
  ALLOCATE(atm_mass(natom_max),STAT=isos)
  IF(isos/=0) CALL stop_memory('wfc_to_molecule','atm_mass')
  ALLOCATE(atm_charge(natom_max),STAT=isos)
  IF(isos/=0) CALL stop_memory('wfc_to_molecule','atm_charge')

  !Zero all the stuff
  DO iatom=1,natom_max
    r(:,iatom) = 0.0_dp
    label_atmname(iatom) = ''
    atm_mol_name(iatom) = ''
    atm_mass(iatom) = 0.0_dp
    atm_charge(iatom) = 0.0_dp
    distance(:) = 1.E10_dp
  END DO

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  counter = 0
  DO imol=1,nmol
    ll_atoms => molecule ( imol ) % ll_atoms
    DO iatom=1,molecule(imol)%natoms_mol
      counter = counter + 1
      r ( :, counter ) = ll_atoms % part % p % r ( : )
      CALL get_atomic_kind (ll_atoms % part % p % atomic_kind, &
                            element_symbol = label_atmname ( counter ) ) 
      atm_mol_name ( counter ) = molecule ( imol ) % mol_name
      ll_atoms => ll_atoms % next
    END DO
  END DO
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  DO istate=1,nstate
    distance(:) = 1.E10_dp
    DO iatom=1,natom
      !WRITE(*,*) "I'm here"
      dr(1) = r(1,iatom) - center(istate,1)
      dr(2) = r(2,iatom) - center(istate,2)
      dr(3) = r(3,iatom) - center(istate,3)
      r_pbc = pbc ( dr, cell )
      distance(iatom) = SQRT ( DOT_PRODUCT ( r_pbc, r_pbc ) )
    END DO

    !combine distance() from all procs

    local_location = MINLOC(distance,DIM=1)
    
    mydist(1) = distance(local_location)
    mydist(2) = para_env%mepos

    CALL mp_minloc(mydist,para_env%group)
    
    IF(mydist(2)==para_env%mepos) THEN
      wfc_to_atom_map(istate)  =  local_location
    ELSE
      wfc_to_atom_map(istate)  =  0
    END IF
  END DO
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  !If info is stored initially, then this is redundant
  IF (natom /= 0) THEN
    ALLOCATE(map_mol_typ(natom),STAT=isos)
    IF(isos/=0) CALL stop_memory('wfc_to_molecule','map_mol_typ')
    ALLOCATE(map_mol_num(natom),STAT=isos)
    IF(isos/=0) CALL stop_memory('wfc_to_molecule','map_mol_num')
  
    CALL topology_generate_bond(cell,r,label_atmname,bond_a,bond_b,&
                                    atm_mass,atm_charge)
  
    nbond = SIZE(bond_a)
  
    CALL topology_generate_molecule(natom,atm_mol_name,nbond,bond_a,bond_b,&
                                        map_mol_typ,map_mol_num)
  END IF
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  IF(natom/=0) THEN
    DO istate = 1,nstate
      iatom = wfc_to_atom_map(istate)
      IF(iatom/=0) THEN
        imoltyp = map_mol_typ(iatom)
        imolnum = map_mol_num(iatom)
     
        imol_now = 0
        oldtyp = 0
        oldnum = 0
        DO i=1,iatom
          IF(.NOT.((map_mol_typ(i)==oldtyp).AND.(map_mol_num(i)==oldnum))) THEN
            imol_now = imol_now + 1
            oldtyp = map_mol_typ(i)
            oldnum = map_mol_num(i)
          END IF
        END DO
        
        i = molecule(imol_now)%nstates
        i = i + 1
        molecule(imol_now)%nstates = i
  
        CALL reallocate(molecule(imol_now)%states,1,i)
        molecule(imol_now)%states(i) = istate
      END IF
    END DO
  END IF

  !DO imol=1,nmol
  !  WRITE(*,*)"imol= ",imol,molecule(imol)%states
  !END DO
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  DEALLOCATE(r,STAT=isos)
  IF(isos/=0) CALL stop_memory('wfc_to_molecule','r')
  DEALLOCATE(label_atmname,STAT=isos)
  IF(isos/=0) CALL stop_memory('wfc_to_molecule','label_atmname')
  DEALLOCATE(atm_mol_name,STAT=isos)
  IF(isos/=0) CALL stop_memory('wfc_to_molecule','atm_mol_name')
  DEALLOCATE(distance,STAT=isos)
  IF(isos/=0) CALL stop_memory('wfc_to_molecule','distance')

  IF(natom/=0) THEN
    DEALLOCATE(map_mol_typ,STAT=isos)
    IF(isos/=0) CALL stop_memory('wfc_to_molecule','map_mol_typ')
    DEALLOCATE(map_mol_num,STAT=isos)
    IF(isos/=0) CALL stop_memory('wfc_to_molecule','map_mol_num')
    DEALLOCATE(atm_mass,STAT=isos)
    IF(isos/=0) CALL stop_memory('wfc_to_molecule','atm_mass')
    DEALLOCATE(bond_b,STAT=isos)
    IF(isos/=0) CALL stop_memory('wfc_to_molecule','bond_b')
    DEALLOCATE(bond_a,STAT=isos)
    IF(isos/=0) CALL stop_memory('wfc_to_molecule','bond_a')
  ENDIF

  DEALLOCATE(atm_charge,STAT=isos)
  IF(isos/=0) CALL stop_memory('wfc_to_molecule','atm_charge')
  DEALLOCATE(wfc_to_atom_map,STAT=isos)
  IF(isos/=0) CALL stop_memory('wfc_to_molecule','wfc_to_atom_map')
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  
END SUBROUTINE wfc_to_molecule
!------------------------------------------------------------------------------

END MODULE qs_localization

