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
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE global_types,                    ONLY: WANNIER_CRAZY,&
                                             WANNIER_DIRECT,&
                                             WANNIER_JACOBI,&
                                             WANNIER_NONE
  USE kinds,                           ONLY: dp,default_string_length
  USE mathconstants,                   ONLY: twopi
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_max,&
                                             mp_minloc
  USE molecule_kind_types,             ONLY: atom_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
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
  USE topology_types,                  ONLY: deallocate_topology,&
                                             init_topology,&
                                             topology_parameters_type
  USE topology_util,                   ONLY: find_boundary,&
                                             topology_generate_bond,&
                                             topology_generate_bend,&
                                             topology_generate_dihe,&
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
       TYPE(distribution_1d_type), POINTER         :: local_molecules

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

  CALL timeset("qs_loc_init","I",'',handle)
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
                   local_molecules=qs_loc_env%local_molecules,&
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
SUBROUTINE qs_loc_compute_properties(qs_loc_env, molecule_set, vectors, center, ionode, iunit)

  TYPE(qs_loc_env_type), POINTER :: qs_loc_env
  TYPE ( molecule_type ), POINTER :: molecule_set ( : )
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

  CALL timeset("qs_loc_compute_properties","I",'',handle)
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
    CALL wfc_to_molecule ( qs_loc_env, center, molecule_set )
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
SUBROUTINE wfc_to_molecule ( qs_loc_env, center, molecule_set ) 
  IMPLICIT NONE
! arguments
  TYPE ( qs_loc_env_type ), INTENT ( IN ) :: qs_loc_env
  TYPE ( molecule_type ), POINTER :: molecule_set ( : )
  TYPE (molecule_kind_type), POINTER :: molecule_kind
  REAL ( KIND = dp ), INTENT ( IN ) :: center ( :, : )

! locals
  TYPE ( particle_type ), POINTER :: particle_set ( : )
  TYPE ( cp_para_env_type ), POINTER :: para_env
  INTEGER :: i, j , k, isos, iatom, istate, imoltyp, imolnum,first_atom
  INTEGER :: counter, imol, natom, nbond, nmol, nstate, first, last, imol_now
  INTEGER :: oldtyp,oldnum
  REAL ( KIND = dp ) :: r_pbc ( 3 ), dr (3)
  REAL ( KIND = dp ), POINTER :: distance ( : )
  REAL ( KIND = dp ), POINTER :: r ( :, : )
  INTEGER, POINTER :: wfc_to_atom_map ( : )

  INTEGER :: local_location
  INTEGER :: natom_max,natom_loc,nkind,ikind
  REAL ( KIND = dp ) :: mydist(2)
  TYPE ( topology_parameters_type ) :: topology
  TYPE(atom_type), DIMENSION(:), POINTER            :: atom_list
  CHARACTER (LEN=default_string_length)  :: mol_name
  TYPE(distribution_1d_type), POINTER :: local_molecules


!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  CALL init_topology(topology)

  topology%cell => qs_loc_env % cell
  particle_set => qs_loc_env % particle_set
  para_env => qs_loc_env % para_env
  local_molecules => qs_loc_env % local_molecules
  nstate = SIZE ( center, 1 )

  ALLOCATE(wfc_to_atom_map(nstate),STAT=isos)
  IF(isos/=0) CALL stop_memory('wfc_to_molecule','wfc_to_atom_map')
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  nkind = SIZE (local_molecules%n_el)
  natom = 0
  natom_max = 0
  DO ikind=1,nkind
    nmol = SIZE ( local_molecules%list(ikind)%array )
    !WRITE(*,*) "IKIND, NMOL = ",ikind,nmol,local_molecules%list(ikind)%array
    DO imol=1,nmol
      i = local_molecules%list(ikind)%array(imol)
      molecule_kind => molecule_set(i)%molecule_kind
      CALL get_molecule_kind(molecule_kind=molecule_kind,&
                             natom=natom)
      natom_max = natom_max + natom
    END DO
  END DO
  natom_loc = natom_max
  natom = natom_max

  topology%natoms = natom_max

  CALL mp_max(natom_max,para_env%group)

  !WRITE(*,*) "natom_loc  = ",natom_loc
  !WRITE(*,*) "natom_max  = ",natom_max
  !WRITE(*,*) "nstate = ",nstate

  ALLOCATE(topology%r(3,natom_max),STAT=isos)
  IF(isos/=0) CALL stop_memory('centers_and_spreads','r')
  ALLOCATE(topology%label_atmname(natom_max),STAT=isos)
  IF(isos/=0) CALL stop_memory('centers_and_spreads','label_atmname')
  ALLOCATE(topology%label_resname(natom_max),STAT=isos)
  IF(isos/=0) CALL stop_memory('centers_and_spreads','label_resname')

  ALLOCATE(distance(natom_max),STAT=isos)
  IF(isos/=0) CALL stop_memory('wfc_to_molecule','distance')
  ALLOCATE(topology%atm_mass(natom_max),STAT=isos)
  IF(isos/=0) CALL stop_memory('wfc_to_molecule','atm_mass')
  ALLOCATE(topology%atm_charge(natom_max),STAT=isos)
  IF(isos/=0) CALL stop_memory('wfc_to_molecule','atm_charge')

  !Zero all the stuff
  topology%r(:,:) = 0.0_dp
  topology%label_atmname(:) = ''
  topology%label_resname(:) = ''
  topology%atm_mass(:) = 0.0_dp
  topology%atm_charge(:) = 0.0_dp
  distance(:) = 1.E10_dp

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  counter = 0
  DO ikind=1,nkind
    nmol = SIZE ( local_molecules%list(ikind)%array )
    DO imol=1,nmol
      i = local_molecules%list(ikind)%array(imol)
      molecule_kind => molecule_set(i)%molecule_kind
      first_atom = molecule_set(i)%first_atom
      CALL get_molecule_kind(molecule_kind=molecule_kind,&
                             natom=natom,&
                             name=mol_name,&
                             atom_list=atom_list)
      DO iatom=1,natom
        counter = counter + 1
        topology%r ( :, counter ) = particle_set(first_atom+iatom-1) % r ( : )
        CALL get_atomic_kind (atom_list(iatom)%atomic_kind, &
                              element_symbol=topology%label_atmname(counter) ) 
        topology%label_resname(counter) = mol_name
      END DO
    END DO
  END DO

! DO iatom=1,48
!   WRITE(*,*) iatom,topology%r(:,iatom)
! END DO

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  DO istate=1,nstate
    distance(:) = 1.E10_dp
    DO iatom=1,natom_loc
      dr(1) = topology%r(1,iatom) - center(istate,1)
      dr(2) = topology%r(2,iatom) - center(istate,2)
      dr(3) = topology%r(3,iatom) - center(istate,3)
      r_pbc = pbc ( dr, topology%cell )
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
  IF (natom_loc /= 0) THEN
    ALLOCATE(topology%map_mol_typ(natom_loc),STAT=isos)
    IF(isos/=0) CALL stop_memory('wfc_to_molecule','map_mol_typ')
    ALLOCATE(topology%map_mol_num(natom_loc),STAT=isos)
    IF(isos/=0) CALL stop_memory('wfc_to_molecule','map_mol_num')
  
    CALL topology_generate_bond(topology)
    CALL topology_generate_bend(topology)
    CALL topology_generate_dihe(topology)
    CALL topology_generate_molecule(topology)
  END IF
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  IF(natom_loc/=0) THEN
    DO istate = 1,nstate
      iatom = wfc_to_atom_map(istate)
      IF(iatom/=0) THEN
        imoltyp = topology%map_mol_typ(iatom)
        imolnum = topology%map_mol_num(iatom)
     
        imol_now = 0
        oldtyp = 0
        oldnum = 0
        DO i=1,iatom
          IF(.NOT.((topology%map_mol_typ(i)==oldtyp).AND.&
            (topology%map_mol_num(i)==oldnum))) THEN
            imol_now = imol_now + 1
            oldtyp = topology%map_mol_typ(i)
            oldnum = topology%map_mol_num(i)
          END IF
        END DO
        
        !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        !------------------------------------------------
        !Need to generalize for more than one kinds!!!---
        imol_now = local_molecules%list(1)%array(imol_now)
        !------------------------------------------------
        !------------------------------------------------
        !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

        i = molecule_set(imol_now)%lmi%nstates
        i = i + 1
        molecule_set(imol_now)%lmi%nstates = i
  
        !WRITE(*,*) "YY",para_env%mepos,istate,imol_now,i
        CALL reallocate(molecule_set(imol_now)%lmi%states,1,i)
        molecule_set(imol_now)%lmi%states(i) = istate
      END IF
    END DO
  END IF

! DO ikind=1,nkind
!   nmol = SIZE ( local_molecules%list(ikind)%array )
!   DO imol=1,nmol
!     i = local_molecules%list(ikind)%array(imol)
!     WRITE(*,*)"XX",para_env%mepos,imol,molecule_set(i)%lmi%states
!     !WRITE(*,*)"XX",para_env%mepos,imol,i
!   END DO
! END DO
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  DEALLOCATE(distance,STAT=isos)
  IF(isos/=0) CALL stop_memory('wfc_to_molecule','distance')

  DEALLOCATE(wfc_to_atom_map,STAT=isos)
  IF(isos/=0) CALL stop_memory('wfc_to_molecule','wfc_to_atom_map')
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  
  CALL deallocate_topology(topology)
  !CALL stop_program("wfc_to_atom","checkpoint")

END SUBROUTINE wfc_to_molecule
!------------------------------------------------------------------------------

END MODULE qs_localization

