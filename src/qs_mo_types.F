!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_mo_types [1.0] *
!!
!!   NAME
!!     qs_mo_types
!!
!!   FUNCTION
!!     Definition and initialisation of the mo data type.
!!
!!   AUTHOR
!!     Matthias Krack (09.05.2001,MK)
!!
!!   MODIFICATION HISTORY
!!     - adapted to the new QS environment data structure (02.04.2002,MK)
!!     - set_mo_occupation added (17.04.02,MK)
!!     - correct_mo_eigenvalues added (18.04.02,MK)
!!     - calculate_density_matrix moved from qs_scf to here (22.04.02,MK)
!!     - mo_set_p_type added (23.04.02,MK)
!!     - PRIVATE attribute set for TYPE mo_set_type (23.04.02,MK)
!!     - calculate_w_matrix added for the calculation of the Pulay forces
!!       (23.05.02,MK)
!!
!!   SOURCE
!******************************************************************************

MODULE qs_mo_types

! *****************************************************************************

  USE kinds, ONLY: wp => dp,&
                   wp_size => dp_size

  USE atomic_kind_types,    ONLY: atomic_kind_type,&
                                  get_atomic_kind,&
                                  get_atomic_kind_set
  USE basis_set_types,      ONLY: get_gto_basis_set,&
                                  gto_basis_set_type
  USE qs_blacs,             ONLY: allocate_blacs_matrix,&
                                  blacs_matrix_type,&
                                  blacs_set_all,&
                                  blacs_syrk, &
                                  blacs_gemm, &
                                  copy_blacs_to_full_matrix,&
                                  deallocate_blacs_matrix,&
                                  get_blacs_info,&
                                  read_blacs_matrix,&
                                  write_blacs_matrix,&
                                  sparse_plus_blacs_blacst, &
                                  blacs_replicated_copy, &
                                  sparse_times_blacs, &
                                  blacs_scale_matrix,&
                                  replicate_blacs_matrix, &
                                  get_blacs_matrix_info

  USE global_types,         ONLY: global_environment_type
  USE particle_types,       ONLY: particle_type
  USE qs_parser,            ONLY: close_file,&
                                  open_file
  USE sparse_matrix_types,  ONLY: add_matrices,&
                                  allocate_matrix,&
                                  get_block_node,&
                                  real_matrix_type,&
                                  scale_matrix,&
                                  set_matrix

  USE message_passing,      ONLY: mp_sum
  USE orbital_pointers,     ONLY: nso
  USE termination,          ONLY: stop_memory,&
                                  stop_program
  USE timings,              ONLY: timeset,&
                                  timestop

  IMPLICIT NONE

  PRIVATE

  TYPE mo_set_type
    PRIVATE
    REAL(wp)                         :: maxocc
    INTEGER                          :: homo,lfomo,nao,nelectron,nmo
    REAL(wp), DIMENSION(:), POINTER  :: eigenvalues,occupation_numbers
    TYPE(blacs_matrix_type), POINTER :: eigenvectors
  END TYPE mo_set_type

  TYPE mo_set_p_type
    TYPE(mo_set_type), POINTER :: mo_set
  END TYPE mo_set_p_type

! *** Public data types ***

  PUBLIC :: mo_set_p_type,&
            mo_set_type

! *** Public subroutines ***

  PUBLIC :: allocate_mo_set,&
            calculate_density_matrix,&
            calculate_w_matrix,&
            correct_mo_eigenvalues,&
            deallocate_mo_set,&
            get_mo_set,&
            read_mo_set,&
            set_mo_occupation,&
            write_mo_set

! *****************************************************************************

  INTERFACE write_mo_set
    MODULE PROCEDURE write_mo_set_to_output_unit,write_mo_set_to_restart_unit, &
                     write_mo_set_to_restart_new
  END INTERFACE
  INTERFACE read_mo_set
    MODULE PROCEDURE read_mo_set_old, &
                     read_mo_set_from_restart_new
  END INTERFACE

  INTERFACE calculate_density_matrix
    MODULE PROCEDURE calculate_density_matrix_blacs,&
                     calculate_density_matrix_sparse
  END INTERFACE

  INTERFACE calculate_w_matrix
    MODULE PROCEDURE calculate_w_matrix_sparse
  END INTERFACE

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_mo_set(mo_set,nao,nmo,nelectron,maxocc,nrow_block,&
                             ncol_block,globenv)

!   Purpose: Allocate a wavefunction data structure.

!   History: - Creation (15.05.2001, MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_type), POINTER                :: mo_set
    REAL(wp), INTENT(IN)                      :: maxocc
    INTEGER, INTENT(IN)                       :: nao,ncol_block,nelectron,nmo,&
                                                 nrow_block

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE allocate_mo_set (MODULE qs_mo_types)"

!   *** Local variables ***

    INTEGER :: istat

!   ---------------------------------------------------------------------------

    ALLOCATE (mo_set,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"mo_set",0)

    IF (maxocc == 2.0_wp) THEN
      mo_set%maxocc = maxocc
    ELSE
      CALL stop_program(routine,"LSD not implemented yet")
    END IF

    IF (MODULO(nelectron,2) == 0) THEN
      mo_set%nelectron = nelectron
    ELSE
      CALL stop_program(routine,"Odd number of electrons")
    END IF

    mo_set%nao = nao
    mo_set%nmo = nmo

    mo_set%homo = 0
    mo_set%lfomo = 0

    ALLOCATE (mo_set%eigenvalues(nmo),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine,"mo_set%eigenvalues",nmo*wp_size)
    END IF
    mo_set%eigenvalues(:) = 0.0_wp

    ALLOCATE (mo_set%occupation_numbers(nmo),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine,"mo_set%occupation_numbers",nmo*wp_size)
    END IF
    CALL set_mo_occupation(mo_set=mo_set,smear=0.0_wp)

    CALL allocate_blacs_matrix(new_matrix=mo_set%eigenvectors,&
                               nrow_global=nao,&
                               ncol_global=nmo,&
                               nrow_block=nrow_block,&
                               ncol_block=ncol_block,&
                               name="MO EIGENVECTORS",&
                               globenv=globenv)

  END SUBROUTINE allocate_mo_set

! *****************************************************************************

  SUBROUTINE calculate_density_matrix_blacs(mo_set,density_matrix,globenv)

!   Purpose: Calculate the density matrix from the MO eigenvectors and the
!            MO occupation numbers.

!   History: - Creation (01.05.2001,MK)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: density_matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_type), POINTER                :: mo_set

!   *** Local variables ***

    INTEGER :: handle,imo

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_density_matrix","I","",handle)

    CALL blacs_syrk("U","N",mo_set%lfomo-1,2.0_wp,mo_set%eigenvectors,1,1,&
                    0.0_wp,density_matrix,globenv)

!   *** Add the density contribution from the fractional occupied MOs ***

    DO imo=mo_set%lfomo,mo_set%homo
      CALL blacs_syrk("U","N",1,mo_set%occupation_numbers(imo),&
                      mo_set%eigenvectors,1,imo,1.0_wp,density_matrix,globenv)
    END DO

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE calculate_density_matrix_blacs

! *****************************************************************************

  SUBROUTINE calculate_density_matrix_sparse(mo_set,density_matrix,globenv)

!   Purpose: Calculate the density matrix 

!   as yet only for fully occupied orbitals !!!!

!   History: - Creation (july 2002 Joost VandeVondele)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER           :: density_matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_type), POINTER                :: mo_set

!   *** Local variables ***

    INTEGER :: handle,imo

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_density_matrix","I","",handle)

    CALL set_matrix(density_matrix,0.0_wp)

    CALL sparse_plus_blacs_blacst(sparse_matrix=density_matrix,&
                                  matrix_v=mo_set%eigenvectors,&
                                  first_col=1,&
                                  ncol=mo_set%lfomo-1,&
                                  alpha=2.0_wp,&
                                  globenv=globenv)

    DO imo=mo_set%lfomo,mo_set%homo
      CALL sparse_plus_blacs_blacst(sparse_matrix=density_matrix,&
                                    matrix_v=mo_set%eigenvectors,&
                                    first_col=imo,&
                                    ncol=imo,&
                                    alpha=mo_set%occupation_numbers(imo),&
                                    globenv=globenv)
    END DO

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE calculate_density_matrix_sparse

! *****************************************************************************

  SUBROUTINE calculate_w_matrix_sparse(mo_set,ks_matrix,w_matrix,globenv)

!   Purpose: Calculate the W matrix from the MO eigenvectors, MO eigenvalues,
!            and the MO occupation numbers.
!   History: - Creation (17.07.2002,MK)
!   ***************************************************************************

    TYPE(real_matrix_type), POINTER           :: w_matrix,ks_matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_type), POINTER                :: mo_set
    TYPE(blacs_matrix_type), POINTER          :: weighted_vectors,h_block

!   *** Local variables ***

    INTEGER  :: handle,imo,ncol_global,nrow_global,nrow_block,ncol_block

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_w_matrix","I","",handle)

    NULLIFY(weighted_vectors,h_block)
    CALL get_blacs_matrix_info(matrix=mo_set%eigenvectors, &
                                 ncol_global=ncol_global, &
                                 nrow_global=nrow_global, &
                                 nrow_block=nrow_block, &
                                 ncol_block=ncol_block,globenv=globenv)

    CALL allocate_blacs_matrix(new_matrix=weighted_vectors,&
                                 nrow_global=nrow_global, &
                                 ncol_global=ncol_global, &
                                 nrow_block=nrow_block, &
                                 ncol_block=ncol_block, &
                                 name="h block", & 
                                 globenv=globenv)

    CALL allocate_blacs_matrix(new_matrix=h_block,&
                                 nrow_global=ncol_global, &
                                 ncol_global=ncol_global, &
                                 nrow_block=nrow_block, &
                                 ncol_block=ncol_block, &
                                 name="h block", &
                                 globenv=globenv)

    !CALL blacs_scale_matrix(weighted_vectors,mo_set%eigenvalues, &
    !                                                    globenv,mo_set%homo)

    if (mo_set%homo .ne. ncol_global) CALL stop_program("calculate_w_matrix","NYI")

    CALL sparse_times_blacs(ks_matrix,mo_set%eigenvectors,weighted_vectors, &
                                                  mo_set%homo,globenv)
    CALL blacs_gemm('T','N',ncol_global,ncol_global,nrow_global,1.0_wp, &
                    mo_set%eigenvectors,weighted_vectors,0.0_wp,h_block,globenv)
    CALL blacs_gemm('N','N',nrow_global,ncol_global,ncol_global,1.0_wp, &
                    mo_set%eigenvectors,h_block,0.0_wp,weighted_vectors,globenv)

    ! scale every column with the occupation
    CALL blacs_scale_matrix(weighted_vectors,mo_set%occupation_numbers, &
                                                        globenv,mo_set%homo)

    CALL set_matrix(w_matrix,0.0_wp)
    CALL sparse_plus_blacs_blacst(sparse_matrix=w_matrix,&
                                  matrix_v=mo_set%eigenvectors,&
                                  matrix_g=weighted_vectors,&
                                  first_col=1,&
                                  ncol=mo_set%homo,&
                                  globenv=globenv)

    CALL deallocate_blacs_matrix(weighted_vectors)
    CALL deallocate_blacs_matrix(h_block)
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE calculate_w_matrix_sparse

! *****************************************************************************

  SUBROUTINE correct_mo_eigenvalues(mo_set,level_shift)

!   Purpose: Correct MO eigenvalues after MO level shifting.

!   History: - Creation (19.04.2002, MK)

!   ***************************************************************************

    TYPE(mo_set_type), POINTER :: mo_set
    REAL(wp), INTENT(IN)       :: level_shift

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE correct_mo_eigenvalues (MODULE qs_mo_types)"

!   *** Local variables ***

    INTEGER :: imo

!   ---------------------------------------------------------------------------

    IF (level_shift == 0.0_wp) RETURN

    DO imo=mo_set%homo+1,mo_set%nmo
      mo_set%eigenvalues(imo) = mo_set%eigenvalues(imo) - level_shift
    END DO

  END SUBROUTINE correct_mo_eigenvalues

! *****************************************************************************

  SUBROUTINE deallocate_mo_set(mo_set)

!   Purpose: Deallocate a wavefunction data structure.

!   History: - Creation (15.05.2001,MK)

!   ***************************************************************************

    TYPE(mo_set_type), POINTER :: mo_set

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE deallocate_mo (MODULE qs_mo_types)"

!   *** Local variables ***

    INTEGER :: istat

!   ---------------------------------------------------------------------------

    DEALLOCATE (mo_set%eigenvalues,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"mo_set%eigenvalues")

    DEALLOCATE (mo_set%occupation_numbers,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"mo_set%occupation_numbers")

    CALL deallocate_blacs_matrix(mo_set%eigenvectors)

    DEALLOCATE (mo_set,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"mo_set")

  END SUBROUTINE deallocate_mo_set

! *****************************************************************************

  SUBROUTINE get_mo_set(mo_set,maxocc,homo,lfomo,nao,nelectron,nmo,&
                        eigenvalues,occupation_numbers,eigenvectors)

!   Purpose: Get the components of a MO set data structure.

!   History: - Creation (22.04.2002, MK)

!   ***************************************************************************

    TYPE(mo_set_type), POINTER                 :: mo_set
    REAL(wp), OPTIONAL, INTENT(OUT)            :: maxocc
    INTEGER, OPTIONAL, INTENT(OUT)             :: homo,lfomo,nao,nelectron,nmo
    REAL(wp), OPTIONAL, DIMENSION(:), POINTER  :: eigenvalues,&
                                                  occupation_numbers
    TYPE(blacs_matrix_type), OPTIONAL, POINTER :: eigenvectors

!   *** Local variables ***

    INTEGER :: istat

!   ---------------------------------------------------------------------------

    IF (PRESENT(maxocc)) maxocc = mo_set%maxocc
    IF (PRESENT(homo)) homo = mo_set%homo
    IF (PRESENT(lfomo)) lfomo = mo_set%lfomo
    IF (PRESENT(nao)) nao = mo_set%nao
    IF (PRESENT(nelectron)) nelectron = mo_set%nelectron
    IF (PRESENT(nmo)) nmo = mo_set%nmo
    IF (PRESENT(eigenvalues)) eigenvalues => mo_set%eigenvalues
    IF (PRESENT(occupation_numbers)) THEN
      occupation_numbers => mo_set%occupation_numbers
    END IF
    IF (PRESENT(eigenvectors)) eigenvectors => mo_set%eigenvectors

  END SUBROUTINE get_mo_set

! *****************************************************************************

  SUBROUTINE read_mo_set_old(mo_set,globenv)

!   Purpose: Read the MO eigenvectors from the restart file.

!   History: - Creation (15.05.2001,MK)
!            - Parallel input (19.05.2001, MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_type), POINTER                :: mo_set

!   *** Local variables ***

    CHARACTER(LEN=6)   :: extension
    CHARACTER(LEN=200) :: file_name
    INTEGER            :: mype,npe,restart_unit

!   ---------------------------------------------------------------------------

    CALL get_blacs_info(globenv=globenv,&
                        my_process_number=mype,&
                        number_of_processes=npe)

    IF (npe > 1) THEN
      WRITE (UNIT=extension,FMT="(I6)") mype
      file_name = TRIM(globenv%restart_file_name)//"."//ADJUSTL(extension)
    ELSE
      file_name = globenv%restart_file_name
    END IF

    CALL open_file(file_name=file_name,&
                   file_action="READ",&
                   file_form="UNFORMATTED",&
                   file_status="OLD",&
                   unit_number=restart_unit)

    CALL read_blacs_matrix(matrix=mo_set%eigenvectors,&
                           lunit=restart_unit,&
                           globenv=globenv)

    CALL close_file(unit_number=restart_unit)

  END SUBROUTINE read_mo_set_old

! *****************************************************************************

  SUBROUTINE set_mo_occupation(mo_set,smear)

!   Purpose: Smearing of the MO occupation due to the given smear value.

!   History: - Creation (17.04.2002, MK)

!   ***************************************************************************

    TYPE(mo_set_type), POINTER :: mo_set
    REAL(wp), INTENT(IN)       :: smear

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE set_mo_occupation (MODULE qs_mo_types)"

!   *** Local variables ***

    REAL(wp) :: e1,e2,edelta,edist,nelec
    INTEGER  :: imo,nomo

!   ---------------------------------------------------------------------------

    nomo = mo_set%nelectron/2

    mo_set%occupation_numbers(:) = 0.0_wp
    mo_set%occupation_numbers(1:nomo) = mo_set%maxocc

    mo_set%homo = nomo
    mo_set%lfomo = nomo + 1

    IF (smear /= 0.0_wp) THEN

!     *** Define the energy window for the eigenvalues ***

      e1 = mo_set%eigenvalues(nomo) - 0.5_wp*smear
      e2 = mo_set%eigenvalues(nomo) + 0.5_wp*smear

!     *** Find the lowest fractional occupied MO ***

      DO imo=1,nomo
        IF (mo_set%eigenvalues(imo) > e1) THEN
          mo_set%lfomo = imo
          EXIT
        END IF
      END DO

!     *** Find the highest fractional occupied MO which will be the HOMO ***

      DO imo=mo_set%nmo,nomo,-1
        IF (mo_set%eigenvalues(imo) < e2) THEN
          mo_set%homo = imo
          EXIT
        END IF
      END DO

!     *** Get the number of electrons to be smeared ***

      edist = 0.0_wp
      nelec = 0.0_wp

      DO imo=mo_set%lfomo,mo_set%homo
        nelec = nelec + mo_set%occupation_numbers(imo)
        edist = edist + ABS(e2 - mo_set%eigenvalues(imo))
      END DO

!     *** Smear electrons inside the window ***

      DO imo=mo_set%lfomo,mo_set%homo
        edelta = ABS(e2 - mo_set%eigenvalues(imo))
        mo_set%occupation_numbers(imo) = MIN(mo_set%maxocc,nelec*edelta/edist)
        nelec = nelec - mo_set%occupation_numbers(imo)
        edist = edist - edelta
      END DO

!     *** Check, if the smearing involves more than one MO ***

      IF (mo_set%lfomo == mo_set%homo) THEN
        mo_set%homo = nomo
        mo_set%lfomo = nomo + 1
      END IF

    END IF

  END SUBROUTINE set_mo_occupation

! *****************************************************************************

  SUBROUTINE write_mo_set_to_restart_unit(mo_set,globenv)

!   Purpose: Write the MO eigenvectors to the restart file.

!   History: - Creation (15.05.2001,MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_type), POINTER                :: mo_set

!   *** Local variables ***

    CHARACTER(LEN=6)   :: extension
    CHARACTER(LEN=200) :: file_name
    INTEGER            :: mype,npe,restart_unit

!   ---------------------------------------------------------------------------

    CALL get_blacs_info(globenv=globenv,&
                        my_process_number=mype,&
                        number_of_processes=npe)

    IF (npe > 1) THEN
      WRITE (UNIT=extension,FMT="(I6)") mype
      file_name = TRIM(globenv%restart_file_name)//"."//ADJUSTL(extension)
    ELSE
      file_name = globenv%restart_file_name
    END IF

    CALL open_file(file_name=file_name,&
                   file_action="WRITE",&
                   file_form="UNFORMATTED",&
                   file_status="REPLACE",&
                   unit_number=restart_unit)

    CALL write_blacs_matrix(matrix=mo_set%eigenvectors,&
                            lunit=restart_unit,&
                            globenv=globenv)

    CALL close_file(unit_number=restart_unit)

  END SUBROUTINE write_mo_set_to_restart_unit

! *****************************************************************************
  SUBROUTINE write_mo_set_to_restart_new(mo_set,atomic_kind_set,particle_set,globenv)

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(global_environment_type), INTENT(IN)     :: globenv
    TYPE(mo_set_type), POINTER                    :: mo_set
    TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set

!   *** Local variables ***
    CHARACTER(LEN=6)   :: extension
    CHARACTER(LEN=200) :: file_name
    INTEGER            :: mype,npe,restart_unit
    TYPE(gto_basis_set_type), POINTER :: orb_basis_set
    INTEGER           :: iatom,icol,irow,iset,isgf,iso,ishell,istat,jcol,&
                         left,lshell,natom,ncol,nset,nsgf,output_unit,right, &
                         nao,nmo,i,nset_max,nshell_max
                         
    INTEGER, DIMENSION(:), POINTER          :: nshell
    REAL(wp), DIMENSION(:,:), POINTER       :: matrix
    INTEGER, DIMENSION(:,:), POINTER        :: l

    INTEGER, DIMENSION(:), POINTER :: nset_info
    INTEGER, DIMENSION(:,:), POINTER :: nshell_info
    INTEGER, DIMENSION(:,:,:), POINTER :: nso_info

    REAL(wp) , DIMENSION(:,:), POINTER :: vecbuffer

    nmo=mo_set%nmo
    nao=mo_set%nao

    if (globenv%ionode) then

      file_name = globenv%restart_file_name
      CALL open_file(file_name=file_name,&
                     file_action="WRITE",&
                     file_form="UNFORMATTED",&
                     file_status="REPLACE",&
                     unit_number=restart_unit)


      ! create some info about the basis set first

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               natom=natom,&
                               nsgf=nsgf)


      nset_max=0
      nshell_max=0

      DO iatom=1,natom
         CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                              orb_basis_set=orb_basis_set)
         CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                nset=nset,&
                                nshell=nshell,&
                                l=l)
  
         nset_max=MAX(nset_max,nset)
         DO iset=1,nset
           nshell_max=MAX(nshell_max,nshell(iset))
         END DO
      END DO


      allocate(nso_info   (nshell_max,nset_max,natom))
      allocate(nshell_info(nset_max,natom))
      allocate(nset_info  (natom))

      nset_info(:)=0
      nshell_info(:,:)=0
      nso_info(:,:,:)=0

      DO iatom=1,natom
         CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                              orb_basis_set=orb_basis_set)
         CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                nset=nset,&
                                nshell=nshell,&
                                l=l)
         nset_info(iatom)=nset
         DO iset=1,nset
           nshell_info(iset,iatom)=nshell(iset)
           DO ishell=1,nshell(iset)
             lshell = l(ishell,iset)
             nso_info(ishell,iset,iatom)= nso(lshell)
           ENDDO
         END DO
      END DO

      ! write the bs information

      write(restart_unit) natom,nmo,nao,nset_max,nshell_max
      write(restart_unit) nset_info
      write(restart_unit) nshell_info
      write(restart_unit) nso_info

      deallocate(nset_info,nshell_info,nso_info)

    endif


    allocate(vecbuffer(1,nao))
    do i=1,nmo
       call blacs_replicated_copy(mo_set%eigenvectors,vecbuffer,i,1,1,"B2R",globenv)
       if (globenv%ionode) write(restart_unit) vecbuffer
    enddo
    deallocate(vecbuffer)

    if (globenv%ionode) CALL close_file(unit_number=restart_unit)

  END SUBROUTINE write_mo_set_to_restart_new

! *****************************************************************************

! *****************************************************************************
  SUBROUTINE read_mo_set_from_restart_new(mo_set,atomic_kind_set,particle_set,globenv)

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(global_environment_type), INTENT(IN)     :: globenv
    TYPE(mo_set_type), POINTER                    :: mo_set
    TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set

!   *** Local variables ***
    CHARACTER(LEN=6)   :: extension
    CHARACTER(LEN=200) :: file_name
    INTEGER            :: mype,npe,restart_unit
    TYPE(gto_basis_set_type), POINTER :: orb_basis_set
    INTEGER           :: iatom,icol,irow,iset,isgf,iso,ishell,istat,jcol,&
                         left,lshell,natom,ncol,nset,nsgf,output_unit,right, &
                         nao,nmo,i,nset_max,nshell_max,natom_read,nao_read,nmo_read
                         
    INTEGER, DIMENSION(:), POINTER          :: nshell
    REAL(wp), DIMENSION(:,:), POINTER       :: matrix
    INTEGER, DIMENSION(:,:), POINTER        :: l

    INTEGER, DIMENSION(:), POINTER :: nset_info
    INTEGER, DIMENSION(:,:), POINTER :: nshell_info
    INTEGER, DIMENSION(:,:,:), POINTER :: nso_info,offset_info

    REAL(wp) , DIMENSION(:,:), POINTER :: vecbuffer,vecbuffer_read

    integer :: irow_read,nset_read,ishell_read,iset_read,offset_read
    logical :: use_this,did_change

    nmo=mo_set%nmo
    nao=mo_set%nao

    if (globenv%ionode) then

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               natom=natom,&
                               nsgf=nsgf)

      file_name = globenv%restart_file_name
      CALL open_file(file_name=file_name,&
                     file_action="READ",&
                     file_form="UNFORMATTED",&
                     file_status="OLD",&
                     unit_number=restart_unit)

      read(restart_unit) natom_read,nmo_read,nao_read,nset_max,nshell_max

      if (nmo_read.ne.nmo) &
         call stop_program("read_mo_set_from_restart","incorrect nmo")
      if (natom_read.ne.natom) &
         call stop_program("read_mo_set_from_restart","incorrect natom")

      allocate(nso_info   (nshell_max,nset_max,natom_read))
      allocate(offset_info   (nshell_max,nset_max,natom_read))
      allocate(nshell_info(nset_max,natom_read))
      allocate(nset_info  (natom_read))

      if (nao_read.ne.nao) then
         write(globenv%scr,*) "READ RESTART : WARNING : DIFFERENT # AOs ",nao,nao_read
         write(globenv%scr,*) "writing out mapping"
         did_change=.true.
      endif

      read(restart_unit) nset_info
      read(restart_unit) nshell_info
      read(restart_unit) nso_info

      i=1
      do iatom=1,natom
       do iset=1,nset_info(iatom)
        do ishell=1,nshell_info(iset,iatom)
           offset_info(ishell,iset,iatom)=i
           i=i+nso_info(ishell,iset,iatom)
        enddo
       enddo
      enddo

      allocate(vecbuffer_read(1,nao_read))

    endif


    allocate(vecbuffer(1,nao))
    do i=1,nmo
       if (globenv%ionode) then 
            read(restart_unit) vecbuffer_read
            ! now, try to assign the read to the real vector
            ! in case the basis set changed this involves some guessing
            irow=1
            DO iatom=1,natom
               CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                                    orb_basis_set=orb_basis_set)
               CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                      nset=nset,&
                                      nshell=nshell,&
                                      l=l)

               use_this=.true.
               iset_read=1
               DO iset=1,nset
                  ishell_read=1
                  DO ishell=1,nshell(iset)
                     lshell = l(ishell,iset)
                     if (iset_read.gt.nset_info(iatom)) use_this=.false.
                     if (use_this .and. nso(lshell).eq.nso_info(ishell_read,iset_read,iatom)) then
                         offset_read=offset_info(ishell_read,iset_read,iatom)
                         ishell_read=ishell_read+1
                         if (ishell_read.gt.nshell_info(iset,iatom)) then
                            ishell_read=1
                            iset_read=iset_read+1
                         endif
                     else
                         use_this=.false.
                     endif
                     DO iso=1,nso(lshell)
                        if (use_this) then
                            if (offset_read-1+iso.lt.1 .or. offset_read-1+iso.gt.nao_read) then
                               ! this indicates a bug
                               if (did_change) write(globenv%scr,*) "-?-",irow,offset_read-1+iso
                               vecbuffer(1,irow)=0.0_wp
                            else
                               vecbuffer(1,irow)=vecbuffer_read(1,offset_read-1+iso)
                               if (did_change) write(globenv%scr,*) "   ",irow,offset_read-1+iso
                            endif
                        else
                            vecbuffer(1,irow)=0.0_wp
                        endif
                        irow=irow+1
                     ENDDO
                     use_this=.true.
                  ENDDO
               ENDDO
            ENDDO

       else
            vecbuffer(1,:)=0.0_wp
       endif
       CALL mp_sum(vecbuffer,globenv%group)  ! the usual ugly way
       call blacs_replicated_copy(mo_set%eigenvectors,vecbuffer,i,1,1,"R2B",globenv)
       did_change=.false.
    enddo
    deallocate(vecbuffer)

    if (globenv%ionode) then
       deallocate(nset_info,nshell_info,nso_info,vecbuffer_read,offset_info)
       CALL close_file(unit_number=restart_unit)
    endif

  END SUBROUTINE read_mo_set_from_restart_new

! *****************************************************************************

  SUBROUTINE write_mo_set_to_output_unit(mo_set,atomic_kind_set,particle_set,&
                                         before,after,globenv)

!   Purpose: Write the MO eigenvalues, MO occupation numbers and
!            MO eigenvectors.

!   History: - Creation (15.05.2001,MK)

!   ***************************************************************************

!   List of variables:

!   after : Number of digits after point.
!   before: Number of digits before point.

!   ---------------------------------------------------------------------------

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(global_environment_type), INTENT(IN)     :: globenv
    TYPE(mo_set_type), POINTER                    :: mo_set
    TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set
    INTEGER, INTENT(IN)                           :: after,before

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE write_mo_set_to_output_unit (MODULE qs_mo_types)"

!   *** Local variables ***

    TYPE(gto_basis_set_type), POINTER :: orb_basis_set

    CHARACTER(LEN=60) :: name
    CHARACTER(LEN=25) :: fmtstr1
    CHARACTER(LEN=22) :: fmtstr2
    CHARACTER(LEN=35) :: fmtstr3
    CHARACTER(LEN=2)  :: element_symbol
    INTEGER           :: from,iatom,icol,irow,iset,isgf,iso,ishell,istat,jcol,&
                         left,lshell,natom,ncol,nset,nsgf,output_unit,right,&
                         to,width

    CHARACTER(LEN=6), DIMENSION(:), POINTER :: sgf_symbol
    INTEGER, DIMENSION(:), POINTER          :: nshell
    REAL(wp), DIMENSION(:,:), POINTER       :: matrix
    INTEGER, DIMENSION(:,:), POINTER        :: l

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%print%mo_eigenvalues) RETURN

    output_unit = globenv%scr

    IF (globenv%print%mo_eigenvectors) THEN
      NULLIFY (matrix)
      CALL copy_blacs_to_full_matrix(mo_set%eigenvectors,matrix,globenv)
    END IF

    IF (.NOT.globenv%ionode) RETURN

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             natom=natom,&
                             nsgf=nsgf)

!   *** Definition of the variable formats ***

    fmtstr1 = "(/,T2,23X,  (  X,I5,  X))"
    fmtstr2 = "(T2,23X,  (1X,F  .  ))"
    fmtstr3 = "(T2,2I5,2X,A2,1X,A8,   (1X,F  .  ))"

!   *** Write headline ***

    IF (globenv%print%mo_eigenvectors) THEN
      name = "MO EIGENVALUES, MO OCCUPATION NUMBERS AND MO EIGENVECTORS"
    ELSE IF (globenv%print%mo_occupation_numbers) THEN
      name = "MO EIGENVALUES AND MO OCCUPATION NUMBERS"
    ELSE
      name = "MO EIGENVALUES"
    END IF

    WRITE (UNIT=output_unit,FMT="(/,/,T2,A)") TRIM(name)

!   *** Write the variable format strings ***

    width = before + after + 3
    ncol = INT(56/width)

    right = MAX((after-2),1)
    left =  width - right - 5

    WRITE (UNIT=fmtstr1(11:12),FMT="(I2)") ncol
    WRITE (UNIT=fmtstr1(14:15),FMT="(I2)") left
    WRITE (UNIT=fmtstr1(21:22),FMT="(I2)") right

    WRITE (UNIT=fmtstr2(9:10),FMT="(I2)") ncol
    WRITE (UNIT=fmtstr2(16:17),FMT="(I2)") width - 1
    WRITE (UNIT=fmtstr2(19:20),FMT="(I2)") after

    WRITE (UNIT=fmtstr3(22:23),FMT="(I2)") ncol
    WRITE (UNIT=fmtstr3(29:30),FMT="(I2)") width - 1
    WRITE (UNIT=fmtstr3(32:33),FMT="(I2)") after

!   *** Write the matrix in the selected format ***

    DO icol=1,nsgf,ncol

      from = icol
      to = MIN((from+ncol-1),nsgf)

      WRITE (UNIT=output_unit,FMT=fmtstr1)&
        (jcol,jcol=from,to)
      WRITE (UNIT=output_unit,FMT=fmtstr2)&
        (mo_set%eigenvalues(jcol),jcol=from,to)
      WRITE (UNIT=output_unit,FMT="(A)") ""

      IF (globenv%print%mo_occupation_numbers) THEN
        WRITE (UNIT=output_unit,FMT=fmtstr2)&
          (mo_set%occupation_numbers(jcol),jcol=from,to)
        WRITE (UNIT=output_unit,FMT="(A)") ""
      END IF

      IF (globenv%print%mo_eigenvectors) THEN

        irow = 1

        DO iatom=1,natom

          IF (iatom /= 1) WRITE (UNIT=output_unit,FMT="(A)") ""

          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                               element_symbol=element_symbol,&
                               orb_basis_set=orb_basis_set)
          CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                 nset=nset,&
                                 nshell=nshell,&
                                 l=l,&
                                 sgf_symbol=sgf_symbol)

          isgf = 1

          DO iset=1,nset
            DO ishell=1,nshell(iset)
              lshell = l(ishell,iset)
              DO iso=1,nso(lshell)
                WRITE (UNIT=output_unit,FMT=fmtstr3)&
                  irow,iatom,element_symbol,sgf_symbol(isgf),&
                  (matrix(irow,jcol),jcol=from,to)
                isgf = isgf + 1
                irow = irow + 1
              END DO
            END DO
          END DO

        END DO

      END IF

    END DO

    WRITE (UNIT=output_unit,FMT="(/)")

!   *** Release work storage ***

    IF (globenv%print%mo_eigenvectors) THEN
      DEALLOCATE (matrix,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"matrix")
    END IF

  END SUBROUTINE write_mo_set_to_output_unit

! *****************************************************************************

END MODULE qs_mo_types
