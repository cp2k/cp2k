!------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations          !
!   Copyright (C) 2000  CP2K developers group                                  !
!------------------------------------------------------------------------------!
!!****** cp2k/tbmd_input [1.0] *
!!
!!   NAME
!!     tbmd_input
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
    MODULE tbmd_input

  USE global_types,                    ONLY: global_environment_type
  USE input_types,                     ONLY: setup_parameters_type
  USE kinds,                           ONLY: dbl
  USE parser,                          ONLY: get_next,&
                                             parser_end,&
                                             parser_init,&
                                             read_line,&
                                             search_label,&
                                             stop_parser,&
                                             test_next
  USE string_utilities,                ONLY: str_search,&
                                             uppercase,&
                                             xstring
  USE tbmd_types,                      ONLY: tb_hopping_type,&
                                             tbatom_param_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

      IMPLICIT NONE
!
      PRIVATE
      PUBLIC :: read_tb_hamiltonian, read_tb_hopping_elements
!
!!*****
!------------------------------------------------------------------------------!
!
    CONTAINS
!
!------------------------------------------------------------------------------!
!!****** tbmd_input/read_tb_hamiltonian [1.0] *
!!
!!   NAME
!!     read_tb_hamiltonian
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE read_tb_hamiltonian ( setup, atom_names, tbatom, tbmdpar )


    TYPE(setup_parameters_type), &
      INTENT(INOUT)                          :: setup
    CHARACTER(LEN=*), INTENT(IN)             :: atom_names( : )
    TYPE(tbatom_param_type), DIMENSION(:), &
      POINTER                                :: tbatom
    TYPE(global_environment_type), &
      INTENT(IN)                             :: tbmdpar

    CHARACTER(len=1), PARAMETER :: lq(0:3) = (/ 'S', 'P', 'D', 'F'/)

    CHARACTER(len=10)                        :: basis
    CHARACTER(len=20)                        :: at, string
    CHARACTER(len=5)                         :: label
    INTEGER                                  :: i, i1, iat, ierror, ilen, &
                                                ios, iw, j, k, nt

  iat = 0
  nt = setup%natom_type
  ALLOCATE (tbatom(nt),STAT=ios)
  IF (ios/=0) CALL stop_memory('read_tb_hamiltonian','tbatom',nt)
!..parse the input section &HAMILTONIAN
  label = '&HAMILTONIAN'
  CALL parser_init(setup%set_file_name,tbmdpar)
  CALL search_label(label,ierror,ignore_case=.TRUE.)
  IF (ierror/=0) THEN
    IF (tbmdpar%ionode) &
      WRITE (tbmdpar%scr,'(a)') &
           ' read_tb_hamiltonian| No input section &HAMILTONIAN found '
  ELSE
    CALL read_line
    DO WHILE (test_next()/='X')
      ilen = 6
      CALL get_next(string,ilen)
      CALL uppercase ( string )
      SELECT CASE (string)
      CASE ('TBATOM')
        iat = iat + 1
        IF (iat>setup%natom_type) THEN
          CALL stop_parser('read_tb_hamiltonian','Too many atoms')
        END IF
        ilen = 20
        CALL get_next(at,ilen)
        i1 = str_search(atom_names,nt,at)
        IF (i1==0) CALL stop_parser('read_tb_hamiltonian','atom not found')
        tbatom(i1) %atom_name = at
        CALL get_next ( tbatom(i1) %atomic_number )
        CALL get_next ( tbatom(i1) %valence_charge )
        ilen = 10
        CALL get_next(basis,ilen)
        CALL set_basis(basis,tbatom(i1)%bsf)
        DO i = 0, 3
          DO j = 1, tbatom(iat) %bsf(i)
            CALL get_next ( tbatom(iat) %orbital_energy(i,j) )
          END DO
        END DO
      END SELECT
      CALL read_line
    END DO
  END IF
  CALL parser_end
!..write some information to output
  IF ( tbmdpar % ionode) THEN
    iw = tbmdpar % scr
    IF (tbmdpar%print_level>=0) THEN
      WRITE (iw,'(A,T71,I10)') ' TBMD| Number of atom types ', &
           setup%natom_type
      WRITE (iw,'(A,T40,A,T55,A,T65,A,T75,A)') ' TBMD|', ' Atom type', &
           ' Atomic# ', ' Charge ', ' Basis'
      DO i = 1, setup%natom_type
        WRITE (iw,'(A,T40,A,T57,I6,T67,F5.1,T71,A)') ' TBMD| ', &
             adjustr(tbatom(i)%atom_name), tbatom(i) %atomic_number, &
             tbatom(i) %valence_charge, adjustr(basis)
      END DO
      WRITE (iw,'(A,T40,A,T55,A,T66,A)') ' TBMD|', ' Atom type', &
           ' Orbital ', ' Orbital energy'
      DO i = 1, setup%natom_type
        DO j = 0, 3
          DO k = 1, tbatom(i) %bsf(j)
            WRITE (iw,'(A,T40,A,T62,A,T68,F10.6,A)') ' TBMD| ', &
                 adjustr(tbatom(i)%atom_name), lq(j), &
                 tbatom(i) %orbital_energy(j,k), ' eV'
          END DO
        END DO
      END DO
      WRITE (iw,*)
    END IF
  END IF

END SUBROUTINE read_tb_hamiltonian

!------------------------------------------------------------------------------!

SUBROUTINE set_basis(basis,bsf)
    CHARACTER(len=*), INTENT(INOUT)          :: basis
    INTEGER, INTENT(OUT)                     :: bsf(0:3)

    INTEGER                                  :: i, ia, ie, lmax, lnew

  lmax = -1
  bsf = 0
  CALL xstring(basis,ia,ie)
  CALL uppercase ( basis(ia:ie) )
  DO i = ia, ie
    SELECT CASE (basis(i:i))
    CASE DEFAULT
      CALL stop_parser('set_basis','angular momentum label not valid')
    CASE ('S')
      bsf(0) = bsf(0) + 1
      lnew = 0
    CASE ('P')
      bsf(1) = bsf(1) + 1
      lnew = 1
    CASE ('D')
      bsf(2) = bsf(2) + 1
      lnew = 2
    CASE ('F')
      bsf(3) = bsf(3) + 1
      lnew = 3
    END SELECT
    IF ( lnew < lmax ) THEN
      CALL stop_parser('set_basis','shells have to be ordered')
    ELSE
      lmax = lnew
    END IF
  END DO
END SUBROUTINE set_basis

!!*****
!------------------------------------------------------------------------------!
!!****** tbmd_input/read_tb_hopping_elements [1.0] *
!!
!!   NAME
!!     read_tb_hopping_elements
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE read_tb_hopping_elements ( setup, atom_names, tbatom, tbhop, tbmdpar )


    TYPE(setup_parameters_type), &
      INTENT(INOUT)                          :: setup
    CHARACTER(LEN=*), INTENT(IN)             :: atom_names( : )
    TYPE(tbatom_param_type), DIMENSION(:)    :: tbatom
    TYPE(tb_hopping_type), DIMENSION(:, :), &
      POINTER                                :: tbhop
    TYPE(global_environment_type), &
      INTENT(INOUT)                          :: tbmdpar

    CHARACTER(len=10)                        :: at1, at2, itype
    CHARACTER(len=20)                        :: string
    CHARACTER(len=5)                         :: label
    INTEGER                                  :: i, iat, ierror, ilen, ios, &
                                                iw, j, k, k1, k2, kk, ku, l, &
                                                l1, l2, lu, n, nk
    REAL(dbl)                                :: rc

  iat = setup%natom_type
  ALLOCATE (tbhop(iat,iat),STAT=ios)
  IF (ios/=0) CALL stop_memory('read_tb_hopping_elements', &
          'allocation of array tbhop failed',iat*iat)
!..calculate number of non-zero integrals
  DO i = 1, iat
    DO j = 1, iat
      n = 0
      IF (i==j) THEN
        DO k = 0, 3
          DO l = k, 3
            n = n + tbatom(i) %bsf(k)*tbatom(j) %bsf(l)*(k+1)
          END DO
        END DO
      ELSE
        DO k = 0, 3
          DO l = 0, 3
            n = n + tbatom(i) %bsf(k)*tbatom(j) %bsf(l)*(min(k,l)+1)
          END DO
        END DO
      END IF
      tbhop(i,j) %nint = n
      ALLOCATE (tbhop(i,j)%lls(3,tbhop(i,j)%nint),STAT=ios)
      IF (ios/=0) CALL stop_memory('read_tb_hopping_elements', &
        'allocation of array tbhop%lls failed',0)
      tbhop(i,j) %rcutsq = 0._dbl
    END DO
  END DO
  iw = tbmdpar%scr
!..parse the input section &HAMILTONIAN
  label = '&HAMILTONIAN'
  CALL parser_init(setup%set_file_name,tbmdpar)
  CALL search_label(label,ierror,ignore_case=.TRUE.)
  IF (ierror/=0) THEN
    IF (tbmdpar%ionode) WRITE (iw,'(a)') &
      ' input_inter| No input section &HAMILTONIAN found '
  ELSE
    CALL read_line
    DO WHILE (test_next()/='X')
      ilen = 5
      CALL get_next(string,ilen)
      CALL uppercase ( string )
      SELECT CASE (string)
      CASE ('TBHAM')
        ilen = 10
        CALL get_next(itype,ilen)
        CALL uppercase ( itype )
        IF (itype=='GOODWIN') THEN
          CALL get_next(at1,ilen)
          CALL get_next(at2,ilen)
          i = str_search(atom_names,iat,at1)
          j = str_search(atom_names,iat,at2)
          IF (i==0 .OR. j==0) CALL stop_parser('read_tb_hopping_elements', &
            'atom type not known')
          tbhop(i,j) %htype = itype
          ALLOCATE (tbhop(i,j)%gsp(tbhop(i,j)%nint),STAT=ios)
          IF (ios/=0) CALL stop_memory('read_tb_hopping_elements', &
                    'allocation of array tbhop%gsp failed')
          DO k = 1, tbhop(i,j) %nint
            CALL read_line
            CALL get_next(itype,ilen)
            CALL set_int(itype,tbhop(i,j)%lls(1:3,k))
            CALL get_next ( tbhop(i,j) %gsp(k) %er0 )
            CALL get_next ( tbhop(i,j) %gsp(k) %r0 )
            CALL get_next ( tbhop(i,j) %gsp(k) %rt )
            CALL get_next ( tbhop(i,j) %gsp(k) %na )
            CALL get_next ( tbhop(i,j) %gsp(k) %nb )
            tbhop(j,i) %gsp(k) = tbhop(i,j) %gsp(k)
          END DO
        ELSE IF (itype=='CUTOFF') THEN
          ilen = 10
          CALL get_next(at1,ilen)
          CALL get_next(at2,ilen)
          CALL get_next ( rc )
          i = str_search(atom_names,iat,at1)
          j = str_search(atom_names,iat,at2)
          IF (i==0 .OR. j==0) CALL stop_parser('read_tb_hopping_elements', &
            'atom type not known')
          tbhop(i,j) %rcutsq = rc*rc
        ELSE
          CALL stop_parser('read_tb_hopping_elements','interaction type not known')
        END IF
      END SELECT
      CALL read_line
    END DO
  END IF
  CALL parser_end
!..check the ordering of the hamiltonian matrix eleements
  DO i = 1, setup%natom_type
    DO j = i, setup%natom_type
      nk = 0
      DO l1 = 0, 3
        DO k1 = 1, tbatom(i) %bsf(l1)
          lu = 0
          IF (i==j) lu = l1
          DO l2 = lu, 3
            ku = 1
            IF (i==j .AND. l1==l2) ku = k1
            DO k2 = ku, tbatom(j) %bsf(l2)
              DO kk = 0, min(l1,l2)
                nk = nk + 1
                IF (tbhop(i,j)%lls(1,nk)/=l1) &
                  CALL stop_program('read_tb_hopping_elements', &
                  'ordering of H matrix elements wrong')
                IF (tbhop(i,j)%lls(2,nk)/=l2) &
                  CALL stop_program('read_tb_hopping_elements', &
                  'ordering of H matrix elements wrong')
                IF (tbhop(i,j)%lls(3,nk)/=kk) &
                  CALL stop_program('read_tb_hopping_elements', &
                  'ordering of H matrix elements wrong')
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
  END DO
!..write some information to output
  IF (tbmdpar%ionode) THEN
    IF (tbmdpar%print_level>=0) THEN
      WRITE (iw,'(A)') ' TBMD| Hopping matrix element parameters'
      DO i = 1, setup%natom_type
        DO j = i, setup%natom_type
          WRITE (iw,'(A,T40,A,T50,A,T71,A)') ' TBMD| ', &
            adjustr(tbatom(i)%atom_name), adjustr(tbatom(j)%atom_name), &
            adjustr(tbhop(i,j)%htype)
          DO k = 1, tbhop(i,j) %nint
            IF (tbhop(i,j)%htype=='GOODWIN') THEN
              WRITE (iw, &
                '(A,T20,3I1,T25,A,F10.5,T44,A,F10.5,T63,A,F11.5)') &
                ' TBMD| ', tbhop(i,j) %lls(:,k), ' T0[eV]=', &
                tbhop(i,j) %gsp(k) %er0, ' r0[A]=', &
                tbhop(i,j) %gsp(k) %r0, ' rt[A]=', &
                tbhop(i,j) %gsp(k) %rt
              WRITE (iw,'(A,T47,A,F10.5,T66,A,F11.5)') ' TBMD| ', &
                ' NA=', tbhop(i,j) %gsp(k) %na, ' NB=', &
                tbhop(i,j) %gsp(k) %nb
            ELSE
              CALL stop_program('read_tb_hopping_elements','unknown type')
            END IF
          END DO
        END DO
      END DO
      DO i = 1, setup%natom_type
        DO j = i, setup%natom_type
          WRITE (iw,'(A,T40,A,T50,A,T71,F10.3)') &
            ' TBMD| Cutoff distance [A] ', adjustr(tbatom(i)%atom_name), &
            adjustr(tbatom(j)%atom_name), sqrt(tbhop(i,j)%rcutsq)
        END DO
      END DO
      WRITE (iw,*)
    END IF
  END IF

END SUBROUTINE read_tb_hopping_elements

!------------------------------------------------------------------------------!

SUBROUTINE set_int(itype,int)
    CHARACTER(len=*), INTENT(IN)             :: itype
    INTEGER, INTENT(OUT)                     :: int(1:3)

    INTEGER                                  :: i, ia, ie, j

  int = 0
  CALL xstring(itype,ia,ie)
  IF (ie-ia+1/=3) call stop_parser('set_int','non valid integral label')
  DO i = ia, ie
    j = i - ia + 1
    SELECT CASE (itype(i:i))
    CASE DEFAULT
    CASE ('S')
      int(j) = 0
    CASE ('P')
      int(j) = 1
    CASE ('D')
      int(j) = 2
    CASE ('F')
      int(j) = 3
    END SELECT
  END DO
END SUBROUTINE set_int

!!*****
!------------------------------------------------------------------------------!
END MODULE tbmd_input
!------------------------------------------------------------------------------!
