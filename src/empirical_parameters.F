!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!***** cp2k/empirical_parameters [1.0] *
!!
!!   NAME
!!     empirical_parameters
!!
!!   FUNCTION
!!   defines  the empirical parameters for the polarization code
!!   (first order hohemberg-kohn kernel and hardness kernel parameters)
!!   reads them from the set file
!!
!!   AUTHOR
!!     gloria
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE empirical_parameters

  USE atomic_kinds,                    ONLY: kind_info_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE parser,                          ONLY: get_next,&
                                             p_error,&
                                             parser_end,&
                                             parser_init,&
                                             read_line,&
                                             search_label,&
                                             stop_parser,&
                                             test_next
  USE string_utilities,                ONLY: str_search,&
                                             uppercase
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_empirical_parameters, empirical_parameter_type

  TYPE empirical_parameter_type
    CHARACTER (len = 4) :: aname
    INTEGER :: nset
    INTEGER, DIMENSION (:), pointer :: l
    REAL(KIND=dp), DIMENSION (:), POINTER :: hk_param
    REAL(KIND=dp), DIMENSION (:), POINTER :: hardness_param
  END TYPE empirical_parameter_type

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

 SUBROUTINE read_empirical_parameters (  ki, set_fn, atom_names, empparm, globenv )


    TYPE(kind_info_type), dimension(:), &
      intent(in)                             :: ki
    CHARACTER(LEN=*), INTENT(IN)             :: set_fn
    CHARACTER(len=*), DIMENSION(:), &
      INTENT(IN)                             :: atom_names
    TYPE(empirical_parameter_type), &
      DIMENSION(:), POINTER                  :: empparm
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=20)                        :: string, string2
    CHARACTER(LEN=6)                         :: aname, label
    CHARACTER(LEN=6), DIMENSION(:), POINTER  :: element_symbol
    INTEGER                                  :: i, iat, icount, ierror, ilen, &
                                                ios, iset, iw, natom_types, &
                                                nkind, nset
    TYPE(empirical_parameter_type), &
      DIMENSION(:), POINTER                  :: read_ep

  natom_types = SIZE ( atom_names )
  iw = globenv % scr
  icount = 0
  nkind = 0

  ALLOCATE (read_ep(natom_types), stat = ios)
  IF (ios /= 0 ) CALL stop_memory &
       ( 'read_empirical_parameters', 'read_ep', natom_types )
  ALLOCATE (element_symbol(natom_types), stat = ios)
  IF (ios /= 0 ) CALL stop_memory &
       ( 'read_emporical_parameters', 'elementsymbol', natom_types )
  ALLOCATE (empparm(size(ki)), stat = ios)
  IF (ios /= 0 ) CALL stop_memory &
       ( 'read_empirical_parameters', 'empparm', nkind )
  

!..parse the input section
  label = '&EMPIRICAL'
  CALL parser_init ( set_fn, globenv )
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF ( ierror /= 0 ) THEN
     IF ( globenv % ionode ) THEN
        WRITE ( iw, '( A )' ) ' No input section &EMPIRICAL found on file '
        WRITE ( iw, '( T2, A )' ) set_fn
     END IF
     CALL stop_parser ( 'read_empirical_parameters', '&EMPIRICAL' )
  ELSE
     CALL read_line
     DO WHILE (test_next()/='X')
        ilen = 8
        CALL get_next(string,ilen)
        iat = str_search(atom_names,natom_types,string)
        nkind = nkind + 1

        DO
           CALL read_line
           ilen = 8
           CALL get_next(string2,ilen)
           CALL uppercase ( string2 )

           SELECT CASE (string2)
           CASE DEFAULT
              CALL p_error()
              CALL stop_parser ( 'read_empirical_parameter','unknown option')

           CASE ( 'NSET')
              CALL get_next ( nset )
              read_ep (iat) % nset =  nset     

              IF ( read_ep (iat) % nset  > 0) THEN
                 ALLOCATE (read_ep(iat)%l(nset), STAT = ios)
                 ALLOCATE (read_ep(iat)%hk_param(nset), stat = ios)
                 ALLOCATE (read_ep(iat)%hardness_param(nset), stat = ios)
                 IF (ios /= 0 ) CALL stop_memory &
                      ( 'read_emporical_parameters', 'read_ep%l', nset )
                 element_symbol (iat)  = string
              END IF                                              

           CASE ( 'EMPAR')
               CALL get_next ( i )
               CALL get_next ( read_ep(iat) % l (i) )
               CALL get_next ( read_ep(iat) % hk_param (i) )
               CALL get_next ( read_ep(iat) % hardness_param (i) )
      
           CASE ( 'END')
                 ilen = 0
                 CALL get_next(string2,ilen)
                 EXIT
           END SELECT

        END DO
   
        CALL read_line

     END DO
      
  END IF

  CALL parser_end                                                    

! map the empirical parameter to the corresponding structure
! and checks that the mapping is the same as in kind_info

  IF (nkind /= size(ki)) CALL stop_program ('read_empirical_parameters', &
     'inconsistent number of polarizable atomic kinds')

  DO i = 1, nkind
  
     aname = ki(i)% element_symbol
     iat = str_search(element_symbol,natom_types,aname)

     empparm (i) % aname = aname
     nset = read_ep(iat)% nset
     empparm (i) % nset = read_ep(iat)% nset

     IF (empparm(i)% nset /= ki(i) % orb_basis_set % nset) CALL stop_program &
       ('read_empirical_parameters', 'number of polarizable sets not matching')
       
     ALLOCATE ( empparm(i) % l (nset), STAT = ios )
     ALLOCATE ( empparm(i) % hk_param (nset), stat = ios )
     ALLOCATE ( empparm(i) % hardness_param (nset), stat = ios )

     IF (ios /= 0 ) CALL stop_memory &
                ( 'read_molecule_section', 'empparm%l', nset )

     DO iset = 1, nset 
      
       empparm(i) % l (iset) = read_ep ( iat ) % l (iset)
       empparm(i) % hk_param (iset) = read_ep ( iat ) % hk_param (iset)
       empparm(i) % hardness_param (iset) = read_ep ( iat ) % hardness_param (iset)

     END DO
  END DO
  
! ..write some information to output

  IF (globenv%ionode) THEN

     IF (globenv%print_level>=0) THEN
     
        DO i = 1, nkind

           WRITE ( iw, '( A,T61,A )' ) ' EMPIRICAL_PARAMETERS| atom type ', &
                ADJUSTR ( empparm(i) % aname )      
           WRITE ( iw, '( A,T71,I10 )' ) ' EMPIRICAL_PARAMETERS| Number of sets ', &
                empparm(i) % nset

           DO iset = 1, empparm (i) % nset     
   
              WRITE ( iw, '( A,T71,I10 )' ) ' EMPIRICAL| angular momentum number ', &
                   empparm(i) % l (iset)

              WRITE ( iw, '( A,T71,F12.6 )' ) ' EMPIRICAL| Hohemberg-Kohn parameter ', &
                   empparm(i) % hk_param (iset)
   
              WRITE ( iw, '( A,T71,F12.6 )' ) ' EMPIRICAL| Hardness parameter ', &
                   empparm(i) % hardness_param (iset)
           END DO                                                       
 
           WRITE ( iw,'( )' )

        END DO

     END IF                        

  END IF                        
  
! deallocates 
  
  DO i = 1, nkind

     IF (ASSOCIATED (read_ep (i) %l)) DEALLOCATE (read_ep (i) %l, stat = ios)
        IF (ios /= 0 ) CALL stop_memory &
                ( 'read_molecule_section', ' deallocate read_ep%l' )

     IF (ASSOCIATED(read_ep (i) %hk_param)) DEALLOCATE (read_ep (i) %hk_param, stat = ios)
        IF (ios /= 0 ) CALL stop_memory &
                ( 'read_molecule_section', ' deallocate read_ep%hk_param' )

     IF (ASSOCIATED(read_ep (i) %hardness_param )) &
              DEALLOCATE (read_ep (i) % hardness_param, stat = ios)
        IF (ios /= 0 ) CALL stop_memory &
                ( 'read_molecule_section', ' deallocate read_ep%l' )

  END DO

  IF (ASSOCIATED(read_ep)) DEALLOCATE (read_ep, stat = ios)
     IF (ios /= 0 ) CALL stop_memory &
             ( 'read_molecule_section', ' deallocate read_ep' )

END SUBROUTINE read_empirical_parameters
 
!******************************************************************************              

END MODULE empirical_parameters

!******************************************************************************
