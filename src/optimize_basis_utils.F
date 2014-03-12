!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!
MODULE optimize_basis_utils
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_methods,               ONLY: parser_get_object,&
                                             parser_search_string
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE f77_blas
  USE input_constants,                 ONLY: do_opt_all,&
                                             do_opt_coeff,&
                                             do_opt_exps,&
                                             do_opt_none
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE machine,                         ONLY: m_getcwd
  USE optimize_basis_types,            ONLY: basis_optimization_type,&
                                             derived_basis_info,&
                                             flex_basis_type,&
                                             subset_type
  USE powell,                          ONLY: opt_state_type
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: stop_program
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'optimize_basis_utils'

  PUBLIC :: optimize_basis_init_read_input, get_set_and_basis_id,&
            update_derived_basis_sets

CONTAINS

! *****************************************************************************
!> \brief initialize all parts of the optimization type and read input settings
!> \author Florian Schiffmann
! *****************************************************************************
   
   SUBROUTINE optimize_basis_init_read_input(opt_bas,root_section,para_env,error)
    TYPE(basis_optimization_type)            :: opt_bas
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'optimize_basis_init_read_input', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: main_dir
    INTEGER                                  :: iset, iweight, nrep
    TYPE(section_vals_type), POINTER         :: kind_section, optbas_section, &
                                                powell_section, train_section

      optbas_section => section_vals_get_subs_vals(root_section,"OPTIMIZE_BASIS",error=error)
      powell_section => section_vals_get_subs_vals(optbas_section,"OPTIMIZATION",error=error)
      train_section => section_vals_get_subs_vals(optbas_section,"TRAINING_FILES",error=error)
      kind_section => section_vals_get_subs_vals(optbas_section,"FIT_KIND",error=error)

      CALL section_vals_val_get(optbas_section,"BASIS_TEMPLATE_FILE",c_val=opt_bas%template_basis_file,error=error)
      CALL section_vals_val_get(optbas_section,"BASIS_WORK_FILE",c_val=opt_bas%work_basis_file,error=error)
      CALL section_vals_val_get(optbas_section,"BASIS_OUTPUT_FILE",c_val=opt_bas%output_basis_file,error=error)
      CALL m_getcwd(main_dir)
      opt_bas%work_basis_file=TRIM(ADJUSTL(main_dir))//"/"//TRIM(ADJUSTL(opt_bas%work_basis_file))

      CALL section_vals_val_get(optbas_section,"WRITE_FREQUENCY",i_val=opt_bas%write_frequency,error=error)
      CALL section_vals_val_get(optbas_section,"USE_CONDITION_NUMBER",l_val=opt_bas%use_condition_number,error=error)

      CALL generate_initial_basis(kind_section,opt_bas,para_env,error)
      
      CALL section_vals_get(train_section,n_repetition=opt_bas%ntraining_sets,error=error)
      CALL cp_assert(opt_bas%ntraining_sets/=0,cp_fatal_level,cp_assertion_failed,&
                      routineP,"No training set was specified in the Input",&
                      only_ionode=.TRUE.)
       
      ALLOCATE(opt_bas%training_input(opt_bas%ntraining_sets))
      ALLOCATE(opt_bas%training_dir(opt_bas%ntraining_sets))
      DO iset=1,opt_bas%ntraining_sets
          CALL section_vals_val_get(train_section,"DIRECTORY",c_val=opt_bas%training_dir(iset),&
                                    i_rep_section=iset,error=error)
          CALL section_vals_val_get(train_section,"INPUT_FILE_NAME",c_val=opt_bas%training_input(iset),&
                                    i_rep_section=iset,error=error)
      END DO

      CALL init_powell_var(opt_bas%powell_param,powell_section,error)
      opt_bas%powell_param%nvar=SIZE(opt_bas%x_opt)

      CALL generate_derived_basis_sets(opt_bas,.TRUE.,opt_bas%work_basis_file,para_env,error)
 
      CALL generate_basis_combinations(opt_bas,optbas_section,error)

      CALL section_vals_val_get(optbas_section,"RESIDUUM_WEIGHT",n_rep_val=nrep,error=error)
      ALLOCATE(opt_bas%fval_weight(0:opt_bas%ncombinations))
      opt_bas%fval_weight=1.0_dp
      DO iweight=1,nrep
          CALL section_vals_val_get(optbas_section,"RESIDUUM_WEIGHT",r_val=opt_bas%fval_weight(iweight-1),&
                                    i_rep_val=iweight,error=error)
      END DO

      CALL section_vals_val_get(optbas_section,"CONDITION_WEIGHT",n_rep_val=nrep,error=error)
      ALLOCATE(opt_bas%condition_weight(0:opt_bas%ncombinations))
      opt_bas%condition_weight=1.0_dp
      DO iweight=1,nrep
          CALL section_vals_val_get(optbas_section,"CONDITION_WEIGHT",r_val=opt_bas%condition_weight(iweight-1),&
                                    i_rep_val=iweight,error=error)
      END DO

      CALL generate_computation_groups(opt_bas,optbas_section,para_env,error)            
     
      CALL print_opt_info(opt_bas,error) 
              
   END SUBROUTINE optimize_basis_init_read_input

   SUBROUTINE print_opt_info(opt_bas,error)
    TYPE(basis_optimization_type)            :: opt_bas
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'print_opt_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: icomb, ikind, unit_nr
    TYPE(cp_logger_type), POINTER            :: logger

    logger => cp_error_get_logger(error)
    unit_nr=-1
    IF(logger%para_env%ionode)&
       unit_nr=cp_logger_get_default_unit_nr(logger)

    IF(unit_nr>0)THEN
       WRITE(unit_nr,'(1X,A,A)')"BASOPT| Total number of calcualtions ",&
                            TRIM(cp_to_string(opt_bas%ncombinations*opt_bas%ntraining_sets))
       WRITE(unit_nr,'(A)')""
       DO icomb=1,opt_bas%ncombinations
          WRITE(unit_nr,'(1X,A,A)')"BASOPT| Content of basis combination ",TRIM(cp_to_string(icomb))
          DO ikind=1,opt_bas%nkind
             WRITE(unit_nr,'(1X,A,A4,4X,A,3X,A)')"BASOPT|     Element: ",TRIM(opt_bas%kind_basis(ikind)%element),&
             "Basis set: ",TRIM(opt_bas%kind_basis(ikind)%flex_basis(opt_bas%combination(icomb,ikind))%basis_name) 
          END DO
          WRITE(unit_nr,'(A)')""
       END DO
    END IF
   END SUBROUTINE print_opt_info
          

! *****************************************************************************
!> \brief Generation of the requested basis set combinations if multiple kinds
!>        are fitted at the same time (if not specified create all possible)
!> \author Florian Schiffmann
! *****************************************************************************
   SUBROUTINE generate_basis_combinations(opt_bas,optbas_section,error)
    TYPE(basis_optimization_type)            :: opt_bas
    TYPE(section_vals_type), POINTER         :: optbas_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'generate_basis_combinations', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ikind, j, n_rep
    INTEGER, DIMENSION(:), POINTER           :: i_vals, tmp_i, tmp_i2
    LOGICAL                                  :: explicit, raise

!setup the basis combinations to optimize

      CALL section_vals_val_get(optbas_section,"BASIS_COMBINATIONS",explicit=explicit,n_rep_val=n_rep,error=error)
      IF(.NOT.explicit)THEN
         opt_bas%ncombinations=1
         ALLOCATE(tmp_i(opt_bas%nkind))
         ALLOCATE(tmp_i2(opt_bas%nkind))
         DO ikind=1,opt_bas%nkind
            opt_bas%ncombinations=opt_bas%ncombinations*SIZE(opt_bas%kind_basis(ikind)%flex_basis)
            tmp_i(ikind)=opt_bas%kind_basis(ikind)%nbasis_deriv
         END DO
         ALLOCATE(opt_bas%combination(opt_bas%ncombinations,opt_bas%nkind))
         tmp_i2=0
         DO i=1,opt_bas%ncombinations
             DO j=1,opt_bas%nkind
                opt_bas%combination(i,j)=tmp_i2(j)
             END DO
             tmp_i2(opt_bas%nkind)=tmp_i2(opt_bas%nkind)+1
             raise=.FALSE.
             DO j=opt_bas%nkind,1,-1
                IF(raise)tmp_i2(j)=tmp_i2(j)+1
                IF( tmp_i2(j).gt.tmp_i(j))THEN
                    tmp_i2(j)=0
                    raise=.TRUE.
                END IF
             END DO
          END DO
          DEALLOCATE(tmp_i)
          DEALLOCATE(tmp_i2)
       ELSE
          opt_bas%ncombinations=n_rep
          ALLOCATE(opt_bas%combination(opt_bas%ncombinations,opt_bas%nkind))
          DO i=1,n_rep
             CALL section_vals_val_get(optbas_section,"BASIS_COMBINATIONS",i_vals=i_vals,i_rep_val=i,error=error)
             opt_bas%combination(i,:)=i_vals(:)
          END DO
       END IF

   END SUBROUTINE generate_basis_combinations

! *****************************************************************************
!> \brief returns a mapping from the calcualtion id to the trainings set id and
!>        basis combination id
!> \author Florian Schiffmann
! *****************************************************************************

   SUBROUTINE get_set_and_basis_id(calc_id,opt_bas,set_id,bas_id)
              
    INTEGER                                  :: calc_id
    TYPE(basis_optimization_type)            :: opt_bas
    INTEGER                                  :: set_id, bas_id

    CHARACTER(len=*), PARAMETER :: routineN = 'get_set_and_basis_id', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ncom, nset

      ncom=opt_bas%ncombinations
      nset=opt_bas%ntraining_sets
 
      set_id=(calc_id)/ncom+1
      bas_id=MOD(calc_id,ncom)+1

   END SUBROUTINE
 
! *****************************************************************************
!> \brief Pack calcualtions in groups. If less MPI tasks than systems are used
!>        multiple systems will be assigned to a single MPI task
!> \author Florian Schiffmann
! *****************************************************************************

   SUBROUTINE  generate_computation_groups(opt_bas,optbas_section,para_env,error)
    TYPE(basis_optimization_type)            :: opt_bas
    TYPE(section_vals_type), POINTER         :: optbas_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'generate_computation_groups', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iadd1, iadd2, icount, igroup, &
                                                isize, j, ncalc, nproc, nptot
    INTEGER, DIMENSION(:), POINTER           :: i_vals
    LOGICAL                                  :: explicit

      nproc=para_env%num_pe
      ncalc=opt_bas%ncombinations*opt_bas%ntraining_sets
      CALL section_vals_val_get(optbas_section,"GROUP_PARTITION",explicit=explicit,error=error)

      ! No input information available, try to equally distribute 
      IF(.NOT.explicit)THEN
         IF(nproc.ge.ncalc)THEN
             iadd1=nproc/ncalc
             iadd2=MOD(nproc,ncalc)
             ALLOCATE(opt_bas%comp_group(ncalc))
             ALLOCATE(opt_bas%group_partition(0:ncalc-1))
             DO igroup=0,ncalc-1
                ALLOCATE(opt_bas%comp_group(igroup+1)%member_list(1))
                opt_bas%comp_group(igroup+1)%member_list(1)=igroup
                opt_bas%group_partition(igroup)=iadd1
                IF(igroup.lt.iadd2)opt_bas%group_partition(igroup)=opt_bas%group_partition(igroup)+1
             END DO
          ELSE
             iadd1=ncalc/nproc
             iadd2=MOD(ncalc,nproc)
             ALLOCATE(opt_bas%comp_group(nproc))
             ALLOCATE(opt_bas%group_partition(0:nproc-1))
             icount=0
             DO igroup=0,nproc-1
                opt_bas%group_partition(igroup)=1
                isize=iadd1
                IF(igroup.lt.iadd2)isize=isize+1
                ALLOCATE(opt_bas%comp_group(igroup+1)%member_list(isize))
                DO j=1,isize
                   opt_bas%comp_group(igroup+1)%member_list(j)=icount
                   icount=icount+1
                END DO
             END DO
          END IF
       ELSE

      ! Group partition from input. see if all systems can be assigned. If not add to existing group 
          CALL section_vals_val_get(optbas_section,"GROUP_PARTITION",i_vals=i_vals,error=error)
          isize=SIZE(i_vals)
          nptot=SUM(i_vals)
          CALL cp_assert(nptot==nproc,cp_fatal_level,cp_assertion_failed,&
                       routineP,"Number of processors in group distribution does not match number of MPI tasks."//&
                       " Please change input.",only_ionode=.TRUE.)
          CALL cp_assert(isize.le.ncalc,cp_fatal_level,cp_assertion_failed,&
                       routineP,"Number of Groups larger than number of calculations"//&
                       " Please change input.",only_ionode=.TRUE.)
          IF(nptot.NE.nproc)STOP
          ALLOCATE(opt_bas%comp_group(isize))
          ALLOCATE(opt_bas%group_partition(0:isize-1))
          IF(isize.lt.ncalc)THEN
             iadd1=ncalc/isize
             iadd2=MOD(ncalc,isize)
             icount=0
             DO igroup=0,isize-1
                opt_bas%group_partition(igroup)=i_vals(igroup+1)
                isize=iadd1
                IF(igroup.lt.iadd2)isize=isize+1
                ALLOCATE(opt_bas%comp_group(igroup+1)%member_list(isize))
                DO j=1,isize
                   opt_bas%comp_group(igroup+1)%member_list(j)=icount
                   icount=icount+1
                END DO
             END DO
           ELSE
              DO igroup=0,isize-1      
                opt_bas%group_partition(igroup)=i_vals(igroup+1)
                ALLOCATE(opt_bas%comp_group(igroup+1)%member_list(1))
                opt_bas%comp_group(igroup+1)%member_list(1)=igroup 
              END DO
           END IF
        END IF
    END SUBROUTINE generate_computation_groups

! *****************************************************************************
!> \brief Regenerate the basis sets from reference 0 after an update from the
!>        optimizer to reference was performed, and print basis file if required
!> \author Florian Schiffmann
! *****************************************************************************
   SUBROUTINE update_derived_basis_sets(opt_bas,write_it,output_file,para_env,error)
    TYPE(basis_optimization_type)            :: opt_bas
    LOGICAL                                  :: write_it
    CHARACTER(LEN=default_path_length)       :: output_file
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'update_derived_basis_sets', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ibasis, ikind, unit_nr

    DO ikind=1,opt_bas%nkind
       DO ibasis=1,opt_bas%kind_basis(ikind)%nbasis_deriv
           CALL update_used_parts(opt_bas%kind_basis(ikind)%deriv_info(ibasis),&
                                  opt_bas%kind_basis(ikind)%flex_basis(0),&
                                  opt_bas%kind_basis(ikind)%flex_basis(ibasis),error)
        END DO
    END DO

    IF(write_it)THEN
       IF(para_env%ionode)THEN
            CALL open_file(file_name=output_file,file_status="UNKNOWN",&
                     file_action="WRITE", unit_number=unit_nr)
       ELSE
            unit_nr=-999
       END IF
       DO ikind=1,opt_bas%nkind
          DO ibasis=0,opt_bas%kind_basis(ikind)%nbasis_deriv
              CALL write_basis(opt_bas%kind_basis(ikind)%flex_basis(ibasis),opt_bas%kind_basis(ikind)%element,&
                               unit_nr)
          END DO
       END DO
       IF (para_env%ionode) CALL close_file(unit_number=unit_nr)
     END IF
   
   END SUBROUTINE update_derived_basis_sets

! *****************************************************************************
!> \brief Update the the information in a given basis set
!> \author Florian Schiffmann
! *****************************************************************************

   SUBROUTINE update_used_parts(info_new,basis,basis_new,error)
    TYPE(derived_basis_info)                 :: info_new
    TYPE(flex_basis_type)                    :: basis, basis_new
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'update_used_parts', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: icont, iset, jcont, jset

      jset=0
      DO iset=1,basis%nsets
         IF(info_new%in_use_set(iset))THEN
            jset=jset+1
            basis_new%subset(jset)%exps=basis%subset(iset)%exps
            jcont=0
            DO icont=1,basis%subset(iset)%ncon_tot
               IF(info_new%use_contr(iset)%in_use(icont))THEN
                  jcont=jcont+1
                  basis_new%subset(jset)%coeff(:,jcont)=basis%subset(iset)%coeff(:,icont)
               END IF
            END DO
          END IF
      END DO

   END SUBROUTINE update_used_parts      

! *****************************************************************************
!> \brief Initial generation of the basis set from the file and DERIVED_SET
!> \author Florian Schiffmann
! *****************************************************************************

   SUBROUTINE generate_derived_basis_sets(opt_bas,write_it,output_file,para_env,error)
    TYPE(basis_optimization_type)            :: opt_bas
    LOGICAL                                  :: write_it
    CHARACTER(LEN=default_path_length)       :: output_file
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'generate_derived_basis_sets', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ibasis, ikind, iref, jbasis, &
                                                unit_nr

      DO ikind=1,opt_bas%nkind
         CALL init_deriv_info_ref(opt_bas%kind_basis(ikind)%deriv_info(0),opt_bas%kind_basis(ikind)%flex_basis(0),error)
         ! initialize the reference set used as template for the rest
         DO ibasis=1,opt_bas%kind_basis(ikind)%nbasis_deriv
             iref=opt_bas%kind_basis(ikind)%deriv_info(ibasis)%reference_set
             DO jbasis=0,opt_bas%kind_basis(ikind)%nbasis_deriv
                IF(iref==jbasis)CALL setup_used_parts_init_basis(opt_bas%kind_basis(ikind)%deriv_info(ibasis),&
                                                      opt_bas%kind_basis(ikind)%deriv_info(iref),&
                                                      opt_bas%kind_basis(ikind)%flex_basis(0),&
                                                      opt_bas%kind_basis(ikind)%flex_basis(ibasis),error)
             END DO
          END DO
      END DO
   
      IF(para_env%ionode)THEN
           CALL open_file(file_name=opt_bas%work_basis_file,file_status="UNKNOWN",&
                    file_action="WRITE", unit_number=unit_nr)
      ELSE 
           unit_nr=-999
      END IF
      DO ikind=1,opt_bas%nkind
         DO ibasis=0,opt_bas%kind_basis(ikind)%nbasis_deriv
             opt_bas%kind_basis(ikind)%flex_basis(ibasis)%basis_name=&
                  TRIM(ADJUSTL(opt_bas%kind_basis(ikind)%basis_name))//"-DERIVED_SET-"//ADJUSTL(cp_to_string(ibasis))
             CALL write_basis(opt_bas%kind_basis(ikind)%flex_basis(ibasis),opt_bas%kind_basis(ikind)%element,&
                              unit_nr) 
         END DO
      END DO
      IF (para_env%ionode) CALL close_file(unit_number=unit_nr)
            
   END SUBROUTINE generate_derived_basis_sets

! *****************************************************************************
!> \brief Write a basis set file which can be used from CP2K
!> \author Florian Schiffmann
! *****************************************************************************

   SUBROUTINE write_basis(basis,element,unit_nr)
    TYPE(flex_basis_type)                    :: basis
    CHARACTER(LEN=default_string_length)     :: element
    INTEGER                                  :: unit_nr

    CHARACTER(len=*), PARAMETER :: routineN = 'write_basis', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iexp, iset

      IF (unit_nr > 0) THEN
        WRITE (UNIT=unit_nr,FMT="(A)") TRIM(ADJUSTL(element))//" "//TRIM(ADJUSTL(basis%basis_name))
        WRITE (UNIT=unit_nr,FMT="(1X,I0)") basis%nsets
        DO iset=1,basis%nsets
           WRITE (UNIT=unit_nr,FMT="(30(1X,I0))") basis%subset(iset)%n,basis%subset(iset)%lmin,basis%subset(iset)%lmax,&
                                                  basis%subset(iset)%nexp,basis%subset(iset)%l
           DO iexp=1,basis%subset(iset)%nexp
              WRITE (UNIT=unit_nr,FMT="(30(1X,ES24.16))") basis%subset(iset)%exps(iexp),basis%subset(iset)%coeff(iexp,:)
           END DO
        END DO
      END IF

   END SUBROUTINE write_basis 
 
! *****************************************************************************
!> \brief Initialize the derived basis sets and the vectors containing their
!>        assembly information ehich is used for regeneration of the sets.
!> \author Florian Schiffmann
! *****************************************************************************

   SUBROUTINE setup_used_parts_init_basis(info_new,info_ref,basis,basis_new,error)
    TYPE(derived_basis_info)                 :: info_new, info_ref
    TYPE(flex_basis_type)                    :: basis, basis_new
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_used_parts_init_basis', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, jset, lind, nsets

! copy the reference information on the new set

      ALLOCATE(info_new%in_use_set(SIZE(info_ref%in_use_set)))
      info_new%in_use_set=info_ref%in_use_set
      ALLOCATE(info_new%use_contr(SIZE(info_ref%in_use_set)))
      DO i=1,SIZE(info_ref%in_use_set)
         ALLOCATE(info_new%use_contr(i)%in_use(SIZE(info_ref%use_contr(i)%in_use)))
         info_new%use_contr(i)%in_use=info_ref%use_contr(i)%in_use
      END DO
      DO i=1,info_new%nsets
         info_new%in_use_set(info_new%remove_set(i))=.FALSE.
      END DO
      DO i=1,info_new%ncontr
         lind=convert_l_contr_to_entry(basis%subset(info_new%remove_contr(i,1))%lmin,&
                                       basis%subset(info_new%remove_contr(i,1))%l,&
                                       info_new%remove_contr(i,3),info_new%remove_contr(i,2))
        
         info_new%use_contr(info_new%remove_contr(i,1))%in_use(lind)=.FALSE.
      END DO
     
      nsets=0
      DO i=1,basis%nsets
         IF(info_new%in_use_set(i))nsets=nsets+1
      END DO
      basis_new%nsets=nsets
      ALLOCATE(basis_new%subset(nsets))
      jset=0
      DO i=1,basis%nsets
         IF(info_new%in_use_set(i))THEN
            jset=jset+1
            CALL create_new_subset(basis%subset(i),basis_new%subset(jset),info_new%use_contr(jset)%in_use,error)
         END IF
      END DO
  
   END SUBROUTINE setup_used_parts_init_basis        

! *****************************************************************************
!> \brief Fill the low level information of the derived basis set.
!> \author Florian Schiffmann
! *****************************************************************************

   SUBROUTINE create_new_subset(subset,subset_new,in_use,error)
    TYPE(subset_type)                        :: subset, subset_new
    LOGICAL, DIMENSION(:)                    :: in_use
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_new_subset', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: icon, iind, il
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: tmp_l

      ALLOCATE(tmp_l(SIZE(subset%l)))
      tmp_l=subset%l
      subset_new%lmin=subset%lmin
      subset_new%lmax=subset%lmin-1
      subset_new%nexp=subset%nexp
      subset_new%n=subset%n
      DO il=1,SIZE(subset%l)
         DO icon=1,subset%l(il)
            iind=convert_l_contr_to_entry(subset%lmin,subset%l,icon,subset%lmin+il-1)
            IF(.NOT.in_use(iind))tmp_l(il)=tmp_l(il)-1
         END DO
         IF(tmp_l(il).gt.0)subset_new%lmax=subset_new%lmax+1
      END DO
      subset_new%nl=subset_new%lmax-subset_new%lmin+1
      subset_new%ncon_tot=SUM(tmp_l)
      ALLOCATE(subset_new%l(subset_new%nl))
      ALLOCATE(subset_new%coeff(subset_new%nexp,subset_new%ncon_tot))
      ALLOCATE(subset_new%exps(subset_new%nexp))
      subset_new%exps=subset%exps
      DO il=1,SIZE(subset%l)
        IF(tmp_l(il)==0)EXIT
        subset_new%l(il)=tmp_l(il)
      END DO
      DEALLOCATE(tmp_l)
      iind=0
      DO icon=1,subset%ncon_tot
         IF(in_use(icon))THEN
            iind=iind+1
            subset_new%coeff(:,iind)=subset%coeff(:,icon)
         END IF
      END DO

   END SUBROUTINE create_new_subset
          
! *****************************************************************************
!> \brief for completeness generate the derived info for set 0(reference from file)
!> \author Florian Schiffmann
! *****************************************************************************

   SUBROUTINE init_deriv_info_ref(info,basis,error)
    TYPE(derived_basis_info)                 :: info
    TYPE(flex_basis_type)                    :: basis
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_deriv_info_ref', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

       ALLOCATE(info%in_use_set(basis%nsets))
       info%in_use_set=.TRUE.
       ALLOCATE(info%use_contr(basis%nsets))
       DO i=1,basis%nsets
          ALLOCATE(info%use_contr(i)%in_use(basis%subset(i)%ncon_tot))
          info%use_contr(i)%in_use=.TRUE.
       END DO
   END SUBROUTINE init_deriv_info_ref

! *****************************************************************************
!> \brief get the general information for the basis sets. E.g. what to optimize,
!>        Basis set name, constraints upon optimization and read the reference basis
!> \author Florian Schiffmann
! *****************************************************************************

   SUBROUTINE generate_initial_basis(kind_section,opt_bas,para_env,error)
    TYPE(section_vals_type), POINTER         :: kind_section
    TYPE(basis_optimization_type)            :: opt_bas
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'generate_initial_basis', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, variable_counter
    LOGICAL                                  :: explicit, failure
    TYPE(section_vals_type), POINTER         :: set_section

      failure=.FALSE.
      CALL section_vals_get(kind_section,n_repetition=opt_bas%nkind,error=error)
      ALLOCATE(opt_bas%kind_basis(opt_bas%nkind))

      ! counter to get the number of free variables in optimization
      variable_counter=0
      DO ikind=1,opt_bas%nkind
         CALL section_vals_val_get(kind_section,"_SECTION_PARAMETERS_",c_val=opt_bas%kind_basis(ikind)%element,&
                                   i_rep_section=ikind,error=error)
         CALL section_vals_val_get(kind_section,"BASIS_SET",c_val=opt_bas%kind_basis(ikind)%basis_name,&
                                   i_rep_section=ikind,error=error)
         set_section=>section_vals_get_subs_vals(kind_section,"DERIVED_BASIS_SETS",&
                                     i_rep_section=ikind,error=error)
         CALL section_vals_get(set_section,n_repetition=opt_bas%kind_basis(ikind)%nbasis_deriv,explicit=explicit,error=error)
         IF(.NOT.explicit)opt_bas%kind_basis(ikind)%nbasis_deriv=0
         ALLOCATE(opt_bas%kind_basis(ikind)%flex_basis(0:opt_bas%kind_basis(ikind)%nbasis_deriv))
         ALLOCATE(opt_bas%kind_basis(ikind)%deriv_info(0:opt_bas%kind_basis(ikind)%nbasis_deriv))

         CALL fill_basis_template(kind_section, opt_bas%kind_basis(ikind)%flex_basis(0), opt_bas%template_basis_file,&
                                 opt_bas%kind_basis(ikind)%element,opt_bas%kind_basis(ikind)%basis_name,para_env,ikind,error)

         CALL setup_exp_constraints(kind_section, opt_bas%kind_basis(ikind)%flex_basis(0),error)

         CALL parse_derived_basis(kind_section,opt_bas%kind_basis(ikind)%deriv_info,ikind,error)

         variable_counter= variable_counter+opt_bas%kind_basis(ikind)%flex_basis(0)%nopt 
      END DO       
      
      ALLOCATE(opt_bas%x_opt(variable_counter))

      variable_counter=0
      DO ikind=1,opt_bas%nkind
          CALL assign_x_to_basis(opt_bas%x_opt, opt_bas%kind_basis(ikind)%flex_basis(0),variable_counter,error)
      END DO

      CPPostcondition(variable_counter==SIZE(opt_bas%x_opt),cp_failure_level,routineP,error,failure)
      
   END SUBROUTINE generate_initial_basis

! *****************************************************************************
!> \brief get low level information about how to construc new basis sets from reference
!> \author Florian Schiffmann
! *****************************************************************************

   SUBROUTINE parse_derived_basis(kind_section,deriv_info,ikind,error)
    TYPE(section_vals_type), POINTER         :: kind_section
    TYPE(derived_basis_info), DIMENSION(:)   :: deriv_info
    INTEGER                                  :: ikind
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'parse_derived_basis', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i_rep, iset, jset, n_rep, &
                                                nsets
    INTEGER, DIMENSION(:), POINTER           :: i_vals
    LOGICAL                                  :: explicit
    TYPE(section_vals_type), POINTER         :: set1_section

      nsets=SIZE(deriv_info)-1
      set1_section=>section_vals_get_subs_vals(kind_section,"DERIVED_BASIS_SETS",&
                        i_rep_section=ikind,error=error)
      DO jset=1,nsets
         ! stracnge but as derive info is allcated from 0 to n the count over here has to be shifted
         iset=jset+1
         CALL section_vals_val_get(set1_section,"REFERENCE_SET",i_vals=i_vals, i_rep_section=jset,error=error)
         deriv_info(iset)%reference_set=i_vals(1)
         CALL section_vals_val_get(set1_section,"REMOVE_CONTRACTION",explicit=explicit,n_rep_val=n_rep,&
                                    i_rep_section=jset,error=error)
         deriv_info(iset)%ncontr=n_rep
         IF(explicit)THEN
            ALLOCATE(deriv_info(iset)%remove_contr(n_rep,3))
            DO i_rep=1,n_rep
               CALL section_vals_val_get(set1_section,"REMOVE_CONTRACTION",i_rep_val=i_rep,i_vals=i_vals,&
                                          i_rep_section=jset,error=error)
               deriv_info(iset)%remove_contr(i_rep,:)=i_vals(:)
            END DO
         END IF
         CALL section_vals_val_get(set1_section,"REMOVE_SET",explicit=explicit,n_rep_val=n_rep,&
                                    i_rep_section=jset,error=error)
         deriv_info(iset)%nsets=n_rep
         IF(explicit)THEN
            ALLOCATE(deriv_info(iset)%remove_set(n_rep))
            DO i_rep=1,n_rep
               CALL section_vals_val_get(set1_section,"REMOVE_SET",i_rep_val=i_rep,i_vals=i_vals,&
                                          i_rep_section=jset,error=error)
               deriv_info(iset)%remove_set(i_rep)=i_vals(1)
            END DO
         END IF
      END DO

   END SUBROUTINE parse_derived_basis
                  
! *****************************************************************************
!> \brief get low level information about constraint on exponents
!> \author Florian Schiffmann
! *****************************************************************************

   SUBROUTINE setup_exp_constraints(kind1_section,flex_basis,error)
    TYPE(section_vals_type), POINTER         :: kind1_section
    TYPE(flex_basis_type)                    :: flex_basis
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_exp_constraints', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ipgf, irep, iset, nrep
    INTEGER, DIMENSION(:), POINTER           :: def_exp
    LOGICAL                                  :: is_bound, is_varlim
    TYPE(section_vals_type), POINTER         :: const_section

      const_section=>section_vals_get_subs_vals(kind1_section,"CONSTRAIN_EXPONENTS",error=error)
      CALL section_vals_get(const_section,n_repetition=nrep,error=error)
      DO irep=1,nrep
         CALL section_vals_val_get(const_section,"USE_EXP",i_vals=def_exp,i_rep_section=irep,error=error)
         CALL section_vals_val_get(const_section,"BOUNDARIES",explicit=is_bound,i_rep_section=irep,error=error)
         CALL section_vals_val_get(const_section,"MAX_VAR_FRACTION",explicit=is_varlim,i_rep_section=irep,error=error)
         CALL cp_assert(.NOT.(is_bound.AND.is_varlim),cp_fatal_level,cp_assertion_failed,&
                        routineP,"Exponent has two constraints. This is not possible at the moment."//&
                        " Please change input.",only_ionode=.TRUE.) 
         CALL cp_assert(is_bound.OR.is_varlim,cp_fatal_level,cp_assertion_failed,&
                        routineP,"Exponent is declared to be constraint but none is given"//&
                        " Please change input.",only_ionode=.TRUE.)
         IF(def_exp(1)==-1)THEN
            DO iset=1,flex_basis%nsets
               IF(def_exp(2)==-1)THEN
                  DO ipgf=1,flex_basis%subset(iset)%nexp
                     CALL set_constraint(flex_basis,iset,ipgf,const_section,is_bound,is_varlim,irep,error)
                  END DO
               ELSE 
                  CALL cp_assert(def_exp(2).gt.flex_basis%subset(iset)%nexp,cp_fatal_level,cp_assertion_failed,&
                       routineP,"Exponent declared in constraint is larger than number of exponents in the set"//&
                       " Please change input.",only_ionode=.TRUE.)
                  CALL set_constraint(flex_basis,iset,def_exp(2),const_section,is_bound,is_varlim,irep,error)
               END IF
            END DO
         ELSE
            CALL cp_assert(def_exp(1).le.flex_basis%nsets,cp_fatal_level,cp_assertion_failed,&
                   routineP,"Set number of constraint is larger than number of sets in the template basis set."//&
                  " Please change input.",only_ionode=.TRUE.)
            IF(def_exp(2)==-1)THEN
               DO ipgf=1,flex_basis%subset(iset)%nexp
                  CALL set_constraint(flex_basis,def_exp(1),ipgf,const_section,is_bound,is_varlim,irep,error)
               END DO
            ELSE
               CALL cp_assert(def_exp(2).le.flex_basis%subset(def_exp(1))%nexp,cp_fatal_level,cp_assertion_failed,&
                       routineP,"Exponent declared in constraint is larger than number of exponents in the set"//&
                       " Please change input.",only_ionode=.TRUE.)
                CALL set_constraint(flex_basis,def_exp(1),def_exp(2),const_section,is_bound,is_varlim,irep,error)
            END IF
         END IF
      END DO
                  
   END SUBROUTINE setup_exp_constraints    
                      
! *****************************************************************************
!> \brief put the constraint information in type and process if requires
!>        BOUNDARIES constraint gets transformed into MAX_VAR_FRACTION constraint.
!> \author Florian Schiffmann
! *****************************************************************************
         
   SUBROUTINE set_constraint(flex_basis,iset,ipgf,const_section,is_bound,is_varlim,irep,error)
    TYPE(flex_basis_type)                    :: flex_basis
    INTEGER                                  :: iset, ipgf
    TYPE(section_vals_type), POINTER         :: const_section
    LOGICAL                                  :: is_bound, is_varlim
    INTEGER                                  :: irep
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_constraint', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: r_val
    REAL(KIND=dp), DIMENSION(:), POINTER     :: r_vals

      CALL cp_assert(.NOT.flex_basis%subset(iset)%exp_has_const(ipgf),cp_fatal_level,cp_assertion_failed,&
            routineP,"Multiple constraints due to collision in CONSTRAIN_EXPONENTS."//&
            " Please change input.",only_ionode=.TRUE.)
      flex_basis%subset(iset)%exp_has_const(ipgf)=.TRUE.
      IF(is_bound)THEN
         flex_basis%subset(iset)%exp_const(ipgf)%const_type=0
         CALL section_vals_val_get(const_section,"BOUNDARIES",r_vals=r_vals,i_rep_section=irep,error=error)
         flex_basis%subset(iset)%exp_const(ipgf)%llim=MINVAL(r_vals)
         flex_basis%subset(iset)%exp_const(ipgf)%ulim=MAXVAL(r_vals)
         r_val=flex_basis%subset(iset)%exps(ipgf)
         IF(flex_basis%subset(iset)%exps(ipgf).gt.MAXVAL(r_vals).OR.flex_basis%subset(iset)%exps(ipgf).lt.MINVAL(r_vals))&
         CALL cp_assert(.FALSE.,cp_fatal_level,cp_assertion_failed,&
                       routineP,"Exponent "//cp_to_string(r_val)//&
                       " declared in constraint is out of bounds of constraint"//cp_to_string(MINVAL(r_vals))//&
                       " to"// cp_to_string(MAXVAL(r_vals))//&
                       " Please change input.",only_ionode=.TRUE.)
         flex_basis%subset(iset)%exp_const(ipgf)%init=SUM(r_vals)/2.0_dp
         flex_basis%subset(iset)%exp_const(ipgf)%var_fac=MAXVAL(r_vals)/flex_basis%subset(iset)%exp_const(ipgf)%init-1.0_dp
      END IF
      IF(is_varlim)THEN
         flex_basis%subset(iset)%exp_const(ipgf)%const_type=1
         CALL section_vals_val_get(const_section,"MAX_VAR_FRACTION",r_vals=r_vals,i_rep_section=irep,error=error)
         flex_basis%subset(iset)%exp_const(ipgf)%var_fac=r_vals(1)
         flex_basis%subset(iset)%exp_const(ipgf)%init=flex_basis%subset(iset)%exps(ipgf)
      END IF
 
   END SUBROUTINE set_constraint

! *****************************************************************************
!> \brief Initialize the optimization vector with the values from the refernece sets
!> \author Florian Schiffmann
! *****************************************************************************
                  
   SUBROUTINE assign_x_to_basis(x,basis,x_ind,error)
    REAL(KIND=dp), DIMENSION(:)              :: x
    TYPE(flex_basis_type)                    :: basis
    INTEGER                                  :: x_ind
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'assign_x_to_basis', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: icont, ipgf, iset

      DO iset=1,basis%nsets
        DO ipgf=1,basis%subset(iset)%nexp
           IF(basis%subset(iset)%opt_exps(ipgf))THEN
              x_ind=x_ind+1
              basis%subset(iset)%exp_x_ind(ipgf)=x_ind
              x(x_ind)= basis%subset(iset)%exps(ipgf)
           END IF
           DO icont=1,basis%subset(iset)%ncon_tot
              IF(basis%subset(iset)%opt_coeff(ipgf,icont))THEN
                 x_ind=x_ind+1
                 basis%subset(iset)%coeff_x_ind(ipgf,icont)=x_ind
                 x(x_ind)= basis%subset(iset)%coeff(ipgf,icont)
              END IF
           END DO
         END DO
      END DO
   
   END SUBROUTINE assign_x_to_basis      

! *****************************************************************************
!> \brief Fill the reference set and get the free varialbles from input
!> \author Florian Schiffmann
! *****************************************************************************

   SUBROUTINE fill_basis_template(kind1_section, flex_basis, template_basis_file, element, basis_name, para_env, ikind,error)
    TYPE(section_vals_type), POINTER         :: kind1_section
    TYPE(flex_basis_type)                    :: flex_basis
    CHARACTER(LEN=default_path_length)       :: template_basis_file
    CHARACTER(LEN=default_string_length)     :: element, basis_name
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: ikind
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fill_basis_template', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: icont, idof, ipgf, irep, &
                                                iset, nrep
    INTEGER, DIMENSION(:), POINTER           :: switch

      CALL parse_basis( flex_basis, template_basis_file, element, basis_name, para_env, error )

      ! get the optimizable parameters. Many way to modify them but in the end only logical matrix
      ! is either set or values get flipped according to the input
      CALL section_vals_val_get(kind1_section,"INITIAL_DEGREES_OF_FREEDOM", i_val=idof,&
                                i_rep_section=ikind, error=error)
      DO iset=1,flex_basis%nsets
         SELECT CASE(idof)
         CASE(do_opt_none)
            ! initialization in parse subset did the job
         CASE(do_opt_all)
            flex_basis%subset(iset)%opt_coeff=.TRUE.
            flex_basis%subset(iset)%opt_exps=.TRUE.
         CASE(do_opt_coeff)
            flex_basis%subset(iset)%opt_coeff=.TRUE.
         CASE(do_opt_exps)
            flex_basis%subset(iset)%opt_exps=.TRUE.
         CASE DEFAULT
            CALL cp_assert(.FALSE.,cp_fatal_level,cp_assertion_failed,&
                      routineP,"No initialization available?????",&
                      only_ionode=.TRUE.)
         END SELECT
      END DO

      CALL section_vals_val_get(kind1_section,"SWITCH_CONTRACTION_STATE",n_rep_val=nrep,i_rep_section=ikind,error=error)
      DO irep=1, nrep
         CALL section_vals_val_get(kind1_section,"SWITCH_CONTRACTION_STATE", i_rep_val=irep,&
                                   i_rep_section=ikind,i_vals=switch,error=error)
         icont=convert_l_contr_to_entry(flex_basis%subset(switch(1))%lmin,flex_basis%subset(switch(1))%l,switch(3),switch(2))
         DO ipgf=1, flex_basis%subset(switch(1))%nexp
            flex_basis%subset(switch(1))%opt_coeff(ipgf,icont)=.NOT.flex_basis%subset(switch(1))%opt_coeff(ipgf,icont)
         END DO
      END DO
 
      CALL section_vals_val_get(kind1_section,"SWITCH_COEFF_STATE",n_rep_val=nrep,i_rep_section=ikind,error=error)
      DO irep=1, nrep
         CALL section_vals_val_get(kind1_section,"SWITCH_COEFF_STATE", i_rep_val=irep,&
                                   i_rep_section=ikind,i_vals=switch,error=error)
         icont=convert_l_contr_to_entry(flex_basis%subset(switch(1))%lmin,flex_basis%subset(switch(1))%l,switch(3),switch(2))
         flex_basis%subset(switch(1))%opt_coeff(switch(4),icont)=&
                     .NOT.flex_basis%subset(switch(1))%opt_coeff(switch(4),icont)
      END DO

      CALL section_vals_val_get(kind1_section,"SWITCH_EXP_STATE",n_rep_val=nrep,i_rep_section=ikind,error=error)
      DO irep=1, nrep
         CALL section_vals_val_get(kind1_section,"SWITCH_EXP_STATE", i_rep_val=irep,&
                                   i_rep_section=ikind,i_vals=switch,error=error)
         flex_basis%subset(switch(1))%opt_exps(switch(2))=.NOT.flex_basis%subset(switch(1))%opt_exps(switch(2))
      END DO

      CALL section_vals_val_get(kind1_section,"SWITCH_SET_STATE",n_rep_val=nrep,i_rep_section=ikind,error=error)
      DO irep=1, nrep
          CALL section_vals_val_get(kind1_section,"SWITCH_SET_STATE", i_rep_val=irep,&
                                    i_rep_section=ikind,i_vals=switch,error=error)
         DO ipgf=1, flex_basis%subset(switch(2))%nexp
            SELECT CASE(switch(1))
            CASE(0) ! switch all states in the set
               DO icont=1,flex_basis%subset(switch(2))%ncon_tot 
                  flex_basis%subset(switch(2))%opt_coeff(ipgf,icont)=&
                       .NOT.flex_basis%subset(switch(2))%opt_coeff(ipgf,icont)
               END DO
               flex_basis%subset(switch(2))%opt_exps(ipgf)=.NOT.flex_basis%subset(switch(2))%opt_exps(ipgf)
            CASE(1) ! switch only exp
               flex_basis%subset(switch(2))%opt_exps(ipgf)=.NOT.flex_basis%subset(switch(2))%opt_exps(ipgf)
            CASE(2) ! switch only coeff
               DO icont=1,flex_basis%subset(switch(2))%ncon_tot
                  flex_basis%subset(switch(2))%opt_coeff(ipgf,icont)=&
                       .NOT.flex_basis%subset(switch(2))%opt_coeff(ipgf,icont)
               END DO
            CASE DEFAULT
               CALL cp_assert(.FALSE.,cp_fatal_level,cp_assertion_failed,&
                         routineP,"Invalid option in SWITCH_SET_STATE, 1st value has to be 0, 1 or 2",&
                         only_ionode=.TRUE.)
            END SELECT
         END DO
      END DO
    
      ! perform a final modification. If basis set is uncontracted coefficient will never have to be optimized
      DO irep=1, flex_basis%nsets
         IF(flex_basis%subset(irep)%nexp==1)THEN
             DO ipgf=1,flex_basis%subset(irep)%nexp
                flex_basis%subset(irep)%opt_coeff(ipgf,1)=.FALSE.
             END DO
         END IF
      END DO

      ! finally count the total number of free parameters
      flex_basis%nopt=0
      DO irep=1, flex_basis%nsets
         DO ipgf=1,flex_basis%subset(irep)%nexp
            DO icont=1,flex_basis%subset(irep)%ncon_tot
               IF(flex_basis%subset(irep)%opt_coeff(ipgf,icont)) flex_basis%nopt= flex_basis%nopt+1
            END DO
            IF(flex_basis%subset(irep)%opt_exps(ipgf)) flex_basis%nopt= flex_basis%nopt+1
         END DO
      END DO

   END SUBROUTINE fill_basis_template

! *****************************************************************************
!> \brief Helper function to parse input. Converts l and index position of
!>        a contraction to index in the contraction array of the set using lmin and nl
!> \author Florian Schiffmann
! *****************************************************************************

   FUNCTION convert_l_contr_to_entry(lmin,nl,icontr,l)    RESULT(ientry)
    INTEGER                                  :: lmin
    INTEGER, DIMENSION(:)                    :: nl
    INTEGER                                  :: icontr, l, ientry

    INTEGER                                  :: i, icon2l, iwork

      iwork=l-lmin
      icon2l=0
      DO i=1,iwork
         icon2l=icon2l+nl(i)
      END DO
      ientry=icon2l+icontr

   END FUNCTION convert_l_contr_to_entry
      

! *****************************************************************************
!> \brief Read the reference basis sets from the template basis file
!> \author Florian Schiffmann
! *****************************************************************************

   SUBROUTINE parse_basis(flex_basis, template_basis_file,element, basis_name, para_env ,error)
    TYPE(flex_basis_type)                    :: flex_basis
    CHARACTER(LEN=default_path_length)       :: template_basis_file
    CHARACTER(LEN=default_string_length)     :: element, basis_name
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'parse_basis', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=240)                       :: line
    CHARACTER(LEN=242)                       :: line2
    CHARACTER(LEN=LEN(basis_name)+2)         :: basis_name2
    CHARACTER(LEN=LEN(element)+2)            :: element2
    INTEGER                                  :: iset, strlen1, strlen2
    LOGICAL                                  :: basis_found, found, match
    TYPE(cp_parser_type), POINTER            :: parser

      basis_found = .FALSE.
      CALL uppercase(element)
      CALL uppercase(basis_name)
      NULLIFY(parser)
      CALL parser_create(parser,template_basis_file,para_env=para_env,error=error)

      search_loop: DO 
         CALL parser_search_string(parser,TRIM(basis_name),.TRUE.,found,line,error=error)
         IF(found)THEN
            match = .FALSE.
            CALL uppercase(line)
            ! Check both the element symbol and the basis set name
            line2 = " "//line//" "
            element2 = " "//TRIM(element)//" "
            basis_name2 = " "//TRIM(basis_name)//" "
            strlen1 = LEN_TRIM(element2) + 1
            strlen2 = LEN_TRIM(basis_name2) + 1
            IF ( (INDEX(line2,element2(:strlen1)) > 0).AND.&
                         (INDEX(line2,basis_name2(:strlen2)) > 0) ) match = .TRUE.
            IF (match) THEN
               CALL parser_get_object(parser,flex_basis%nsets,newline=.TRUE.,error=error)
               ALLOCATE(flex_basis%subset(flex_basis%nsets))
               DO iset=1,flex_basis%nsets
                  CALL parse_subset(parser,flex_basis%subset(iset),error)
               END DO
               basis_found = .TRUE.
               EXIT
            END IF
         ELSE
            EXIT search_loop
         END IF     
      END DO search_loop
      CALL parser_release(parser,error=error)

      IF( .NOT. basis_found ) CALL stop_program(routineN,moduleN,__LINE__,&
                              "The requested basis set <"//TRIM(basis_name)//&
                              "> for element <"//TRIM(element)//"> was not "//&
                              "found in the template basis set file ",para_env)

   END SUBROUTINE parse_basis

! *****************************************************************************
!> \brief Read the subset information from the template basis file
!> \author Florian Schiffmann
! *****************************************************************************

   SUBROUTINE parse_subset(parser,subset,error)
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(subset_type)                        :: subset
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'parse_subset', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=20*default_string_length)  :: line_att
    INTEGER                                  :: icon1, icon2, il, ipgf, &
                                                ishell, istart
    REAL(KIND=dp)                            :: gs_scale
    REAL(KIND=dp), POINTER                   :: r_val

      line_att=""
      CALL parser_get_object(parser,subset%n,newline=.TRUE.,error=error)
      CALL parser_get_object(parser,subset%lmin,error=error)
      CALL parser_get_object(parser,subset%lmax,error=error)
      CALL parser_get_object(parser,subset%nexp,error=error)
      subset%nl=subset%lmax-subset%lmin+1
      ALLOCATE(r_val)
      ALLOCATE(subset%l(subset%nl))
      ALLOCATE(subset%exps(subset%nexp))
      ALLOCATE(subset%exp_has_const(subset%nexp))
      subset%exp_has_const=.FALSE.
      ALLOCATE(subset%opt_exps(subset%nexp))
      subset%opt_exps=.FALSE.
      ALLOCATE(subset%exp_const(subset%nexp))
      ALLOCATE(subset%exp_x_ind(subset%nexp))
      DO ishell=1,subset%nl
         CALL parser_get_object(parser,subset%l(ishell),error=error)
      END DO
      subset%ncon_tot=SUM(subset%l)
      ALLOCATE(subset%coeff(subset%nexp,subset%ncon_tot))
      ALLOCATE(subset%opt_coeff(subset%nexp,subset%ncon_tot))
      subset%opt_coeff=.FALSE.
      ALLOCATE(subset%coeff_x_ind(subset%nexp,subset%ncon_tot))
      DO ipgf=1,subset%nexp
         CALL parser_get_object(parser,r_val,newline=.TRUE.,error=error)
         subset%exps(ipgf)=r_val
         DO ishell=1,subset%ncon_tot
            CALL parser_get_object(parser,r_val,error=error)
            subset%coeff(ipgf,ishell)=r_val
         END DO
      END DO

      ! orthonormalize contraction coefficients using gram schmidt
      istart=1
      DO il=1,subset%nl
         DO icon1=istart,istart+subset%l(il)-2
            DO icon2=icon1+1,istart+subset%l(il)-1
               gs_scale=DOT_PRODUCT(subset%coeff(:,icon2),subset%coeff(:,icon1))/&
                        DOT_PRODUCT(subset%coeff(:,icon1),subset%coeff(:,icon1))
               subset%coeff(:,icon2)=subset%coeff(:,icon2)-gs_scale*subset%coeff(:,icon1)
            END DO
         END DO
         istart=istart+subset%l(il)
      END DO

      ! just to get an understandable basis normalize coefficients
      DO icon1=1,subset%ncon_tot
         subset%coeff(:,icon1)=subset%coeff(:,icon1)/&
               SQRT(DOT_PRODUCT(subset%coeff(:,icon1),subset%coeff(:,icon1)))
      END DO
      DEALLOCATE(r_val)
             
   END  SUBROUTINE parse_subset   
      
! *****************************************************************************
!> \brief Initialize the variables for the powell optimizer
!> \author Florian Schiffmann
! *****************************************************************************

   SUBROUTINE init_powell_var(p_param,powell_section,error)
      TYPE(opt_state_type)  :: p_param
      TYPE(cp_error_type), INTENT(inout)       :: error

      CHARACTER(len=*), PARAMETER :: routineN = 'init_powell_var',&
        routineP = moduleN//':'//routineN

      TYPE(section_vals_type), POINTER         :: optbas_section,powell_section

      p_param%state=0
      p_param%nvar=0
      p_param%iprint=0
      p_param%unit=6
      CALL section_vals_val_get(powell_section,"ACCURACY",r_val=p_param%rhoend,error=error)
      CALL section_vals_val_get(powell_section,"STEP_SIZE",r_val=p_param%rhobeg,error=error)
      CALL section_vals_val_get(powell_section,"MAX_FUN",i_val=p_param%maxfun,error=error)

   END SUBROUTINE init_powell_var
 
END MODULE optimize_basis_utils
