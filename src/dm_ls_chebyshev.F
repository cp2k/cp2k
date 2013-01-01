!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines using linear scaling chebyshev methods 
!> \par History
!>       2012.10 created [Jinwoong Cha]
!> \author Jinwoong Cha
! *****************************************************************************
MODULE dm_ls_chebyshev
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_add_on_diag, cp_dbcsr_copy, cp_dbcsr_create, &
       cp_dbcsr_filter, cp_dbcsr_frobenius_norm, cp_dbcsr_get_info, &
       cp_dbcsr_get_occupation, cp_dbcsr_init, cp_dbcsr_multiply, &
       cp_dbcsr_release, cp_dbcsr_scale, cp_dbcsr_set, cp_dbcsr_trace
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_type
  USE cp_dbcsr_util,                   ONLY: lanczos_alg_serial
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE dbcsr_types,                     ONLY: dbcsr_type_no_symmetry
  USE dm_ls_scf_types,                 ONLY: ls_scf_env_type
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_walltime
  USE mathconstants,                   ONLY: pi
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dm_ls_chebyshev'

  PUBLIC :: compute_chebyshev

CONTAINS

! *****************************************************************************
!> \brief compute chebyshev polynomials up to order n for a given value of x
!> \par History
!>       2012.11 created [Jinwoong Cha]
!> \author Jinwoong Cha
! *****************************************************************************
  SUBROUTINE chebyshev_poly(value,x,n)
    REAL(KIND=dp), INTENT(OUT)               :: value
    REAL(KIND=dp), INTENT(IN)                :: x
    INTEGER, INTENT(IN)                      :: n

!polynomial values 
!number of chev polynomials

    value = COS((n-1)*ACOS(x))     

  END SUBROUTINE chebyshev_poly

! *****************************************************************************
!> \brief kernel for chebyshev polynomials expansion (Jackson kernel)
!> \par History
!>       2012.11 created [Jinwoong Cha]
!> \author Jinwoong Cha
! *****************************************************************************
  SUBROUTINE kernel(value,n,nc)
    REAL(KIND=dp), INTENT(OUT)               :: value
    INTEGER, INTENT(IN)                      :: n, nc

!kernel at n
!n-1 order of chebyshev polynomials
!number of total chebyshev polynomials
!Kernel define

     value = 1.0_dp/(nc+1.0_dp)*((nc-(n-1)+1.0_dp)*&
                        COS(pi*(n-1)/(nc+1.0_dp))+SIN(pi*(n-1)/(nc+1.0_dp))*1.0_dp/TAN(pi/(nc+1.0_dp)))

  END SUBROUTINE kernel

! *****************************************************************************
!> \brief compute properties based on chebyshev expansion
!> \par History
!>       2012.10 created [Jinwoong Cha]
!> \author Jinwoong Cha
! *****************************************************************************
  SUBROUTINE compute_chebyshev(ls_scf_env,error)
    TYPE(ls_scf_env_type)                    :: ls_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'compute_chebyshev', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: lanczos_threshold = 1.0E-4_dp,&
                                                scale_evals = 1.01_dp

    INTEGER                                  :: handle, icheb, igrid, iinte, &
                                                ispin, n_gridpoint_dos, &
                                                ncheb, ninte, Nrows, &
                                                unit_dos, unit_nr
    REAL(KIND=dp) :: chev_T, chev_T_dos, dummy1, dummy2, ev1, ev2, final, &
      frob_matrix, initial, interval_a, interval_b, max_ev, min_ev, occ, &
      sev1, sev2, summa, t1, t2, trace_dm
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: aitchev_T, chev_E, &
                                                chev_Es_dos, dos, E_inte, &
                                                gdensity, kernel_g, mu, &
                                                sqrt_vec
    TYPE(cp_dbcsr_type)                      :: matrix_dummy1, matrix_dummy2, &
                                                matrix_F, matrix_tmp1, &
                                                matrix_tmp2, matrix_tmp3
    TYPE(cp_logger_type), POINTER            :: logger

    IF (.NOT.ls_scf_env%chebyshev%compute_chebyshev) RETURN

    CALL timeset(routineN,handle)

    ncheb=ls_scf_env%chebyshev%n_chebyshev
    ninte=2*ncheb
    n_gridpoint_dos=ls_scf_env%chebyshev%n_gridpoint_dos
    ev1 = ls_scf_env%chebyshev%min_energy
    ev2 = ls_scf_env%chebyshev%max_energy

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    IF (unit_nr>0) THEN
       WRITE(unit_nr,'()') 
       WRITE(unit_nr,'(T2,A)') "STARTING CHEBYSHEV CALCULATION"
    ENDIF
     
    ! create 3 temporary matrices
    CALL cp_dbcsr_init(matrix_tmp1,error=error) !initialization
    CALL cp_dbcsr_create(matrix_tmp1,template=ls_scf_env%matrix_s,matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_init(matrix_tmp2,error=error)
    CALL cp_dbcsr_create(matrix_tmp2,template=ls_scf_env%matrix_s,matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_init(matrix_tmp3,error=error)
    CALL cp_dbcsr_create(matrix_tmp3,template=ls_scf_env%matrix_s,matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_init(matrix_F,error=error)
    CALL cp_dbcsr_create(matrix_F,template=ls_scf_env%matrix_s,matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_init(matrix_dummy1,error=error)
    CALL cp_dbcsr_create(matrix_dummy1,template=ls_scf_env%matrix_s,matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_init(matrix_dummy2,error=error) ! density matrix
    CALL cp_dbcsr_create(matrix_dummy2,template=ls_scf_env%matrix_s,matrix_type=dbcsr_type_no_symmetry,error=error)
 
    DO ispin=1,SIZE(ls_scf_env%matrix_ks)

       ! create matrix_F=inv(sqrt(S))*H*inv(sqrt(S))
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s_sqrt_inv, ls_scf_env%matrix_ks(ispin), &
                              0.0_dp, matrix_tmp1, filter_eps=ls_scf_env%eps_filter,error=error)
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, ls_scf_env%matrix_s_sqrt_inv, &
                              0.0_dp, matrix_F, filter_eps=ls_scf_env%eps_filter,error=error)

       ! find largest and smallest eigenvalues
       CALL lanczos_alg_serial(matrix_F, max_ev, min_ev , lanczos_threshold, error=error) !Lanczos algorithm to calculate eigenvalue
       IF (unit_nr>0) WRITE(unit_nr,'(T2,A,2F16.8)') "smallest largest eigenvalue", min_ev, max_ev
       IF (unit_nr>0) WRITE(unit_nr,'(T2,A,2F16.8)') "requested interval", ev1, ev2
       interval_a=(max_ev-min_ev)*scale_evals/2
       interval_b=(max_ev+min_ev)/2
       
       sev1 = (ev1-interval_b)/interval_a !scaled ev1
       sev2 = (ev2-interval_b)/interval_a !scaled ev2
      
    
       !chebyshev domain,pi*sqrt(1-x^2) vector construction and chebyshev polynomials for integration (for g(E))
       ALLOCATE (E_inte(1:ninte+1))
       ALLOCATE (sqrt_vec(1:ninte+1))
       DO iinte = 1,ninte+1
          E_inte(iinte) = sev1 + ((sev2-sev1)/ninte)*(iinte-1)
          sqrt_vec(iinte) = pi*SQRT(1.0_dp -E_inte(iinte)*E_inte(iinte))
       END DO        
       
       !integral.. (identical to the coefficient for g(E))
       ALLOCATE (aitchev_T(1:ncheb)) !after intergral. =>ainte
               
       DO icheb = 1,ncheb
          CALL chebyshev_poly(initial,E_inte(1),icheb)
          CALL chebyshev_poly(final,E_inte(1),icheb)
          summa = (sev2-sev1)/(2.0_dp*ninte)*(initial/sqrt_vec(1)+final/sqrt_vec(ninte+1))
          DO iinte = 2,ninte
             CALL chebyshev_poly(chev_T,E_inte(iinte),icheb) 
             summa = summa + ((sev2-sev1)/ninte)*(chev_T/sqrt_vec(iinte))
          END DO
          aitchev_T(icheb) = summa
          summa = 0
       END DO

       ! scale the matrix to get evals in the interval -1,1
       CALL cp_dbcsr_add_on_diag(matrix_F,-interval_b,error=error)
       CALL cp_dbcsr_scale(matrix_F,1/interval_a,error=error)

       ! compute chebyshev matrix recursion
       CALL cp_dbcsr_get_info(matrix=matrix_F,nfullrows_total = Nrows) !get information about a matrix
       CALL cp_dbcsr_set(matrix_dummy1,0.0_dp,error=error) !empty matrix creation(for density matrix)
       CALL cp_dbcsr_set(matrix_dummy2,0.0_dp,error=error) !empty matrix creation(for density matrix)

       ALLOCATE(mu(1:ncheb))
       ALLOCATE (kernel_g(1:ncheb))
       CALL kernel(kernel_g(1), 1, ncheb)
       CALL kernel(kernel_g(2), 2, ncheb)       

       CALL cp_dbcsr_set(matrix_tmp1,0.0_dp,error=error) !matrix creation
       CALL cp_dbcsr_add_on_diag(matrix_tmp1,1.0_dp,error=error) !add a only number to diagonal elements
       CALL cp_dbcsr_trace(matrix_tmp1,trace=mu(1),error=error)
       CALL cp_dbcsr_copy(matrix_dummy1,matrix_tmp1,error=error)
       CALL cp_dbcsr_scale(matrix_dummy1,kernel_g(1)*aitchev_T(1),error=error) !first term of chebyshev poly(matrix)
       CALL cp_dbcsr_copy(matrix_tmp2,matrix_F,error=error) !make matrix_tmp2 = matrix_F
       CALL cp_dbcsr_trace(matrix_tmp2,trace=mu(2),error=error) 
       CALL cp_dbcsr_copy(matrix_dummy2,matrix_tmp2,error=error) !matrix_dummy2=
       CALL cp_dbcsr_scale(matrix_dummy2,2.0_dp*kernel_g(2)*aitchev_T(2),error=error) !second term of chebyshev poly(matrix)
     
       CALL cp_dbcsr_add(matrix_dummy2, matrix_dummy1, 1.0_dp, 1.0_dp, error=error)
       DO icheb=2,ncheb-1
          t1 = m_walltime()
          CALL kernel(kernel_g(icheb+1), icheb+1, ncheb)
          CALL cp_dbcsr_multiply("N", "N", 2.0_dp, matrix_F, matrix_tmp2, &
                                  -1.0_dp, matrix_tmp1, filter_eps=ls_scf_env%eps_filter,error=error) !matrix pultiplication(Recursion)
          CALL cp_dbcsr_filter(matrix_tmp1,ls_scf_env%eps_filter,error=error) !filtering
          CALL cp_dbcsr_copy(matrix_tmp3,matrix_tmp1,error=error)
          CALL cp_dbcsr_copy(matrix_tmp1,matrix_tmp2,error=error)
          CALL cp_dbcsr_copy(matrix_tmp2,matrix_tmp3,error=error)
          CALL cp_dbcsr_trace(matrix_tmp2,trace=mu(icheb+1),error=error) !icheb+1 th coefficient
          CALL cp_dbcsr_copy(matrix_dummy1,matrix_tmp2,error=error)
          CALL cp_dbcsr_scale(matrix_dummy1,2.0_dp*kernel_g(icheb+1)*aitchev_T(icheb+1),error=error) !second term of chebyshev poly(matrix)
          CALL cp_dbcsr_add(matrix_dummy2, matrix_dummy1, 1.0_dp, 1.0_dp, error=error)
          CALL cp_dbcsr_trace(matrix_dummy2,trace=trace_dm,error=error) !icheb+1 th coefficient
          occ=cp_dbcsr_get_occupation(matrix_tmp2)
          t2 = m_walltime()
          IF (unit_nr>0 .AND. MOD(icheb,20)==0) THEN
             WRITE(unit_nr,'(T2,A,I5,1X,A,1X,F8.3,1X,A,1X,F8.6,1X,A,1X,F16.8)') &
                  "Iter.",icheb,"time=",t2-t1,"occ=",occ,"trace=",trace_dm
          ENDIF
       ENDDO
 

       IF (.FALSE.) THEN
          ! we can check the difference between the density matrix obtained during the SCF and from the Chebyshev expansion
          ! create matrix_P=inv(sqrt(S))*G(H)*inv(sqrt(S))
          CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s_sqrt_inv, matrix_dummy2 , &
                                 0.0_dp, matrix_tmp1, filter_eps=ls_scf_env%eps_filter,error=error)
          CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, ls_scf_env%matrix_s_sqrt_inv, &
                                 0.0_dp, matrix_tmp2, filter_eps=ls_scf_env%eps_filter,error=error)

          ! look at the difference with the density matrix from the ls routines
          CALL cp_dbcsr_add(matrix_tmp2, ls_scf_env%matrix_p(ispin) , 2.0_dp, -1.0_dp, error=error)
          frob_matrix=cp_dbcsr_frobenius_norm(matrix_tmp2)
          IF (unit_nr>0) WRITE(unit_nr,*) "Difference between Chebyshev DM and LS DM",frob_matrix
       ENDIF

    ENDDO

    ! Chebyshev expansion with calculated coefficient
    ! grid construction and rescaling (by J)
    unit_dos=cp_print_key_unit_nr(logger,ls_scf_env%chebyshev%print_key_dos,"",extension=".xy",&
                                    middle_name="DOS",log_filename=.FALSE.,error=error)
    IF (unit_dos>0) THEN
      ALLOCATE (dos(1:n_gridpoint_dos))
      ALLOCATE (gdensity(1:n_gridpoint_dos))
      ALLOCATE (chev_E(1:n_gridpoint_dos))
      ALLOCATE (chev_Es_dos(1:n_gridpoint_dos))        
      ! are these the right bounds, max_ev, min_ev or interval_a,interval_b ?
      DO igrid = 1,n_gridpoint_dos 
         chev_E(igrid) = interval_b-interval_a+(igrid-1)*2*interval_a/(n_gridpoint_dos-1)
         chev_Es_dos(igrid) = (chev_E(igrid) - interval_b)/interval_a
      END DO
      DO igrid= 1,n_gridpoint_dos
         dummy1 = 0.0_dp !summation of polynomials
         dummy2 = 0.0_dp !summation of polynomials
         DO icheb = 2, ncheb
            CALL chebyshev_poly(chev_T_dos,chev_Es_dos(igrid),icheb)
            dummy1 = dummy1 + kernel_g(icheb)*mu(icheb)*chev_T_dos
            dummy2 = dummy2 + kernel_g(icheb)*aitchev_T(icheb)*chev_T_dos
         END DO
         dos(igrid) = 1.0_dp/(interval_a*Nrows*&
                          (pi*SQRT(1.0_dp-chev_Es_dos(igrid)*chev_Es_dos(igrid))))*(kernel_g(1)*mu(1) +2.0_dp*dummy1)
         gdensity(igrid) =kernel_g(1)*aitchev_T(1) +2.0_dp*dummy2
         WRITE(unit_dos,'(3F16.8)') chev_E(igrid), dos(igrid), gdensity(igrid)
      END DO
      DEALLOCATE (chev_Es_dos, chev_E, dos, gdensity)
    ENDIF
    CALL cp_print_key_finished_output(unit_dos,logger,ls_scf_env%chebyshev%print_key_dos,"",error=error)
     
    ! free the matrices
    CALL cp_dbcsr_release(matrix_tmp1,error=error)
    CALL cp_dbcsr_release(matrix_tmp2,error=error)
    CALL cp_dbcsr_release(matrix_tmp3,error=error)
    CALL cp_dbcsr_release(matrix_F,error=error)
    CALL cp_dbcsr_release(matrix_dummy1,error=error)
    CALL cp_dbcsr_release(matrix_dummy2,error=error)
     
    !Need deallocation 
    DEALLOCATE (mu,kernel_g, aitchev_T, E_inte,sqrt_vec)
   
    IF (unit_nr>0) WRITE(unit_nr,'(T2,A)') "ENDING CHEBYSHEV CALCULATION"

    CALL timestop(handle)

  END SUBROUTINE compute_chebyshev

END MODULE dm_ls_chebyshev
