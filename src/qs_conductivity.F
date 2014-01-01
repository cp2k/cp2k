!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief should be able to compute the optical conductivity
!> \author Joost VandeVondele (04.2003)
! *****************************************************************************
MODULE qs_conductivity
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE message_passing,                 ONLY: mp_sum
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_conductivity'

! *** Public ***
  PUBLIC :: optical_conductivity

CONTAINS

! *****************************************************************************
SUBROUTINE optical_conductivity(matrix_s,occupied_orbs,occupied_evals, &
                                unoccupied_orbs,unoccupied_evals,volume,scr,error)
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_fm_type), POINTER                :: occupied_orbs
    REAL(KIND=dp), DIMENSION(:), POINTER     :: occupied_evals
    TYPE(cp_fm_type), POINTER                :: unoccupied_orbs
    REAL(KIND=dp), DIMENSION(:), POINTER     :: unoccupied_evals
    REAL(KIND=dp)                            :: volume
    INTEGER                                  :: scr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'optical_conductivity', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nbands = 4, npoints = 1000

    INTEGER                                  :: component, handle, i, iband, &
                                                index, j, nao, ncol_local, &
                                                nocc, nrow_local, nunocc
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    REAL(KIND=dp)                            :: dE_bin, emax, emin, energy, &
                                                occ_max, occ_min, strength, &
                                                unocc_max, unocc_min
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: RESULT
    TYPE(cp_fm_struct_type), POINTER         :: me_struct, template
    TYPE(cp_fm_type), POINTER                :: buf, matrix_elements

! bands can be used to split the spectrum divinding the occupied evals in nbands
! CALL stop_program(routineN,moduleN,__LINE__,"NYI")
! very much *not* tested and most likely to be wrong
! the core is present

   CALL timeset(routineN,handle)

   IF (SIZE(matrix_s).lt.4) CALL stop_program(routineN,moduleN,__LINE__,&
                                              "Where is the derivative?")

   NULLIFY(buf,me_struct,template,RESULT,row_indices,col_indices)
   ALLOCATE(RESULT(0:npoints-1,nbands))
   !
   ! we need the matrix elements of the momentum operator.
   ! apart from minus signs and factors of i and a couple of units,
   ! these are the derivatives of S with respect to the atomic coordinates
   !
   CALL cp_fm_get_info(occupied_orbs,ncol_global=nocc,nrow_global=nao,error=error)
   CALL cp_fm_get_info(unoccupied_orbs,ncol_global=nunocc,error=error)

   template => unoccupied_orbs%matrix_struct
   CALL cp_fm_create(buf,template,"buffer",error=error)
   CALL cp_fm_struct_create(me_struct,para_env=template%para_env,context=template%context, &
                            nrow_global=nocc,ncol_global=nunocc,error=error)
   CALL cp_fm_create(matrix_elements,me_struct,"matrix_elements",error=error)
   CALL cp_fm_get_info(matrix_elements,nrow_local=nrow_local,ncol_local=ncol_local, &
                       row_indices=row_indices,col_indices=col_indices,error=error)
   CALL cp_fm_struct_release(me_struct,error=error)

   occ_min=MINVAL(occupied_evals)
   occ_max=MAXVAL(occupied_evals)
   unocc_min=MINVAL(unoccupied_evals)
   unocc_max=MAXVAL(unoccupied_evals)
   emin=(unocc_min-occ_max-1.0E-3_dp)
   ! old strategy ....
   ! we retain only these transitions that have an energy in the range
   ! where all transitions are computed (i.e. as soon as we start neglecting HOMO -> unoccupied,
   ! we stop recording other transitions)
   emax=MAX((unocc_max-occ_max),emin+1.0E-3_dp)
   ! the other option is to compute the full spectrum anyway.
   ! the advantage of this option is that it allows for alligning the spectrum with
   ! experiment, using features from low occupied bands.
   emax=MAX((unocc_max-occ_min),emin)+1.0E-3_dp
   dE_bin=(emax-emin)/npoints
   DO component=1,3
      RESULT=0.0_dp
      CALL cp_dbcsr_sm_fm_multiply(matrix_s(component+1)%matrix,unoccupied_orbs,buf,nunocc,error=error)
      CALL cp_fm_gemm("T","N",nocc,nunocc,nao,1.0_dp,occupied_orbs,buf,0.0_dp,matrix_elements,error=error)
      DO j=1,ncol_local
        DO i=1,nrow_local

           ! the energy of the transition
           energy=-occupied_evals(row_indices(i))+unoccupied_evals(col_indices(j))

           ! this should be the strength of the transition.
           strength=Pi*(matrix_elements%local_data(i,j))**2 / (volume*energy)
           ! write(6,*) i,j,matrix_elements%local_data(i,j)

           ! put the transitions in the correct bin
           index=FLOOR(((energy-emin)/(emax-emin))*npoints)
           IF (index.GE.0 .AND. index.LT.npoints) THEN
              iband=1+(row_indices(i)*nbands-1)/SIZE(occupied_evals)
              RESULT(index,iband)=RESULT(index,iband)+strength/dE_bin
           ENDIF

        ENDDO
      ENDDO
      CALL mp_sum(RESULT,template%para_env%group)
      IF (scr>0) THEN
         DO iband=1,nbands
            DO i=0,npoints-1
               WRITE(scr,'(2E20.10)') emin+(REAL(i,dp)/npoints)*(emax-emin),RESULT(i,iband)
            ENDDO
         ENDDO
      ENDIF
   ENDDO

   CALL cp_fm_release(buf,error=error)
   CALL cp_fm_release(matrix_elements,error=error)
   DEALLOCATE(RESULT)

   CALL timestop(handle)

END SUBROUTINE

END MODULE qs_conductivity
