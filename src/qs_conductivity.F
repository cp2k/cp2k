!   CP2K: A general program to perform molecular dynamics simulations         !
!   copyright 2003 the cp2k developers group !

!!****** cp2k/qs_conductivity [1.0] *
!!
!!   NAME
!!     qs_conductivity
!!   CVS id
!!   $Id$
!!
!!   FUNCTION
!! 
!!     should be able to compute the optical conductivity
!!
!!   AUTHORS
!!     Joost VandeVondele (04.2003)
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE qs_conductivity
! *****************************************************************************
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_p_type,&
                                             cp_fm_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE message_passing,                 ONLY: mp_sum
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE

! *** Public ***
  PUBLIC :: optical_conductivity 
! *****************************************************************************

CONTAINS

SUBROUTINE optical_conductivity(matrix_s,occupied_orbs,occupied_evals, &
                                unoccupied_orbs,unoccupied_evals,volume,ionode,scr)
   TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: matrix_s
   TYPE(cp_fm_type), POINTER              :: occupied_orbs, unoccupied_orbs
   REAL(KIND = dp), DIMENSION(:), POINTER                 :: occupied_evals, unoccupied_evals
   REAL(KIND = dp)                                        :: volume
   LOGICAL                                         :: ionode
   INTEGER                                         :: scr


   TYPE(cp_fm_type), POINTER              :: buf,matrix_elements
   TYPE(cp_fm_struct_type), POINTER                :: me_struct,template
   INTEGER                                         :: component,nocc,nunocc,nao,i,j,nrow_local, &
                                                      ncol_local,index,handle

   ! bands can be used to split the spectrum divinding the occupied evals in nbands
   INTEGER, PARAMETER                              :: npoints=1000,nbands=4
   REAL(KIND = dp), DIMENSION(:,:), POINTER                 :: result
   REAL(KIND = dp)                                        :: emin,emax,strength,occ_min,occ_max, &
                                                      unocc_min,unocc_max,energy,dE_bin
   INTEGER, DIMENSION(:), POINTER                  :: row_indices,col_indices
   INTEGER                                         :: iband

   ! CALL stop_program("optical_conductivity","NYI")
   ! very much *not* tested and most likely to be wrong
   ! the core is present 
   CALL timeset("optical_conductivity","I","",handle)

   IF (SIZE(matrix_s).lt.4) CALL stop_program("optical_conductivity","where is the derivative ?")

   NULLIFY(buf,me_struct,template,result,row_indices,col_indices)
   ALLOCATE(result(0:npoints-1,nbands))
   ! 
   ! we need the matrix elements of the momentum operator.
   ! apart from minus signs and factors of i and a couple of units,
   ! these are the derivatives of S with respect to the atomic coordinates
   !
   CALL cp_fm_get_info(occupied_orbs,ncol_global=nocc,nrow_global=nao)
   CALL cp_fm_get_info(unoccupied_orbs,ncol_global=nunocc)

   template => unoccupied_orbs%matrix_struct
   CALL cp_fm_create(buf,template,"buffer")
   CALL cp_fm_struct_create(me_struct,para_env=template%para_env,context=template%context, &
                            nrow_global=nocc,ncol_global=nunocc)
   CALL cp_fm_create(matrix_elements,me_struct,"matrix_elements")
   CALL cp_fm_get_info(matrix_elements,nrow_local=nrow_local,ncol_local=ncol_local, &
                       row_indices=row_indices,col_indices=col_indices)
   CALL cp_fm_struct_release(me_struct)

   occ_min=MINVAL(occupied_evals)
   occ_max=MAXVAL(occupied_evals)
   unocc_min=MINVAL(unoccupied_evals)
   unocc_max=MAXVAL(unoccupied_evals)
   emin=(unocc_min-occ_max-1.0E-3)
   ! old strategy ....
   ! we retain only these transitions that have an energy in the range
   ! where all transitions are computed (i.e. as soon as we start neglecting HOMO -> unoccupied,
   ! we stop recording other transitions)
   emax=MAX((unocc_max-occ_max),emin+1.0E-3)
   ! the other option is to compute the full spectrum anyway. 
   ! the advantage of this option is that it allows for alligning the spectrum with
   ! experiment, using features from low occupied bands.
   emax=MAX((unocc_max-occ_min),emin)+1.0E-3
   dE_bin=(emax-emin)/npoints
   DO component=1,3
      result=0.0_dp
      call cp_sm_fm_multiply(matrix_s(component+1)%matrix,unoccupied_orbs,buf,nunocc,template%para_env)
      call cp_fm_gemm("T","N",nocc,nunocc,nao,1.0_dp,occupied_orbs,buf,0.0_dp,matrix_elements)
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
              result(index,iband)=result(index,iband)+strength/dE_bin
           ENDIF

        ENDDO
      ENDDO
      CALL mp_sum(result,template%para_env%group)
      IF (ionode) THEN
      DO iband=1,nbands
      DO i=0,npoints-1
         write(scr,'(2E20.10)') emin+(REAL(i,dp)/npoints)*(emax-emin),result(i,iband)
      ENDDO
      ENDDO
      ENDIF
   ENDDO

   CALL cp_fm_release(buf)
   CALL cp_fm_release(matrix_elements)
   DEALLOCATE(result)

   CALL timestop(0.0_dp,handle)

END SUBROUTINE

END MODULE qs_conductivity
