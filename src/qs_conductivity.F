!   CP2K: A general program to perform molecular dynamics simulations         !
!   copyright 2003 the cp2k developers group !

!!****** cp2k/qs_conductivity [1.0] *
!!
!!   NAME
!!     qs_conductivity
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
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_sm_fm_multiply
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_element,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_full_matrix_p_type,&
                                             cp_full_matrix_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
  USE kinds,                           ONLY: wp => dp
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
   TYPE(cp_full_matrix_type), POINTER              :: occupied_orbs, unoccupied_orbs
   REAL(wp), DIMENSION(:), POINTER                 :: occupied_evals, unoccupied_evals
   REAL(wp)                                        :: volume
   LOGICAL                                         :: ionode
   INTEGER                                         :: scr


   TYPE(cp_full_matrix_type), POINTER              :: buf,matrix_elements
   TYPE(cp_fm_struct_type), POINTER                :: me_struct,template
   INTEGER                                         :: component,nocc,nunocc,nao,i,j,nrow_local, &
                                                      ncol_local,index,handle

   INTEGER, PARAMETER                              :: npoints=1000
   REAL(wp), DIMENSION(:), POINTER                 :: result
   REAL(wp)                                        :: emin,emax,strength,occ_max,unocc_min,unocc_max,energy
   INTEGER, DIMENSION(:), POINTER                  :: row_indices,col_indices

   ! CALL stop_program("optical_conductivity","NYI")
   ! very much *not* tested and most likely to be wrong
   ! the core is present 
   CALL timeset("optical_conductivity","I","",handle)

   NULLIFY(buf,me_struct,template,result,row_indices,col_indices)
   ALLOCATE(result(0:npoints-1))
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

   occ_max=MAXVAL(occupied_evals)
   unocc_max=MAXVAL(unoccupied_evals)
   unocc_min=MINVAL(unoccupied_evals)
   ! we retain only these transitions that have an energy in the range
   ! where all transitions are computed (i.e. as soon as we start neglecting HOMO -> unoccupied,
   ! we stop recording other transitions)
   emin=(unocc_min-occ_max)
   emax=MAX((unocc_max-occ_max),emin+1.0E-16)
   DO component=1,3
      result=0.0_wp
      call cp_sm_fm_multiply(matrix_s(component+1)%matrix,unoccupied_orbs,buf,nunocc,template%para_env)
      call cp_fm_gemm("T","N",nocc,nunocc,nao,1.0_wp,occupied_orbs,buf,0.0_wp,matrix_elements)
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
              result(index)=result(index)+strength
           ENDIF

        ENDDO
      ENDDO
      CALL mp_sum(result,template%para_env%group)
      IF (ionode) THEN
      DO i=0,npoints-1
         write(scr,'(2E20.10)') emin+(REAL(i,wp)/npoints)*(emax-emin),result(i)
      ENDDO
      ENDIF
   ENDDO

   CALL cp_fm_release(buf)
   CALL cp_fm_release(matrix_elements)
   DEALLOCATE(result)

   CALL timestop(0.0_wp,handle)

END SUBROUTINE

END MODULE qs_conductivity
