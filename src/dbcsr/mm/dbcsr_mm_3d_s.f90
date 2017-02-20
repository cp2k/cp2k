! **************************************************************************************************
!> \brief Calculates norms of each image with minimal overhead.
!> \param buffer ...
!> \param norms ...
! **************************************************************************************************
  SUBROUTINE calc_image_norms_s(images,norms,uf,ul)
  TYPE(dbcsr_1d_array_type), INTENT(IN)    :: images
  REAL(kind=sp), DIMENSION(:, :), INTENT(INOUT) :: norms
  INTEGER, INTENT(IN)                      :: uf, ul

  INTEGER, DIMENSION(:), POINTER    :: row, col, bps, rbs, cbs, &
                                       local_rows, local_cols
  REAL(kind=real_4), DIMENSION(:), POINTER    :: data
  INTEGER                           :: ui, blk, bpe

  !$omp parallel default(none) &
  !$omp private(ui,row,col,bps,blk,bpe,data,&
  !$omp         rbs,cbs,local_rows,local_cols) &
  !$omp shared(norms,images,uf,ul)
  DO ui=uf,ul
     IF (images%mats(ui)%nblks.EQ.0) CYCLE
     row => images%mats(ui)%coo_l(1::3)
     col => images%mats(ui)%coo_l(2::3)
     bps => images%mats(ui)%coo_l(3::3)
     rbs => array_data(images%mats(ui)%row_blk_size)
     cbs => array_data(images%mats(ui)%col_blk_size)
     local_rows => array_data(images%mats(ui)%local_rows)
     local_cols => array_data(images%mats(ui)%local_cols)
     data => dbcsr_get_data_p_s (images%mats(ui)%data_area)
     !$omp do
     DO blk = 1, images%mats(ui)%nblks
        IF (bps(blk).NE.0) THEN
           bpe = bps(blk) + rbs(local_rows(row(blk))) * cbs(local_cols(col(blk))) - 1
           norms(blk,ui) = SQRT (REAL (SUM(ABS(data(bps(blk):bpe))**2), KIND=sp))
        ELSE
           norms(blk,ui) = 0.0_sp
        ENDIF
     ENDDO
     !$omp end do
  ENDDO
  !$omp end parallel

END SUBROUTINE calc_image_norms_s

