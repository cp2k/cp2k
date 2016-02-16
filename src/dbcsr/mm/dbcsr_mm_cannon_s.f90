! *****************************************************************************
!> \brief Calculates norms of the entire matrix with minimal overhead.
!> \param norms ...
!> \param nrows ...
!> \param row_p ...
!> \param col_i ...
!> \param blk_p ...
!> \param rbs ...
!> \param cbs ...
!> \param DATA ...
!> \param local ...
!> \param local2global ...
! *****************************************************************************
  SUBROUTINE calc_norms_s(norms, nrows,&
       row_p, col_i, blk_p, rbs, cbs, DATA, local, local2global)
    REAL(kind=sp), DIMENSION(:), INTENT(OUT) :: norms
    INTEGER, INTENT(IN)                      :: nrows
    INTEGER, DIMENSION(1:nrows+1), &
      INTENT(IN)                             :: row_p
    INTEGER, DIMENSION(*), INTENT(IN)        :: col_i, blk_p, rbs, cbs
    REAL(kind=real_4), DIMENSION(*), &
      INTENT(IN)                             :: DATA
    LOGICAL, INTENT(IN)                      :: local
    INTEGER, DIMENSION(*), INTENT(IN)        :: local2global

    INTEGER                                  :: blk, bp, bpe, row, row_i, &
                                                row_size
    REAL(kind=sp)                            :: val

!   ---------------------------------------------------------------------------

    !$omp parallel default(none) &
    !$omp          private (row_i, row, row_size, blk, bp, bpe, val) &
    !$omp          shared (nrows, local) &
    !$omp          shared (local2global, rbs, cbs, row_p, col_i, blk_p, &
    !$omp                  data, norms)
    !$omp do
    DO row_i = 1, nrows
       IF (local) THEN
          row = local2global(row_i)
       ELSE
          row = row_i
       ENDIF
       row_size = rbs(row)
       DO blk = row_p(row_i)+1, row_p(row_i+1)
          IF (blk_p(blk) .NE. 0) THEN
             bp = ABS(blk_p(blk))
             bpe = bp + row_size * cbs(col_i(blk)) - 1
             val = SQRT (REAL (SUM(ABS(DATA(bp:bpe))**2), KIND=sp))
          ELSE
             val = 0.0_sp
          ENDIF
          norms(blk) = val
       ENDDO
    ENDDO
    !$omp end do
    !$omp end parallel
  END SUBROUTINE calc_norms_s

! *****************************************************************************
!> \brief Calculates norms of the entire matrix with minimal overhead.
!> \param norms ...
!> \param nblks ...
!> \param blki ...
!> \param rbs ...
!> \param cbs ...
!> \param DATA ...
!> \param local ...
!> \param local2global_rows ...
!> \param local2global_cols ...
! *****************************************************************************
  SUBROUTINE calc_norms_list_s(norms, nblks,&
       blki, rbs, cbs, DATA, local, local2global_rows, local2global_cols)
    REAL(kind=sp), DIMENSION(:), INTENT(OUT) :: norms
    INTEGER, INTENT(IN)                      :: nblks
    INTEGER, DIMENSION(3,nblks), INTENT(IN)  :: blki
    INTEGER, DIMENSION(:), INTENT(IN)        :: rbs, cbs
    REAL(kind=real_4), DIMENSION(:), &
      INTENT(IN)                             :: DATA
    LOGICAL, INTENT(IN)                      :: local
    INTEGER, DIMENSION(:), INTENT(IN)        :: local2global_rows
    INTEGER, DIMENSION(:), INTENT(IN)        :: local2global_cols

    INTEGER                                  :: blk, bp, bpe, row, col
    REAL(kind=sp)                            :: val

!   ---------------------------------------------------------------------------

    !$omp parallel default(none) &
    !$omp          private (row, col, blk, bp, bpe, val) &
    !$omp          shared (local, nblks) &
    !$omp          shared (rbs, cbs, blki, &
    !$omp                  data, norms, local2global_rows, local2global_cols)
    !$omp do
    DO blk = 1, nblks
       IF (blki(3,blk) .NE. 0) THEN
          bp = blki(3,blk)
          IF (local) THEN
             row = local2global_rows(blki(1,blk))
             col = local2global_cols(blki(2,blk))
          ELSE
             row = blki(1,blk)
             col = blki(2,blk)
          ENDIF
          bpe = bp + rbs(row) * cbs(col) - 1
          val = SQRT (REAL (SUM(ABS(DATA(bp:bpe))**2), KIND=sp))
       ELSE
          val = 0.0_sp
       ENDIF
       norms(blk) = val
    ENDDO
    !$omp end do
    !$omp end parallel
  END SUBROUTINE calc_norms_list_s

! *****************************************************************************
!> \brief Calculates max norms of each cluster with minimal overhead.
!> \param meta ...
!> \param data ...
!> \param refs ...
!> \param img_map_row ...
!> \param img_map_col ...
!> \param img_offset_rows ...
!> \param img_offset_cols ...
!> \param row_blk_size ...
!> \param col_blk_size ...
!> \param local_row ...
!> \param local_col ...
!> \param max_norms ...
!> \param is_left ...
! *****************************************************************************
  SUBROUTINE calc_max_image_norms_s(meta,DATA,&
     refs_size,refs_displ,&
     img_map_row,img_map_col,&
     img_offset_rows,img_offset_cols,&
     row_blk_size,col_blk_size,&
     local_row,local_col,&
     max_norms, is_left)
  INTEGER, DIMENSION(:), TARGET, INTENT(IN) :: meta
  REAL(kind=real_4), DIMENSION(:), &
       INTENT(IN)                         :: DATA
  INTEGER, DIMENSION(:, :, :), INTENT(IN), POINTER :: refs_size, refs_displ
  INTEGER, DIMENSION(:), INTENT(IN)       :: img_map_row, img_map_col,&
                                             img_offset_rows, img_offset_cols, &
                                             row_blk_size, col_blk_size, &
                                             local_row, local_col
  REAL(kind=sp), DIMENSION(:, :), INTENT(INOUT) :: max_norms
  LOGICAL                                 :: is_left

  INTEGER, DIMENSION(:), POINTER    :: row, col, bps
  INTEGER                           :: nblks, blk, bpe
  INTEGER, TARGET                   :: mi, ui
  INTEGER, POINTER                  :: rowi, coli

  !$omp parallel default(none) &
  !$omp          private (mi, ui, nblks, blk, row, col, bps, bpe,&
  !$omp                   rowi, coli) &
  !$omp          shared (max_norms, data, meta, refs_size, img_offset_rows,&
  !$omp                  img_map_row, img_map_col, refs_displ,&
  !$omp                  img_offset_cols, row_blk_size, col_blk_size,&
  !$omp                  local_row, local_col, is_left)
  IF (is_left) THEN
     rowi => mi
  ELSE
     coli => mi
  ENDIF
  !$omp do schedule(dynamic)
  DO ui = 1, SIZE(max_norms,2)
     IF (is_left) THEN
        coli => ui
     ELSE
        rowi => ui
     ENDIF
     DO mi = 1, SIZE(max_norms,1)
        max_norms(mi,ui) = 0.0_sp
        IF (refs_size(imeta_size,mi,ui).NE.0) THEN
           nblks = meta(refs_displ(imeta_displ,mi,ui)+dbcsr_slot_nblks)
           row => meta(meta(refs_displ(imeta_displ,mi,ui)+dbcsr_slot_coo_l)+&
                refs_displ(imeta_displ,mi,ui):&
                meta(refs_displ(imeta_displ,mi,ui)+dbcsr_num_slots)+&
                refs_displ(imeta_displ,mi,ui):3)
           col => meta(meta(refs_displ(imeta_displ,mi,ui)+dbcsr_slot_coo_l)+&
                refs_displ(imeta_displ,mi,ui)+1:&
                meta(refs_displ(imeta_displ,mi,ui)+dbcsr_num_slots)+&
                refs_displ(imeta_displ,mi,ui):3)
           bps => meta(meta(refs_displ(imeta_displ,mi,ui)+dbcsr_slot_coo_l)+&
                refs_displ(imeta_displ,mi,ui)+2:&
                meta(refs_displ(imeta_displ,mi,ui)+dbcsr_num_slots)+&
                refs_displ(imeta_displ,mi,ui):3)
           DO blk = 1, nblks
              IF (bps(blk).NE.0) THEN
                 bpe = bps(blk) + row_blk_size(local_row(img_map_row(row(blk)+img_offset_rows(rowi))))*&
                      col_blk_size(local_col(img_map_col(col(blk)+img_offset_cols(coli)))) - 1
                 max_norms(mi,ui) = MAX(max_norms(mi,ui),&
                      SQRT (REAL(SUM(ABS(DATA(bps(blk)+refs_displ(idata_displ,mi,ui):&
                                              bpe+refs_displ(idata_displ,mi,ui)))**2), KIND=sp)))
              ENDIF
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  !$omp end do
  !$omp end parallel
END SUBROUTINE calc_max_image_norms_s

! *****************************************************************************
!> \brief Calculates norms of each cluster with minimal overhead.
!> \param buffer ...
!> \param norms ..
! *****************************************************************************
  SUBROUTINE calc_image_norms_s(images,norms)
  TYPE(dbcsr_1d_array_type), INTENT(IN)    :: images
  REAL(kind=sp), DIMENSION(:, :), INTENT(INOUT) :: norms

  INTEGER, DIMENSION(:), POINTER    :: row, col, bps, rbs, cbs, &
                                       local_rows, local_cols
  REAL(kind=real_4), DIMENSION(:), POINTER    :: DATA
  INTEGER                           :: ui, blk, bpe

  !$omp parallel default(none) &
  !$omp private(ui,row,col,bps,blk,bpe,data,&
  !$omp         rbs,cbs,local_rows,local_cols) &
  !$omp shared(norms,images)
  !$omp do schedule(dynamic)
  DO ui=1,SIZE(norms,2)
     IF (images%mats(ui)%m%nblks.EQ.0) CYCLE
     row => images%mats(ui)%m%coo_l(1::3)
     col => images%mats(ui)%m%coo_l(2::3)
     bps => images%mats(ui)%m%coo_l(3::3)
     rbs => array_data(images%mats(ui)%m%row_blk_size)
     cbs => array_data(images%mats(ui)%m%col_blk_size)
     local_rows => array_data(images%mats(ui)%m%local_rows)
     local_cols => array_data(images%mats(ui)%m%local_cols)
     DATA => dbcsr_get_data_p_s (images%mats(ui)%m%data_area)
     DO blk = 1, images%mats(ui)%m%nblks
        IF (bps(blk).NE.0) THEN
           bpe = bps(blk) + rbs(local_rows(row(blk))) * cbs(local_cols(col(blk))) - 1
           norms(blk,ui) = SQRT (REAL (SUM(ABS(DATA(bps(blk):bpe))**2), KIND=sp))
        ELSE
           norms(blk,ui) = 0.0_sp
        ENDIF
     ENDDO
  ENDDO
  !$omp end do
  !$omp end parallel

END SUBROUTINE calc_image_norms_s

