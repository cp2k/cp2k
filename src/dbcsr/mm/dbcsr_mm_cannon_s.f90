! **************************************************************************************************
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
! **************************************************************************************************
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

! **************************************************************************************************
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
! **************************************************************************************************
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

! **************************************************************************************************
!> \brief Calculates max norms of each cluster with minimal overhead.
!> \param meta ...
!> \param data ...
!> \param refs ...
!> \param img_map ...
!> \param img_offset ...
!> \param row_blk_size ...
!> \param col_blk_size ...
!> \param local_row ...
!> \param local_col ...
!> \param max_norms ...
!> \param is_left ...
!> \param is_diagonal ...
! **************************************************************************************************
  SUBROUTINE calc_max_image_norms_s(meta,data,&
     refs_size,refs_displ,&
     img_map,&
     img_offset,&
     row_blk_size,col_blk_size,&
     local_row,local_col,&
     slot_coo_l,&
     max_norms, is_left, off_diagonal)
  INTEGER, DIMENSION(:), TARGET, INTENT(IN) :: meta
  REAL(kind=real_4), DIMENSION(:), &
       INTENT(IN)                         :: data
  INTEGER, DIMENSION(:, :), &
       INTENT(IN)                         :: refs_size
  INTEGER, DIMENSION(:, :, :), &
       POINTER, INTENT(IN)                :: refs_displ
  INTEGER, DIMENSION(:), INTENT(IN)       :: img_map, img_offset, &
                                             row_blk_size, col_blk_size, &
                                             local_row, local_col
  INTEGER, INTENT(IN)                     :: slot_coo_l
  REAL(kind=sp), DIMENSION(:), INTENT(INOUT) :: max_norms
  LOGICAL, INTENT(IN)                        :: is_left, off_diagonal

  INTEGER, DIMENSION(:), POINTER    :: row, col, bps
  INTEGER                           :: nblks, blk, bpe, icluster
  INTEGER, TARGET                   :: mi, ui
  INTEGER, POINTER                  :: rowi, coli
  REAL(kind=sp)                     :: max_norm

  icluster = 1
  max_norms(:) = 0.0_sp
  !
  !$omp parallel default(none) &
  !$omp          private (mi, ui, nblks, blk, row, col, bps, bpe,&
  !$omp                   rowi, coli, max_norm) &
  !$omp          shared (max_norms, data, meta, refs_size, img_offset,&
  !$omp                  img_map, refs_displ, slot_coo_l,&
  !$omp                  row_blk_size, col_blk_size, icluster,&
  !$omp                  local_row, local_col, is_left, off_diagonal)
  IF (is_left) THEN
     rowi => mi
     coli => ui
  ELSE
     coli => mi
     rowi => ui
  ENDIF
  DO ui = 1, SIZE(refs_size,2)
     DO mi = 1, SIZE(refs_size,1)
        IF (refs_size(mi,ui).NE.0) THEN
           IF (off_diagonal.OR.ui.NE.mi) THEN
              max_norm = 0.0_sp
              nblks = meta(refs_displ(imeta_displ,mi,ui)+dbcsr_slot_nblks)
              row => meta(refs_displ(imeta_displ,mi,ui)+slot_coo_l:&
                   meta(refs_displ(imeta_displ,mi,ui)+dbcsr_slot_size)+&
                   refs_displ(imeta_displ,mi,ui):3)
              col => meta(refs_displ(imeta_displ,mi,ui)+slot_coo_l+1:&
                   meta(refs_displ(imeta_displ,mi,ui)+dbcsr_slot_size)+&
                   refs_displ(imeta_displ,mi,ui):3)
              bps => meta(refs_displ(imeta_displ,mi,ui)+slot_coo_l+2:&
                   meta(refs_displ(imeta_displ,mi,ui)+dbcsr_slot_size)+&
                   refs_displ(imeta_displ,mi,ui):3)
              !$omp do
              DO blk = 1, nblks
                 IF (bps(blk).NE.0) THEN
                    IF (is_left) THEN
                       bpe = bps(blk) + row_blk_size(local_row(row(blk)))*&
                            col_blk_size(local_col(img_map(col(blk)+img_offset(coli)))) - 1
                    ELSE
                       bpe = bps(blk) + row_blk_size(local_row(img_map(row(blk)+img_offset(rowi))))*&
                            col_blk_size(local_col(col(blk))) - 1
                    ENDIF
                    max_norm = MAX(max_norm,&
                         SQRT (REAL(SUM(ABS(data(bps(blk)+refs_displ(idata_displ,mi,ui):&
                                                 bpe+refs_displ(idata_displ,mi,ui)))**2), KIND=sp)))
                 ENDIF
              ENDDO
              !$omp end do
              !$omp critical(cannon_max_norm)
              max_norms(icluster) = MAX(max_norms(icluster),max_norm)
              !$omp end critical(cannon_max_norm)
           ENDIF
           !$omp barrier
           !$omp master
           icluster = icluster + 1
           !$omp end master
        ENDIF
     ENDDO
  ENDDO
  !$omp end parallel
END SUBROUTINE calc_max_image_norms_s

! **************************************************************************************************
!> \brief Calculates norms of each cluster with minimal overhead.
!> \param buffer ...
!> \param norms ..
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
     IF (images%mats(ui)%m%nblks.EQ.0) CYCLE
     row => images%mats(ui)%m%coo_l(1::3)
     col => images%mats(ui)%m%coo_l(2::3)
     bps => images%mats(ui)%m%coo_l(3::3)
     rbs => array_data(images%mats(ui)%m%row_blk_size)
     cbs => array_data(images%mats(ui)%m%col_blk_size)
     local_rows => array_data(images%mats(ui)%m%local_rows)
     local_cols => array_data(images%mats(ui)%m%local_cols)
     data => dbcsr_get_data_p_s (images%mats(ui)%m%data_area)
     !$omp do
     DO blk = 1, images%mats(ui)%m%nblks
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

