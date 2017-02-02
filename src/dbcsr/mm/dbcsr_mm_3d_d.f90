! **************************************************************************************************
!> \brief Calculates max norms of each image with minimal overhead.
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
! **************************************************************************************************
  SUBROUTINE calc_max_image_norms_d(meta,data,&
     refs_size,refs_displ,&
     img_map,&
     img_offset,&
     row_blk_size,col_blk_size,&
     local_row,local_col,&
     slot_coo_l,&
     max_norms, is_left)
  INTEGER, DIMENSION(:), TARGET, INTENT(IN) :: meta
  REAL(kind=real_8), DIMENSION(:), TARGET, &
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
  LOGICAL, INTENT(IN)                        :: is_left

  INTEGER, DIMENSION(:), POINTER    :: row, col, bps
  INTEGER                           :: nblks, blk, bpe, iimage
  INTEGER, TARGET                   :: mi, ui
  INTEGER, POINTER                  :: rowi, coli
  REAL(kind=sp)                     :: max_norm

  iimage = 1
  max_norms(:) = 0.0_sp
  !
  !$omp parallel default(none) &
  !$omp          private (mi, ui, nblks, blk, row, col, bps, bpe,&
  !$omp                   rowi, coli, max_norm) &
  !$omp          shared (max_norms, data, meta, refs_size, img_offset,&
  !$omp                  img_map, refs_displ, slot_coo_l,&
  !$omp                  row_blk_size, col_blk_size, iimage,&
  !$omp                  local_row, local_col, is_left)
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
           max_norm = 0.0_sp
           nblks = meta(refs_displ(imeta,mi,ui)+dbcsr_slot_nblks)
           row => meta(refs_displ(imeta,mi,ui)+slot_coo_l:&
                       meta(refs_displ(imeta,mi,ui)+dbcsr_slot_size)+&
                            refs_displ(imeta,mi,ui):3)
           col => meta(refs_displ(imeta,mi,ui)+slot_coo_l+1:&
                       meta(refs_displ(imeta,mi,ui)+dbcsr_slot_size)+&
                            refs_displ(imeta,mi,ui):3)
           bps => meta(refs_displ(imeta,mi,ui)+slot_coo_l+2:&
                       meta(refs_displ(imeta,mi,ui)+dbcsr_slot_size)+&
                            refs_displ(imeta,mi,ui):3)
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
                      SQRT (REAL(SUM(ABS(data(bps(blk)+refs_displ(idata,mi,ui):&
                                              bpe+refs_displ(idata,mi,ui)))**2), KIND=sp)))
              ENDIF
           ENDDO
           !$omp end do
           !$omp critical(cannon_max_norm)
           max_norms(iimage) = MAX(max_norms(iimage),max_norm)
           !$omp end critical(cannon_max_norm)
           !$omp barrier
           !$omp master
           iimage = iimage + 1
           !$omp end master
        ELSE IF (SIZE(refs_size,1) .EQ. 1) THEN
           !$omp barrier
           !$omp master
           max_norms(iimage) = 0.0_sp
           iimage = iimage + 1
           !$omp end master
        ENDIF
     ENDDO
  ENDDO
  !$omp end parallel
END SUBROUTINE calc_max_image_norms_d

! **************************************************************************************************
!> \brief Calculates norms of each image with minimal overhead.
!> \param buffer ...
!> \param norms ...
! **************************************************************************************************
  SUBROUTINE calc_image_norms_d(images,norms,uf,ul)
  TYPE(dbcsr_1d_array_type), INTENT(IN)    :: images
  REAL(kind=sp), DIMENSION(:, :), INTENT(INOUT) :: norms
  INTEGER, INTENT(IN)                      :: uf, ul

  INTEGER, DIMENSION(:), POINTER    :: row, col, bps, rbs, cbs, &
                                       local_rows, local_cols
  REAL(kind=real_8), DIMENSION(:), POINTER    :: data
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
     data => dbcsr_get_data_p_d (images%mats(ui)%data_area)
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

END SUBROUTINE calc_image_norms_d

