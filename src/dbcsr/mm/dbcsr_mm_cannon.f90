#:include '../data/dbcsr.fypp'
#:for n, nametype1, base1, prec1, kind1, type1, dkind1, normname1 in inst_params_float

! **************************************************************************************************
!> \brief Prepare buffers for multiplications
! **************************************************************************************************
  SUBROUTINE prepare_buffers_${nametype1}$(negate_real, negate_imaginary, &
                                           iter, row, col, blk, blk_p, bp, &
                                           row_size, col_size, nze, nsymmetries, symmetry_i, &
                                           stored_row, stored_col, tr_row_size, tr_col_size, tr, &
                                           row_img, col_img, nrow_images, ncol_images, &
                                           row_img_dist, col_img_dist, predist_type_fwd, blacs2mpi, &
                                           target_imgdist, prow, pcol, rowi, coli, &
                                           row_dist, col_dist, dst_p, sm_pos, myt_smp, metalen, &
                                           sd_pos, myt_sdp, send_meta, sd_disp, &
                                           data_area, send_data_area, scale_neg_one, scale_value)
    LOGICAL, INTENT(IN)                                     :: negate_real, negate_imaginary
    TYPE(dbcsr_iterator), INTENT(INOUT)                     :: iter
    INTEGER, INTENT(INOUT)                                  :: row, col, blk, blk_p, row_size, col_size, &
                                                               nze, bp, symmetry_i, &
                                                               stored_row, stored_col, tr_row_size, tr_col_size, &
                                                               row_img, col_img, prow, pcol, rowi, coli, &
                                                               dst_p, sm_pos, sd_pos
    INTEGER, INTENT(IN)                                     :: nsymmetries, nrow_images, ncol_images, metalen
    LOGICAL, INTENT(INOUT)                                  :: tr
    INTEGER, DIMENSION(:), INTENT(IN), POINTER              :: row_img_dist, col_img_dist, row_dist, col_dist
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(IN)          :: sd_disp
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT)       :: myt_smp, myt_sdp, send_meta
    TYPE(dbcsr_imagedistribution_obj), INTENT(IN)           :: target_imgdist
    INTEGER, DIMENSION(:, :), INTENT(IN), POINTER           :: blacs2mpi
    CHARACTER, INTENT(IN)                                   :: predist_type_fwd
    ${type1}$, DIMENSION(:), INTENT(IN)                     :: data_area
    TYPE(dbcsr_data_obj), INTENT(INOUT)                     :: send_data_area
    TYPE(dbcsr_scalar_type), INTENT(IN)                     :: scale_neg_one
    TYPE(dbcsr_scalar_type), INTENT(IN), OPTIONAL           :: scale_value

    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, blk, blk_p=blk_p, &
                                      row_size=row_size, col_size=col_size)
       nze = row_size*col_size
       IF (nze .EQ. 0) CYCLE
       bp = ABS(blk_p)
       DO symmetry_i = 1, nsymmetries
          IF (symmetry_i .EQ. 1) THEN
             stored_row = row; stored_col = col; tr = blk_p .LT. 0
             tr_row_size = col_size; tr_col_size = row_size
          ELSE
             IF (row .EQ. col) CYCLE
             stored_row = col; stored_col = row; tr = blk_p .GT. 0
             tr_row_size = row_size; tr_col_size = col_size
          ENDIF
          ! Where do we send this block?
          row_img = 1
          IF (nrow_images .GT. 1) row_img = row_img_dist(stored_row)
          col_img = 1
          IF (ncol_images .GT. 1) col_img = col_img_dist(stored_col)
          CALL image_calculator(target_imgdist, &
                                prow=prow, rowi=rowi, &
                                pcol=pcol, coli=coli, &
                                myprow=row_dist(stored_row), myrowi=row_img, &
                                mypcol=col_dist(stored_col), mycoli=col_img, &
                                shifting=predist_type_fwd)
          dst_p = blacs2mpi(prow, pcol)
          sm_pos = myt_smp(dst_p)
          myt_smp(dst_p) = myt_smp(dst_p)+metalen
          sd_pos = myt_sdp(dst_p)
          myt_sdp(dst_p) = myt_sdp(dst_p)+nze
          IF (tr) THEN
             IF (PRESENT(scale_value)) THEN
                CALL dbcsr_block_transpose(send_data_area%d%${base1}$_${prec1}$(sd_pos:myt_sdp(dst_p)-1), &
                                           data_area(bp:bp+nze-1)*scale_value%${base1}$_${prec1}$, &
                                           tr_row_size, tr_col_size)
             ELSE
                CALL dbcsr_block_transpose(send_data_area%d%${base1}$_${prec1}$(sd_pos:myt_sdp(dst_p)-1), &
                                           data_area(bp:bp+nze-1), &
                                           tr_row_size, tr_col_size)
             ENDIF
             IF (negate_real .AND. negate_imaginary) THEN
                send_data_area%d%${base1}$_${prec1}$(sd_pos:myt_sdp(dst_p)-1) = & 
                     send_data_area%d%${base1}$_${prec1}$(sd_pos:myt_sdp(dst_p)-1)*scale_neg_one%${base1}$_${prec1}$
             ELSEIF (negate_real) THEN
#:if nametype1=="s" or nametype1=="d"
                send_data_area%d%${base1}$_${prec1}$(sd_pos:myt_sdp(dst_p)-1) = &
                     -send_data_area%d%${base1}$_${prec1}$(sd_pos:myt_sdp(dst_p)-1)
#:else
                send_data_area%d%${base1}$_${prec1}$(sd_pos:myt_sdp(dst_p)-1) = &
                     CMPLX( &
                     -REAL(send_data_area%d%${base1}$_${prec1}$(sd_pos:myt_sdp(dst_p)-1), KIND=${kind1}$), &
                     AIMAG(send_data_area%d%${base1}$_${prec1}$(sd_pos:myt_sdp(dst_p)-1)), &
                     KIND=${kind1}$)
#:endif
             ELSEIF (negate_imaginary) THEN
#:if nametype1=="c" or nametype1=="z"
                send_data_area%d%${base1}$_${prec1}$(sd_pos:myt_sdp(dst_p)-1) = &
                     CONJG(send_data_area%d%${base1}$_${prec1}$(sd_pos:myt_sdp(dst_p)-1))
#:endif
             ENDIF
          ELSE
             ! Copy the block
             IF (PRESENT(scale_value)) THEN
                send_data_area%d%${base1}$_${prec1}$(sd_pos:myt_sdp(dst_p)-1) = &
                     data_area(bp:bp+nze-1)*scale_value%${base1}$_${prec1}$
             ELSE
                send_data_area%d%${base1}$_${prec1}$(sd_pos:myt_sdp(dst_p)-1) = data_area(bp:bp+nze-1)
             ENDIF
          END IF

          send_meta(sm_pos) = stored_row
          send_meta(sm_pos+1) = stored_col
          send_meta(sm_pos+2) = sd_pos-sd_disp(dst_p)+1
          send_meta(sm_pos+3) = rowi
          send_meta(sm_pos+4) = coli
       ENDDO
    ENDDO
  END SUBROUTINE prepare_buffers_${nametype1}$
#:endfor
