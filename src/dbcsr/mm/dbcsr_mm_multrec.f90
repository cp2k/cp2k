#:include '../data/dbcsr.fypp'
#:for n, nametype1, base1, prec1, kind1, type1, dkind1, normname1 in inst_params_float

! **************************************************************************************************
!> \brief Applying in-place filtering on the workspace.
!> \brief Use Frobenius norm
! **************************************************************************************************
  SUBROUTINE multrec_filtering_${nametype1}$(filter_eps, nblks, rowi, coli, blkp, &
                                             rbs, cbs, nze, DATA)
    REAL(kind=real_8), INTENT(IN)              :: filter_eps
    INTEGER, INTENT(INOUT)                     :: nblks, nze
    INTEGER, DIMENSION(1:nblks), INTENT(INOUT) :: rowi, coli, blkp
    INTEGER, DIMENSION(:), INTENT(IN)          :: rbs, cbs
    ${type1}$, DIMENSION(:), &
      INTENT(INOUT)                            :: DATA

    INTEGER                                    :: blk, lastblk, blk_nze, blk_p
    REAL(kind=real_8)                          :: nrm

    REAL(KIND=real_8), EXTERNAL                :: DZNRM2, DDOT
#if defined (__ACCELERATE)
    REAL(KIND=real_8), EXTERNAL                :: SCNRM2, SDOT
#else
    REAL(KIND=real_4), EXTERNAL                :: SCNRM2, SDOT
#endif

    lastblk = 0
    nze = 0
    !
    DO blk = 1, nblks
       blk_p = blkp(blk)
       IF (blk_p .EQ. 0) CYCLE
       blk_nze = rbs(rowi(blk)) * cbs(coli(blk))
       IF (blk_nze .EQ. 0) CYCLE ! Skip empty blocks
       nrm = REAL(${normname1}$(blk_nze, data(blk_p), 1, data(blk_p), 1)), KIND=real_8)
       IF (nrm .GE. filter_eps) THEN
          ! Keep block
          lastblk = lastblk+1
          IF (lastblk .LT. blk) THEN
             rowi(lastblk) = rowi(blk)
             coli(lastblk) = coli(blk)
             blkp(lastblk) = blkp(blk)
          ENDIF
          nze = nze + blk_nze
       ENDIF
    ENDDO
    ! 
    nblks = lastblk

  END SUBROUTINE multrec_filtering_${nametype1}$
#:endfor
