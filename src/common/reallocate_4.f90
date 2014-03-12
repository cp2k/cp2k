    INTEGER :: istat, lb1, lb1_old, lb2, lb2_old, lb3, lb3_old, lb4, lb4_old, &
      ub1, ub1_old, ub2, ub2_old, ub3, ub3_old, ub4, ub4_old

    IF (ASSOCIATED(p)) THEN
       lb1_old = LBOUND(p,1)
       ub1_old = UBOUND(p,1)
       lb2_old = LBOUND(p,2)
       ub2_old = UBOUND(p,2)
       lb3_old = LBOUND(p,3)
       ub3_old = UBOUND(p,3)
       lb4_old = LBOUND(p,4)
       ub4_old = UBOUND(p,4)
       lb1 = MAX(lb1_new,lb1_old)
       ub1 = MIN(ub1_new,ub1_old)
       lb2 = MAX(lb2_new,lb2_old)
       ub2 = MIN(ub2_new,ub2_old)
       lb3 = MAX(lb3_new,lb3_old)
       ub3 = MIN(ub3_new,ub3_old)
       lb4 = MAX(lb4_new,lb4_old)
       ub4 = MIN(ub4_new,ub4_old)
       ALLOCATE (work(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4),STAT=istat)
       IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,&
                           "work",t_size*(ub1-lb1+1)*&
                                         (ub2-lb2+1)*&
                                         (ub3-lb3+1)*&
                                         (ub4-lb4+1))
       END IF
       work(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4) = p(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4)
       DEALLOCATE (p,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"p")
    END IF

    ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new,lb3_new:ub3_new,lb4_new:ub4_new),&
              STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,__LINE__,&
                        "p",t_size*(ub1_new-lb1_new+1)*&
                                   (ub2_new-lb2_new+1)*&
                                   (ub3_new-lb3_new+1)*&
                                   (ub4_new-lb4_new+1))
    END IF
    p(:,:,:,:) = zero

    IF (ASSOCIATED(p).AND.ALLOCATED(work)) THEN
       p(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4) = work(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4)
       DEALLOCATE (work,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"work")
    END IF
