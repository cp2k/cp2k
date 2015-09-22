    INTEGER                                  :: istat, lb1, lb1_old, &
                                                ub1,  ub1_old
    IF (ASSOCIATED(p)) THEN
       lb1_old = LBOUND(p,1)
       ub1_old = UBOUND(p,1)
       lb1 = MAX(lb1_new,lb1_old)
       ub1 = MIN(ub1_new,ub1_old)
       ALLOCATE (work(lb1:ub1),STAT=istat)
       IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,&
                           "work",t_size*(ub1-lb1+1))
       END IF
       work(lb1:ub1) = p(lb1:ub1)
       DEALLOCATE (p)
    END IF

    ALLOCATE (p(lb1_new:ub1_new),STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,__LINE__,&
                        "p",t_size*(ub1_new-lb1_new+1))
    END IF
    p(:) = zero

    IF (ASSOCIATED(p).AND.ALLOCATED(work)) THEN
       p(lb1:ub1) = work(lb1:ub1)
       DEALLOCATE (work)
    END IF
