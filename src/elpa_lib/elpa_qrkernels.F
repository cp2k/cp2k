! calculates A = A - Y*T'*Z (rev=0)
! calculates A = A - Y*T*Z (rev=1)
! T upper triangle matrix
! assuming zero entries in matrix in upper kxk block
SUBROUTINE tum_pdlarfb_kernel_local(m,n,k,a,lda,v,ldv,t,ldt,z,ldz)
    IMPLICIT NONE
 
    ! input variables (local)
    INTEGER lda,ldv,ldt,ldz
    DOUBLE PRECISION a(lda,*),v(ldv,*),t(ldt,*),z(ldz,*)

    ! input variables (global)
    INTEGER m,n,k

    ! local variables
    DOUBLE PRECISION t11
    DOUBLE PRECISION t12,t22,sum1,sum2
    DOUBLE PRECISION t13,t23,t33,sum3
    DOUBLE PRECISION sum4,t44
    DOUBLE PRECISION y1,y2,y3,y4
    DOUBLE PRECISION a1
    INTEGER icol,irow,v1col,v2col,v3col
  
    ! reference implementation
    IF (k .EQ. 1) THEN
        t11 = t(1,1)
        DO icol=1,n
            sum1 = z(1,icol)
            a(1:m,icol) = a(1:m,icol) - t11*sum1*v(1:m,1)
        ENDDO
        RETURN
    ELSE IF (k .EQ. 2) THEN
            v1col = 2
            v2col = 1
            t22 = t(1,1)
            t12 = t(1,2)
            t11 = t(2,2)

        DO icol=1,n
            sum1 = t11 * z(v1col,icol)
            sum2 = t12 * z(v1col,icol) + t22 * z(v2col,icol)

            DO irow=1,m
                a(irow,icol) = a(irow,icol) - v(irow,v1col) * sum1 - v(irow,v2col) * sum2
            END DO
        END DO
    ELSE IF (k .EQ. 3) THEN
            v1col = 3
            v2col = 2
            v3col = 1

            t33 = t(1,1)

            t23 = t(1,2)
            t22 = t(2,2)

            t13 = t(1,3)
            t12 = t(2,3)
            t11 = t(3,3)
 
        DO icol=1,n
            ! misusing variables for fetch of z parts
            y1 = z(v1col,icol)
            y2 = z(v2col,icol)
            y3 = z(v3col,icol)

            sum1 = t11 * y1!+ 0   * y2!+ 0   * y3
            sum2 = t12 * y1 + t22 * y2!+ 0   * y3
            sum3 = t13 * y1 + t23 * y2 + t33 * y3

            DO irow=1,m
                a(irow,icol) = a(irow,icol) - v(irow,v1col) * sum1 - v(irow,v2col) * sum2 - v(irow,v3col) * sum3
            END DO
        END DO
    ELSE IF (k .EQ. 4) THEN
            DO icol=1,n
                ! misusing variables for fetch of z parts
                y1 = z(1,icol)
                y2 = z(2,icol)
                y3 = z(3,icol)
                y4 = z(4,icol)

                ! dtrmv like - starting from main diagonal and working
                ! upwards
                t11 = t(1,1)
                t22 = t(2,2)
                t33 = t(3,3)
                t44 = t(4,4)
                
                sum1 = t11 * y1
                sum2 = t22 * y2
                sum3 = t33 * y3
                sum4 = t44 * y4
 
                t11 = t(1,2)
                t22 = t(2,3)
                t33 = t(3,4)
 
                sum1 = sum1 + t11 * y2
                sum2 = sum2 + t22 * y3
                sum3 = sum3 + t33 * y4
  
                t11 = t(1,3)
                t22 = t(2,4)
 
                sum1 = sum1 + t11 * y3
                sum2 = sum2 + t22 * y4
  
                t11 = t(1,4)
                sum1 = sum1 + t11 * y4
 
                ! one column of V is calculated 
                ! time to calculate A - Y * V
                DO irow=1,m ! TODO: loop unrolling
                    y1 = v(irow,1)
                    y2 = v(irow,2)
                    y3 = v(irow,3)
                    y4 = v(irow,4)

                    a1 = a(irow,icol)

                    a1 = a1 - y1*sum1
                    a1 = a1 - y2*sum2 
                    a1 = a1 - y3*sum3
                    a1 = a1 - y4*sum4

                    a(irow,icol) = a1
                END DO
            END DO
    ELSE
        ! reference implementation
            ! V' = T * Z'
            CALL dtrmm("Left","Upper","Notrans","Nonunit",k,n,1.0d0,t,ldt,z,ldz)
            ! A = A - Y * V'
            CALL dgemm("Notrans","Notrans",m,n,k,-1.0d0,v,ldv,z,ldz,1.0d0,a,lda)
    END IF

END SUBROUTINE
SUBROUTINE tum_pdlarft_merge_kernel_local(oldk,k,t,ldt,yty,ldy)
    IMPLICIT NONE

    ! input variables (local)
    INTEGER ldt,ldy
    DOUBLE PRECISION t(ldt,*),yty(ldy,*)

    ! input variables (global)
    INTEGER k,oldk
 
    ! output variables (global)

    ! local scalars
    INTEGER icol,leftk,rightk

    ! local scalars for optimized versions
    INTEGER irow
    DOUBLE PRECISION t11
    DOUBLE PRECISION yty1,yty2,yty3,yty4,yty5,yty6,yty7,yty8
    DOUBLE PRECISION reg01,reg02,reg03,reg04,reg05,reg06,reg07,reg08
    DOUBLE PRECISION final01,final02,final03,final04,final05,final06,final07,final08

    IF (oldk .EQ. 0) RETURN ! nothing to be done

        leftk = k
        rightk = oldk
     
    ! optimized implementations:
    IF (leftk .EQ. 1) THEN
        DO icol=1,rightk
            ! multiply inner products with right t matrix
            ! (dtrmv like)
            yty1 = yty(1,1)
            t11 = t(leftk+1,leftk+icol)

            reg01 = yty1 * t11

            DO irow=2,icol
                yty1 = yty(1,irow)
                t11 = t(leftk+irow,leftk+icol)

                reg01 = reg01 + yty1 * t11
            END DO

            ! multiply intermediate results with left t matrix and store in final t
            ! matrix
            t11 = -t(1,1)
            final01 = t11 * reg01
            t(1,leftk+icol) = final01
        END DO

        !print *,'efficient tmerge - leftk=1'
    ELSE IF (leftk .EQ. 2) THEN
        DO icol=1,rightk
            ! multiply inner products with right t matrix
            ! (dtrmv like)
            yty1 = yty(1,1)
            yty2 = yty(2,1)

            t11  = t(leftk+1,leftk+icol)

            reg01 = yty1 * t11
            reg02 = yty2 * t11

            DO irow=2,icol
                yty1 = yty(1,irow)
                yty2 = yty(2,irow)
                t11 = t(leftk+irow,leftk+icol)

                reg01 = reg01 + yty1 * t11
                reg02 = reg02 + yty2 * t11
            END DO

            ! multiply intermediate results with left t matrix and store in final t
            ! matrix
            yty1 = -t(1,1)
            yty2 = -t(1,2)
            yty3 = -t(2,2)

            final01 = reg02 * yty2
            final02 = reg02 * yty3

            final01 = final01 + reg01 * yty1

            t(1,leftk+icol) = final01
            t(2,leftk+icol) = final02
        END DO
 
        !print *,'efficient tmerge - leftk=2'
    ELSE IF (leftk .EQ. 4) THEN
        DO icol=1,rightk
            ! multiply inner products with right t matrix
            ! (dtrmv like)
            yty1 = yty(1,1)
            yty2 = yty(2,1)
            yty3 = yty(3,1)
            yty4 = yty(4,1)

            t11  = t(leftk+1,leftk+icol)

            reg01 = yty1 * t11
            reg02 = yty2 * t11
            reg03 = yty3 * t11
            reg04 = yty4 * t11

            DO irow=2,icol
                yty1 = yty(1,irow)
                yty2 = yty(2,irow)
                yty3 = yty(3,irow)
                yty4 = yty(4,irow)

                t11 = t(leftk+irow,leftk+icol)

                reg01 = reg01 + yty1 * t11
                reg02 = reg02 + yty2 * t11
                reg03 = reg03 + yty3 * t11
                reg04 = reg04 + yty4 * t11
            END DO

            ! multiply intermediate results with left t matrix and store in final t
            ! matrix (start from diagonal and move upwards)
            yty1 = -t(1,1)
            yty2 = -t(2,2)
            yty3 = -t(3,3)
            yty4 = -t(4,4)

            ! main diagonal
            final01 = reg01 * yty1
            final02 = reg02 * yty2
            final03 = reg03 * yty3
            final04 = reg04 * yty4

            ! above main diagonal
            yty1 = -t(1,2)
            yty2 = -t(2,3)
            yty3 = -t(3,4)

            final01 = final01 + reg02 * yty1
            final02 = final02 + reg03 * yty2
            final03 = final03 + reg04 * yty3

            ! above first side diagonal
            yty1 = -t(1,3)
            yty2 = -t(2,4)

            final01 = final01 + reg03 * yty1
            final02 = final02 + reg04 * yty2

            ! above second side diagonal
            yty1 = -t(1,4)

            final01 = final01 + reg04 * yty1

            ! write back to final matrix
            t(1,leftk+icol) = final01
            t(2,leftk+icol) = final02
            t(3,leftk+icol) = final03
            t(4,leftk+icol) = final04
        END DO
 
        !print *,'efficient tmerge - leftk=4'
    ELSE IF (leftk .EQ. 8) THEN
        DO icol=1,rightk
            ! multiply inner products with right t matrix
            ! (dtrmv like)
            yty1 = yty(1,1)
            yty2 = yty(2,1)
            yty3 = yty(3,1)
            yty4 = yty(4,1)
            yty5 = yty(5,1)
            yty6 = yty(6,1)
            yty7 = yty(7,1)
            yty8 = yty(8,1)

            t11  = t(leftk+1,leftk+icol)

            reg01 = yty1 * t11
            reg02 = yty2 * t11
            reg03 = yty3 * t11
            reg04 = yty4 * t11
            reg05 = yty5 * t11
            reg06 = yty6 * t11
            reg07 = yty7 * t11
            reg08 = yty8 * t11

            DO irow=2,icol
                yty1 = yty(1,irow)
                yty2 = yty(2,irow)
                yty3 = yty(3,irow)
                yty4 = yty(4,irow)
                yty5 = yty(5,irow)
                yty6 = yty(6,irow)
                yty7 = yty(7,irow)
                yty8 = yty(8,irow)

                t11 = t(leftk+irow,leftk+icol)

                reg01 = reg01 + yty1 * t11
                reg02 = reg02 + yty2 * t11
                reg03 = reg03 + yty3 * t11
                reg04 = reg04 + yty4 * t11
                reg05 = reg05 + yty5 * t11
                reg06 = reg06 + yty6 * t11
                reg07 = reg07 + yty7 * t11
                reg08 = reg08 + yty8 * t11
            END DO

            ! multiply intermediate results with left t matrix and store in final t
            ! matrix (start from diagonal and move upwards)
            yty1 = -t(1,1)
            yty2 = -t(2,2)
            yty3 = -t(3,3)
            yty4 = -t(4,4)
            yty5 = -t(5,5)
            yty6 = -t(6,6)
            yty7 = -t(7,7)
            yty8 = -t(8,8)

            ! main diagonal
            final01 = reg01 * yty1
            final02 = reg02 * yty2
            final03 = reg03 * yty3
            final04 = reg04 * yty4
            final05 = reg05 * yty5
            final06 = reg06 * yty6
            final07 = reg07 * yty7
            final08 = reg08 * yty8

            ! above main diagonal
            yty1 = -t(1,2)
            yty2 = -t(2,3)
            yty3 = -t(3,4)
            yty4 = -t(4,5)
            yty5 = -t(5,6)
            yty6 = -t(6,7)
            yty7 = -t(7,8)

            final01 = final01 + reg02 * yty1
            final02 = final02 + reg03 * yty2
            final03 = final03 + reg04 * yty3
            final04 = final04 + reg05 * yty4
            final05 = final05 + reg06 * yty5
            final06 = final06 + reg07 * yty6
            final07 = final07 + reg08 * yty7

            ! above first side diagonal
            yty1 = -t(1,3)
            yty2 = -t(2,4)
            yty3 = -t(3,5)
            yty4 = -t(4,6)
            yty5 = -t(5,7)
            yty6 = -t(6,8)

            final01 = final01 + reg03 * yty1
            final02 = final02 + reg04 * yty2
            final03 = final03 + reg05 * yty3
            final04 = final04 + reg06 * yty4
            final05 = final05 + reg07 * yty5
            final06 = final06 + reg08 * yty6

            !above second side diagonal

            yty1 = -t(1,4)
            yty2 = -t(2,5)
            yty3 = -t(3,6)
            yty4 = -t(4,7)
            yty5 = -t(5,8)

            final01 = final01 + reg04 * yty1
            final02 = final02 + reg05 * yty2
            final03 = final03 + reg06 * yty3
            final04 = final04 + reg07 * yty4
            final05 = final05 + reg08 * yty5

            ! i think you got the idea by now
 
            yty1 = -t(1,5)
            yty2 = -t(2,6)
            yty3 = -t(3,7)
            yty4 = -t(4,8)

            final01 = final01 + reg05 * yty1
            final02 = final02 + reg06 * yty2
            final03 = final03 + reg07 * yty3
            final04 = final04 + reg08 * yty4

            ! .....

            yty1 = -t(1,6)
            yty2 = -t(2,7)
            yty3 = -t(3,8)

            final01 = final01 + reg06 * yty1
            final02 = final02 + reg07 * yty2
            final03 = final03 + reg08 * yty3

            ! .....

            yty1 = -t(1,7)
            yty2 = -t(2,8)

            final01 = final01 + reg07 * yty1
            final02 = final02 + reg08 * yty2
 
            ! .....

            yty1 = -t(1,8)

            final01 = final01 + reg08 * yty1

            ! write back to final matrix
            t(1,leftk+icol) = final01
            t(2,leftk+icol) = final02
            t(3,leftk+icol) = final03
            t(4,leftk+icol) = final04
            t(5,leftk+icol) = final05
            t(6,leftk+icol) = final06
            t(7,leftk+icol) = final07
            t(8,leftk+icol) = final08
        END DO

        !print *,'efficient tmerge - leftk=8'
    ELSE
        ! reference implementation
        DO icol=1,rightk
            t(1:leftk,leftk+icol) = yty(1:leftk,icol)
        END DO
            
        ! -T1 * Y1'*Y2
        CALL dtrmm("Left","Upper","Notrans","Nonunit",leftk,rightk,-1.0d0,t(1,1),ldt,t(1,leftk+1),ldt)
        ! (-T1 * Y1'*Y2) * T2
        CALL dtrmm("Right","Upper","Notrans","Nonunit",leftk,rightk,1.0d0,t(leftk+1,leftk+1),ldt,t(1,leftk+1),ldt)
    END IF

END SUBROUTINE
! yty structure
! Y1'*Y2   Y1'*Y3  Y1'*Y4 ...
!    0     Y2'*Y3  Y2'*Y4 ...
!    0        0    Y3'*Y4 ...
!    0        0       0   ...
SUBROUTINE tum_tmerge_set_kernel(k,blocksize,t,ldt,yty,ldy)
    IMPLICIT NONE
 
    ! input variables (local)
    INTEGER ldt,ldy
    DOUBLE PRECISION t(ldt,*),yty(ldy,*)

    ! input variables (global)
    INTEGER k,blocksize
 
    ! output variables (global)

    ! local scalars
    INTEGER nr_blocks,current_block
    INTEGER remainder,oldk
    INTEGER yty_column,toffset
  
    IF (k .LE. blocksize) RETURN ! nothing to merge

    nr_blocks = k / blocksize
    remainder = k - nr_blocks*blocksize

        ! work in "negative" direction:
        ! start with latest T matrix part and add older ones
        toffset = 1
        yty_column = 1
 
        IF (remainder .GT. 0) THEN
            CALL tum_pdlarft_merge_kernel_local(blocksize,remainder,t(toffset,toffset),ldt,yty(1,yty_column),ldy)
            current_block = 1
            oldk = remainder+blocksize
            yty_column =  yty_column + blocksize
        ELSE
            CALL tum_pdlarft_merge_kernel_local(blocksize,blocksize,t(toffset,toffset),ldt,yty(1,yty_column),ldy)
            current_block = 2
            oldk = 2*blocksize
            yty_column = yty_column + blocksize
        END IF
 
        DO WHILE (current_block .LT. nr_blocks)
            CALL tum_pdlarft_merge_kernel_local(blocksize,oldk,t(toffset,toffset),ldt,yty(toffset,yty_column),ldy)

            current_block = current_block + 1
            oldk = oldk + blocksize
            yty_column = yty_column + blocksize
        END DO

END SUBROUTINE
! yty structure
! Y1'*Y2   Y1'*Y3  Y1'*Y4 ...
!    0     Y2'*Y3  Y2'*Y4 ...
!    0        0    Y3'*Y4 ...
!    0        0       0   ...

SUBROUTINE tum_tmerge_tree_kernel(k,blocksize,treeorder,t,ldt,yty,ldy)
    IMPLICIT NONE
 
    ! input variables (local)
    INTEGER ldt,ldy
    DOUBLE PRECISION t(ldt,*),yty(ldy,*)

    ! input variables (global)
    INTEGER k,blocksize,treeorder
 
    ! output variables (global)

    ! local scalars
    INTEGER temp_blocksize,nr_sets,current_set,setsize,nr_blocks
    INTEGER remainder,max_treeorder,remaining_size
    INTEGER toffset,yty_column
    INTEGER toffset_start,yty_column_start
    INTEGER yty_end,total_remainder,yty_remainder

    IF (treeorder .EQ. 0) RETURN ! no merging

    IF (treeorder .EQ. 1) THEN
        CALL tum_tmerge_set_kernel(k,blocksize,t,ldt,yty,ldy)
        RETURN
    END IF
  
    nr_blocks = k / blocksize
    max_treeorder = MIN(nr_blocks,treeorder)

    IF (max_treeorder .EQ. 1) THEN
        CALL tum_tmerge_set_kernel(k,blocksize,t,ldt,yty,ldy)
        RETURN
    END IF
 
        ! work in "negative" direction: from latest set to oldest set
        ! implementation differs from rev=0 version due to issues with
        ! calculating the remainder parts
        ! compared to the rev=0 version we split remainder parts directly from
        ! parts which can be easily merged in a recursive way

        yty_end = (k / blocksize) * blocksize
        IF (yty_end .EQ. k) THEN
            yty_end = yty_end - blocksize
        END IF

        !print *,'tree',yty_end,k,blocksize

        yty_column_start = 1
        toffset_start = 1

        ! is there a remainder block?
        nr_blocks = k / blocksize
        remainder = k - nr_blocks * blocksize
        IF (remainder .EQ. 0) THEN
            !print *,'no initial remainder'

            ! set offsets to the very beginning as there is no remainder part
            yty_column_start = 1
            toffset_start = 1
            total_remainder = 0
            remaining_size = k
            yty_remainder = 0
        ELSE
            !print *,'starting with initial remainder'
            ! select submatrix and make remainder block public
            yty_column_start = 1 + blocksize
            toffset_start = 1 + remainder
            total_remainder = remainder
            remaining_size = k - remainder
            yty_remainder = 1
        END IF
 
        ! from now on it is a clean set of blocks with sizes of multiple of
        ! blocksize

        temp_blocksize = blocksize

        !-------------------------------
        DO WHILE (remaining_size .GT. 0)
            nr_blocks = remaining_size / temp_blocksize
            max_treeorder = MIN(nr_blocks,treeorder)

            IF (max_treeorder .EQ. 1) THEN
                remainder = 0
                nr_sets = 0
                setsize = 0

                IF (yty_remainder .GT. 0) THEN
                    yty_column = yty_remainder
                    !print *,'final merging with remainder',temp_blocksize,k,remaining_size,yty_column
                    CALL tum_tmerge_set_kernel(k,temp_blocksize,t,ldt,yty(1,yty_column),ldy)
                ELSE
                    !print *,'no remainder - no merging needed',temp_blocksize,k,remaining_size
                ENDIF
  
                remaining_size = 0
             
                RETURN ! done
            ELSE
                nr_sets = nr_blocks / max_treeorder
                setsize = max_treeorder*temp_blocksize
                remainder = remaining_size - nr_sets*setsize
            END IF
  
            IF (remainder .GT. 0) THEN
                IF (remainder .GT. temp_blocksize) THEN
                    toffset = toffset_start
                    yty_column = yty_column_start
 
                    !print *,'set merging', toffset, yty_column,remainder
                    CALL tum_tmerge_set_kernel(remainder,temp_blocksize,t(toffset,toffset),ldt,yty(toffset,yty_column),ldy)

                    IF (total_remainder .GT. 0) THEN
                        ! merge with existing global remainder part
                        !print *,'single+set merging',yty_remainder,total_remainder,remainder

                        CALL tum_pdlarft_merge_kernel_local(remainder,total_remainder,t(1,1),ldt,yty(1,yty_remainder),ldy)
      
                        yty_remainder = yty_remainder + remainder
                        toffset_start = toffset_start + remainder

                        !print *,'single+set merging (new offsets)',yty_remainder,yty_column_start,toffset_start

                        yty_column_start = yty_column_start + remainder
                    ELSE
                        ! create new remainder part
                        !print *,'new remainder+set',yty_remainder
                        yty_remainder = yty_column_start + remainder - temp_blocksize 
                        yty_column_start = yty_column_start + remainder
                        toffset_start = toffset_start + remainder
                        !print *,'new remainder+set (new offsets)',yty_remainder,yty_column_start,toffset_start
                    END IF

                ELSE
                    IF (total_remainder .GT. 0) THEN
                        ! merge with existing global remainder part
                        !print *,'single merging',yty_remainder,total_remainder,remainder

                        CALL tum_pdlarft_merge_kernel_local(remainder,total_remainder,t(1,1),ldt,yty(1,yty_remainder),ldy)
      
                        yty_remainder = yty_remainder + remainder
                        toffset_start = toffset_start + remainder

                        !print *,'single merging (new offsets)',yty_remainder,yty_column_start,toffset_start

                        yty_column_start = yty_column_start + remainder
                    ELSE
                        ! create new remainder part
                        !print *,'new remainder',yty_remainder
                        yty_remainder = yty_column_start
                        yty_column_start = yty_column_start + temp_blocksize
                        toffset_start = toffset_start + remainder
                        !print *,'new remainder (new offsets)',yty_remainder,yty_column_start,toffset_start
                    END IF
                END IF
 
                total_remainder = total_remainder + remainder
                remaining_size = remaining_size - remainder
            END IF

            current_set = 0
            DO WHILE (current_set .LT. nr_sets)
                toffset = toffset_start + current_set * setsize
                yty_column = yty_column_start + current_set * setsize

                !print *,'recursive merging', toffset, yty_column,setsize

                CALL tum_tmerge_set_kernel(setsize,temp_blocksize,t(toffset,toffset),ldt,yty(toffset,yty_column),ldy)
                
                current_set = current_set +  1
            END DO

            !print *,'increasing blocksize', temp_blocksize, setsize
            yty_column_start = yty_column_start + (setsize - temp_blocksize)
            temp_blocksize = setsize
        END DO
END SUBROUTINE
! yty should not contain the inner products vi'*vi
SUBROUTINE tum_dlarft_kernel(n,tau,yty,ldy,t,ldt)
    IMPLICIT NONE

    ! input variables
    INTEGER n,ldy,ldt
    DOUBLE PRECISION tau(*),yty(ldy,*)
    
    ! output variables
    DOUBLE PRECISION t(ldt,*)

    ! local variables
    INTEGER icol
 
    ! DEBUG: clear buffer first
    !t(1:n,1:n) = 0.0d0

        ! T1 = tau1
        ! | tauk  Tk-1' * (-tauk * Y(:,1,k+1:n) * Y(:,k))' |
        ! | 0           Tk-1                           |
        t(n,n) = tau(n)
        DO icol=n-1,1,-1
            t(icol,icol+1:n) = -tau(icol)*yty(icol,icol:n-1)
            CALL dtrmv("Upper","Trans","Nonunit",n-icol,t(icol+1,icol+1),ldt,t(icol,icol+1),ldt)
            t(icol,icol) = tau(icol)
        END DO
END SUBROUTINE
