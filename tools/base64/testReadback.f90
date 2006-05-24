PROGRAM testReadback
integer, parameter :: wp=selected_real_kind(14,200) ! SELECTED_REAL_KIND(5,30)
real(wp) :: r1(6)=(/1._wp/3._wp,234.345e1_wp,tiny(1._wp),tiny(1._wp)/2._wp,tiny(1._wp)/4._wp,tiny(1._wp)/8._wp/),&
     r2(6),diff1,diff2,diffAtt
character(len=100) :: str
integer :: i
! tests if a compiler manages to recreate exactly the same number
! from a formatted output
str=" "
diff1=0._wp
do i=1,size(r1)
   write(str,"(es35.20e5)")r1(i)
   read(str,"(es35.20e5)") r2(i)
   diffAtt=(abs(r1(i)-r2(i)))/max(abs(r1(i)),tiny(1._wp))
   diff1=MAX(diff1,diffAtt)
   PRINT *," diff",diffAtt," r1:",r1(i)," r2:",r2(i)
end do
PRINT *
do i=1,size(r1)
   WRITE(str,*)r1(i)
   READ(str,*) r2(i)
   diffAtt=(abs(r1(i)-r2(i)))/max(abs(r1(i)),tiny(1._wp))
   diff2=MAX(diff2,diffAtt)
   PRINT *," diff",diffAtt," r1:",r1(i)," r2:",r2(i)
end do
PRINT *
PRINT *,"diff es35.10e5 format: ",diff1," diff '*' format: ",diff2
END PROGRAM testReadback

