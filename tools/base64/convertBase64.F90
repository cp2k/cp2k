module machine
#ifdef __NAG
  USE f90_unix
  USE f90_unix_dir
  USE f90_unix_proc
#endif
  implicit none
contains
FUNCTION m_iargc() RESULT (ic)
    INTEGER                                  :: ic
    ic = iargc()
END FUNCTION m_iargc
SUBROUTINE m_getarg(i,arg)
    INTEGER, INTENT(IN)                      :: i
    CHARACTER(len=*), INTENT(OUT)            :: arg

    CALL getarg(i,arg)
END SUBROUTINE m_getarg

end module machine

module base64_types

#ifdef __NAG
#define __IEEE
#define __IEEE_EXCEPTIONS
use IEEE_ARITHMETIC
use IEEE_EXCEPTIONS
#endif
#ifdef __XLF
#define __IEEE
#endif

implicit none

integer, parameter :: wp= SELECTED_REAL_KIND ( 14, 200 )! SELECTED_REAL_KIND(5,30) ! 
integer, parameter :: mt=selected_int_kind(precision(1._wp))
character(len=64) :: base64="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"

character(len=250) :: refString,ref32String=&
     "lI1+BFAAABAAAgAABAAgAAAAACw///3fAAAA0AAAA+3////fAAAg/AAAAAQ"//&
     "JStfwBAAQACAAQAYAAAIgAAAgA+///9PAAAAtAAAg/7t183DAAAEA",&
     ref64String=&
     "Y4jImSka/AEAAAAAAAAA4AAAAAAAAA9NAAAAAAAAQgDAAAAAAAAEA8/////"//&
     "///+fAAAAAAAAwyDAAAAAAAA8/9/////////fAAAAAAAAw/DAAAAAAAAAAg"//&
     "hPiYKRq9DwAAAAAAAAAgLAAAAAAAA43CAAAAAAAABuAAAAAAAAQA4//////"//&
     "//7/DAAAAAAAALvAAAAAAAAw/vcajPubr5v/AAAAAAAAgAAAAAAAAAAEA43"//&
     "P14luJBAA"

real(wp) :: infinity,nan,small_val

integer, private, save :: bit_mantissa=0, bit_exp=0
integer(mt), private, save :: pad_exp, pad_mantissa,pad6
logical, private, save :: initialized_module=.false.

type base64c_type
   integer(mt)::pad,pad_exp_in
   integer :: len_pad,len_str,pos_str, bit_mantissa_in, bit_exp_in, unit_nr
   character(len=380) :: str
   logical :: use_fast_conversion
end type base64c_type

type failures_type
   logical :: support_normal,support_denormal,support_infinity,support_nan
end type failures_type

contains

function get_bit_mantissa()
  integer :: get_bit_mantissa
  call init_module()
  get_bit_mantissa=bit_mantissa
end function get_bit_mantissa

function get_bit_exp()
  integer :: get_bit_exp
  call init_module()
  get_bit_exp=bit_exp
end function get_bit_exp

subroutine print_bits(unit_nr,bits,nbits)
  integer, intent(in) :: unit_nr
  integer(mt), intent(in) :: bits
  integer, intent(in), optional :: nbits

  integer :: my_nbits,i
  my_nbits=bit_size(bits)
  if (present(nbits)) my_nbits=min(my_nbits,nbits)
  do i=0,my_nbits-1
     write (unit_nr,"(l1)",advance="NO") btest(bits,i)
  end do
  write (unit_nr,"(' ',i20)") bits
end subroutine print_bits
   
function mantissa_m1(r1) result(res)
  real(wp),intent(in) :: r1
  integer(mt) :: res
  integer(mt), dimension(1) :: m
  m=transfer(r1,res,1)
  res=m(1)
end function mantissa_m1

function mantissa_m2(r1) result(res)
  real(wp),intent(in) :: r1
  integer(mt) :: res,frac
  integer :: shift
  
  frac=int(scale(fraction(abs(r1)),digits(r1)),mt)  
  shift=int(minexponent(r1),mt)-min(int(minexponent(r1),mt),my_exponent(r1))
  res=my_shift(frac,-shift)
end function mantissa_m2

function mantissa_m3(r1) result(res)
  real(wp),intent(in) :: r1
  integer(mt) :: res
  real(wp) :: ffrac
  integer :: i
  
  if (my_is_denormal(r1)) then
     ffrac=fraction(abs(r1)+small_val)
  else
     ffrac=fraction(abs(r1))
  end if
  if (ffrac>=1._wp) STOP "s1"
  res=0
  do i=1,digits(r1)
     if (ffrac==0._wp) exit
     !ffrac=ffrac*2._wp
     ffrac=scale(ffrac,1)!wrong if radix/=2
     if (ffrac>=1._wp) then
        ffrac=ffrac-1._wp
        res=ibset(res,digits(r1)-i)
     end if
  end do
end function mantissa_m3

function from_mantissa_m3(mantissa,expo,bit_exp_in,pad_e) result(res)
  integer(mt), intent(in) :: mantissa,expo
  integer, intent(in) :: bit_exp_in
  integer(mt), intent(in) :: pad_e
  real(wp) :: res

  integer :: expon,pow2,i
  real(wp) :: ffrac
  
  ffrac=0._wp
  do i=0,bit_mantissa-1
     if (btest(mantissa,i)) then
        ffrac=ffrac+1._wp
     end if
     ffrac=scale(ffrac,-1)
  end do
  if (expo/=0) then
     ffrac=ffrac+1._wp
     if (expo==pad_e) then
        if (mantissa==0_mt) then
           !infinity
           ffrac=infinity
        else
           ! nan
           ffrac=nan
        end if
     else
        expon=0
        pow2=1
        do i=0,bit_exp_in-1
           if (btest(expo,i)) expon=expon+pow2
           pow2=pow2*2
        end do
        expon=expon-pow2/2+1
        if (expon> maxexponent(ffrac)) then
           ffrac=infinity
        else if (expon< minexponent(ffrac)) then
           if (expon< minexponent(ffrac)-digits(ffrac)+1) then
              ffrac=0._wp
           else
              ! tryies to represent it as denormalized number
              ffrac=scale(ffrac,minexponent(ffrac))
              ffrac=ffrac*0.5_wp**(minexponent(ffrac)-expon)
           end if
        else
           ffrac=scale(ffrac,expon)
        end if
     end if
  else
     if (bit_exp_in>bit_exp) then
        ffrac=0._wp
     else
        ffrac=scale(ffrac,-2**(bit_exp_in-1)+digits(ffrac)+1)
        ffrac=ffrac*0.5_wp**(digits(ffrac)-1)
     end if
  end if
  res=ffrac
end function from_mantissa_m3

function exponent_m1(r1) result(res)
  real(wp),intent(in) :: r1
  integer(mt) :: res
  integer(mt), dimension(1) :: m
  m=transfer(r1,res,1)
  res=my_shift(m(1),-digits(r1)+1)
end function exponent_m1

function exponent_m2(r1) result(res)
  real(wp),intent(in) :: r1
  integer(mt) :: res,e
  e=my_exponent(r1)-minexponent(r1)+1
  if (my_is_denormal(r1).or.r1==0._wp) e=0
  res=e
end function exponent_m2

function my_shift(val,sh) result(res)
  integer(mt),intent(in) :: val
  integer,intent(in) :: sh
  integer(mt) :: res
  res=ishft(val,sh)
end function my_shift

subroutine base64c_init(b64,bit_mantissa_in,bit_exp_in,log_unit)
  type(base64c_type), intent(out) :: b64
  integer, intent(in), optional :: bit_mantissa_in, bit_exp_in, log_unit

  integer :: unit_nr,i

  call init_module()
  b64%pad=0
  b64%len_pad=0
  b64%len_str=0
  b64%pos_str=1
  b64%str=" "
  b64%bit_mantissa_in=bit_mantissa
  if (present(bit_mantissa_in)) b64% bit_mantissa_in=bit_mantissa_in
  b64%bit_exp_in=bit_exp
  if (present(bit_exp_in)) b64%bit_exp_in=bit_exp_in
  if (bit_size(0_mt)-1< b64%bit_exp_in) stop "bitsize too big"
  b64%pad_exp_in=0
  do i=0,b64%bit_exp_in-1
     b64%pad_exp_in=ibset(b64%pad_exp_in,i)
  end do
  b64%use_fast_conversion=.true.
  unit_nr=-1
  if (present(log_unit)) unit_nr=log_unit
  call base64c_do_tests(b64,unit_nr)
end subroutine base64c_init

subroutine base64c_dealloc_ref(b64)
  type(base64c_type), intent(inout) :: b64
  
  ! no cleanup to do at the moment
end subroutine base64c_dealloc_ref

  subroutine base64c_do_tests(b64,unit_nr)
    type(base64c_type), intent(inout) :: b64
    integer, intent(in) :: unit_nr
    logical :: ref_string_match,use_fast_conversion
    type(failures_type) :: incoerences,ref_failures,fast_failures
    real(wp) :: r1
    real(wp), dimension(25) :: test2,test_r=&
         (/ 31.41511_wp,5.8774717541114375e-39_wp,&
         2.9387358770557188e-39_wp/4._wp,1.1754943508222875e-38_wp,&
         tiny(r1),huge(r1),epsilon(r1),0._wp,0._wp,1._wp,0._wp,&
         -31.41511_wp,-5.8774717541114375e-39_wp,&
         -2.9387358770557188e-39_wp/2._wp,-1.1754943508222875e-38_wp,&
         -tiny(r1),-huge(r1),-epsilon(r1),-1._wp,0.1234567_wp,tiny(r1)/2._wp,&
         -tiny(r1)/4._wp,1.0045_wp*tiny(r1)-tiny(r1),0._wp,0._wp /)
    integer :: ii
    integer (mt) :: expo(2),mantissa(3)

    use_fast_conversion=.true.
    call failures_init(incoerences)
    call failures_init(ref_failures)
    call failures_init(fast_failures)
    test_r(8)=infinity
    test_r(9)=nan
    test_r(19)=-infinity
    call random_number(test_r(24))
    call random_number(test_r(25))
    test_r(25)=(test_r(25)-0.5_wp)/test_r(24)
    r1=test_r(5)
    if (unit_nr>0) write (unit_nr,"(a)") "*** data on number type ***"
    if (unit_nr>0) then
       write (unit_nr,"(a,i6,a,i20,a,i20,a,i6)") " digits ",digits(r1),&
            " maxexp ",maxexponent(r1)," minexp ",minexponent(r1)," radix ",radix(r1)
       write (unit_nr,"(a,i6,a,i6)") "precision ",precision(r1)," range ",range(r1)
       write (unit_nr,"(a,i6,a,i6)") "bit_mant ",bit_mantissa," bit_exp ",bit_exp
       !ceiling(log(real(radix(r1),wp))/log(2._wp)*real(digits(r1)+1,wp))
    end if

    if (bit_exp/=ceiling(log(real(maxexponent(1._wp),wp))/log(2._wp))+1) then
       if (unit_nr>0) then
          write (unit_nr,"(a)") "***** WARNING *****, strange exponent range"
       end if
       use_fast_conversion=.false.
    end if
    if (unit_nr>0) then
       write (unit_nr,"(a,i6)") "total bitsize of integer container of mantissa",bit_size(mantissa(1))
       write (unit_nr,"(a,i6)") "guaranteed bitsize",bit_mantissa+1

       write (unit_nr,"(/,a,es35.20e5,a)") "** data on r1=",r1," **"
       write (unit_nr,"(a,l1)") "r1==0:",r1==0._wp
       write (unit_nr,"(a,es35.20e5,a,es35.20e5)")"rrspacing ",rrspacing(r1)," spacing ",spacing(r1)
       write (unit_nr,"(a,es35.20e5,' ',es35.20e5)")"nearest ",nearest(r1, 1._wp),nearest(r1, -1._wp)

       expo=my_exponent(r1)
       write (unit_nr,"(a,i20,a,es35.20e5)")"my_exponent(r1)",expo(1)," fraction(r1) ",fraction(r1)
       write (unit_nr,"(a,l1,a,l1,a,l1)") "my_is_nan:",my_is_nan(r1),&
            " my_is_finite:",my_is_finite(r1),&
            " my_is_denormal:",my_is_denormal(r1)
    end if

    mantissa(1)=mantissa_m1(r1)
    mantissa(2)=mantissa_m2(r1)
    mantissa(3)=mantissa_m3(r1)

    if (unit_nr>0) then
       do ii=1,3
          write (unit_nr,"('mantissa',i3,'=')",advance="NO") ii
          call print_bits(unit_nr,mantissa(ii),digits(r1)-1)
       end do
    end if

    if (unit_nr>0) then
       do ii=1,2
          if (iand(pad_mantissa,mantissa(ii))/=iand(pad_mantissa,mantissa(3))) then
             write (unit_nr,"(a,i2,a)") "*** WARNING *** mantissa(",ii,") /= mantissa(3)"
          end if
       end do
    end if

    expo(1)=exponent_m1(r1)
    expo(2)=exponent_m2(r1)

    if (unit_nr>0) then
       do ii=1,2
          write (unit_nr,"('expo',i3,'=')",advance="NO") ii
          call print_bits(unit_nr,expo(ii),bit_exp)
       end do
    end if

    if (unit_nr>0) then
       do ii=1,1
          if (iand(pad_exp,expo(ii))/=iand(pad_exp,expo(2))) then
             write (unit_nr,"(a,i2,a)") "*** WARNING *** expo(",ii,") /= expo(2)"
          end if
       end do

       write (unit_nr,"(/)")
    end if

    ! test ref in
    do ii=1,size(test_r)
       r1=test_r(ii)
       call base64c_real_ref(b64,r1,err_rep=incoerences,unit_nr=unit_nr)
    end do
    call base64c_flush(b64)

    if (bit_mantissa/=b64%bit_mantissa_in.or.bit_exp/=b64%bit_exp_in) then
       use_fast_conversion=.false.
    else
       ! test ref out
       do ii=1,size(test_r)
          test2(ii)=base64c_decode_real_ref(b64)
          if (test2(ii)/=test_r(ii).and..not.(&
               my_is_nan(test_r(ii)).and.my_is_nan(test2(ii))) ) then
             call add_failed_nr(ref_failures,test_r(ii))
             if (unit_nr>0) then
                write (unit_nr,"(a,i3,a,es35.20e5,a,es35.20e5)") "ii",ii," test_r ",test_r(ii)," test2 ",test2(ii)
             end if
          end if
       end do
    end if
    ref_string_match=(b64%str(1:len_trim(refString))==refString)
    b64%pos_str=1
    if (unit_nr>0) then
       write (unit_nr,"('base64 reference string')")
       call base64c_write(b64,unit_nr)
    end if

    if (bit_size(mantissa(1))/=bit_mantissa+bit_exp+1) use_fast_conversion=.false.
    ! test fast conversion
    if (use_fast_conversion) then
       if (unit_nr>0) write (unit_nr,"('TEST Fast conversion')")
       call base64c_clear(b64)
       do ii=1,size(test_r)
          r1=test_r(ii)
          call base64c_real_fast(b64,r1)
       end do
       call base64c_flush(b64)
       do ii=1,size(test_r)
          test2(ii)=base64c_decode_real_fast(b64)
          if (test2(ii)/=test_r(ii).and..not.(&
               my_is_nan(test_r(ii)).and.my_is_nan(test2(ii))) ) then
             call add_failed_nr(fast_failures,test_r(ii))
             if (unit_nr>0) then
                write (unit_nr,"(a,i2,a,es35.20e5,a,es35.20e5)") "ii",ii,"test_r",test_r(ii),"test2",test2(ii)
             end if
          end if
       end do
       b64%pos_str=1
       if (unit_nr>0) then
          write(unit_nr,"('base64 fast string')")
          call base64c_write(b64,unit_nr)
          write (unit_nr,"('** fast conversions **')")
          call write_failures(fast_failures,unit_nr)
       end if
    end if

    ! add a test that decodes the reference string?
    call base64c_clear(b64)

    if (unit_nr>0) then
       write(unit_nr,"('** ref conversions **')")
       call write_failures(ref_failures,unit_nr)
       write(unit_nr,"('** incoerences **')")
       call write_failures(incoerences,unit_nr)
       if (len_trim(refString)>0) then   
          if (ref_string_match) then
             write(unit_nr,"('refString matches')")
          else
             write(unit_nr,"('ERROR refString DOES NOT MATCH')")
             write(unit_nr,"('refString:')")
             write(unit_nr,"(a)") refString
          end if
       else
          write(unit_nr,"('no refString available')")
       end if
       write(unit_nr,"('**** end tests ****')")
    end if

    if (.not.ref_failures%support_normal) STOP "reference implementation broken"
    if (.not.incoerences%support_normal) use_fast_conversion=.false.
    if (use_fast_conversion) then
       if (.not.fast_failures%support_normal) use_fast_conversion=.false.
    end if

    if (.not.use_fast_conversion) b64%use_fast_conversion=.false.
  end subroutine base64c_do_tests

subroutine base64c_real_ref(b64,r1,err_rep,unit_nr)
  type(base64c_type), intent(inout) :: b64
  real(wp),intent(in) :: r1
  type(failures_type), intent(inout), optional :: err_rep
  integer, intent(in), optional :: unit_nr

  integer(mt) :: mantissa,pad_m,expo,pad_e,segno,comb_fast,comb_ref
  integer :: bit_m
  
  bit_m=bit_mantissa
  pad_m=pad_mantissa
  pad_e=pad_exp
  
  segno=0
  if (.not.my_is_finite(r1)) then
     expo=pad_e
     if (my_is_nan(r1)) then
        mantissa=pad_m !quiet
        ! mantissa=ibclr(pad_m,bit_mantissa-1) !signalling
     else
        if (r1<0._wp) segno=1
        mantissa=0
     end if
  else
     if (r1<0._wp) segno=1
     mantissa=mantissa_m3(r1)
     if (my_is_denormal(r1)) then
        expo=0
     else
        expo=exponent_m2(r1)
     end if
  end if
  comb_fast=mantissa_m1(r1)
  comb_ref=ior(ior(ibits(mantissa,0,bit_m),my_shift(ibits(expo,0,bit_exp),bit_m)),&
       my_shift(segno,bit_m+bit_exp))
  if (comb_fast/=comb_ref.and..not.my_is_nan(r1)) then
     if (present(err_rep)) call add_failed_nr(err_rep,r1)
     if (present(unit_nr)) then
        if (unit_nr>0) then
           if (my_is_denormal(r1)) then
              write(unit_nr,"('denormal r1=',es35.20e5)") r1
           else
              write (unit_nr,"('r1=',es35.20e5)")r1
           end if
           write (unit_nr,"('comb_ref =')",advance="NO")
           call print_bits(unit_nr,comb_ref)
           write (unit_nr,"('comb_fast=')",advance="NO")
           call print_bits(unit_nr,comb_fast)
        end if
     end if
  end if
  
  call base64c_feed(b64,mantissa,bit_m)
  call base64c_feed(b64,expo,bit_exp)
  call base64c_feed(b64,segno,1)
end subroutine base64c_real_ref

subroutine base64c_real_fast(b64,r1)
  type(base64c_type), intent(inout) :: b64
  real(wp),intent(in) :: r1
  
  integer(mt),dimension(1) :: nr
  
  if (bit_size(nr)/=bit_exp+bit_mantissa+1) &
       stop "no fast_conversion possible"
  nr=transfer(r1,nr(1),1)
  call base64c_feed(b64,nr(1),bit_mantissa+bit_exp+1)
end subroutine base64c_real_fast

subroutine base64c_real(b64,r1)
  type(base64c_type), intent(inout) :: b64
  real(wp),intent(in) :: r1
  
  if (b64%use_fast_conversion) then
     call base64c_real_fast(b64,r1)
  else
     call base64c_real_ref(b64,r1)
  end if
end subroutine base64c_real

subroutine base64c_feed(b64,bits,nbits)
  type(base64c_type), intent(inout) :: b64
  integer(mt),intent(in) :: bits
  integer, intent(in) :: nbits

  integer(mt) :: indice,my_bits
  integer :: my_nbits

  my_bits=bits
  my_nbits=nbits
  do while (my_nbits>0)
!     print *,"my_bits"
!     call print_bits(my_bits,my_nbits)
     indice=iand(ior(my_shift(my_bits,b64%len_pad),b64%pad),pad6)
!     print *, "pad,indice"
!     call print_bits(b64%pad,b64%len_pad)
!     call print_bits(indice,6)
     if (b64%len_pad+my_nbits>=6) then
        if (b64%len_str>=len(b64%str)) stop "maxLine"
        if (indice>=64.or.indice<0) stop "indice"
        b64%len_str=b64%len_str+1
        b64%str(b64%len_str:b64%len_str)=base64(indice+1:indice+1)
!        print *,"encoded ",base64(indice+1:indice+1)
        my_bits=my_shift(my_bits,-6+b64%len_pad)
        my_nbits=my_nbits-6+b64%len_pad
        b64%len_pad=0
        b64%pad=0
     else
        b64%len_pad=b64%len_pad+my_nbits
        b64%pad=iand(indice,my_shift(pad6,-6+b64%len_pad))
        my_nbits=0
        my_bits=0
     end if
  end do
end subroutine base64c_feed

subroutine base64c_flush(b64)
  type(base64c_type), intent(inout) :: b64
  integer(mt) :: indice
  
  if (b64%len_pad>0) then
     if (b64%len_str>=len(b64%str)) stop "maxLine2"
     indice=b64%pad
     if (indice>=64.or.indice<0) stop "indice2"
     b64%len_str=b64%len_str+1
     b64%str(b64%len_str:b64%len_str)=base64(indice+1:indice+1)
     b64%pad=0
     b64%len_pad=0
  end if
end subroutine base64c_flush

subroutine base64c_clear(b64,clear_pad)
  type(base64c_type), intent(inout) :: b64
  logical, intent(in), optional :: clear_pad
  
  logical :: my_clear_pad

  b64%str=" "
  b64%len_str=0
  my_clear_pad=.true.
  if (present(clear_pad)) my_clear_pad=clear_pad
  if (my_clear_pad) then
     b64%pad=0
     b64%len_pad=0
  end if
end subroutine base64c_clear

subroutine base64c_write(b64,unit_nr)
  type(base64c_type), intent(inout) :: b64
  integer, intent(in) :: unit_nr

  write(unit_nr,"(a)") b64%str(1:b64%len_str)
  call base64c_clear(b64,.false.)
end subroutine base64c_write

subroutine base64c_getbits(b64,bits,nbits)
  type(base64c_type), intent(inout) :: b64
  integer(mt),intent(out) :: bits
  integer, intent(in) :: nbits

  integer(mt), parameter :: pad6=63
  integer(mt) :: indice,my_bits
  integer :: my_nbits

  if (bit_size(indice)<nbits) stop "nbits too big"
  my_bits=0
  my_nbits=0
  do while (my_nbits<nbits)
     if (b64%len_pad==0) then
        if (b64%len_str<b64%pos_str) stop "no string to decode"
        b64%pad=base64_char_to_num(b64%str(b64%pos_str:b64%pos_str))
        b64%pos_str=b64%pos_str+1
        b64%len_pad=6
     end if
     if (nbits-my_nbits>b64%len_pad) then
        my_bits=ior(my_shift(b64%pad,my_nbits),my_bits)
!        print *,"my_bits,pad"
!        call print_bits(my_bits,my_nbits+b64%len_pad)
!        call print_bits(b64%pad,b64%len_pad)
        b64%pad=0
        my_nbits=my_nbits+b64%len_pad
        b64%len_pad=0
     else
        my_bits=ior(my_shift(iand(my_shift(pad6,-6+nbits-my_nbits),&
             b64%pad),my_nbits),my_bits)
!        print *,"my_bits,pad"
!        call print_bits(my_bits,nbits)
!        call print_bits(b64%pad,b64%len_pad)
        b64%pad=my_shift(b64%pad,-nbits+my_nbits)
        b64%len_pad=b64%len_pad-nbits+my_nbits
        my_nbits=nbits
     end if
  end do
  bits=my_bits
!  print *, "bits"
!  call print_bits(bits,nbits)
end subroutine base64c_getbits

subroutine base64c_capacity_info(b64,chars_contained,chars_free,&
     bits_contained,bits_free,reals_contained,reals_free)
  type(base64c_type), intent(in) :: b64
  integer, intent(out), optional :: chars_contained,chars_free,&
       bits_contained,bits_free,reals_contained,reals_free

  integer :: my_bits_free,my_bits_contained
  
  if (present(chars_contained))chars_contained=b64%len_str-b64%pos_str+1
  if (present(chars_free))chars_free=len(b64%str)-chars_contained
  my_bits_contained=max(0,(b64%len_str-b64%pos_str+1)*6)
  my_bits_free=6*len(b64%str)-my_bits_contained
  my_bits_contained=my_bits_contained+b64%len_pad
  if (present(bits_free))bits_free=my_bits_free
  if (present(bits_contained))bits_contained=my_bits_contained
  if (present(reals_free)) reals_free=my_bits_free/&
       (b64%bit_mantissa_in+b64%bit_exp_in+1)
  if (present(reals_contained)) reals_contained=my_bits_contained/&
       (b64%bit_mantissa_in+b64%bit_exp_in+1)
end subroutine base64c_capacity_info

function base64c_add_string(b64,str,ignore_garbage,skip_space) result(parsed_chars)
  type(base64c_type), intent(inout) :: b64
  character(len=*), intent(in) :: str
  logical, intent(in), optional :: ignore_garbage, skip_space
  integer :: parsed_chars
  ! ignore_garbage: wether to ignore garbage or stop at it, defaults to false
  ! skip_space: if space should be ignored (defaults to true)

  logical :: my_ignore_garbage,my_skip_space
  character :: c
  integer :: n,i

  my_ignore_garbage=.false.
  if (present(ignore_garbage)) my_ignore_garbage=ignore_garbage
  my_skip_space=.true.
  if (present(skip_space)) my_skip_space=skip_space
  b64%len_str=b64%len_str-b64%pos_str+1
  if (b64%pos_str/=1) then
     b64%str(1:b64%len_str)=b64%str(b64%pos_str:b64%pos_str+b64%len_str-1)
  end if
  b64%pos_str=1
  parsed_chars=0
  do i=1,len_trim(str)
     c=str(i:i)
     n=base64_char_to_num(c)
     if (n==-1) then
        if (c==' ') then
           if (.not.my_skip_space) exit
        else if (.not.my_ignore_garbage) then
           exit
        end if
     else
        if (b64%len_str>=len(b64%str)) STOP "buffer too small to add the whole string"
        b64%len_str=b64%len_str+1
        b64%str(b64%len_str:b64%len_str)=c
     end if
     parsed_chars=i
  end do
end function base64c_add_string

function base64_char_to_num(c) result(res)
  !! result: a number from 0 to 63
  character :: c
  integer :: res

  integer :: aval

  aval=iachar(c)
  if (aval>=iachar('A') .and. aval<=iachar('Z')) then
     res=aval-iachar('A')
  else if (aval>=iachar('a').and. aval<=iachar('z')) then
     res=aval-iachar('a')+26
  else if (aval>=iachar('1').and.aval<=iachar('9')) then
     res=aval-iachar('1')+53
  else if (c=='0') then
     res=52
  else if (c=='+') then
     res=62
  else if (c=='/') then
     res=63
  else
     res=-1
  end if
end function base64_char_to_num

function base64c_decode_real_ref(b64,err_rep,log_unit) result(res)
  type(base64c_type), intent(inout) :: b64
  type(failures_type), intent(inout), optional :: err_rep
  integer, intent(in), optional :: log_unit
  real(wp) :: res

  integer(mt) :: mantissa,pad_m,expo,pad_e,segno,comb_fast,comb_ref
  integer :: bit_m
  real(wp) :: r1
  
  bit_m=bit_mantissa
  pad_m=pad_mantissa
  pad_e=pad_exp
  
  ! skip over precision
  call base64c_getbits(b64,mantissa,max(0,b64%bit_mantissa_in-bit_mantissa))
  mantissa=0
  call base64c_getbits(b64,mantissa,min(bit_mantissa,b64%bit_mantissa_in))
  ! align if under precision
  mantissa=my_shift(mantissa,max(0,bit_mantissa-b64%bit_mantissa_in))
  if (b64%bit_exp_in>bit_size(expo)-1) stop "exponent too big"
  call base64c_getbits(b64,expo,b64%bit_exp_in)
  call base64c_getbits(b64,segno,1)
  r1=from_mantissa_m3(mantissa,expo,b64%bit_exp_in,b64%pad_exp_in)
  if (segno==1) r1=-r1
  
  comb_fast=mantissa_m1(r1)
  comb_ref=ior(ior(ibits(mantissa,0,bit_m),my_shift(ibits(expo,0,bit_exp),bit_m)),&
       my_shift(segno,bit_m+bit_exp))
  if (comb_fast/=comb_ref.and..not.my_is_nan(r1)) then
     if (present(err_rep)) call add_failed_nr(err_rep,r1)
     if (present(log_unit)) then
        if (log_unit>0) then
           write (log_unit,"('mantissa =')",advance="NO")
           call print_bits(log_unit,mantissa,digits(r1)-1)
           if (my_is_denormal(r1)) then
              write(log_unit,"('denormal r=',es35.20e5)")r1
           else
              write(log_unit,"('r=',es35.20e5)")r1
           end if
           write (log_unit,"('comb_ref =')",advance="NO")
           call print_bits(log_unit,comb_ref)
           write (log_unit,"('comb_fast=')",advance="NO")
           call print_bits(log_unit,comb_fast)
        end if
     end if
  end if
  res=r1
end function base64c_decode_real_ref

function base64c_decode_real_fast(b64) result(res)
  type(base64c_type), intent(inout) :: b64
  real(wp) :: res

  integer(mt) :: nr
  real(wp),dimension(1) :: r1
  
  if (bit_size(nr)/=bit_exp+bit_mantissa+1 .or.&
       bit_exp/=b64%bit_exp_in .or.&
       bit_mantissa/=b64%bit_mantissa_in) then
     stop "fast bit conversion not possible"
  end if

  call base64c_getbits(b64,nr,bit_mantissa+bit_exp+1)
  r1=transfer(nr,res,1)
  res=r1(1)
end function base64c_decode_real_fast

function base64c_decode_real(b64) result(res)
  type(base64c_type), intent(inout) :: b64
  real(wp) :: res

  if (b64%use_fast_conversion.and.&
       bit_exp==b64%bit_exp_in .or.&
       bit_mantissa==b64%bit_mantissa_in) then
     res=base64c_decode_real_fast(b64)
  else
     res=base64c_decode_real_ref(b64)
  end if
end function base64c_decode_real

subroutine failures_init(fail)
  type(failures_type), intent(out) :: fail

  fail%support_normal=.true.
  fail%support_denormal=.true.
  fail%support_infinity=.true.
  fail%support_nan=.true.
end subroutine failures_init

subroutine add_failed_nr(fail,nr)
  type(failures_type), intent(inout) :: fail
  real(wp) :: nr
  
  if (my_is_finite(nr)) then
     if (my_is_denormal(nr)) then
        fail%support_denormal=.false.
     else
        fail%support_normal=.false.
     end if
  else if (my_is_nan(nr)) then
     fail%support_nan=.false.
  else
     fail%support_infinity=.false.
  end if
end subroutine add_failed_nr

subroutine write_failures(fail,unit_nr)
  type(failures_type), intent(inout) :: fail
  integer :: unit_nr

  if (.not.fail%support_normal)   write(unit_nr,"('NORMAL_NUMBERS FAILED')")
  if (.not.fail%support_denormal) write(unit_nr,"('DENORMAL_NUMBERS FAILED')")
  if (.not.fail%support_infinity) write(unit_nr,"('INFINITY FAILED')")
  if (.not.fail%support_nan)      write(unit_nr,"('NotANumber FAILED')")
end subroutine write_failures

subroutine init_module()
  real(wp) :: one
  integer :: i
#ifdef __IEEE_EXCEPTIONS
     type(ieee_status_type) :: stat
     logical :: flag,flag2
#endif

  if (.not.initialized_module) then
#ifdef __IEEE_EXCEPTIONS
     call ieee_get_halting_mode(ieee_divide_by_zero,flag)
     call ieee_get_halting_mode(ieee_invalid,flag2)
     call ieee_get_status(stat)
     call ieee_set_halting_mode(ieee_divide_by_zero,.false.)
     call ieee_set_halting_mode(ieee_invalid,.false.)
#endif
     call random_number(one)
     one=one+1._wp
     nan=(one-one)/(one-one)
     infinity=1._wp/(one-one)
     small_val=scale(0.5_wp,minexponent(small_val))  
#ifdef __IEEE_EXCEPTIONS
     call ieee_set_status(stat)
     call ieee_set_halting_mode(ieee_divide_by_zero,flag)
     call ieee_set_halting_mode(ieee_invalid,flag2)
#endif

     if (radix(1._wp)/=2) then
        bit_mantissa=ceiling(log(real(radix(1._wp),wp))/log(10._wp)*real(digits(1._wp),wp))
     else
        bit_mantissa=digits(1._wp)-1
     end if
     bit_exp=ceiling(log(real(maxexponent(1._wp),wp))/log(2._wp))+1
     if (.not.(2**(bit_exp-1)>=maxexponent(1._wp).and.2**(bit_exp-1)-3>=-minexponent(1._wp))) then
        bit_exp=bit_exp+1
     end if
     pad_exp=0
     do i=0,bit_exp-1
        pad_exp=ibset(pad_exp,i)
     end do
     pad_mantissa=0
     do i=0,bit_mantissa-1
        pad_mantissa=ibset(pad_mantissa,i)
     end do
     pad6=0
     do i=0,5
        pad6=ibset(pad6,i)
     end do

     select case(bit_size(1_mt))
     case(32)
        refString=ref32String
     case(64)
        refString=ref64String
     case default
        refString=" "
     end select
     initialized_module=.true.
  end if
end subroutine init_module

function my_is_nan(r1) result(res)
  real(wp),intent(in) :: r1
  logical :: res
#ifdef __IEEE
  res=ieee_is_nan(r1)
#else
  res=isnan(r1)
#endif
end function my_is_nan

function my_is_finite(r1) result(res)
  real(wp),intent(in) :: r1
  logical :: res
#ifdef __IEEE
  res=ieee_is_finite(r1)
#else
  res=.not.my_is_nan(r1)
  res=res.and.(huge(r1)>=abs(r1))
#endif
end function my_is_finite

function my_is_denormal(r1) result(res)
  real(wp),intent(in) :: r1
  logical :: res
#ifdef __IEEE
  res=my_is_finite(r1).and..not.ieee_is_normal(r1).and.r1/=0._wp
#else
  res=exponent(r1)<minexponent(r1)
  if (my_is_finite(r1)) then
     if (r1/=scale(fraction(r1),exponent(r1)).or.&
          (fraction(abs(r1))<1._wp/real(radix(r1),wp).and.r1/=0._wp)) res=.true.
  end if
#endif
end function my_is_denormal

function my_exponent(r1) result(res)
  real(wp), intent(in) :: r1
  integer(mt) :: res
  res=exponent(r1)
end function my_exponent

! start/end tags

subroutine parse_start_tag(pos,str,bit_mantissa,bit_exp,dims,failure)
  integer, intent(inout) :: pos
  character(len=*), intent(in) :: str
  integer, intent(out) :: bit_mantissa,bit_exp
  integer, dimension(:), pointer :: dims
  logical, intent(inout) :: failure

  integer :: i,j,k,l,m,ii,ij,ndims,iostat,stat

  if (associated(dims)) stop "dims should not be associated"
  i=pos
  bit_mantissa=-1
  bit_exp=-1
  call skip_white(i,str)
  if (i>len(str)) failure=.true.
  if (.not.failure) then
     if (str(i:i)/="<") failure=.true.
  end if
  call skip_white(i,str)
  if (i+6>len(str)) failure=.true.
  if (.not.failure) then
     i=i+1
     if (str(i:i+5)=="reals ".or.str(i:i+5)=="REALS ") then
        i=i+6
     else
        failure=.true.
     end if
  end if
  if (.not.failure) then
!     print *, "parsed '<reals', starting attributes"
     bit_mantissa=-1
     bit_exp=-1
     do
        ! parse attributes

        call skip_white(i,str)
        if (str(i:i)==">") then
           i=i+1
!           print *, "exiting start tag parsing, parsed '",str(pos:i-1),"'"
           exit
        end if
        j=i
        call skip_alnum(j,str)
        if (j==i) then
           print *,"non alnum attr"
           failure=.true.
           exit
        end if
        k=j
        call skip_white(k,str)
        if (k>len(str)) then
           failure=.true.
           exit
        end if
        if (str(k:k)/="=") then
           print *,"no = after attribute"
           failure=.true.
           exit
        else
           k=k+1
        end if
        call skip_white(k,str)
        if (k+1>len(str)) then
           failure=.true.
           exit
        end if
        if (str(k:k)=='"') then
           k=k+1
           l=index(str(k:len(str)),'"')
           if (l<1) then
              print *,"no closing quote"
              failure=.true.
              exit
           end if
           l=l+k-1
           m=l+1
        else
           l=k
           call skip_num(l,str)
           m=l
        end if
!        print *,"attr:",str(i:j-1),"=",str(k:l-1)
        select case(str(i:j-1))
        case("bit_mantissa")
           read(str(k:l-1),*,iostat=iostat)bit_mantissa
           if (iostat/=0) then
              failure=.true.
              exit
           end if
        case("bit_exp","bit_exponent")
           read(str(k:l-1),*,iostat=iostat) bit_exp
           if (iostat/=0) then
              failure=.true.
              exit
           end if
        case("dims")
           ndims=0
           ii=k
           call skip_white(ii,str)
           do
              if (ii>=l) exit
              ij=ii
              call skip_num(ij,str)
              if (ij==ii) then
                 failure=.true.
                 exit
              end if
              ndims=ndims+1
              if (ij>=l) exit
              ii=ij
              call skip_white(ii,str)
              if (ii==ij) then
                 failure=.true.
                 exit
              end if
           end do
           if (failure) then
              print *,"dims parsing failed"
              exit
           end if
           allocate(dims(ndims),stat=stat)
           if (stat/=0) stop "alloc dims"
           ii=k
           ndims=0
           do
              call skip_white(ii,str)
              if (ii>=l) exit
              ij=ii
              call skip_num(ij,str)
              ndims=ndims+1
              read(str(ii:ij-1),*,iostat=iostat)dims(ndims)
              if (iostat/=0) then
                 failure=.true.
                 exit
              end if
              ii=ij
           end do
        end select
        i=m
        if (failure) exit
     end do
  end if
  pos=i
end subroutine parse_start_tag

subroutine skip_end_tag(pos,str,failure)
  character(len=*), intent(in) :: str
  integer, intent(inout) :: pos
  logical, intent(inout) :: failure
  
  call skip_white(pos,str)
  if (pos+1>len(str)) failure=.true.
  if (.not.failure) then
     failure=str(pos:pos+1)=="</"
  end if
  if (.not.failure) then
     call skip_white(pos,str)
     if (pos+5>len(str)) failure=.true.
  end if
  if (.not.failure) then
     if (str(pos:pos+4)=="reals".or.str(pos:pos+4)=="REALS") then
        pos=pos+5
     else
        failure=.true.
     end if
  end if
  call skip_white(pos,str)
  if (pos>len(str)) failure=.true.
  if (.not.failure) then
     if (str(pos:pos)==">") then
        pos=pos+1
     else
        failure=.true.
     end if
  end if
end subroutine skip_end_tag

subroutine skip_white(pos,str)
  integer, intent(inout) :: pos
  character(len=*), intent(in) :: str
  do while (pos<=len(str))
     if (str(pos:pos)/=" ") exit
     pos=pos+1
  end do
end subroutine skip_white

subroutine skip_alnum(j,str)
  integer, intent(inout) :: j
  character(len=*), intent(in) :: str
  integer :: ival
  
  do while(j<=len(str))
     ival=iachar(str(j:j))
     if (.not.(ival>=iachar("A").and.ival<=iachar("Z").or.&
          ival>=iachar("a").and.ival<=iachar("z").or.&
          ival>=iachar("0").and.ival<=iachar("9").or.&
          str(j:j)=="_".or.str(j:j)=="-")) exit
     j=j+1
  end do
end subroutine skip_alnum

subroutine skip_num(j,str)
  integer, intent(inout) :: j
  character(len=*), intent(in) :: str
  integer :: ival
  
  if (j<=len(str)) then
     if (str(j:j)=="-") j=j+1
  end if
  do while(j<=len(str))
     ival=iachar(str(j:j))
     if (.not.(ival>=iachar("0").and.ival<=iachar("9"))) exit
     j=j+1
  end do
end subroutine skip_num

subroutine skip_float(j,str)
  integer, intent(inout) :: j
  character(len=*), intent(in) :: str
  integer :: ival

  if (j<=len(str)) then
     if (str(j:j)=="-".or.str(j:j)=="+") j=j+1
  end if
  if (index("INFTYAinftyaq",str(j:j))>0) then
     do while(j<=len(str))
        if (index("INFTYAQinftyaq",str(j:j))<=0) exit
        j=j+1
     end do
  else
     do while(j<=len(str))
        ival=iachar(str(j:j))
        if (.not.(ival>=iachar("0").and.ival<=iachar("9"))) exit
        j=j+1
     end do
     if (j<len(str)) then
        if (str(j:j)==".") j=j+1
     end if
     do while(j<=len(str))
        ival=iachar(str(j:j))
        if (.not.(ival>=iachar("0").and.ival<=iachar("9"))) exit
        j=j+1
     end do
     if (j<len(str)) then
        if (str(j:j)=="e".or.str(j:j)=="E".or.str(j:j)=="d".or.str(j:j)=="D") j=j+1
     end if
     if (j<=len(str)) then
        if (str(j:j)=="-".or.str(j:j)=="+") j=j+1
     end if
     do while(j<=len(str))
        ival=iachar(str(j:j))
        if (.not.(ival>=iachar("0").and.ival<=iachar("9"))) exit
        j=j+1
     end do
  end if
end subroutine skip_float

subroutine lowercase(str)
  character(len=*), intent(inout) :: str
  integer :: i,ival
  do i=1,len(str)
     ival=iachar(str(i:i))
     if (ival>=iachar("A").and.ival<=iachar("Z")) then
        str(i:i)=achar(ival-iachar("A")+iachar("a"))
     end if
  end do
end subroutine lowercase

end module base64_types

program t
use base64_types
USE machine,           ONLY : m_iargc,m_getarg
implicit none

type(base64c_type) :: b64
integer :: n_arg,i,direction,log_unit,out_unit,in_unit,nFile,nread,&
     ireals,nreals,nproc,bit_mantissa_in,bit_exp_in,iostat,pos,stat,&
     end_str
character(len=200) :: arg_att,infile,outfile,my_format,cmndName,str
logical :: error,exists,ignore_garbage,failure,echo,at_end
real(wp) :: r1
integer, dimension(:), pointer :: dims

failure=.false.
n_arg=m_iargc()
log_unit=6
CALL m_getarg(0,cmndName)

direction=0
nFile=0
infile="in"
outfile=" "
error=.false.
my_format=" "
ignore_garbage=.false.

i=1
echo=.false.
do while(i<=n_arg)
   call m_getarg(i,arg_att)
   select case(arg_att)
   case ("--from-base64")
      direction=-1
   case ("--to-base64")
      direction=1
   case ("--echo")
      echo=.true.
   case ("--format")
      i=i+1
      if (i<=n_arg) then
         call m_getarg(i,my_format)
         print *,"my_format",my_format
      else
         error=.true.
         write (log_unit,"(a)") "ERROR --format should be followed by the format string"
      end if
   case("--ignore-garbage")
      ignore_garbage=.true.
   case default
      select case(nFile)
      case(0)
         infile=arg_att
         nFile=1
      case(1)
         outfile=arg_att
         nFile=2
      case default
         write (log_unit,"(a,a)") "ERROR unexpected extra argument",trim(arg_att)
         error=.true.
      end select
   end select
   i=i+1
end do
if (direction/=0) then
   inquire(file=infile,exist=exists)
   if (.not.exists) then
      write (log_unit,"(a,a,a)") "ERROR infile '",trim(infile),"' does not exist"
      error=.true.
   end if
end if
if (error.or.n_arg<1) then
   WRITE (log_unit,"(a,a)") trim(cmndName),&
        " [--to-base64] [--from-base64] [--format f90OutputFormat]"
   write (log_unit,"(a,/)") "[--ignore-garbage] [--echo] [infile] [outfile]"
   WRITE (log_unit,"(a)") "converts to/from base64 format."
   WRITE (log_unit,"(a,/,a,/,a,/,a,/)")&
        "Base64 format begins with <reals bit_mantissa=xx bit_expo=xx>:",&
        " xx is bit size of the mantissa and the bit size of the exponent",&
        " (i.e. 52 11 for double precision and 23 8 for single precision).",&
        " The bit pattern should be stored in the IEEE 754 format in big endian sequence."
end if
if (.not.error) then
   if (direction/=0) then
      in_unit=10
      out_unit=11
      open(unit=in_unit,file=trim(infile),form="FORMATTED",action="READ")
      if (outfile==" ") then
         out_unit=6
      else
         open(unit=out_unit,file=trim(outfile),form="FORMATTED",&
              status="UNKNOWN",action="WRITE")
      end if
   end if
   ireals=0
   if (direction==-1) then
      read(in_unit,"(a100)",advance="no",iostat=iostat,size=nread) str
      if (iostat>0) then
         write(log_unit,"(a,/,a,/)")&
              'ERROR', 'could not read header <reals ...>'
         failure=.true.
      end if
      if (.not.failure) then
         pos=1
         nullify(dims)
         call parse_start_tag(pos,str,bit_mantissa_in,bit_exp_in,dims,failure)
         if (failure.or.bit_mantissa_in<1.or.bit_exp_in<1) then
            write(log_unit,"(a,/,a,/,a,/,a,/)")&
                 'ERROR',&
                 ' could not parse header <reals ...> and find the number of bits of mantissa and exponent',&
                 ' (double:<reals bit_mantissa="52" bit_expo="11">,',&
                 '  single:<reals bit_mantissa="23" bit_expo="8">)'
            failure=.true.
         end if
      end if
      if (.not.failure) then
         if (associated(dims)) then
            write(out_unit,"('# dims ')",advance="no")
            do i=1,size(dims)
               write(out_unit,"(' ',i6)",advance="no") dims(i)
            end do
            write(out_unit,"()")
            deallocate(dims,stat=stat)
            if (stat/=0) stop "dealloc dims"
         end if
         call base64c_init(b64,bit_mantissa_in=bit_mantissa_in,&
              bit_exp_in=bit_exp_in)
            nproc=base64c_add_string(b64,str(pos:nread),ignore_garbage=ignore_garbage)

         do
            read(in_unit,"(a100)",advance="no",eor=10,end=120,size=nread)str
10          continue
            nproc=base64c_add_string(b64,str(1:nread),ignore_garbage=ignore_garbage)
            if (nproc<len_trim(str(1:nread))) exit
            do
               call base64c_capacity_info(b64,reals_contained=nreals)
               if (nreals<1) exit
               r1=base64c_decode_real(b64)
               ireals=ireals+1
               if (my_format/=" ") then
                  write(out_unit,my_format) r1
               else
                  write(out_unit,"(es35.20e5)") r1
               end if
            end do
         end do
120      continue
         write(log_unit,"('decoded ',i20,' reals')")ireals
      end if
   else if (direction==1) then
      call base64c_init(b64)
      write(out_unit,"(a,i4,a,i4,a)",advance="no")'<reals bit_mantissa="',&
           get_bit_mantissa(),'" bit_exp="',get_bit_exp(),'"'
      read (in_unit,"(a100)",iostat=iostat)str
      if (iostat>0) then
         stop "error reading input file"
      end if
      pos=1
      call skip_white(pos,str)
      if (pos<=len(str)) then
         if (str(pos:pos)=="#") then
            pos=pos+1
            call skip_white(pos,str)
            if (pos+4<=len(str)) then
               if (str(pos:pos+4)=="dims ") then
                  pos=pos+4
                  write (out_unit,"(a,a,a)",advance="no")&
                       ' dims="',str(pos:len_trim(str)),'"'
               end if
            end if
         else
            backspace(in_unit)
         end if
      end if
      write(out_unit,"('>')")
      ireals=0
      pos=1
      end_str=1
      at_end=.false.
      do
         if (my_format/=" ") then
            read(in_unit,my_format,end=121,iostat=iostat)r1
         else
!            print *,"1,str(",pos,":",end_str-1,")=",str(pos:end_str-1)
            call skip_white(pos,str)
            do while(end_str-pos<50.and..not.at_end)
!               print *,"2,pre-filling,str(",pos,":",end_str-1,")=",str(pos:end_str-1)
               if (pos>=end_str)then
                  pos=1
                  end_str=1
                  str=" "
               else if (pos>1) then
                  str(1:end_str-pos)=str(pos:end_str-1)
                  str(end_str:len(str))=" "
                  end_str=end_str-pos+1
                  pos=1
               end if
               
               at_end=.true.
               read(in_unit,"(a100)",eor=122,end=123,advance="no",size=nread)&
                    str(end_str:end_str+100)
               end_str=end_str+nread
               goto 124
122            end_str=end_str+nread+1
               str(end_str:len(str))=" "
124            continue
               at_end=.false.
               call skip_white(pos,str)
!               print *,"2,post-filling,str(",pos,":",end_str-1,")=",str(pos:end_str-1)
            end do
123         continue
!            print *,"3,post-filling,str(",pos,":",end_str-1,")=",str(pos:end_str-1)
            if (pos>=end_str) exit
            i=pos
            call skip_float(i,str)
            if (i==pos) then
               print *,"non float at",i," ",str(i:i)
               exit
            end if
!            print *,"pre-lowercase",str(pos:i-1)
            call lowercase(str(pos:i-1))
!            print *,"post-lowercase",str(pos:i-1)
            select case(str(pos:i-1))
            case("inf","+inf","infinity","+infinity")
               r1=infinity
            case("-inf","-infinity")
               r1=-infinity
            case("nan","nanq")
               r1=nan
            case default
               read (str(pos:i-1),*,iostat=iostat)r1
               if (iostat/=0) then
                  print *,"ERROR trying to interpret '",str(pos:i),"' as real"
                  stop "IO error"
               end if
            end select
            pos=i
         end if
         ireals=ireals+1
         if (echo) print *,ireals,r1
         call base64c_real(b64,r1)
         if (modulo(ireals,10)==0) call base64c_write(b64,out_unit)
      end do
121   continue
      call base64c_flush(b64)
      call base64c_write(b64,out_unit)
      write(out_unit,"('</reals>')")
      write(log_unit,"('encoded ',i20,' reals')")ireals
   else
      call base64c_init(b64,log_unit=log_unit)
   end if
   call base64c_dealloc_ref(b64)
end if

end program t


