!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Utilities for string manipulations
!> \par History
!>      Adapted compress and uppercase for use in CP2K (JGH)
!>      string_to_integer and integer_to_string added (06.02.2001, MK)
!>      Cleaned (04.01.2004,MK)
!> \author MK & JGH
! *****************************************************************************
MODULE string_utilities

  USE f77_blas
  USE glob_matching,                   ONLY: pattern_match=>string_match
  USE kinds,                           ONLY: default_blank_character

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=1), PARAMETER :: backslash = '\\'
  CHARACTER(LEN=1), PARAMETER :: star      = '*'
  CHARACTER(LEN=1), PARAMETER :: question  = '?'
  CHARACTER(LEN=1), PARAMETER :: newline   = ACHAR(10)

  INTEGER, PARAMETER :: maxlen_entity_name = 10
  CHARACTER(LEN=maxlen_entity_name), DIMENSION(252), PARAMETER :: html_entity_table =&
   (/"&quot;    ","&#34;     ","&amp;     ","&#38;#38; ","&apos;    ","&#39;     ","&lt;      ","&#38;#60; ",&
     "&gt;      ","&#62;     ","&nbsp;    ","&#160;    ","&copy;    ","&#169;    ","&deg;     ","&#176;    ",&
     "&plusmn;  ","&#177;    ","&sup2;    ","&#178;    ","&sup3;    ","&#179;    ","&middot;  ","&#183;    ",&
     "&sup1;    ","&#185;    ","&Auml;    ","&#196;    ","&Aring;   ","&#197;    ","&Ccedil;  ","&#199;    ",&
     "&Egrave;  ","&#200;    ","&Eacute;  ","&#201;    ","&Ecirc;   ","&#202;    ","&Ograve;  ","&#210;    ",&
     "&Oacute;  ","&#211;    ","&Ocirc;   ","&#212;    ","&Ouml;    ","&#214;    ","&times;   ","&#215;    ",&
     "&Uuml;    ","&#220;    ","&agrave;  ","&#224;    ","&aacute;  ","&#225;    ","&acirc;   ","&#226;    ",&
     "&auml;    ","&#228;    ","&aring;   ","&#229;    ","&ccedil;  ","&#231;    ","&egrave;  ","&#232;    ",&
     "&eacute;  ","&#233;    ","&ecirc;   ","&#234;    ","&ograve;  ","&#242;    ","&oacute;  ","&#243;    ",&
     "&ocirc;   ","&#244;    ","&ouml;    ","&#246;    ","&ugrave;  ","&#249;    ","&uacute;  ","&#250;    ",&
     "&ucirc;   ","&#251;    ","&uuml;    ","&#252;    ","&Gamma;   ","&#915;    ","&Delta;   ","&#916;    ",&
     "&Theta;   ","&#920;    ","&Lambda;  ","&#923;    ","&Xi;      ","&#926;    ","&Pi;      ","&#928;    ",&
     "&Sigma;   ","&#931;    ","&Phi;     ","&#934;    ","&Psi;     ","&#936;    ","&Omega;   ","&#937;    ",&
     "&alpha;   ","&#945;    ","&beta;    ","&#946;    ","&gamma;   ","&#947;    ","&delta;   ","&#948;    ",&
     "&epsilon; ","&#949;    ","&zeta;    ","&#950;    ","&eta;     ","&#951;    ","&theta;   ","&#952;    ",&
     "&iota;    ","&#953;    ","&kappa;   ","&#954;    ","&lambda;  ","&#955;    ","&mu;      ","&#956;    ",&
     "&nu;      ","&#957;    ","&xi;      ","&#958;    ","&omicron; ","&#959;    ","&pi;      ","&#960;    ",&
     "&rho;     ","&#961;    ","&sigmaf;  ","&#962;    ","&sigma;   ","&#963;    ","&tau;     ","&#964;    ",&
     "&upsilon; ","&#965;    ","&phi;     ","&#966;    ","&chi;     ","&#967;    ","&psi;     ","&#968;    ",&
     "&omega;   ","&#969;    ","&thetasym;","&#977;    ","&ndash;   ","&#8211;   ","&mdash;   ","&#8212;   ",&
     "&dagger;  ","&#8224;   ","&Dagger;  ","&#8225;   ","&bull;    ","&#8230;   ","&hellip;  ","&#8230;   ",&
     "&permil;  ","&#8240;   ","&larr;    ","&#8592;   ","&uarr;    ","&#8593;   ","&rarr;    ","&#8594;   ",&
     "&darr;    ","&#8595;   ","&harr;    ","&#8596;   ","&lArr;    ","&#8656;   ","&uArr;    ","&#8657;   ",&
     "&rArr;    ","&#8658;   ","&dArr;    ","&#8659;   ","&hArr;    ","&#8660;   ","&forall;  ","&#8704;   ",&
     "&part;    ","&#8706;   ","&nabla;   ","&#8711;   ","&isin;    ","&#8712;   ","&notin;   ","&#8713;   ",&
     "&ni;      ","&#8715;   ","&prop;    ","&#8733;   ","&infin;   ","&#8734;   ","&ang;     ","&#8736;   ",&
     "&and;     ","&#8743;   ","&or;      ","&#8744;   ","&cap;     ","&#8745;   ","&cup;     ","&#8746;   ",&
     "&sim;     ","&#8764;   ","&cong;    ","&#8773;   ","&asymp;   ","&#8776;   ","&ne;      ","&#8800;   ",&
     "&equiv;   ","&#8801;   ","&le;      ","&#8804;   ","&ge;      ","&#8805;   ","&sub;     ","&#8834;   ",&
     "&sup;     ","&#8835;   ","&nsub;    ","&#8836;   ","&sube;    ","&#8838;   ","&supe;    ","&#8839;   ",&
     "&oplus;   ","&#8853;   ","&otimes;  ","&#8855;   ","&perp;    ","&#8869;   ","&sdot;    ","&#8901;   ",&
     "&lang;    ","&#9001;   ","&rang;    ","&#9002;   "/)

  INTEGER, PARAMETER :: maxlen_tag_name = 10
  CHARACTER(LEN=maxlen_tag_name), DIMENSION(38), PARAMETER :: html_tag_table =&
   (/"a         ","b         ","big       ","blockquote","br/       ","code      ","dd        ","del       ",&
     "div       ","dl        ","dt        ","em        ","h1        ","h2        ","h3        ","h4        ",&
     "h5        ","h6        ","hr        ","i         ","ins       ","li        ","ol        ","p         ",&
     "span      ","sub       ","sup       ","table     ","tbody     ","td        ","tfoot     ","th        ",&
     "thead     ","tr        ","tt        ","u         ","ul        ","pre       "/)

  PUBLIC :: ascii_to_string,&
            compress,&
            integer_to_string,&
            is_whitespace,&
            remove_word,&
            s2a, a2s,&
            str_comp,&
            string_to_ascii,&
            substitute_special_xml_tokens,&
            typo_match,&
            uppercase,&
            write_html_tables,&
            xstring

  PUBLIC :: html_entity_table,&
            newline

  INTERFACE s2a
     MODULE PROCEDURE s2a_1,s2a_2,s2a_3, s2a_4,  s2a_5,  s2a_6,  s2a_7,  s2a_8,  s2a_9, &
          s2a_10,s2a_11,s2a_12, s2a_13, s2a_14, s2a_15, s2a_16, s2a_17, s2a_18, s2a_19, &
          s2a_20,s2a_21,s2a_22, s2a_23, s2a_24, s2a_25, s2a_26, s2a_27, s2a_28, s2a_29, &
          s2a_30,s2a_31,s2a_32, s2a_33, s2a_34, s2a_35, s2a_36, s2a_37, s2a_38, s2a_39 ! should be clear how to add more
  END INTERFACE

CONTAINS

! *****************************************************************************
!> \brief returns a non-zero positive value if typo_string equals string apart from a few typos.
!>     It is case sensitive, apart from typos.
!> \note
!>     could maybe be made a bit smarter
!> \par History
!>      02.2006 created [Joost VandeVondele]
! *****************************************************************************
  FUNCTION typo_match(string,typo_string) RESULT(match)
    CHARACTER(LEN=*), INTENT(IN)             :: string, typo_string
    INTEGER                                  :: match

    CHARACTER(LEN=1)                         :: kind
    CHARACTER(LEN=LEN(string))               :: tmp2
    CHARACTER(LEN=LEN(typo_string))          :: tmp
    INTEGER                                  :: i, j

    match=0
    IF (LEN_TRIM(typo_string).LE.4) THEN
       kind=question
    ELSE
       kind=star
    ENDIF
    DO i=1,LEN_TRIM(typo_string)
       DO j=i,LEN_TRIM(typo_string)
          tmp=typo_string
          tmp(i:i)=kind
          tmp(j:j)=kind
          IF (i==j .AND. LEN_TRIM(tmp)>2 ) tmp(i:i)=star
          IF (pattern_match(string=string,pattern=tmp)) match=match+1
       ENDDO
    ENDDO
    IF (LEN_TRIM(string).LE.4) THEN
       kind=question
    ELSE
       kind=star
    ENDIF
    DO i=1,LEN_TRIM(string)
       DO j=i,LEN_TRIM(string)
          tmp2=string
          tmp2(i:i)=kind
          tmp2(j:j)=kind
          IF (i==j .AND. LEN_TRIM(tmp2)>2 ) tmp2(i:i)=star
          IF (pattern_match(string=typo_string,pattern=tmp2)) match=match+1
       ENDDO
    ENDDO

  END FUNCTION typo_match


! *****************************************************************************
!> \brief Converts a character-array into a string
!> \par History
!>      10.2013 created [Ole Schuett]
! *****************************************************************************
  PURE FUNCTION a2s(array) RESULT(string)
    CHARACTER, DIMENSION(:), INTENT(IN)      :: array
    CHARACTER(LEN=SIZE(array))               :: string

    INTEGER                                  :: i

    DO i=1,SIZE(array)
      string(i:i) = array(i)
    ENDDO
  END FUNCTION a2s


! *****************************************************************************
!> \brief converts a bunch of strings of different length to an array of
!>        strings of the same length
!> \note
!>     can be used instead of the illegal (/"12","1234"/) generating
!>     s2a("12","1234").EQ.(/"12  ","1234"/)
!> \par History
!>      11.2004 created [Joost VandeVondele ]
! *****************************************************************************
  PURE FUNCTION s2a_1(s1) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1
    CHARACTER(LEN=1000), DIMENSION(1)        :: a

    a(1)=s1
  END FUNCTION s2a_1
  ! *****************************************************************************
  PURE FUNCTION s2a_2(s1,s2) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2
    CHARACTER(LEN=1000), DIMENSION(2)        :: a

    a(1)=s1; a(2)=s2
  END FUNCTION s2a_2
  ! *****************************************************************************
  PURE FUNCTION s2a_3(s1,s2,s3) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3
    CHARACTER(LEN=1000), DIMENSION(3)        :: a

    a(1)=s1; a(2)=s2; a(3)=s3
  END FUNCTION s2a_3
  ! *****************************************************************************
  PURE FUNCTION s2a_4(s1,s2,s3,s4) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4
    CHARACTER(LEN=1000), DIMENSION(4)        :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4
  END FUNCTION s2a_4
  ! *****************************************************************************
  PURE FUNCTION s2a_5(s1,s2,s3,s4,s5) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5
    CHARACTER(LEN=1000), DIMENSION(5)        :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5
  END FUNCTION s2a_5
  ! *****************************************************************************
  PURE FUNCTION s2a_6(s1,s2,s3,s4,s5,s6) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6
    CHARACTER(LEN=1000), DIMENSION(6)        :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6
  END FUNCTION s2a_6
  ! *****************************************************************************
  PURE FUNCTION s2a_7(s1,s2,s3,s4,s5,s6,s7) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7
    CHARACTER(LEN=1000), DIMENSION(7)        :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
  END FUNCTION s2a_7
  ! *****************************************************************************
  PURE FUNCTION s2a_8(s1,s2,s3,s4,s5,s6,s7,s8) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7, s8
    CHARACTER(LEN=1000), DIMENSION(8)        :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8
  END FUNCTION s2a_8
  ! *****************************************************************************
  PURE FUNCTION s2a_9(s1,s2,s3,s4,s5,s6,s7,s8,s9) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7, &
                                                s8, s9
    CHARACTER(LEN=1000), DIMENSION(9)        :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9
  END FUNCTION s2a_9
  ! *****************************************************************************
  PURE FUNCTION s2a_10(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7, &
                                                s8, s9, s10
    CHARACTER(LEN=1000), DIMENSION(10)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10
  END FUNCTION s2a_10
  ! *****************************************************************************
  PURE FUNCTION s2a_11(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7, &
                                                s8, s9, s10, s11
    CHARACTER(LEN=1000), DIMENSION(11)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11
  END FUNCTION s2a_11
  ! *****************************************************************************
  PURE FUNCTION s2a_12(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7, &
                                                s8, s9, s10, s11, s12
    CHARACTER(LEN=1000), DIMENSION(12)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12
  END FUNCTION s2a_12
  ! *****************************************************************************
  PURE FUNCTION s2a_13(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7, &
                                                s8, s9, s10, s11, s12, s13
    CHARACTER(LEN=1000), DIMENSION(13)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12; a(13)=s13
  END FUNCTION s2a_13
  ! *****************************************************************************
  PURE FUNCTION s2a_14(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7, &
                                                s8, s9, s10, s11, s12, s13, &
                                                s14
    CHARACTER(LEN=1000), DIMENSION(14)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12; a(13)=s13; a(14)=s14
  END FUNCTION s2a_14
  ! *****************************************************************************
  PURE FUNCTION s2a_15(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7, &
                                                s8, s9, s10, s11, s12, s13, &
                                                s14, s15
    CHARACTER(LEN=1000), DIMENSION(15)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12; a(13)=s13; a(14)=s14
    a(15)=s15
  END FUNCTION s2a_15
  ! *****************************************************************************
  PURE FUNCTION s2a_16(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15, &
       s16) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7, &
                                                s8, s9, s10, s11, s12, s13, &
                                                s14, s15, s16
    CHARACTER(LEN=1000), DIMENSION(16)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12; a(13)=s13; a(14)=s14
    a(15)=s15; a(16)=s16
  END FUNCTION s2a_16
  ! *****************************************************************************
  PURE FUNCTION s2a_17(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15, &
       s16,s17) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7, &
                                                s8, s9, s10, s11, s12, s13, &
                                                s14, s15, s16, s17
    CHARACTER(LEN=1000), DIMENSION(17)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12; a(13)=s13; a(14)=s14
    a(15)=s15; a(16)=s16; a(17)=s17
  END FUNCTION s2a_17
  ! *****************************************************************************
  PURE FUNCTION s2a_18(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15, &
       s16,s17,s18) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7, &
                                                s8, s9, s10, s11, s12, s13, &
                                                s14, s15, s16, s17, s18
    CHARACTER(LEN=1000), DIMENSION(18)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12; a(13)=s13; a(14)=s14
    a(15)=s15; a(16)=s16; a(17)=s17; a(18)=s18
  END FUNCTION s2a_18
  ! *****************************************************************************
  PURE FUNCTION s2a_19(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15, &
       s16,s17,s18,s19) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7, &
                                                s8, s9, s10, s11, s12, s13, &
                                                s14, s15, s16, s17, s18, s19
    CHARACTER(LEN=1000), DIMENSION(19)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12; a(13)=s13; a(14)=s14
    a(15)=s15; a(16)=s16; a(17)=s17; a(18)=s18; a(19)=s19
  END FUNCTION s2a_19
  ! *****************************************************************************
  PURE FUNCTION s2a_20(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15, &
       s16,s17,s18,s19,s20) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7, &
                                                s8, s9, s10, s11, s12, s13, &
                                                s14, s15, s16, s17, s18, s19, &
                                                s20
    CHARACTER(LEN=1000), DIMENSION(20)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12; a(13)=s13; a(14)=s14
    a(15)=s15; a(16)=s16; a(17)=s17; a(18)=s18; a(19)=s19; a(20)=s20
  END FUNCTION s2a_20
  ! *****************************************************************************
  PURE FUNCTION s2a_21(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15, &
       s16,s17,s18,s19,s20,s21) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7, &
                                                s8, s9, s10, s11, s12, s13, &
                                                s14, s15, s16, s17, s18, s19, &
                                                s20, s21
    CHARACTER(LEN=1000), DIMENSION(21)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12; a(13)=s13; a(14)=s14
    a(15)=s15; a(16)=s16; a(17)=s17; a(18)=s18; a(19)=s19; a(20)=s20
    a(21)=s21
  END FUNCTION s2a_21
  ! *****************************************************************************
  PURE FUNCTION s2a_22(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15, &
       s16,s17,s18,s19,s20,s21,s22) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7, &
                                                s8, s9, s10, s11, s12, s13, &
                                                s14, s15, s16, s17, s18, s19, &
                                                s20, s21, s22
    CHARACTER(LEN=1000), DIMENSION(22)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12; a(13)=s13; a(14)=s14
    a(15)=s15; a(16)=s16; a(17)=s17; a(18)=s18; a(19)=s19; a(20)=s20
    a(21)=s21; a(22)=s22
  END FUNCTION s2a_22
  ! *****************************************************************************
  PURE FUNCTION s2a_23(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15, &
       s16,s17,s18,s19,s20,s21,s22,s23) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7, &
                                                s8, s9, s10, s11, s12, s13, &
                                                s14, s15, s16, s17, s18, s19, &
                                                s20, s21, s22, s23
    CHARACTER(LEN=1000), DIMENSION(23)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12; a(13)=s13; a(14)=s14
    a(15)=s15; a(16)=s16; a(17)=s17; a(18)=s18; a(19)=s19; a(20)=s20
    a(21)=s21; a(22)=s22; a(23)=s23
  END FUNCTION s2a_23
  ! *****************************************************************************
  PURE FUNCTION s2a_24(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15, &
       s16,s17,s18,s19,s20,s21,s22,s23,s24) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7, &
                                                s8, s9, s10, s11, s12, s13, &
                                                s14, s15, s16, s17, s18, s19, &
                                                s20, s21, s22, s23, s24
    CHARACTER(LEN=1000), DIMENSION(24)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12; a(13)=s13; a(14)=s14
    a(15)=s15; a(16)=s16; a(17)=s17; a(18)=s18; a(19)=s19; a(20)=s20
    a(21)=s21; a(22)=s22; a(23)=s23; a(24)=s24
  END FUNCTION s2a_24
  ! *****************************************************************************
  PURE FUNCTION s2a_25(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15, &
       s16,s17,s18,s19,s20,s21,s22,s23,s24,s25) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7, &
                                                s8, s9, s10, s11, s12, s13, &
                                                s14, s15, s16, s17, s18, s19, &
                                                s20, s21, s22, s23, s24, s25
    CHARACTER(LEN=1000), DIMENSION(25)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12; a(13)=s13; a(14)=s14
    a(15)=s15; a(16)=s16; a(17)=s17; a(18)=s18; a(19)=s19; a(20)=s20
    a(21)=s21; a(22)=s22; a(23)=s23; a(24)=s24; a(25)=s25
  END FUNCTION s2a_25
  ! *****************************************************************************
  PURE FUNCTION s2a_26(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15, &
       s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7, &
                                                s8, s9, s10, s11, s12, s13, &
                                                s14, s15, s16, s17, s18, s19, &
                                                s20, s21, s22, s23, s24, s25, &
                                                s26
    CHARACTER(LEN=1000), DIMENSION(26)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12; a(13)=s13; a(14)=s14
    a(15)=s15; a(16)=s16; a(17)=s17; a(18)=s18; a(19)=s19; a(20)=s20
    a(21)=s21; a(22)=s22; a(23)=s23; a(24)=s24; a(25)=s25; a(26)=s26
  END FUNCTION s2a_26
  ! *****************************************************************************
  PURE FUNCTION s2a_27(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15, &
       s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7, &
                                                s8, s9, s10, s11, s12, s13, &
                                                s14, s15, s16, s17, s18, s19, &
                                                s20, s21, s22, s23, s24, s25, &
                                                s26, s27
    CHARACTER(LEN=1000), DIMENSION(27)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12; a(13)=s13; a(14)=s14
    a(15)=s15; a(16)=s16; a(17)=s17; a(18)=s18; a(19)=s19; a(20)=s20
    a(21)=s21; a(22)=s22; a(23)=s23; a(24)=s24; a(25)=s25; a(26)=s26
    a(27)=s27
  END FUNCTION s2a_27
  ! *****************************************************************************
  PURE FUNCTION s2a_28(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15, &
       s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7, &
                                                s8, s9, s10, s11, s12, s13, &
                                                s14, s15, s16, s17, s18, s19, &
                                                s20, s21, s22, s23, s24, s25, &
                                                s26, s27, s28
    CHARACTER(LEN=1000), DIMENSION(28)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12; a(13)=s13; a(14)=s14
    a(15)=s15; a(16)=s16; a(17)=s17; a(18)=s18; a(19)=s19; a(20)=s20
    a(21)=s21; a(22)=s22; a(23)=s23; a(24)=s24; a(25)=s25; a(26)=s26
    a(27)=s27; a(28)=s28
  END FUNCTION s2a_28
  ! *****************************************************************************
  PURE FUNCTION s2a_29(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15, &
       s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN) :: s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, &
      s11, s12, s13, s14, s15, s16, s17, s18, s19, s20, s21, s22, s23, s24, &
      s25, s26, s27, s28, s29
    CHARACTER(LEN=1000), DIMENSION(29)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12; a(13)=s13; a(14)=s14
    a(15)=s15; a(16)=s16; a(17)=s17; a(18)=s18; a(19)=s19; a(20)=s20
    a(21)=s21; a(22)=s22; a(23)=s23; a(24)=s24; a(25)=s25; a(26)=s26
    a(27)=s27; a(28)=s28; a(29)=s29
  END FUNCTION s2a_29
  ! *****************************************************************************
  PURE FUNCTION s2a_30(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15, &
       s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29,s30) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN) :: s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, &
      s11, s12, s13, s14, s15, s16, s17, s18, s19, s20, s21, s22, s23, s24, &
      s25, s26, s27, s28, s29, s30
    CHARACTER(LEN=1000), DIMENSION(30)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12; a(13)=s13; a(14)=s14
    a(15)=s15; a(16)=s16; a(17)=s17; a(18)=s18; a(19)=s19; a(20)=s20
    a(21)=s21; a(22)=s22; a(23)=s23; a(24)=s24; a(25)=s25; a(26)=s26
    a(27)=s27; a(28)=s28; a(29)=s29; a(30)=s30
  END FUNCTION s2a_30
  ! *****************************************************************************
  PURE FUNCTION s2a_31(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15, &
       s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29,s30,s31) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN) :: s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, &
      s11, s12, s13, s14, s15, s16, s17, s18, s19, s20, s21, s22, s23, s24, &
      s25, s26, s27, s28, s29, s30, s31
    CHARACTER(LEN=1000), DIMENSION(31)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12; a(13)=s13; a(14)=s14
    a(15)=s15; a(16)=s16; a(17)=s17; a(18)=s18; a(19)=s19; a(20)=s20
    a(21)=s21; a(22)=s22; a(23)=s23; a(24)=s24; a(25)=s25; a(26)=s26
    a(27)=s27; a(28)=s28; a(29)=s29; a(30)=s30; a(31)=s31
  END FUNCTION s2a_31
  ! *****************************************************************************
  PURE FUNCTION s2a_32(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15, &
       s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29,s30,s31,s32) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN) :: s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, &
      s11, s12, s13, s14, s15, s16, s17, s18, s19, s20, s21, s22, s23, s24, &
      s25, s26, s27, s28, s29, s30, s31, s32
    CHARACTER(LEN=1000), DIMENSION(32)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12; a(13)=s13; a(14)=s14
    a(15)=s15; a(16)=s16; a(17)=s17; a(18)=s18; a(19)=s19; a(20)=s20
    a(21)=s21; a(22)=s22; a(23)=s23; a(24)=s24; a(25)=s25; a(26)=s26
    a(27)=s27; a(28)=s28; a(29)=s29; a(30)=s30; a(31)=s31; a(32)=s32
  END FUNCTION s2a_32
  ! *****************************************************************************
  PURE FUNCTION s2a_33(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15, &
       s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29,s30,s31,s32, &
       s33) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN) :: s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, &
      s11, s12, s13, s14, s15, s16, s17, s18, s19, s20, s21, s22, s23, s24, &
      s25, s26, s27, s28, s29, s30, s31, s32, s33
    CHARACTER(LEN=1000), DIMENSION(33)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12; a(13)=s13; a(14)=s14
    a(15)=s15; a(16)=s16; a(17)=s17; a(18)=s18; a(19)=s19; a(20)=s20
    a(21)=s21; a(22)=s22; a(23)=s23; a(24)=s24; a(25)=s25; a(26)=s26
    a(27)=s27; a(28)=s28; a(29)=s29; a(30)=s30; a(31)=s31; a(32)=s32
    a(33)=s33
  END FUNCTION s2a_33
  ! *****************************************************************************
  PURE FUNCTION s2a_34(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15, &
       s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29,s30,s31,s32, &
       s33,s34) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN) :: s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, &
      s11, s12, s13, s14, s15, s16, s17, s18, s19, s20, s21, s22, s23, s24, &
      s25, s26, s27, s28, s29, s30, s31, s32, s33, s34
    CHARACTER(LEN=1000), DIMENSION(34)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12; a(13)=s13; a(14)=s14
    a(15)=s15; a(16)=s16; a(17)=s17; a(18)=s18; a(19)=s19; a(20)=s20
    a(21)=s21; a(22)=s22; a(23)=s23; a(24)=s24; a(25)=s25; a(26)=s26
    a(27)=s27; a(28)=s28; a(29)=s29; a(30)=s30; a(31)=s31; a(32)=s32
    a(33)=s33; a(34)=s34
  END FUNCTION s2a_34
  ! *****************************************************************************
  PURE FUNCTION s2a_35(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15, &
       s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29,s30,s31,s32, &
       s33,s34,s35) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN) :: s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, &
      s11, s12, s13, s14, s15, s16, s17, s18, s19, s20, s21, s22, s23, s24, &
      s25, s26, s27, s28, s29, s30, s31, s32, s33, s34, s35
    CHARACTER(LEN=1000), DIMENSION(35)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12; a(13)=s13; a(14)=s14
    a(15)=s15; a(16)=s16; a(17)=s17; a(18)=s18; a(19)=s19; a(20)=s20
    a(21)=s21; a(22)=s22; a(23)=s23; a(24)=s24; a(25)=s25; a(26)=s26
    a(27)=s27; a(28)=s28; a(29)=s29; a(30)=s30; a(31)=s31; a(32)=s32
    a(33)=s33; a(34)=s34; a(35)=s35
  END FUNCTION s2a_35
  ! *****************************************************************************
  PURE FUNCTION s2a_36(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15, &
       s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29,s30,s31,s32, &
       s33,s34,s35,s36) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN) :: s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, &
      s11, s12, s13, s14, s15, s16, s17, s18, s19, s20, s21, s22, s23, s24, &
      s25, s26, s27, s28, s29, s30, s31, s32, s33, s34, s35, s36
    CHARACTER(LEN=1000), DIMENSION(36)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12; a(13)=s13; a(14)=s14
    a(15)=s15; a(16)=s16; a(17)=s17; a(18)=s18; a(19)=s19; a(20)=s20
    a(21)=s21; a(22)=s22; a(23)=s23; a(24)=s24; a(25)=s25; a(26)=s26
    a(27)=s27; a(28)=s28; a(29)=s29; a(30)=s30; a(31)=s31; a(32)=s32
    a(33)=s33; a(34)=s34; a(35)=s35; a(36)=s36
  END FUNCTION s2a_36
  ! *****************************************************************************
  PURE FUNCTION s2a_37(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15, &
       s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29,s30,s31,s32, &
       s33,s34,s35,s36,s37) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN) :: s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, &
      s11, s12, s13, s14, s15, s16, s17, s18, s19, s20, s21, s22, s23, s24, &
      s25, s26, s27, s28, s29, s30, s31, s32, s33, s34, s35, s36, s37
    CHARACTER(LEN=1000), DIMENSION(37)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12; a(13)=s13; a(14)=s14
    a(15)=s15; a(16)=s16; a(17)=s17; a(18)=s18; a(19)=s19; a(20)=s20
    a(21)=s21; a(22)=s22; a(23)=s23; a(24)=s24; a(25)=s25; a(26)=s26
    a(27)=s27; a(28)=s28; a(29)=s29; a(30)=s30; a(31)=s31; a(32)=s32
    a(33)=s33; a(34)=s34; a(35)=s35; a(36)=s36; a(37)=s37
  END FUNCTION s2a_37
  ! *****************************************************************************
  PURE FUNCTION s2a_38(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15, &
       s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29,s30,s31,s32, &
       s33,s34,s35,s36,s37,s38) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN) :: s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, &
      s11, s12, s13, s14, s15, s16, s17, s18, s19, s20, s21, s22, s23, s24, &
      s25, s26, s27, s28, s29, s30, s31, s32, s33, s34, s35, s36, s37, s38
    CHARACTER(LEN=1000), DIMENSION(38)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12; a(13)=s13; a(14)=s14
    a(15)=s15; a(16)=s16; a(17)=s17; a(18)=s18; a(19)=s19; a(20)=s20
    a(21)=s21; a(22)=s22; a(23)=s23; a(24)=s24; a(25)=s25; a(26)=s26
    a(27)=s27; a(28)=s28; a(29)=s29; a(30)=s30; a(31)=s31; a(32)=s32
    a(33)=s33; a(34)=s34; a(35)=s35; a(36)=s36; a(37)=s37; a(38)=s38
  END FUNCTION s2a_38
  ! *****************************************************************************
  PURE FUNCTION s2a_39(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15, &
       s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29,s30,s31,s32, &
       s33,s34,s35,s36,s37,s38,s39) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN) :: s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, &
      s11, s12, s13, s14, s15, s16, s17, s18, s19, s20, s21, s22, s23, s24, &
      s25, s26, s27, s28, s29, s30, s31, s32, s33, s34, s35, s36, s37, s38, &
      s39
    CHARACTER(LEN=1000), DIMENSION(39)       :: a

    a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
    a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12; a(13)=s13; a(14)=s14
    a(15)=s15; a(16)=s16; a(17)=s17; a(18)=s18; a(19)=s19; a(20)=s20
    a(21)=s21; a(22)=s22; a(23)=s23; a(24)=s24; a(25)=s25; a(26)=s26
    a(27)=s27; a(28)=s28; a(29)=s29; a(30)=s30; a(31)=s31; a(32)=s32
    a(33)=s33; a(34)=s34; a(35)=s35; a(36)=s36; a(37)=s37; a(38)=s38
    a(39)=s39
  END FUNCTION s2a_39
! *****************************************************************************
!> \brief Convert a sequence of integer numbers (ASCII code) to a string.
!>         Blanks are inserted for invalid ASCII code numbers.
!> \author  MK
!> \date    19.10.2000
!> \version 1.0
! *****************************************************************************
  SUBROUTINE ascii_to_string(nascii,string)

    INTEGER, DIMENSION(:), INTENT(IN)        :: nascii
    CHARACTER(LEN=*), INTENT(OUT)            :: string

    INTEGER                                  :: i

    string = ""

    DO i=1,MIN(LEN(string),SIZE(nascii))
       IF ((nascii(i) >= 0).AND.(nascii(i) <= 127)) THEN
          string(i:i) = CHAR(nascii(i))
       ELSE
          string(i:i) = " "
       END IF
    END DO

  END SUBROUTINE ascii_to_string

! *****************************************************************************
!> \brief   Eliminate multiple space characters in a string.
!>          If full is .TRUE., then all spaces are eliminated.
!> \author  MK
!> \date    23.06.1998
!> \version 1.0
! *****************************************************************************
  SUBROUTINE compress(string,full)

    CHARACTER(LEN=*), INTENT(INOUT)          :: string
    LOGICAL, INTENT(IN), OPTIONAL            :: full

    CHARACTER                                :: tmp
    INTEGER                                  :: i, z
    LOGICAL                                  :: remove_all

    IF (PRESENT(full)) THEN
       remove_all = full
    ELSE
       remove_all = .FALSE.
    END IF

    z = 1

    DO i=1,LEN_TRIM(string)
       IF ((z == 1).OR.remove_all) THEN
          IF (string(i:i) /= " ") THEN
             tmp = string(i:i)
             string(z:z) = tmp
             z = z + 1
          END IF
       ELSE
          IF ((string(i:i) /= " ").OR.(string(z-1:z-1) /= " ")) THEN
             tmp = string(i:i)
             string(z:z) = tmp
             z = z + 1
          END IF
       END IF
    END DO

    string(z:) = ""

  END SUBROUTINE compress

! *****************************************************************************
!> \brief   Converts an integer number to a string.
!>          The WRITE statement will return an error message, if the number of
!>          digits of the integer number is larger the than the length of the
!>          supplied string.
!> \author  MK
!> \date    05.01.2004
!> \version 1.0
! *****************************************************************************
  SUBROUTINE integer_to_string(inumber,string)

    INTEGER, INTENT(IN)                      :: inumber
    CHARACTER(LEN=*), INTENT(OUT)            :: string

    WRITE (UNIT=string,FMT='(I0)') inumber
  END SUBROUTINE integer_to_string

! *****************************************************************************
!> \brief   Convert a string to sequence of integer numbers.
!> \author  MK
!> \date    19.10.2000
!> \version 1.0
! *****************************************************************************
  SUBROUTINE string_to_ascii(string,nascii)

    CHARACTER(LEN=*), INTENT(IN)             :: string
    INTEGER, DIMENSION(:), INTENT(OUT)       :: nascii

    INTEGER                                  :: i

    nascii(:) = 0

    DO i=1,MIN(LEN(string),SIZE(nascii))
       nascii(i) = ICHAR(string(i:i))
    END DO

  END SUBROUTINE string_to_ascii

! *****************************************************************************
!> \brief   remove a word from a string (words are separated by white spaces)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE remove_word(string)
    CHARACTER(LEN=*), INTENT(INOUT)          :: string

    INTEGER                                  :: i

    i = 1
    ! possibly clean white spaces
    DO WHILE(string(i:i)==" ")
       i = i + 1
    END DO
    ! now remove the word
    DO WHILE(string(i:i)/=" ")
       i = i + 1
    END DO
    string = string(i:)

  END SUBROUTINE remove_word

! *****************************************************************************
!> \brief  Substitute special XML tokens like "<" or ">" in inp_string.
!>         Optionally convert also all lowercase characters to uppercase, if
!>         ltu is true.
!> \author  Matthias Krack (MK)
!> \date    10.03.2005
!> \version 1.1
!> \par History
!>      - Enable the use of HTML entity names (06.03.13,MK)
! *****************************************************************************
  FUNCTION substitute_special_xml_tokens(inp_string,ltu) RESULT(out_string)

    CHARACTER(LEN=*), INTENT(IN)             :: inp_string
    LOGICAL, INTENT(IN), OPTIONAL            :: ltu
    CHARACTER(LEN=2*LEN(inp_string))         :: out_string

    CHARACTER(LEN=LEN(inp_string))           :: string
    CHARACTER(LEN=maxlen_entity_name)        :: entity_name
    CHARACTER(LEN=maxlen_tag_name)           :: tag_name
    INTEGER                                  :: i, ientry, ilen, j, k

    string = inp_string
    out_string = ""

    IF (PRESENT(ltu)) THEN
       IF (ltu) CALL uppercase(string)
    END IF

    i = 0
    j = 1
    string_loop: DO
       i = i + 1
       IF (i > LEN_TRIM(string)) EXIT string_loop
       IF (string(i:i) == "<") THEN
          ! Detect valid HTML tags and keep them
          ientry = 0
          ilen = INDEX(string(i:),">")
          IF ((ilen > 2).AND.(ilen <= maxlen_tag_name + 3)) THEN
             IF (string(i+1:i+1) == "/") THEN
                tag_name(1:) = string(i+2:i+ilen-2)
             ELSE
                tag_name(1:) = string(i+1:i+ilen-2)
             END IF
             CALL lowercase(tag_name)
             tag_loop: DO k=1,SIZE(html_tag_table)
                IF (tag_name == html_tag_table(k)) THEN
                   ientry = k
                   EXIT tag_loop
                END IF
             END DO tag_loop
             IF (ientry > 0) THEN ! HTML tag found in table
                IF (string(i+1:i+1) == "/") THEN
                   out_string(j:j+ilen+7) = "&#60;/"//TRIM(tag_name)//"&#62;"
                 ELSE
                   out_string(j:j+ilen+7) = "&#60;"//TRIM(tag_name)//"&#62;"
                END IF
                i = i + ilen - 1
                j = j + ilen + 8
             END IF
          END IF
          ! HTML tag not found in table
          IF (ientry == 0) THEN
             out_string(j:j+4) = "&#60;"
             j = j + 5
          END IF
       ELSE IF (string(i:i) == ">") THEN
          out_string(j:j+4) = "&#62;"
          j = j + 5
       ELSE IF (string(i:i) == "&") THEN
          ! Substitute HTML entity names by the corresponding entity number
          ientry = 0
          ilen = INDEX(string(i:),";")
          IF ((ilen > 2).AND.(ilen <= maxlen_entity_name)) THEN
             entity_name(1:) = string(i:i+ilen-1)
             entity_loop: DO k=1,SIZE(html_entity_table),2
                IF (entity_name == html_entity_table(k)) THEN
                   ientry = k + 1
                   EXIT entity_loop
                END IF
             END DO entity_loop
             i = i + ilen - 1
             IF (ientry > 0) THEN ! HTML entity found in table
                ilen = LEN_TRIM(html_entity_table(ientry))
                out_string(j:j+ilen-1) = TRIM(ADJUSTL(html_entity_table(ientry)))
                j = j + ilen
             END IF
          END IF
          IF (ientry == 0) THEN
             out_string(j:j+4) = "&#38;"
             j = j + 5
          END IF
       ELSE
          out_string(j:j) = string(i:i)
          j = j + 1
       END IF
    END DO string_loop

  END  FUNCTION substitute_special_xml_tokens

! *****************************************************************************
!> \brief   Write the HTML entity and HTML tag table to unit output_unit
!> \author  Matthias Krack (MK)
!> \date    14.08.2013
!> \version 1.0
! *****************************************************************************
  SUBROUTINE write_html_tables(output_unit)

    INTEGER, INTENT(IN)                      :: output_unit

    INTEGER                                  :: i

    WRITE (UNIT=output_unit,FMT="(T2,A)") "<CP2K_HTML>"

    DO i=1,SIZE(html_tag_table)
       WRITE (UNIT=output_unit,FMT="(T3,A)")&
        "<TAG>",&
        " <NAME>&#60;"//TRIM(html_tag_table(i))//"&#62;</NAME>",&
        "</TAG>"
    END DO
    DO i=1,SIZE(html_entity_table),2
       WRITE (UNIT=output_unit,FMT="(T3,A)")&
        "<ENTITY>",&
        " <NAME>&#38;"//TRIM(html_entity_table(i)(2:))//"</NAME>",&
        " <CODE>&#38;"//TRIM(html_entity_table(i+1)(2:))//"</CODE>",&
        "</ENTITY>"
    END DO

    WRITE (UNIT=output_unit,FMT="(T2,A)") "</CP2K_HTML>"

  END SUBROUTINE write_html_tables

! *****************************************************************************
!> \brief   Convert all upper case characters in a string to lower case.
!> \author  Matthias Krack (MK)
!> \date    14.08.2013
!> \version 1.0
! *****************************************************************************
  SUBROUTINE lowercase(string)
    CHARACTER(LEN=*), INTENT(INOUT)          :: string

    INTEGER                                  :: i, iascii

    DO i=1,LEN_TRIM(string)
       iascii = ICHAR(string(i:i))
       IF ((iascii >= 65).AND.(iascii <= 90)) THEN
          string(i:i) = CHAR(iascii + 32)
       END IF
    END DO

  END SUBROUTINE lowercase

! *****************************************************************************
!> \brief   Convert all lower case characters in a string to upper case.
!> \author  MK
!> \date    22.06.1998
!> \version 1.0
! *****************************************************************************
  SUBROUTINE uppercase(string)
    CHARACTER(LEN=*), INTENT(INOUT)          :: string

    INTEGER                                  :: i, iascii

    DO i=1,LEN_TRIM(string)
       iascii = ICHAR(string(i:i))
       IF ((iascii >= 97).AND.(iascii <= 122)) THEN
          string(i:i) = CHAR(iascii - 32)
       END IF
    END DO

  END SUBROUTINE uppercase

  ! *****************************************************************************
  SUBROUTINE xstring(string,ia,ib)

    CHARACTER(LEN=*), INTENT(IN)             :: string
    INTEGER, INTENT(OUT)                     :: ia, ib

    ia = 1
    ib = LEN_TRIM(string)
    IF (ib>0) THEN
       DO WHILE (string(ia:ia)==' ')
          ia = ia + 1
       END DO
    END IF

  END SUBROUTINE xstring

  ! *****************************************************************************
  SUBROUTINE make_tuple(int,nt,na,tuple)

    INTEGER, INTENT(IN)                      :: INT( :, : ), nt, na
    CHARACTER(LEN=*), INTENT(OUT)            :: tuple( : )

    INTEGER                                  :: i, nm

    nm = MAXVAL(INT(1:nt,1:na))
    SELECT CASE (nt)
    CASE DEFAULT
       STOP 'make_tuple: case not programmed'
    CASE (1)
       IF (nm<100) THEN
          DO i = 1, na
             WRITE (tuple(i),'(A,I2,A )' ) '[', INT(1,i), ']'
          END DO
       ELSE
          DO i = 1, na
             WRITE (tuple(i),'(A,I4,A )' ) '[', INT(1,i), ']'
          END DO
       END IF
    CASE (2)
       IF (nm<100) THEN
          DO i = 1, na
             WRITE (tuple(i),'(A,I2,A,I2,A )' ) '[', INT(1,i), '-', INT(2,i), &
                  ']'
          END DO
       ELSE
          DO i = 1, na
             WRITE (tuple(i),'(A,I4,A,I4,A )' ) '[', INT(1,i), '-', INT(2,i), &
                  ']'
          END DO
       END IF
    CASE (3)
       IF (nm<100) THEN
          DO i = 1, na
             WRITE (tuple(i),'(A,I2,A,I2,A,I2,A )' ) '[', INT(1,i), '-', &
                  INT(2,i), '-', INT(3,i), ']'
          END DO
       ELSE
          DO i = 1, na
             WRITE (tuple(i),'(A,I4,A,I4,A,I4,A )' ) '[', INT(1,i), '-', &
                  INT(2,i), '-', INT(3,i), ']'
          END DO
       END IF
    CASE (4)
       IF (nm<100) THEN
          DO i = 1, na
             WRITE (tuple(i),'(A,I2,A,I2,A,I2,A,I2,A )' ) '[', INT(1,i), '-', &
                  INT(2,i), '-', INT(3,i), '-', INT(4,i), ']'
          END DO
       ELSE
          DO i = 1, na
             WRITE (tuple(i),'(A,I4,A,I4,A,I4,A,I4,A )' ) '[', INT(1,i), '-', &
                  INT(2,i), '-', INT(3,i), '-', INT(4,i), ']'
          END DO
       END IF
    END SELECT

  END SUBROUTINE make_tuple

  ! *****************************************************************************
  FUNCTION str_comp(str1,str2) RESULT (equal)

    CHARACTER(LEN=*), INTENT(IN)             :: str1, str2
    LOGICAL                                  :: equal

    INTEGER                                  :: i1, i2, j1, j2

    i1 = 0
    i2 = 0
    j1 = 0
    j2 = 0
    CALL xstring(str1,i1,i2)
    CALL xstring(str2,j1,j2)
    equal = (str1(i1:i2)==str2(j1:j2))
  END FUNCTION str_comp

  ! *****************************************************************************
  FUNCTION str_search(str1,n,str2) RESULT (pos)
    CHARACTER(LEN=*), INTENT(IN)             :: str1( : )
    INTEGER, INTENT(IN)                      :: n
    CHARACTER(LEN=*), INTENT(IN)             :: str2
    INTEGER                                  :: pos

    INTEGER                                  :: i

    pos = 0
    DO i = 1, n
       IF (str_comp(str1(i),str2)) THEN
          pos = i
          EXIT
       END IF
    END DO
  END FUNCTION str_search

! *****************************************************************************
!> \brief returns .true. if the character passed is a whitespace char.
!> \par History
!>      02.2008 created, AK
! *****************************************************************************
  FUNCTION is_whitespace(testchar) RESULT(resval)
    CHARACTER(LEN=1), INTENT(IN)             :: testchar
    LOGICAL                                  :: resval

    resval=.FALSE.
    IF (ANY(default_blank_character==testchar)) resval = .TRUE.
  END FUNCTION is_whitespace

END MODULE string_utilities
