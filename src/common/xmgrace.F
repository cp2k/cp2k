!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to facilitate writing XMGRACE files
!> \par History
!>      none
!> \author JGH (10.02.2025)
! **************************************************************************************************
MODULE xmgrace

   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_datum
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xmgrace'

   PUBLIC :: xm_write_defaults, xm_write_frameport, xm_write_frame, xm_graph_info, &
             xm_graph_data

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param iw ...
! **************************************************************************************************
   SUBROUTINE xm_write_defaults(iw)
      INTEGER, INTENT(IN)                                :: iw

      CHARACTER(len=20)                                  :: date

      IF (iw > 0) THEN
         WRITE (iw, '(A)') '# CP2K Grace file', '#', '@version 50125', '@page size 792, 612', &
            '@page scroll 5%', '@page inout 5%', '@link page off'
         WRITE (iw, '(A)') '@map font 0 to "Times-Roman", "Times-Roman"', &
            '@map font 1 to "Times-Italic", "Times-Italic"', &
            '@map font 2 to "Times-Bold", "Times-Bold"', &
            '@map font 3 to "Times-BoldItalic", "Times-BoldItalic"', &
            '@map font 4 to "Helvetica", "Helvetica"', &
            '@map font 5 to "Helvetica-Oblique", "Helvetica-Oblique"', &
            '@map font 6 to "Helvetica-Bold", "Helvetica-Bold"', &
            '@map font 7 to "Helvetica-BoldOblique", "Helvetica-BoldOblique"', &
            '@map font 8 to "Courier", "Courier"', &
            '@map font 9 to "Courier-Oblique", "Courier-Oblique"', &
            '@map font 10 to "Courier-Bold", "Courier-Bold"', &
            '@map font 11 to "Courier-BoldOblique", "Courier-BoldOblique"', &
            '@map font 12 to "Symbol", "Symbol"', &
            '@map font 13 to "ZapfDingbats", "ZapfDingbats"'
         WRITE (iw, '(A)') '@map color 0 to (255, 255, 255), "white"', &
            '@map color 1 to (0, 0, 0), "black"', &
            '@map color 2 to (255, 0, 0), "red"', &
            '@map color 3 to (0, 255, 0), "green"', &
            '@map color 4 to (0, 0, 255), "blue"', &
            '@map color 5 to (255, 255, 0), "yellow"', &
            '@map color 6 to (188, 143, 143), "brown"', &
            '@map color 7 to (220, 220, 220), "grey"', &
            '@map color 8 to (148, 0, 211), "violet"', &
            '@map color 9 to (0, 255, 255), "cyan"', &
            '@map color 10 to (255, 0, 255), "magenta"', &
            '@map color 11 to (255, 165, 0), "orange"', &
            '@map color 12 to (114, 33, 188), "indigo"', &
            '@map color 13 to (103, 7, 72), "maroon"', &
            '@map color 14 to (64, 224, 208), "turquoise"', &
            '@map color 15 to (0, 139, 0), "green4"'
         WRITE (iw, '(A)') '@reference date 0', '@date wrap off', '@date wrap year 1950'
         WRITE (iw, '(A)') '@default linewidth 1.0', '@default linestyle 1', &
            '@default color 1', '@default pattern 1', '@default font 0', &
            '@default char size 1.000000', '@default symbol size 1.000000', &
            '@default sformat "%.8g"', '@background color 0', '@page background fill on'
         WRITE (iw, '(A)') '@timestamp off', '@timestamp 0.03, 0.03', '@timestamp color 1', &
            '@timestamp rot 0', '@timestamp font 0', '@timestamp char size 1.000000'
         CALL m_datum(date)

         WRITE (iw, '(A)') '@timestamp def '//TRIM(date)

      END IF

   END SUBROUTINE xm_write_defaults

! **************************************************************************************************
!> \brief ...
!> \param iw ...
! **************************************************************************************************
   SUBROUTINE xm_write_frameport(iw)
      INTEGER, INTENT(IN)                                :: iw

      IF (iw > 0) THEN
         WRITE (iw, '(A)') '@r0 off', '@link r0 to g0', '@r0 type above', '@r0 linestyle 1', &
            '@r0 linewidth 1.0', '@r0 color 1', '@r0 line 0, 0, 0, 0'
         WRITE (iw, '(A)') '@r1 off', '@link r1 to g0', '@r1 type above', '@r1 linestyle 1', &
            '@r1 linewidth 1.0', '@r1 color 1', '@r1 line 0, 0, 0, 0'
         WRITE (iw, '(A)') '@r2 off', '@link r2 to g0', '@r2 type above', '@r2 linestyle 1', &
            '@r2 linewidth 1.0', '@r2 color 1', '@r2 line 0, 0, 0, 0'
         WRITE (iw, '(A)') '@r3 off', '@link r3 to g0', '@r3 type above', '@r3 linestyle 1', &
            '@r3 linewidth 1.0', '@r3 color 1', '@r3 line 0, 0, 0, 0'

         WRITE (iw, '(A)') '@g0 on', '@g0 hidden false', '@g0 type XY', '@g0 stacked false', &
            '@g0 bar hgap 0.000000', '@g0 fixedpoint off', '@g0 fixedpoint type 0', &
            '@g0 fixedpoint xy 0.000000, 0.000000', '@g0 fixedpoint format general general', &
            '@g0 fixedpoint prec 6, 6'
      END IF
   END SUBROUTINE xm_write_frameport

! **************************************************************************************************
!> \brief ...
!> \param iw ...
!> \param wcoord ...
!> \param title ...
!> \param subtitle ...
!> \param xlabel ...
!> \param ylabel ...
! **************************************************************************************************
   SUBROUTINE xm_write_frame(iw, wcoord, title, subtitle, xlabel, ylabel)
      INTEGER, INTENT(IN)                                :: iw
      REAL(KIND=dp), DIMENSION(:)                        :: wcoord
      CHARACTER(len=*)                                   :: title, subtitle, xlabel, ylabel

      REAL(KIND=dp)                                      :: x1, x2, y1, y2

      x1 = wcoord(1)
      y1 = wcoord(2)
      x2 = wcoord(3)
      y2 = wcoord(4)
      IF (iw > 0) THEN
         WRITE (iw, '(A)') '@with g0'
         WRITE (iw, FMT='(A)', ADVANCE='NO') '@    world '
         WRITE (iw, FMT='(4(F8.1,A))') x1, ',', y1, ',', x2, ',', y2
         WRITE (iw, '(A)') '@    stack world 0, 0, 0, 0'
         WRITE (iw, '(A)') '@    znorm 1', &
            '@    view 0.150000, 0.150000, 1.150000, 0.850000'
         WRITE (iw, '(A)') '@    title "'//title//'"'
         WRITE (iw, '(A)') '@    title font 0', &
            '@    title size 1.500000', &
            '@    title color 1'
         WRITE (iw, '(A)') '@    subtitle "'//subtitle//'"'
         WRITE (iw, '(A)') '@    title font 0', &
            '@    title size 1.000000', &
            '@    title color 1'
         !
         WRITE (iw, '(A)') '@    xaxes scale Normal'
         WRITE (iw, '(A)') '@    yaxes scale Normal'
         WRITE (iw, '(A)') '@    xaxes invert off'
         WRITE (iw, '(A)') '@    yaxes invert off'
         ! xaxis
         WRITE (iw, '(A)') '@    xaxis  on', &
            '@    xaxis  type zero false', &
            '@    xaxis  offset 0.000000 , 0.000000', &
            '@    xaxis  bar on', &
            '@    xaxis  bar color 1', &
            '@    xaxis  bar linestyle 1', &
            '@    xaxis  bar linewidth 1.0'
         WRITE (iw, '(A)') '@    xaxis  label "'//xlabel//'"'
         WRITE (iw, '(A)') '@    xaxis  label layout para', &
            '@    xaxis  label place auto', &
            '@    xaxis  label char size 1.480000', &
            '@    xaxis  label font 0', &
            '@    xaxis  label color 1', &
            '@    xaxis  label place normal'
         WRITE (iw, '(A)') '@    xaxis  tick on', '@    xaxis  tick major 2', '@    xaxis  tick minor ticks 1', &
            '@    xaxis  tick default 6', '@    xaxis  tick place rounded true', '@    xaxis  tick in', &
            '@    xaxis  tick major size 1.000000', '@    xaxis  tick major color 1', &
            '@    xaxis  tick major linewidth 3.0', '@    xaxis  tick major linestyle 1', &
            '@    xaxis  tick major grid off', '@    xaxis  tick minor color 1', &
            '@    xaxis  tick minor linewidth 3.0', '@    xaxis  tick minor linestyle 1', &
            '@    xaxis  tick minor grid off', '@    xaxis  tick minor size 0.500000'
         WRITE (iw, '(A)') '@    xaxis  ticklabel on', '@    xaxis  ticklabel format general', &
            '@    xaxis  ticklabel prec 5', '@    xaxis  ticklabel formula ""', '@    xaxis  ticklabel append ""', &
            '@    xaxis  ticklabel prepend ""', '@    xaxis  ticklabel angle 0', '@    xaxis  ticklabel skip 0', &
            '@    xaxis  ticklabel stagger 0', '@    xaxis  ticklabel place normal', &
            '@    xaxis  ticklabel offset auto', &
            '@    xaxis  ticklabel offset 0.000000 , 0.010000', '@    xaxis  ticklabel start type auto', &
            '@    xaxis  ticklabel start 0.000000', '@    xaxis  ticklabel stop type auto', &
            '@    xaxis  ticklabel stop 0.000000', '@    xaxis  ticklabel char size 1.480000', &
            '@    xaxis  ticklabel font 0', '@    xaxis  ticklabel color 1', &
            '@    xaxis  tick place both', '@    xaxis  tick spec type none'
         ! yaxis
         WRITE (iw, '(A)') '@    yaxis  on', &
            '@    yaxis  type zero false', &
            '@    yaxis  offset 0.000000 , 0.000000', &
            '@    yaxis  bar on', &
            '@    yaxis  bar color 1', &
            '@    yaxis  bar linestyle 1', &
            '@    yaxis  bar linewidth 1.0'
         WRITE (iw, '(A)') '@    yaxis  label "'//ylabel//'"'
         WRITE (iw, '(A)') '@    yaxis  label layout para', &
            '@    yaxis  label place auto', &
            '@    yaxis  label char size 1.000000', &
            '@    yaxis  label font 0', &
            '@    yaxis  label color 1', &
            '@    yaxis  label place normal'
         WRITE (iw, '(A)') '@    yaxis  tick on', '@    yaxis  tick major 0.5', '@    yaxis  tick minor ticks 1', &
            '@    yaxis  tick default 6', '@    yaxis  tick place rounded true', '@    yaxis  tick in', &
            '@    yaxis  tick major size 1.480000', '@    yaxis  tick major color 1', &
            '@    yaxis  tick major linewidth 3.0', '@    yaxis  tick major linestyle 1', &
            '@    yaxis  tick major grid off', '@    yaxis  tick minor color 1', &
            '@    yaxis  tick minor linewidth 3.0', '@    yaxis  tick minor linestyle 1', &
            '@    yaxis  tick minor grid off', '@    yaxis  tick minor size 0.500000'
         WRITE (iw, '(A)') '@    yaxis  ticklabel on', '@    yaxis  ticklabel format general', &
            '@    yaxis  ticklabel prec 5', '@    yaxis  ticklabel formula ""', '@    yaxis  ticklabel append ""', &
            '@    yaxis  ticklabel prepend ""', '@    yaxis  ticklabel angle 0', '@    yaxis  ticklabel skip 0', &
            '@    yaxis  ticklabel stagger 0', '@    yaxis  ticklabel place normal', &
            '@    yaxis  ticklabel offset auto', &
            '@    yaxis  ticklabel offset 0.000000 , 0.010000', '@    yaxis  ticklabel start type auto', &
            '@    yaxis  ticklabel start 0.000000', '@    yaxis  ticklabel stop type auto', &
            '@    yaxis  ticklabel stop 0.000000', '@    yaxis  ticklabel char size 1.480000', &
            '@    yaxis  ticklabel font 0', '@    yaxis  ticklabel color 1', &
            '@    yaxis  tick place both', '@    yaxis  tick spec type none'
         WRITE (iw, '(A)') '@    altxaxis  off', '@    altyaxis  off'
         ! Legend
         WRITE (iw, '(A)') '@    legend on', &
            '@    legend loctype view', &
            '@    legend 0.8, 0.4', &
            '@    legend box color 1', &
            '@    legend box pattern 1', &
            '@    legend box linewidth 2.0', &
            '@    legend box linestyle 1', &
            '@    legend box fill color 0', &
            '@    legend box fill pattern 1', &
            '@    legend font 0', &
            '@    legend char size 1.000000', &
            '@    legend color 1', &
            '@    legend length 4', &
            '@    legend vgap 1', &
            '@    legend hgap 1', &
            '@    legend invert false'
         ! Frame
         WRITE (iw, '(A)') '@    frame type 0', &
            '@    frame linestyle 1', &
            '@    frame linewidth 3.0', &
            '@    frame color 1', &
            '@    frame pattern 1', &
            '@    frame background color 0', &
            '@    frame background pattern 0'
      END IF
   END SUBROUTINE xm_write_frame

! **************************************************************************************************
!> \brief ...
!> \param iw ...
!> \param gnum ...
!> \param linewidth ...
!> \param legend ...
! **************************************************************************************************
   SUBROUTINE xm_graph_info(iw, gnum, linewidth, legend)
      INTEGER, INTENT(IN)                                :: iw, gnum
      REAL(KIND=dp), INTENT(IN)                          :: linewidth
      CHARACTER(LEN=*)                                   :: legend

      CHARACTER(LEN=8)                                   :: cin, cnum, cval

      IF (iw > 0) THEN
         WRITE (cnum, '(I2)') gnum
         WRITE (cval, '(F3.1)') linewidth
         cin = "@    s"//TRIM(ADJUSTL(cnum))
         WRITE (cnum, '(I2)') gnum + 1
         WRITE (iw, '(A)') TRIM(cin)//' hidden false'
         WRITE (iw, '(A)') TRIM(cin)//' type xy'
         WRITE (iw, '(A)') TRIM(cin)//' symbol 0 '
         WRITE (iw, '(A)') TRIM(cin)//' symbol size 1.000000'
         WRITE (iw, '(A)') TRIM(cin)//' symbol color '//TRIM(ADJUSTL(cnum))
         WRITE (iw, '(A)') TRIM(cin)//' symbol pattern 1'
         WRITE (iw, '(A)') TRIM(cin)//' symbol fill color 1'
         WRITE (iw, '(A)') TRIM(cin)//' symbol fill pattern 0'
         WRITE (iw, '(A)') TRIM(cin)//' symbol linewidth 1.0'
         WRITE (iw, '(A)') TRIM(cin)//' symbol linestyle 1'
         WRITE (iw, '(A)') TRIM(cin)//' symbol char 65 '
         WRITE (iw, '(A)') TRIM(cin)//' symbol char font 0'
         WRITE (iw, '(A)') TRIM(cin)//' symbol skip 0'
         WRITE (iw, '(A)') TRIM(cin)//' line type 1'
         WRITE (iw, '(A)') TRIM(cin)//' line linestyle 1'
         WRITE (iw, '(A)') TRIM(cin)//' line linewidth '//TRIM(cval)
         WRITE (iw, '(A)') TRIM(cin)//' line color '//TRIM(ADJUSTL(cnum))
         WRITE (iw, '(A)') TRIM(cin)//' line pattern 1'
         WRITE (iw, '(A)') TRIM(cin)//' baseline type 0'
         WRITE (iw, '(A)') TRIM(cin)//' baseline off'
         WRITE (iw, '(A)') TRIM(cin)//' dropline off'
         WRITE (iw, '(A)') TRIM(cin)//' fill type 0'
         WRITE (iw, '(A)') TRIM(cin)//' fill rule 0'
         WRITE (iw, '(A)') TRIM(cin)//' fill color '//TRIM(ADJUSTL(cnum))
         WRITE (iw, '(A)') TRIM(cin)//' fill pattern 1'
         WRITE (iw, '(A)') TRIM(cin)//' avalue off'
         WRITE (iw, '(A)') TRIM(cin)//' avalue type 2'
         WRITE (iw, '(A)') TRIM(cin)//' avalue char size 1.000000'
         WRITE (iw, '(A)') TRIM(cin)//' avalue font 0'
         WRITE (iw, '(A)') TRIM(cin)//' avalue color '//TRIM(ADJUSTL(cnum))
         WRITE (iw, '(A)') TRIM(cin)//' avalue rot 0'
         WRITE (iw, '(A)') TRIM(cin)//' avalue format general'
         WRITE (iw, '(A)') TRIM(cin)//' avalue prec 3'
         WRITE (iw, '(A)') TRIM(cin)//' avalue prepend ""'
         WRITE (iw, '(A)') TRIM(cin)//' avalue append ""'
         WRITE (iw, '(A)') TRIM(cin)//' avalue offset 0.000000 , 0.000000'
         WRITE (iw, '(A)') TRIM(cin)//' errorbar on'
         WRITE (iw, '(A)') TRIM(cin)//' errorbar place both'
         WRITE (iw, '(A)') TRIM(cin)//' errorbar color '//TRIM(ADJUSTL(cnum))
         WRITE (iw, '(A)') TRIM(cin)//' errorbar pattern 1'
         WRITE (iw, '(A)') TRIM(cin)//' errorbar size 1.000000'
         WRITE (iw, '(A)') TRIM(cin)//' errorbar linewidth 1.0'
         WRITE (iw, '(A)') TRIM(cin)//' errorbar linestyle 1'
         WRITE (iw, '(A)') TRIM(cin)//' errorbar riser linewidth 1.0'
         WRITE (iw, '(A)') TRIM(cin)//' errorbar riser linestyle 1'
         WRITE (iw, '(A)') TRIM(cin)//' errorbar riser clip off'
         WRITE (iw, '(A)') TRIM(cin)//' errorbar riser clip length 0.100000'
         WRITE (iw, '(A)') TRIM(cin)//' comment "Cols 1:2"'
         WRITE (iw, '(A)') TRIM(cin)//' legend  "'//TRIM(legend)//'"'
      END IF
   END SUBROUTINE xm_graph_info

! **************************************************************************************************
!> \brief ...
!> \param iw ...
!> \param gnum ...
!> \param gdata ...
! **************************************************************************************************
   SUBROUTINE xm_graph_data(iw, gnum, gdata)
      INTEGER, INTENT(IN)                                :: iw, gnum
      REAL(KIND=dp), DIMENSION(:, :)                     :: gdata

      CHARACTER(LEN=8)                                   :: cin, cnum
      INTEGER                                            :: i, m

      IF (iw > 0) THEN
         WRITE (cnum, '(I2)') gnum
         cin = "@@target G0.S"//TRIM(ADJUSTL(cnum))
         WRITE (iw, '(A)') TRIM(cin)
         WRITE (iw, '(A)') '@type xy'
         m = SIZE(gdata, 1)
         DO i = 1, m
            WRITE (iw, '(2G18.7)') gdata(i, 1), gdata(i, 2)
         END DO
         WRITE (iw, '(A)') '&'
      END IF
   END SUBROUTINE xm_graph_data

END MODULE xmgrace

