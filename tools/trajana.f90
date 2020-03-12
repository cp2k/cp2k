  PROGRAM trajana
! --------------------------------------------------------------------
! reads an TRAJECTORY or xyz file
! give the filename on the command line or it defaults to TRAJECTORY
!   Features:
!    -skips '<<NEW DATA' lines 
!    -skips blank lines 
!    -finds the number of atoms 
!   LIMITATIONS
!    -assumes that the number of atoms is the same throughout the  
!       trajectory
! --------------------------------------------------------------------

  IMPLICIT NONE


  CHARACTER (80) :: line, output_filename, in_filename, out_filename, &
                          todo_filename
  INTEGER :: read_status, natoms, argcount
  INTEGER :: iu_infile = 25
  INTEGER :: iu_outfile = 35
  INTEGER :: iu_todofile = 45
  INTEGER :: nframes 
  DOUBLE PRECISION, PARAMETER :: au2ang = 0.5291772
  INTEGER, DIMENSION (1000,5) :: actiontype
  INTEGER :: nactions
  INTEGER :: iostatus
  INTEGER :: style
  INTEGER :: IARGC


  argcount = COMMAND_ARGUMENT_COUNT()
  if (argcount==0) then
   in_filename   = "TRAJECTORY"
  elseif (argcount==1) then
   CALL GET_COMMAND_ARGUMENT(1,in_filename)
  else
   STOP "TOO MANY ARGUMENTS"
  endif


  out_filename  = "trajana.out" 
  todo_filename = "trajana.in"
 
  iostatus = 0 
  open (iu_infile,FILE=in_filename,IOSTAT=iostatus,STATUS='OLD')
  if (iostatus.NE.0) CALL give_info("unable to open TRAJECTORY","")
  open (iu_todofile,FILE=todo_filename,IOSTAT=iostatus,STATUS='OLD')
  if (iostatus.NE.0) CALL give_info("unable to open trajana.in","")
  open (iu_outfile,FILE=out_filename,IOSTAT=iostatus,STATUS='UNKNOWN')
  if (iostatus.NE.0) CALL give_info("unable to open trajana.out","")

  write(*,*) "Opened file ",in_filename
  call determine_style(iu_infile, style)


 
  call determine_number_of_atoms (iu_infile, natoms, style)
  write(*,*) "atoms per frame : ", natoms
  if (natoms.LE.0) CALL give_info("number of atoms should be > 0","")

  


  call determine_number_of_frames (iu_infile, natoms, nframes,style)
  write(*,*) "number of frames : ",nframes
  if (nframes.LE.0) CALL give_info("number of frames should be > 0","")
  call read_todo(iu_todofile, actiontype, nactions,natoms)
  write(*,*) "number of quantities to be calculated : ",nactions

  call process_data(iu_infile,iu_outfile,natoms,nframes, actiontype, nactions,style)

  close (iu_infile, status="KEEP")
  close (iu_outfile, status="KEEP")
  close (iu_todofile, status="KEEP")
  write(*,'(/)')

CONTAINS

!==============================================================================
!==============================================================================
  SUBROUTINE give_info(S,S2)
	character(*) , INTENT(IN) :: S,S2
	write(*,*) "**** TRAJectory ANAlyzer ****" 
	write(*,*) "Following error occurred : "
	write(*,*) ""
	write(*,*) S 
	write(*,*) S2
	write(*,*) ""
	write(*,*) "**** TRAJectory ANAlyzer ****"
	write(*,*) "****        Manual       ****"
	write(*,*) "*****************************"
	write(*,*) "It can calculate - Distances      in angstrom"
	write(*,*) "                 - Angles         in degrees"
	write(*,*) "                 - Torsion angles in degrees"
	write(*,*) "This program will read from the file TRAJECTORY"
        write(*,*) "or from the file given on the command line"
        write(*,*) "This file can either be in xmol or TRAJECTORY format"
        write(*,*) "                                file trajana.in"
	write(*,*) "             and write to the   file trajana.out"
	write(*,*) ""
	write(*,*) ""
	write(*,*) "the structure of trajana.in is as follows"
	write(*,*) "----------------------------------"
        write(*,*) "D[istance] #1 #1"
	write(*,*) "A[ngle]    #1 #2 #3"
	write(*,*) "T[orsion]  #1 #2 #3 #4"
	write(*,*) "! COMMENT"
	write(*,*) "# COMMENT"
	write(*,*) "C COMMENT"
	write(*,*) "----------------------------------"
	write(*,*) "the program only checks for the first character of the line"
	write(*,*) "so abbrev. may be used"
	write(*,*) "the order of the lines is not important."
	write(*,*) "The number of lines is limited to 1000"
	write(*,*) "#1 #2 #3 #4, are the numbers of the atoms as they appear in"
	write(*,*) " the trajectory / movie / output file with 1 being the first"
	write(*,*) ""
	write(*,*) ""
	write(*,*) "the structure of trajana.out is as follows"
	write(*,*) "----------------------------------"
	write(*,*) "#1 #2 #3 #4 #5 ....."
	write(*,*) " ...."
	write(*,*) "----------------------------------"
	write(*,*) "#1 is the number of the frame analyzed"
	write(*,*) "   #2 #3 ... is the distance,angle,torsion, in the same order "
	write(*,*) "             as it appears in the trajana.in file"
        write(*,*) "the following line will be for the following frame" 
	STOP
  END SUBROUTINE give_info
!==============================================================================
!==============================================================================

  SUBROUTINE determine_style(iu, style)
    INTEGER, INTENT(IN) :: iu
    INTEGER, INTENT(OUT):: style
    CHARACTER (120):: line
    INTEGER:: ls, read_status
    double precision :: rtemp

    rewind (iu)
    call read_line(iu, read_status, line)
    if (read_status .ne. 0) then
       CALL give_info("error reading file","empty")
    endif
    ls = 1
    call csgetr (line, ls, rtemp)
    call csgetr (line, ls, rtemp)
    
    if (ls == 0)then
      style = 1
      write (*,*) "Taking the XMOL format "
    else 
      style = 2
      write (*,*) "Taking the TRAJECTORY format"
    endif
  END SUBROUTINE determine_style

!==============================================================================
!==============================================================================
   


  SUBROUTINE determine_number_of_frames (iu, natoms, nframes, style)

    INTEGER, INTENT(IN) :: iu
    INTEGER, INTENT(IN) :: natoms, style
    INTEGER, INTENT(OUT) :: nframes 

    INTEGER :: read_status, itemp, iatom_type
    INTEGER :: iframes, i
    LOGICAL :: continue_read

    INTEGER :: iatom, ls
    CHARACTER (120) :: line

    rewind (iu)
    iframes = 0 
    continue_read = .true.

    do 
      iframes = iframes + 1
      if (style==1)then   !xmol - file
          do i=1,2  ! number of atoms and comment
          call read_line (iu, read_status, line)
          if (read_status .ne. 0) then
               continue_read = .false.
               exit
          end if
          end do
      endif
      do iatom = 1, natoms
        call read_line (iu, read_status, line)
        if (read_status .ne. 0) then
          continue_read = .false.
          exit
        end if
      end do

      if (.not.continue_read) then
          exit
      endif

    end do

    nframes = iframes - 1
 
  END SUBROUTINE determine_number_of_frames 

!==============================================================================
!==============================================================================

  SUBROUTINE determine_number_of_atoms(iu, natoms, style)

    INTEGER, INTENT(IN) :: iu, style
    INTEGER, INTENT(OUT) :: natoms

    INTEGER :: read_status, itemp, frame
    LOGICAL :: continue_read

    INTEGER :: iatom, ls
    CHARACTER (120) :: line

    rewind (iu)
    natoms = 0
    continue_read = .true.
    frame = -1
    
    if (style==1)then       !xmol - file
        call read_line (iu, read_status, line)
        ls = 1
        CALL csgeti(line,ls,itemp)
        if (ls.EQ.0) then
            CALL give_info("error reading xmol",line)
            endif
        natoms=itemp
    else                     !trajectory - file
    do
        call read_line_2 (iu, read_status, line)
        if (read_status .ne. 0) then
          exit
        end if
        ls = 1
        CALL csgeti(line,ls,itemp)
        if (ls.eq.0) then
	 exit	
	endif
        if (frame.eq.-1) then
	  natoms=natoms+1
          frame=itemp
        else
          if (itemp.EQ.frame) then
	     natoms=natoms+1
          else
	     exit
          endif
	endif
    end do
    endif

 
  END SUBROUTINE determine_number_of_atoms

!==============================================================================

!==============================================================================

  SUBROUTINE read_todo (iu,actiontype,nactions,natoms)

    INTEGER, INTENT(IN) :: iu, natoms
    INTEGEr, INTENT(OUT) :: nactions
    INTEGER, DIMENSION(1000,5), INTENT(OUT) :: actiontype  

    INTEGER :: read_status, itemp, iatom_type
    LOGICAL :: continue_read

    INTEGER :: iatom, ls

    CHARACTER (120) :: line
    CHARACTER :: firstchar

    rewind (iu)
    continue_read = .true.
    nactions = 0

    do
	
        call read_line_2 (iu, read_status, line)
        if (read_status .ne. 0) then
		continue_read=.false.
		exit
        end if
        if (read_status .eq. 0) then
	   read (line,'(A1)') firstchar 
	   select case (firstchar)
		case('D','d')
		   nactions = nactions + 1
                   actiontype(nactions,1) = 1
		   ls = 1
                   CALL csgeti(line,ls,itemp)
		   actiontype(nactions,2)=itemp
                   if ((ls.EQ.0).OR.(itemp.LE.0).OR.(itemp.GT.natoms)) CALL give_info("error in trajana.in",line)
                   CALL csgeti(line,ls,itemp)
		   actiontype(nactions,3)=itemp 
                   if ((ls.EQ.0).OR.(itemp.LE.0).OR.(itemp.GT.natoms)) CALL give_info("error in trajana.in",line)
		case('A','a')
                   nactions = nactions + 1
                   actiontype(nactions,1) = 2
                   ls = 1
                   CALL csgeti(line,ls,itemp)
                   actiontype(nactions,2)=itemp
                   if ((ls.EQ.0).OR.(itemp.LE.0).OR.(itemp.GT.natoms)) CALL give_info("error in trajana.in",line)
                   CALL csgeti(line,ls,itemp)
                   actiontype(nactions,3)=itemp
                   if ((ls.EQ.0).OR.(itemp.LE.0).OR.(itemp.GT.natoms)) CALL give_info("error in trajana.in",line)
                   CALL csgeti(line,ls,itemp)
                   actiontype(nactions,4)=itemp
                   if ((ls.EQ.0).OR.(itemp.LE.0).OR.(itemp.GT.natoms)) CALL give_info("error in trajana.in",line)

		
		case('T','t')
                   nactions = nactions + 1
                   actiontype(nactions,1) = 3
                   ls = 1
                   CALL csgeti(line,ls,itemp)
                   actiontype(nactions,2)=itemp
                   if ((ls.EQ.0).OR.(itemp.LE.0).OR.(itemp.GT.natoms)) CALL give_info("error in trajana.in",line)
                   CALL csgeti(line,ls,itemp)
                   actiontype(nactions,3)=itemp
                   if ((ls.EQ.0).OR.(itemp.LE.0).OR.(itemp.GT.natoms)) CALL give_info("error in trajana.in",line)
                   CALL csgeti(line,ls,itemp)
                   actiontype(nactions,4)=itemp
                   if ((ls.EQ.0).OR.(itemp.LE.0).OR.(itemp.GT.natoms)) CALL give_info("error in trajana.in",line)
                   CALL csgeti(line,ls,itemp)
                   actiontype(nactions,5)=itemp
                   if ((ls.EQ.0).OR.(itemp.LE.0).OR.(itemp.GT.natoms)) CALL give_info("error in trajana.in",line)
		
                case('!','C','#',' ')
		   ! this is just a comment line 
                case default
		   !something wrong in syntax
                   if ((ls.EQ.0).OR.(itemp.LE.0).OR.(itemp.GT.natoms)) CALL give_info("error in trajana.in",line)
	   end select
        end if
    end do
    if (read_status.EQ.2) CALL give_info("Error reading trajana.in","")

  END SUBROUTINE read_todo

!==============================================================================

  SUBROUTINE process_data(iu_infile,iu_outfile,natoms,nframes,actiontype,nactions,style)

    INTEGER, INTENT(IN) :: iu_outfile
    INTEGER, INTENT(IN) :: iu_infile

    INTEGER, INTENT(IN) :: nframes
    INTEGER, INTENT(IN) :: natoms, style
    INTEGER, INTENT(IN)  :: nactions
    INTEGER, DIMENSION(1000,5), INTENT(IN) :: actiontype
    DOUBLE PRECISION, DIMENSION (1000) :: results
    DOUBLE PRECISION, DIMENSION (3,natoms) :: xyz
    DOUBLE PRECISION, PARAMETER :: au2ang = 0.5291772
    DOUBLE PRECISION, PARAMETER :: rad2deg = 57.29578
    DOUBLE PRECISION, DIMENSION (3,8) :: dummyvec
    INTEGER :: iframe, iaction, I, K, iatom, J
    DOUBLE PRECISION :: distance, scal, norm, scal2
    CHARACTER (80) :: LINE_OUT
    INTEGER :: read_status, itemp,  iatom_type, ls
    CHARACTER(120) :: line
    CHARACTER(120) :: element_type


    WRITE (LINE_OUT,'(A,I5,A)') '(I10.7,',nactions,'E15.7)'
    rewind (iu_infile)

    DO iframe=1, nframes ! TO GET THE data 

      if (style==1)then
          call read_line (iu_infile, read_status, line)
          call read_line (iu_infile, read_status, line)
      endif
     
      do iatom = 1, natoms
        call read_line (iu_infile, read_status, line)
        ls = 1
        if (style==1)then
           call csgets(line,ls, element_type)
        else
           call csgeti (line, ls, itemp)
        endif
        call csgetr (line, ls, xyz(1,iatom)) 
        call csgetr (line, ls, xyz(2,iatom)) 
        call csgetr (line, ls, xyz(3,iatom)) 
        if (ls.EQ.0) CALL give_info("error in trajana.in","")
!       call csgetr (line, ls, dtemp) ! to get the velocities
!       call csgetr (line, ls, dtemp) 
!       call csgetr (line, ls, dtemp)
      end do
      
      DO iaction=1, nactions
	select case (actiontype(iaction,1))
    	 case (1)
		CALL addvec( 1.0D0, xyz(1,actiontype(iaction,2)),  &
			    -1.0D0, xyz(1,actiontype(iaction,3)),  &
			     dummyvec(1,1)   )
		CALL getnorm(dummyvec(1,1), distance)
                if (style==1)then
                results(iaction) = distance
                else
		results(iaction) = distance * au2ang
                endif
	 case (2)
		CALL addvec( 1.0D0, xyz(1,actiontype(iaction,2)),  &
			    -1.0D0, xyz(1,actiontype(iaction,3)),  &
			     dummyvec(1,1)   )
		CALL normvec(dummyvec(1,1), dummyvec(1,1))
		CALL addvec( 1.0D0, xyz(1,actiontype(iaction,4)),  &
	  		    -1.0D0, xyz(1,actiontype(iaction,3)),  &
			     dummyvec(1,2)   )
		CALL normvec(dummyvec(1,2), dummyvec(1,2))
		CALL getscal(dummyvec(1,1), dummyvec(1,2), scal)
                results(iaction) = acos ( scal )*rad2deg
	 case (3)
		
                CALL addvec( 1.0D0, xyz(1,actiontype(iaction,2)),  &
                            -1.0D0, xyz(1,actiontype(iaction,3)),  &
                             dummyvec(1,1)   )
                CALL normvec(dummyvec(1,1), dummyvec(1,1))
                CALL addvec( 1.0D0, xyz(1,actiontype(iaction,4)),  &
                            -1.0D0, xyz(1,actiontype(iaction,3)),  &
                             dummyvec(1,2)   )
                CALL normvec(dummyvec(1,2), dummyvec(1,2))
		CALL vecprod(dummyvec(1,1),dummyvec(1,2),dummyvec(1,5)) 
		CALL normvec(dummyvec(1,5),dummyvec(1,5))

                CALL addvec( 1.0D0, xyz(1,actiontype(iaction,3)),  &
                            -1.0D0, xyz(1,actiontype(iaction,4)),  &
                             dummyvec(1,3)   )
                CALL normvec(dummyvec(1,3), dummyvec(1,3))
                CALL addvec( 1.0D0, xyz(1,actiontype(iaction,5)),  &
                            -1.0D0, xyz(1,actiontype(iaction,4)),  &
                             dummyvec(1,4)   )
                CALL normvec(dummyvec(1,4), dummyvec(1,4))
		CALL vecprod(dummyvec(1,3),dummyvec(1,4),dummyvec(1,6)) 
		CALL normvec(dummyvec(1,6),dummyvec(1,6))

		CALL addvec( 1.0D0, xyz(1,actiontype(iaction,4)),  &
			    -1.0D0, xyz(1,actiontype(iaction,3)),  &
			     dummyvec(1,8)   )

		CALL getscal(dummyvec(1,5),dummyvec(1,6),scal)
		CALL vecprod(dummyvec(1,5),dummyvec(1,6),dummyvec(1,7))
		CALL getnorm(dummyvec(1,7),norm)

		CALL getscal(dummyvec(1,7),dummyvec(1,8),scal2)
		
		results(iaction) = ACOS(scal)
		IF (scal2.LT.0) THEN
			results(iaction) = -results(iaction) +&
						  2*3.1415926535
		END IF
		results(iaction)=results(iaction)*rad2deg
		
	end select
      END DO

      WRITE (iu_outfile,LINE_OUT) iframe,(results(J),J=1,nactions)

    END DO

  END SUBROUTINE process_data

!==============================================================================



!==============================================================================

  SUBROUTINE read_line ( iu, status_out, line)

!   ---------------------------------------------------------
!   reads the from input and skips selected kinds of lines
!
!   status_out = 0  everything normal
!                1  end of file reached
!                2  some other kind of error 
!   ---------------------------------------------------------

    INTEGER, INTENT(IN) :: iu
    INTEGER, INTENT(OUT) :: status_out
    CHARACTER(*), INTENT(IN OUT) :: line
    INTEGER :: ios


    status_out = 0
    do
      read (iu, '(a120)', IOSTAT=ios, END=999) line       
      if (ios .ne. 0) then
        status_out = 2
        exit
      end if
      if (line(4:5) == '<<') cycle
      if (line .ne. "") exit
    end do
   
    return

!   ---------------------------------------
!   If end of file is reached reading input
!   ---------------------------------------
    999  continue
 
   status_out = 1
   
  END SUBROUTINE read_line 

!==============================================================================

!==============================================================================

  SUBROUTINE read_line_2 ( iu, status_out, line)

!   ---------------------------------------------------------
!   reads the from input and skips selected kinds of lines
!
!   status_out = 0  everything normal
!                1  end of file reached
!                2  some other kind of error
!   ---------------------------------------------------------

    INTEGER, INTENT(IN) :: iu
    INTEGER, INTENT(OUT) :: status_out
    CHARACTER(*), INTENT(IN OUT) :: line
    INTEGER :: ios


    status_out = 0
    do
      read (iu, '(a120)', IOSTAT=ios, END=9999) line
      if (ios .ne. 0) then
        status_out = 2
        exit
      end if
      if (line .ne. "") exit
    end do
 
    return

!   ---------------------------------------
!   If end of file is reached reading input
!   ---------------------------------------
    9999  continue
 
   status_out = 1
 
  END SUBROUTINE read_line_2

!==============================================================================

END PROGRAM trajana



      subroutine csbegn (string, is)
!
! ======================================================================
!
! purpose: determine the position of the first non-blank character
!          in a string.
!
! input  : string - character string to be analyzed
!
! output : is     - position of first non-blank (0 if not found)
!
! *=====================================================================
!
      implicit   integer (i-n)
      implicit   double precision (a-h,o-z)
!
      parameter (lchars = 80)
      character*(lchars) line, word
      logical    valid
! ----------------------------------------------------------------------
      character*(*) string
!
      parameter (nchop = 8)
      character*(nchop) space8
      parameter (space8 = '        ')
!
! ======================================================================
!
      l2 = len (string)
      l1 = 0
   10 if (l2-l1 .ge.nchop) then
        if (string(l1+1:l1+nchop) .ne.space8) goto 20
        l1 = l1+nchop
        goto 10
      endif
!
   20 l1 = l1+1
      if (l1.gt.l2) then
        l1 = 0
      elseif (string(l1:l1) .eq.' ') then
        goto 20
      endif
!
      is = l1
!
      return
      end
      subroutine csend (string, ls)
!
! ======================================================================
!
! purpose: determine the end of a character string, not counting the
!          end-blanks.
!
! input  : string - string to be analyzed
!
! output : ls     - position of the last non-blank character
!
! *=====================================================================
!
      implicit   integer (i-n)
      implicit   double precision (a-h,o-z)
!
      parameter (lchars = 80)
      character*(lchars) line, word
      logical    valid
! ----------------------------------------------------------------------
      character*(*) string
!
      parameter (nchop = 8)
!
! ======================================================================
!
      ll = len (string)
      if (ll.le.0) then 
        stop
      end if
!
      l2 = 0
   10 if (string(l2+1:) .eq.' ') goto 20
      if (l2+nchop .lt.ll) then
        l2 = l2+nchop
        goto 10
      endif
      l2 = ll
!
   20 if (l2.eq.0) goto 50
      if (string(l2:l2) .ne.' ') goto 50
      l2 = l2-1
      goto 20
!
   50 ls = l2
!
      return
      end
      logical function cseq (strng1, strng2)
!
! ======================================================================
!
! purpose: test whether two strings are equal, neglecting end-blanks
!          and initial blanks. neglect lower/uppercase.
!
! input  : strng1, strngs - two strings
!
! *=====================================================================
!
      implicit   integer (i-n)
      implicit   double precision (a-h,o-z)
!
      parameter (lchars = 80)
      character*(lchars) line, word
      logical    valid
! ----------------------------------------------------------------------
!
      character*(*)  strng1, strng2
      character*(26) lowerc, upperc
      character*(1)  c1,c2
!
      data lowerc / 'abcdefghijklmnopqrstuvwxyz' /,    &
&          upperc / 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' /
!
! ======================================================================
!
      call csend (strng1, ls1)
      call csend (strng2, ls2)
!
      if (ls1.gt.0 .and. ls2.gt.0) then
!
        call csbegn (strng1, is1)
        call csbegn (strng2, is2)
!
        if (ls1-is1 .ne. ls2-is2) then
          cseq = .false.
        elseif (strng1(is1:ls1) .eq. strng2(is2:ls2)) then
          cseq = .true.
!
        else
!
! upper/lowercase
!
          cseq = .false.
!
          k = -1
   10     k = k+1
!
            k1 = is1+k
            c1 = strng1(k1:k1)
            j1 = index (lowerc, c1) + index (upperc, c1)
!
            k2 = is2+k
            c2 = strng2(k2:k2)
            j2 = index (lowerc, c2) + index (upperc, c2)
!
            if (j1.ne.j2) then
              return
            elseif (j1.eq.0) then
              if (c1.ne.c2) return
            endif
!
          if (k1.lt.ls1) goto 10
!
          cseq = .true.
        endif
!
      else
        cseq = (ls1.eq.0) .and. (ls2.eq.0)
      endif
!
      return
      end
      subroutine csgeti (string, is, ival)
!
! ======================================================================
!
! purpose: get the next integer from a string (skip any preceding
!          non-numeric data)
!
! input  : string
!
! in-out : is - in : position to start searching
!                    n.b.: if non-positive, start at 1
!        :      out: (only when positive on input)
!                    0 : not found
!                    other: position after item
!
! output : ival - integer (if found)
!
! *=====================================================================
!
      implicit   integer (i-n)
      implicit   double precision (a-h,o-z)
!
      parameter (lchars = 80)
      character*(lchars) line, word
      logical    valid
! ----------------------------------------------------------------------
      character*(*) digits, string
!
      parameter (digits = '0123456789')
!
! ======================================================================
!
      istart = max (1,is)
      lx = len (string)
      ip = 0
!
! locate the first digit
!
      do 50 i0 = istart, lx
        if (index (digits,string(i0:i0)) .gt.0) then
!
! scan the integer, by transforming all successive digits from character
! into numerical value, and updating the total integer value
!
          m0 = ichar('0')
          ival = 0
!
          ip = i0
   10     ival = 10*ival + ichar (string(ip:ip)) - m0
          ip = ip+1
!
          if (ip.le.lx) then
            if (index (digits,string(ip:ip)) .gt.0) goto 10
          endif
!
! check for a preceding minus sign (if the integer did not start at the
! first character of the input string)
!
          if (i0.gt.istart) then
            if (string(i0-1:i0-1) .eq.'-') ival = -ival
          endif
!
          goto 60
!
        endif
   50 continue
!
   60 if (is.gt.0) is = ip
!
      return
      end
      subroutine csgetr (string, is, rval)
!
! ======================================================================
!
! purpose: get the next real value from a string (skipping preceding
!          non-numeric data)
!
! input  : string
!
! in-out : is - in : position to start searching
!                    n.b.: if non-positive, start at 1
!        :      out: (only when positive on input)
!                    0 : not found
!                    other: position after item
!
! output : rval - real value (if found)
!
! remark * a real number r is written as : r = s i1 . i2 e
!            (no blanks between the parts)
!          s is a sign (optional)
!          i1 and i2 are unsigned integers. one of them may be absent
!          the decimal point is optional if i2 is absent
!          e is the exponential part (optional), written as : e = p s i
!            (no blanks between the parts)
!          p is the power indicator (d or e)
!          s is a sign (optional)
!          i is an unsigned integer
!
! *=====================================================================
!
      implicit   integer (i-n)
      implicit   double precision (a-h,o-z)
!
      parameter (lchars = 80)
      character*(lchars) line, word
      logical    valid
! ----------------------------------------------------------------------
!
      character*(*) digits, string
      character*(1) ch
      logical     cseq
!
      parameter (digits = '0123456789')
!
! ======================================================================
!
! get the starting digit of the first integer part (i1 or i2)
!
      istart = max (1,is)
      lx = len(string)
!
      do 400 i0 = istart, lx
        ch = string(i0:i0)
        if (index (digits,ch) .ne.0) goto 410
  400 continue
      goto 1000
!
! end of this integer part
!
  410 do 420 j = i0+1, lx
        ch = string(j:j)
        if (index (digits,ch) .eq.0) goto 450
  420 continue
      ix = lx+1
      goto 600
!
! integer preceded by a decimal point. if yes then proceed to
! the exponential part (i1 is absent)
!
  450 if (i0.gt.istart) then
        if (string(i0-1:i0-1) .eq.'.') then
          i0 = i0-1
          goto 550
        endif
      endif
!
! integer not followed by decimal point: goto expo part
!
      if (ch.ne.'.') goto 550
!
! end of decimal fraction
!
      k = j+1
      do 500 j = k, lx
        ch = string(j:j)
        if (index (digits,ch) .eq.0) goto 550
  500 continue
      ix = lx+1
      goto 600
!
! check presence of i2 (fractional part)
!
  550 ix = j
      if (cseq (ch,'E') .or. cseq (ch,'D')) then
!
        j = j+1
        if (j.gt.lx) goto 600
        if (string(j:j).eq.'+' .or. string(j:j).eq.'-') then
          j = j+1
          if (j.gt.lx) goto 600
        endif
        if (index (digits,string(j:j)) .eq.0) goto 600
!
        do 560 ix = j+1, lx
          ch = string(ix:ix)
          if (index (digits,ch) .eq.0) goto 600
  560   continue
        ix = lx+1
!
      endif
!
! complete number may be preceded by a sign
!
  600 if (i0.gt.istart) then
        if (string(i0-1:i0-1).eq.'+'.or.string(i0-1:i0-1).eq.'-') i0 = i0-1
      endif
!
      if (ix-i0 .gt. 80) stop
      line = string(i0:ix-1)
      read (line,'(BN,F80.0)') rval
      ip = ix
!
      goto 1050
!
! ----------------------------------------------------------------------
 1000 ip = 0
 1050 if (is.gt.0) is = ip
!
      return
      end
      subroutine csgets (string, is, text)
!
! ======================================================================
!
! purpose: get the next non-empty substring from a string (i.e. the part
!          separated from the remaining part by any from a set of
!          delimiters: comma, blank, ...)
!
! input  : string
!
! in-out : is - in : position to start searching
!                    n.b.: if non-positive, start at 1
!        :      out: (only when positive on input)
!                    0 : not found
!                    other: position after item
!
! output : text - substring (if found)
!
! *=====================================================================
!
      implicit   integer (i-n)
      implicit   double precision (a-h,o-z)
!
      parameter (lchars = 80)
      character*(lchars) line, word
      logical    valid
! ----------------------------------------------------------------------
      character*(*) delim, string, text
!
      parameter (delim = ' ,=')
!
! ======================================================================
!
      istart = max (1,is)
      lx = len (string)
!
      ip = 0
      do 50 i = istart, lx
        if (index (delim, string(i:i)) .le.0) then
!
          ip = i
   10     ip = ip+1
          if (ip.le.lx) then
            if (index (delim, string(ip:ip)) .le.0) goto 10
          endif
!
          text = string(i:ip-1)
          goto 60
        endif
   50 continue
!
   60 if (is.gt.0) is = ip
!
      return
      end

!     ==================================================================
      SUBROUTINE getnorm(vec,dnorm)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION vec(3)
!     ==--------------------------------------------------------------==
      dnorm = sqrt ( vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3) )
!     ==--------------------------------------------------------------==
      RETURN
      END
!     ==================================================================
!     ==================================================================
      SUBROUTINE getscal(vec1,vec2,scal)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION vec1(3), vec2(3)
!     ==--------------------------------------------------------------==
      scal = vec1(1)*vec2(1) + vec1(2)*vec2(2) + vec1(3)*vec2(3)
!     ==--------------------------------------------------------------==
      RETURN
      END
!     ==================================================================

!     ==================================================================
      SUBROUTINE normvec(vec_old,vec_norm)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION vec_old(3), vec_norm(3)
!     ==--------------------------------------------------------------==
      CALL getnorm(vec_old, dnorm)
      vec_norm(1) = vec_old(1) / dnorm
      vec_norm(2) = vec_old(2) / dnorm       
      vec_norm(3) = vec_old(3) / dnorm
!     ==--------------------------------------------------------------==
      RETURN
      END
!     ==================================================================
!     ==================================================================
      SUBROUTINE vecmul(a,vec1,vec2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION vec1(3), vec2(3)
!     ==--------------------------------------------------------------==
      vec2(1) = a * vec1(1)
      vec2(2) = a * vec1(2)
      vec2(3) = a * vec1(3)
!     ==--------------------------------------------------------------==
      RETURN
      END
!     ==================================================================

!     ==================================================================
      SUBROUTINE addvec(a,vec1,b,vec2,vec3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION vec1(3), vec2(3), vec3(3)
!     ==--------------------------------------------------------------==
      vec3(1)=a*vec1(1)+b*vec2(1)
      vec3(2)=a*vec1(2)+b*vec2(2)
      vec3(3)=a*vec1(3)+b*vec2(3)
!     ==--------------------------------------------------------------==
      RETURN
      END
!     ==================================================================
!     ==================================================================
      SUBROUTINE copyvec(vec1,vec2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION vec1(3), vec2(3)
!     ==--------------------------------------------------------------==
      vec2(1) = vec1(1)
      vec2(2) = vec1(2)
      vec2(3) = vec1(3)
!     ==--------------------------------------------------------------==
      RETURN
      END
!     ==================================================================

!     ==================================================================
      SUBROUTINE vecprod(vec1,vec2,vec3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION vec1(3), vec2(3), vec3(3)
!     ==--------------------------------------------------------------==
      vec3(1) = vec1(2) * vec2(3) - vec1(3) * vec2(2)
      vec3(2) = vec1(3) * vec2(1) - vec1(1) * vec2(3)
      vec3(3) = vec1(1) * vec2(2) - vec1(2) * vec2(1)
!     ==--------------------------------------------------------------==
      RETURN
      END
!     ==================================================================

