!-------------------------------------------------------------------------------
! cubecruncher, a tool to handle cubes
!
! Copyright (C) 2004, Joost VandeVondele
! 
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
! 
!-------------------------------------------------------------------------------
! Execute program without arguments to get a little help  
! 
! tool to handle simple Gaussian cube files 
!   - make cubes containing less points
!   - make cubes centered around a point
!   - merge in coordinates 
!   - fold atom coordinates in the box
!
! currently limited to
!    - orthorombic cubes
!    - machines with enough memory
!    - correct inputs
!    - tested only with Quickstep cubes
!-------------------------------------------------------------------------------
MODULE cubecruncher
  IMPLICIT NONE
  REAL, PARAMETER :: a2au = 1.0/0.5291772083 ! multiply angstrom by a2au to get bohr
  TYPE cube_type
     INTEGER :: npoints(3)
     REAL, DIMENSION(:,:,:), POINTER :: grid
     REAL    :: origin(3)
     REAL    :: dr(3)
     CHARACTER(LEN=80) :: title1,title2
     INTEGER :: Natom
     REAL, DIMENSION(:,:), POINTER :: coords
     INTEGER, DIMENSION(:), POINTER  :: Zatom
     REAL, DIMENSION(:), POINTER :: auxfield
  END TYPE cube_type

CONTAINS

!-------------------------------------------------------------------------------
! initializes a cube to zero
!-------------------------------------------------------------------------------
  SUBROUTINE init_cube(cube)
    TYPE(cube_type), POINTER :: cube

    cube%npoints=0
    NULLIFY(cube%grid)
    cube%origin=0.0
    cube%dr=0.0
    cube%title1=""
    cube%title2=""
    cube%Natom=0
    NULLIFY(cube%coords)
    NULLIFY(cube%Zatom)
    NULLIFY(cube%auxfield)
  END SUBROUTINE init_cube

!-------------------------------------------------------------------------------
! deallocate the internal fields of a cube, resetting stuff
!-------------------------------------------------------------------------------
  SUBROUTINE deallocate_cube(cube)
    TYPE(cube_type), POINTER :: cube
    IF (ASSOCIATED(cube%coords))   DEALLOCATE(cube%coords)
    IF (ASSOCIATED(cube%Zatom))    DEALLOCATE(cube%Zatom)
    IF (ASSOCIATED(cube%auxfield)) DEALLOCATE(cube%auxfield)
    IF (ASSOCIATED(cube%grid))     DEALLOCATE(cube%grid)
    cube%origin=0.0
    cube%dr=0.0
    cube%title1=""
    cube%title2=""
    cube%Natom=0
    cube%npoints=0
    DEALLOCATE(cube)
  END SUBROUTINE deallocate_cube

!-------------------------------------------------------------------------------
! writes to an open unit in a cube respecting format
! stride can be used to reduce the number of points that are actually written
! notice that it just truncates the last points if stride does not divide the number
! of actual points
!-------------------------------------------------------------------------------
  SUBROUTINE write_cube(cube,iunit,stride)
    TYPE(cube_type), POINTER :: cube
    INTEGER                  :: iunit
    INTEGER                  :: stride

    INTEGER                  :: I1,I2,I3

    write(iunit,'(A80)') cube%title1
    write(iunit,'(A80)') cube%title2
    WRITE(iunit,'(I5,3f12.6)') cube%Natom,cube%origin
    WRITE(iunit,'(I5,3f12.6)') cube%npoints(1)/stride,cube%dr(1)*stride,0.0,0.0
    WRITE(iunit,'(I5,3f12.6)') cube%npoints(2)/stride,0.0,cube%dr(2)*stride,0.0
    WRITE(iunit,'(I5,3f12.6)') cube%npoints(3)/stride,0.0,0.0,cube%dr(3)*stride
    DO I1=1,cube%Natom
       WRITE(iunit,'(I5,4f12.6)') cube%Zatom(I1),cube%auxfield(I1),cube%coords(:,I1)
    ENDDO
    DO I1=1,cube%npoints(1),stride
       DO I2=1,cube%npoints(2),stride
        WRITE(iunit,'(6E13.5)') (cube%grid(I1,I2,I3),I3=1,cube%npoints(3),stride)
       ENDDO
    ENDDO
  END SUBROUTINE write_cube

!-------------------------------------------------------------------------------
! reads cube from an open unit to a freshly initialised cube type 
!-------------------------------------------------------------------------------
  SUBROUTINE read_cube(cube,iunit)
    TYPE(cube_type), POINTER :: cube
    INTEGER                  :: iunit

    INTEGER                  :: I1,I2,I3
    REAL                     :: dum

    READ(iunit,'(A80)') cube%title1
    READ(iunit,'(A80)') cube%title2
    READ(iunit,'(I5,3f12.6)') cube%Natom,cube%origin
    READ(iunit,'(I5,3f12.6)') cube%npoints(1), cube%dr(1), dum, dum
    READ(iunit,'(I5,3f12.6)') cube%npoints(2), dum, cube%dr(2), dum
    READ(iunit,'(I5,3f12.6)') cube%npoints(3), dum, dum, cube%dr(3)

    ALLOCATE(cube%Zatom(cube%Natom))
    ALLOCATE(cube%auxfield(cube%Natom))
    ALLOCATE(cube%coords(3,cube%Natom))
    ALLOCATE(cube%grid(cube%npoints(1),cube%npoints(2),cube%npoints(3)))

    DO I1=1,cube%Natom
       READ(iunit,'(I5,4f12.6)') cube%Zatom(I1),cube%auxfield(I1),cube%coords(:,I1)
    ENDDO
    DO I1=1,cube%npoints(1)
       DO I2=1,cube%npoints(2)
        READ(iunit,'(6E13.5)') (cube%grid(I1,I2,I3),I3=1,cube%npoints(3))
       ENDDO
    ENDDO
  END SUBROUTINE read_cube

!-------------------------------------------------------------------------------
! reads xyz into an existing cube type, overwriting the exising atomic data 
!-------------------------------------------------------------------------------
  SUBROUTINE read_xyz(cube,iunit)
    USE periodic_table, ONLY          : label2Z
    TYPE(cube_type), POINTER :: cube
    INTEGER                  :: iunit

    INTEGER                  :: I1
    character(LEN=2)         :: label

    IF (ASSOCIATED(cube%Zatom))    DEALLOCATE(cube%Zatom)
    IF (ASSOCIATED(cube%auxfield)) DEALLOCATE(cube%auxfield)
    IF (ASSOCIATED(cube%coords))   DEALLOCATE(cube%coords)
    
    READ(iunit,*) cube%natom
    READ(iunit,*) ! skip title

    ALLOCATE(cube%Zatom(cube%Natom))
    ALLOCATE(cube%auxfield(cube%Natom))
    ALLOCATE(cube%coords(3,cube%Natom))
    cube%auxfield=0.0

    DO I1=1,cube%natom
       READ(iunit,*) label,cube%coords(:,I1) 
       cube%Zatom(I1)=label2Z(label) 
    ENDDO
    cube%coords=cube%coords*a2au  ! cube coords are in atomic units

  END SUBROUTINE read_xyz

!-------------------------------------------------------------------------------
! fold the cube to center around a given center
!-------------------------------------------------------------------------------
  SUBROUTINE center_cube(cube,center)
    TYPE(cube_type), POINTER :: cube
    REAL, DIMENSION(3) :: center

    INTEGER                  :: I1,I2,I3,O1,O2,O3
    INTEGER, DIMENSION(3)    :: C

    REAL, DIMENSION(:,:,:), POINTER :: grid
    ALLOCATE(grid(cube%npoints(1),cube%npoints(2),cube%npoints(3)))

    ! find coordinates of the lower left corner, centered around the atom
    C=ANINT((center-cube%origin)/cube%dr)-cube%npoints/2
    cube%origin=C*cube%dr+cube%origin
    DO I3=1,cube%npoints(3)
       O3=MODULO(C(3)+I3-1, cube%npoints(3))+1
    DO I2=1,cube%npoints(2)
       O2=MODULO(C(2)+I2-1, cube%npoints(2))+1
    DO I1=1,cube%npoints(1)
       O1=MODULO(C(1)+I1-1, cube%npoints(1))+1
       grid(I1,I2,I3)=cube%grid(O1,O2,O3) 
    ENDDO
    ENDDO
    ENDDO
    cube%grid=grid
    DEALLOCATE(grid) 

  END SUBROUTINE

!-------------------------------------------------------------------------------
! folds the atoms in the cube assuming periodic boundary conditions
! smart_hydrogen=.TRUE. means that hydrogen move to the closest heavy atoms
! and are not really folded
!-------------------------------------------------------------------------------
  SUBROUTINE fold_coords(cube, smart_hydrogen)
    USE utils, ONLY : fold
    TYPE(cube_type), POINTER :: cube
    LOGICAL :: smart_hydrogen 

    INTEGER                  :: I1,I2,closest
    REAL, DIMENSION(3)       :: edges,center,displ
    REAL                     :: d2

    edges=cube%dr*cube%npoints
    center=cube%origin+edges/2.0

    DO I1=1,cube%natom
       cube%coords(:,I1)=fold(cube%coords(:,I1),center=center,edges=edges)
    ENDDO

    ! quadratic search for nearest heavy atom and move hydrogen to it
    IF (smart_hydrogen) THEN
       DO I1=1,cube%natom
          IF (cube%Zatom(I1)==1) THEN
             d2=HUGE(d2)
             closest=-1
             DO I2=1,cube%natom
                IF (cube%Zatom(I2)/=1) THEN
                   displ=fold(cube%coords(:,I1)-cube%coords(:,I2),center=(/0.0,0.0,0.0/),edges=edges) 
                   IF (d2>SUM(displ**2)) THEN
                      closest=I2
                      d2=SUM(displ**2)
                   ENDIF
                ENDIF
             ENDDO
             IF (closest/=-1) THEN
                displ=fold(cube%coords(:,I1)-cube%coords(:,closest),center=(/0.0,0.0,0.0/),edges=edges)
                cube%coords(:,I1)=cube%coords(:,closest)+displ
             ENDIF
          ENDIF
      ENDDO
    ENDIF

  END SUBROUTINE fold_coords

END MODULE cubecruncher
!-------------------------------------------------------------------------------
! simple module to deal with the command line args
! uses the non-Fortran getarg / iargc
!-------------------------------------------------------------------------------
MODULE command_line_tools
  TYPE input_type
    CHARACTER(LEN=200) :: cube_name_in
    CHARACTER(LEN=200) :: cube_name_out
    CHARACTER(LEN=200) :: xyz_name
    LOGICAL :: do_xyz, do_center, do_fold, do_foldsmart,do_center_atom,do_center_geo
    LOGICAL :: have_input, have_output, write_help, write_version
    INTEGER :: atom_center, stride
  END TYPE
CONTAINS
!-------------------------------------------------------------------------------
!  initialise the input with some defaults
!-------------------------------------------------------------------------------
  SUBROUTINE init_input(input)
    TYPE(input_type), INTENT(OUT) :: input
    ! some defaults
    input%do_xyz=.FALSE.
    input%do_center=.FALSE.
    input%do_center_atom=.FALSE.
    input%do_center_geo=.FALSE.
    input%do_fold=.FALSE.
    input%do_foldsmart=.FALSE.
    input%write_help=.FALSE.
    input%write_version=.FALSE.
    input%stride=1
    input%have_input=.FALSE.
    input%have_output=.FALSE.
  END SUBROUTINE init_input
!-------------------------------------------------------------------------------
!  just write an informational message to the output
!-------------------------------------------------------------------------------
   SUBROUTINE write_help()
      write(6,'(A80)') ""
      write(6,'(A80)') "Usage: cubecruncher.x -i input.cube -o output.cube [options]                   " 
      write(6,'(A80)') ""
      write(6,'(A80)') "        transforms 'input.cube' into 'output.cube' doing                       "
      write(6,'(A80)') "        transformations based on the specified options.                        " 
      write(6,'(A80)') "        The following options are available and can normally                   "
      write(6,'(A80)') "        be combined:                                                           "
      write(6,'(A80)') ""
      write(6,'(A80)') " -xyz coords.xyz         : merge in a xyz coordinate file                      "
      write(6,'(A80)') "                              (e.g. for VMD usage)                             "
      write(6,'(A80)') " -center {#|geo}         : #=1..Natom center the cube around atom #            "
      write(6,'(A80)') "                           geo        center the cube so that the cell         "
      write(6,'(A80)') "                              boundaries are as far as possible from the       "
      write(6,'(A80)') "                              molecule                                         "
      write(6,'(A80)') " -fold                   : fold atoms inside the box defined by the cube       "
      write(6,'(A80)') " -foldsmart              : idem, but place hydrogens near heavier atoms        " 
      write(6,'(A80)') " -stride #               : reduces resolution writing grids with stride #      " 
      write(6,'(A80)') " -help                   : print this message                                  "
      write(6,'(A80)') " -v                      : print a version string                              "
      write(6,'(A80)') ""
   END SUBROUTINE write_help
!-------------------------------------------------------------------------------
!  get the command line arguments 
!-------------------------------------------------------------------------------
  SUBROUTINE parse_command_line(input,narg)
    TYPE(input_type), INTENT(INOUT) :: input
    INTEGER, INTENT(OUT) :: narg

    INTEGER :: I
    INTEGER :: iargc 
    CHARACTER(LEN=200) :: arg,nextarg

    narg=iargc()
    DO I=1,narg
      CALL getarg(I,arg)
      SELECT CASE(TRIM(arg))
      CASE("-o")
        input%have_output=.TRUE.
        CALL getarg(I+1,input%cube_name_out)
        WRITE(6,*) "Writing to   cube          ",TRIM(input%cube_name_out)
      CASE("-i")
        input%have_input=.TRUE.
        CALL getarg(I+1,input%cube_name_in)
        WRITE(6,*) "Reading from cube          ",TRIM(input%cube_name_in)
      CASE("-xyz")
        input%do_xyz=.TRUE.
        CALL getarg(I+1,input%xyz_name)
        WRITE(6,*) "Reading from xyz           ",TRIM(input%xyz_name)
      CASE("-center")
        input%do_center=.TRUE.
        CALL getarg(I+1,nextarg)
        IF (nextarg.EQ."geo") THEN
            input%do_center_geo=.TRUE.
            WRITE(6,*) "Centering around the geometric center "
        ELSE
            input%do_center_atom=.TRUE.
            READ(nextarg,*) input%atom_center
            WRITE(6,*) "Centering around atom      ",input%atom_center
        ENDIF
      CASE("-stride")
        CALL getarg(I+1,nextarg)
        READ(nextarg,*) input%stride
        WRITE(6,*) "write cube with stride     ",input%stride
      CASE("-fold")
        input%do_fold=.TRUE.
        WRITE(6,*) "Folding atoms in the cube "
      CASE("-foldsmart")
        input%do_foldsmart=.TRUE.
        WRITE(6,*) "Folding atoms in the cube placing hydrogens near heavy atoms"
      CASE("-help")
        input%write_help=.TRUE.
        WRITE(6,*) "Writing help message "
      CASE("-v")
        input%write_version=.TRUE.
        WRITE(6,*) "Writing version number "
      END SELECT
    ENDDO
  END SUBROUTINE parse_command_line
END MODULE

PROGRAM main
   USE cubecruncher
   USE periodic_table
   USE command_line_tools
   IMPLICIT NONE

   TYPE(cube_type), POINTER :: cube

   INTEGER :: iunit_cube_in,iunit_cube_out,iunit_xyz,narg
   LOGICAL :: did_something
   REAL, DIMENSION(3) :: center
   TYPE(input_type) :: input
  
   CALL init_input(input)
   CALL parse_command_line(input,narg)
   did_something=.FALSE.
   !
   ! do simple checking of input and output
   !
   IF (input%write_version) THEN
      write(6,*) "This is cubecruncher version 1.0"
      did_something=.TRUE.
   ENDIF
   IF (input%write_help .OR. narg<1) THEN
      CALL write_help()
      did_something=.TRUE.
   ENDIF
   IF (.NOT. (input%have_input .AND. input%have_output)) THEN
      IF (.NOT. did_something) THEN
        write(6,*) "require the input and output cube file names"
        write(6,*) "use the -help option for more information"
      ENDIF
   ELSE
     ! init periodic table
     CALL init_periodic_table() 

     ALLOCATE(cube)
     CALL init_cube(cube)
     CALL init_cube(cube)
   
     iunit_cube_in=117
     iunit_cube_out=118
     iunit_xyz=119

     ! read cube
     write(6,FMT='(A)',ADVANCE="No") "Reading cube ... "
     OPEN(iunit_cube_in,FILE=TRIM(input%cube_name_in)) 
     CALL read_cube(cube,iunit_cube_in)
     CLOSE(iunit_cube_in)
     write(6,*) "Done"

     IF (input%do_xyz) THEN  
        ! read xyz
        write(6,FMT='(A)',ADVANCE="No") "Reading xyz ... "
        OPEN(iunit_xyz,FILE=TRIM(input%xyz_name))
        CALL read_xyz(cube,iunit_xyz)
        CLOSE(iunit_xyz)
        write(6,*) "Done"
     ENDIF
  
     ! fold grid around a given atom or the geometrical center of the molecule
     IF (input%do_center) THEN
        write(6,FMT='(A)',ADVANCE="No") "Centering cube ... "
        IF (input%do_center_atom) THEN
           IF (input%atom_center .LT. 1 .OR. input%atom_center.GT.cube%Natom) THEN
              STOP "wrong atom center"
           ENDIF
           center=cube%coords(:,input%atom_center)
        ENDIF
        IF (input%do_center_geo) THEN
           center(1)=0.5*(MINVAL(cube%coords(1,:))+MAXVAL(cube%coords(1,:)))
           center(2)=0.5*(MINVAL(cube%coords(2,:))+MAXVAL(cube%coords(2,:)))
           center(3)=0.5*(MINVAL(cube%coords(3,:))+MAXVAL(cube%coords(3,:)))
        ENDIF
        CALL center_cube(cube,center)
        write(6,*) "Done"
     ENDIF

     ! fold atoms in box, bringing H's near their parents
     IF (input%do_fold .OR. input%do_foldsmart) THEN
        write(6,FMT='(A)',ADVANCE="No") "Folding coordinates ... "
        CALL fold_coords(cube,input%do_foldsmart)
        write(6,*) "Done"
     ENDIF

     ! final out
     write(6,FMT='(A)',ADVANCE="No") "Writing cube ... "
     OPEN(iunit_cube_out,FILE=TRIM(input%cube_name_out)) 
     CALL write_cube(cube,iunit_cube_out,input%stride)
     CLOSE(iunit_cube_out)
     write(6,*) "Done"

     ! deallocate stuff
     CALL deallocate_cube(cube)

  ENDIF
END PROGRAM
