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
!    - machines with enough memory
!    - correct inputs
!    - tested only with Quickstep cubes
!-------------------------------------------------------------------------------
MODULE cubecruncher
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND ( 14, 200 )
  REAL(KIND=dp), PARAMETER :: a2au = 1.0/0.5291772083 ! multiply angstrom by a2au to get bohr
  TYPE cube_type
     INTEGER :: npoints(3)
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER :: grid
     REAL(KIND=dp)    :: origin(3)
     REAL(KIND=dp)    :: dh(3,3), h(3,3),  hinv(3,3)
     CHARACTER(LEN=80) :: title1,title2
     INTEGER :: Natom
     REAL(KIND=dp), DIMENSION(:,:), POINTER :: coords
     INTEGER, DIMENSION(:), POINTER  :: Zatom
     REAL(KIND=dp), DIMENSION(:), POINTER :: auxfield
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
    cube%dh=0.0
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
    cube%dh=0.0
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
  SUBROUTINE write_cube(cube,iunit,stride,mult_fact)
    TYPE(cube_type), POINTER :: cube
    INTEGER                  :: iunit
    INTEGER                  :: stride
    REAL(KIND=dp)                     :: mult_fact

    INTEGER                  :: I1,I2,I3

    write(iunit,'(A80)') cube%title1
    write(iunit,'(A80)') cube%title2
    WRITE(iunit,'(I5,3f12.6)') cube%Natom,cube%origin
    WRITE(iunit,'(I5,3f12.6)') 1+(cube%npoints(1)-1)/stride,cube%dh(1,1)*stride,&
                               cube%dh(2,1)*stride,cube%dh(3,1)*stride
    WRITE(iunit,'(I5,3f12.6)') 1+(cube%npoints(2)-1)/stride,cube%dh(1,2)*stride,&
                               cube%dh(2,2)*stride,cube%dh(3,2)*stride
    WRITE(iunit,'(I5,3f12.6)') 1+(cube%npoints(3)-1)/stride,cube%dh(1,3)*stride,&
                               cube%dh(2,3)*stride,cube%dh(3,3)*stride
    DO I1=1,cube%Natom
       WRITE(iunit,'(I5,4f12.6)') cube%Zatom(I1),cube%auxfield(I1),cube%coords(:,I1)
    ENDDO
    DO I1=1,cube%npoints(1),stride
       DO I2=1,cube%npoints(2),stride
        WRITE(iunit,'(6E13.5)') (cube%grid(I1,I2,I3)*mult_fact,I3=1,cube%npoints(3),stride)
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
    REAL(KIND=dp)                     :: dum

    READ(iunit,'(A80)') cube%title1
    READ(iunit,'(A80)') cube%title2
    READ(iunit,'(I5,3f12.6)') cube%Natom,cube%origin
    READ(iunit,'(I5,3f12.6)') cube%npoints(1), cube%dh(1,1), cube%dh(2,1),cube%dh(3,1)
    cube%h(1,1) = cube%npoints(1)*cube%dh(1,1)
    cube%h(2,1) = cube%npoints(1)*cube%dh(2,1)
    cube%h(3,1) = cube%npoints(1)*cube%dh(3,1)
    READ(iunit,'(I5,3f12.6)') cube%npoints(2), cube%dh(1,2), cube%dh(2,2),cube%dh(3,2)
    cube%h(1,2) = cube%npoints(2)*cube%dh(1,2)
    cube%h(2,2) = cube%npoints(2)*cube%dh(2,2)
    cube%h(3,2) = cube%npoints(2)*cube%dh(3,2)
    READ(iunit,'(I5,3f12.6)') cube%npoints(3), cube%dh(1,3), cube%dh(2,3),cube%dh(3,3)
    cube%h(1,3) = cube%npoints(3)*cube%dh(1,3)
    cube%h(2,3) = cube%npoints(3)*cube%dh(2,3)
    cube%h(3,3) = cube%npoints(3)*cube%dh(3,3)

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
! reads xyz into an existing cube type, overwriting the existing atomic data
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

  SUBROUTINE boxmat(cube)

    TYPE(cube_type), POINTER :: cube
    integer :: i
    REAL(KIND=dp) :: uvol, vol


  write(*,*) cube%h(1,1), cube%h(2,1), cube%h(3,1)
  write(*,*) cube%h(1,2), cube%h(2,2), cube%h(3,2)
  write(*,*) cube%h(1,3), cube%h(2,3), cube%h(3,3)

    vol=cube%h(1,1)*(cube%h(2,2)*cube%h(3,3)-cube%h(2,3)*cube%h(3,2)) &
       -cube%h(2,1)*(cube%h(1,2)*cube%h(3,3)-cube%h(3,2)*cube%h(1,3)) &
       +cube%h(3,1)*(cube%h(1,2)*cube%h(2,3)-cube%h(2,2)*cube%h(1,3))

    vol = ABS(vol)
  write(*,*) vol
    IF(vol < 1.0E-10) stop
    uvol=1.d0/vol
    cube%hinv(1,1)= (cube%h(2,2)*cube%h(3,3)-cube%h(3,2)*cube%h(2,3))*uvol
    cube%hinv(2,1)=-(cube%h(2,1)*cube%h(3,3)-cube%h(3,1)*cube%h(2,3))*uvol
    cube%hinv(3,1)= (cube%h(2,1)*cube%h(3,2)-cube%h(3,1)*cube%h(2,2))*uvol
    cube%hinv(2,2)= (cube%h(1,1)*cube%h(3,3)-cube%h(1,3)*cube%h(3,1))*uvol
    cube%hinv(3,2)=-(cube%h(1,1)*cube%h(3,2)-cube%h(1,2)*cube%h(3,1))*uvol
    cube%hinv(3,3)= (cube%h(1,1)*cube%h(2,2)-cube%h(1,2)*cube%h(2,1))*uvol
    cube%hinv(1,2)=-(cube%h(1,2)*cube%h(3,3)-cube%h(1,3)*cube%h(3,2))*uvol
    cube%hinv(1,3)= (cube%h(1,2)*cube%h(2,3)-cube%h(1,3)*cube%h(2,2))*uvol
    cube%hinv(2,3)=-(cube%h(1,1)*cube%h(2,3)-cube%h(2,1)*cube%h(1,3))*uvol

  END SUBROUTINE

function pbc(r,cube) result(r_pbc)
   
    TYPE(cube_type), POINTER :: cube
    REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: r
    REAL(KIND=dp), DIMENSION(3) :: r_pbc
    REAL(KIND=dp), DIMENSION(3) :: s

    s=MATMUL(cube%hinv,r)
    s=s-ANINT(s)
    r_pbc=MATMUL(cube%h,s)

end function pbc

!-------------------------------------------------------------------------------
! fold the cube to center around a given center
!-------------------------------------------------------------------------------
  SUBROUTINE center_cube(cube,center)
    TYPE(cube_type), POINTER :: cube
    REAL(KIND=dp), DIMENSION(3) :: center

    INTEGER                  :: I1,I2,I3,O1,O2,O3
    INTEGER, DIMENSION(3)    :: C
    REAL(KIND=dp) :: h(3,3), hinv(3,3), uvol, vol, dvec(3)
    REAL(KIND=dp), DIMENSION(:,:,:), POINTER :: grid
    ALLOCATE(grid(cube%npoints(1),cube%npoints(2),cube%npoints(3)))

    ! notice h and hinv here correspond to dh and dhinv, do not overwrite cube%xxx fields with them
    h = cube%dh
    vol=h(1,1)*(h(2,2)*h(3,3)-h(2,3)*h(3,2)) &
          -h(2,1)*(h(1,2)*h(3,3)-h(3,2)*h(1,3)) &
          +h(3,1)*(h(1,2)*h(2,3)-h(2,2)*h(1,3))
    uvol=1.d0/vol
    hinv(1,1)= (h(2,2)*h(3,3)-h(3,2)*h(2,3))*uvol
    hinv(2,1)=-(h(2,1)*h(3,3)-h(3,1)*h(2,3))*uvol
    hinv(3,1)= (h(2,1)*h(3,2)-h(3,1)*h(2,2))*uvol
    hinv(2,2)= (h(1,1)*h(3,3)-h(1,3)*h(3,1))*uvol
    hinv(3,2)=-(h(1,1)*h(3,2)-h(1,2)*h(3,1))*uvol
    hinv(3,3)= (h(1,1)*h(2,2)-h(1,2)*h(2,1))*uvol
    hinv(1,2)=-(h(1,2)*h(3,3)-h(1,3)*h(3,2))*uvol
    hinv(1,3)= (h(1,2)*h(2,3)-h(1,3)*h(2,2))*uvol
    hinv(2,3)=-(h(1,1)*h(2,3)-h(2,1)*h(1,3))*uvol

    ! find coordinates of the lower left corner, centered around the atom

    dvec(1) = hinv(1,1)*(center(1)-cube%origin(1)) + &
              hinv(1,2)*(center(2)-cube%origin(2)) + &
              hinv(1,3)*(center(3)-cube%origin(3))
    dvec(2) = hinv(2,1)*(center(1)-cube%origin(1)) + &
              hinv(2,2)*(center(2)-cube%origin(2)) + &
              hinv(2,3)*(center(3)-cube%origin(3))
    dvec(3) = hinv(3,1)*(center(1)-cube%origin(1)) + &
              hinv(3,2)*(center(2)-cube%origin(2)) + &
              hinv(3,3)*(center(3)-cube%origin(3))
    C=ANINT(dvec)-cube%npoints/2
    cube%origin(1)=C(1)*h(1,1)+C(2)*h(1,2)+C(3)*h(1,3)+cube%origin(1)
    cube%origin(2)=C(1)*h(2,1)+C(2)*h(2,2)+C(3)*h(2,3)+cube%origin(2)
    cube%origin(3)=C(1)*h(3,1)+C(2)*h(3,2)+C(3)*h(3,3)+cube%origin(3)
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
! subtract cubes, cube1=cube1-cube2
!-------------------------------------------------------------------------------
  SUBROUTINE subtract_cubes(cube1,cube2)
    TYPE(cube_type), POINTER :: cube1,cube2
    cube1%grid=cube1%grid-cube2%grid
  END SUBROUTINE


!-------------------------------------------------------------------------------
! folds the atoms in the cube assuming periodic boundary conditions
! smart_hydrogen=.TRUE. means that hydrogen move to the closest heavy atoms
! and are not really folded
!-------------------------------------------------------------------------------
  SUBROUTINE fold_coords(cube, smart_hydrogen)
    TYPE(cube_type), POINTER :: cube
    LOGICAL :: smart_hydrogen 

    INTEGER                  :: i, I1,I2,closest
    INTEGER :: n1, n2, n3, no1, no2, no3
    REAL(KIND=dp), DIMENSION(3)       :: edges,center,displ
    REAL(KIND=dp)                     :: d2, x,y,z
    REAL(KIND=dp) :: h11h22,h11h33,h11h32,h11h23,h33h22,h31h22,h13h22,h31h13,h31h23,&
           h12h21,h12h23,h12h31,h12h33,h32h21,h32h23,h32h13,h21h32,h21h33,&
           h21h13,h32h11h23,h12h23h31,h21h32h13,h11h22h33,h31h22h13,h21h12h33,deth
    REAL(KIND=dp), DIMENSION(3)       :: A, O, r, r_pbc

    center(:)=cube%origin(:)+cube%h(:,1)/2.0 + cube%h(:,2)/2.0 + cube%h(:,3)/2.0

    DO I1=1,cube%natom
       r = cube%coords(:,I1)-center
       r_pbc = pbc(r,cube)
       cube%coords(:,I1) = r_pbc+center
    ENDDO

    ! quadratic search for nearest heavy atom and move hydrogen to it
    IF (smart_hydrogen) THEN
       DO I1=1,cube%natom
          IF (cube%Zatom(I1)==1) THEN
             d2=HUGE(d2)
             closest=-1
             DO I2=1,cube%natom
                IF (cube%Zatom(I2)/=1) THEN
                   displ=pbc(cube%coords(:,I1)-cube%coords(:,I2),cube)
                   IF (d2>SUM(displ**2)) THEN
                      closest=I2
                      d2=SUM(displ**2)
                   ENDIF
                ENDIF
             ENDDO
             IF (closest/=-1) THEN
                displ=pbc(cube%coords(:,I1)-cube%coords(:,closest),cube)
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
  USE cubecruncher, ONLY : dp
  TYPE input_type
    CHARACTER(LEN=200) :: cube_name_in
    CHARACTER(LEN=200) :: cube_name_out
    CHARACTER(LEN=200) :: cube_name_subtract
    CHARACTER(LEN=200) :: cube_name_espot
    CHARACTER(LEN=200) :: xyz_name
    LOGICAL :: do_xyz, do_center, do_fold, do_foldsmart,do_center_atom,do_center_geo, do_center_point
    LOGICAL :: have_input, have_output, write_help, write_version, do_subtract
    LOGICAL :: do_iso, do_slice, do_1d_profile, espot_over_iso, iso_current, have_espot_cube
    LOGICAL :: from_top, check_atom_from_down
    INTEGER :: atom_center, stride, cube, index_profile
    REAL(KIND=dp)    :: center_point(3), slice(3), iso_level(3), iso_hmin, iso_hmax, delta_profile,&
                mult_fact, iso_delta_grid, bwidth
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
    input%do_center_point=.FALSE.
    input%center_point=(/0.0,0.0,0.0/)
    input%slice=(/0.0,0.0,0.0/)
    input%do_slice=.FALSE.
    input%iso_level=(/0.0,0.0,0.0/)
    input%do_iso=.FALSE.
    input%iso_hmin=-1.0
    input%iso_hmax=-1.0
    input%iso_delta_grid=0.5D0
    input%iso_current = .FALSE.
    input%bwidth=2.0D0
    input%espot_over_iso = .FALSE.
    input%from_top = .FALSE.
    input%check_atom_from_down = .FALSE.
    input%index_profile=0
    input%delta_profile=0.0D0
    input%do_1d_profile=.FALSE.
    input%do_fold=.FALSE.
    input%do_foldsmart=.FALSE.
    input%do_subtract=.FALSE.
    input%write_help=.FALSE.
    input%write_version=.FALSE.
    input%stride=1
    input%mult_fact=1.0D0
    input%have_input=.FALSE.
    input%have_output=.FALSE.
    input%have_espot_cube=.FALSE.
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
      write(6,'(A80)') " -subtract min.cube      : subtract min.cube from input.cube                   "
      write(6,'(A80)') " -center {#|geo|point}   : #=1..Natom center the cube around atom #            "
      write(6,'(A80)') "                           geo        center the cube so that the cell         "
      write(6,'(A80)') "                              boundaries are as far as possible from the       "
      write(6,'(A80)') "                              molecule                                         "
      write(6,'(A80)') "                           point      center around a fixed point              "
      write(6,'(A80)') " -fold                   : fold atoms inside the box defined by the cube       "
      write(6,'(A80)') " -point x y z            : coordinates of the center point                     "
      write(6,'(A80)') " -foldsmart              : idem, but place hydrogens near heavier atoms        " 
      write(6,'(A80)') " -stride #               : reduces resolution writing grids with stride #      " 
      write(6,'(A80)') " -multiply #             : multiply all values of cube by the given factor     " 
      write(6,'(A80)') " -1d_profile dir delta   : compute the profile along the cartesian axis dir    "
      write(6,'(A80)') "                          and a smoothed profile with smoothing interval delta " 
      write(6,'(A80)') " -slice h1 h2 h3         : extract the values of the cube slice at height h#,  "
      write(6,'(A80)') "                           along the cartesian axis #, if h#/=0,               "
      write(6,'(A80)') "                           reads three reals, takes the first /=0              "
      write(6,'(A80)') " -iso l1 l2 l3           : compute the isosurface at level l#, recording height"
      write(6,'(A80)') "                           along the axis #, if l#/=0,                         "
      write(6,'(A80)') "                           reads three reals, takes the first /=0              "
      write(6,'(A80)') " -iso_delta_grid  delta  : step of the regular grid to output the iso-surface  " 
      write(6,'(A80)') " -isocurrent W           : correct the iso-density surface to approximate the  "
      write(6,'(A80)') "                           iso-current surface, where                          "
      write(6,'(A80)') "                           I(z) \propto rho(z)*exp(alpha*W*sqrt(ES-pot(z)))   "
      write(6,'(A80)') "                           The ES-pot cube must be loaded                      "            
      write(6,'(A80)') " -espot                  : assign ES-pot values to the calculated iso-surface  "            
      write(6,'(A80)') "                           The ES-pot cube must be loaded                      "            
      write(6,'(A80)') " -espot_cube             : name cube file of the electrostatic potential       " 
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
      CASE("-subtract")
        input%do_subtract=.TRUE.
        CALL getarg(I+1,input%cube_name_subtract)
        WRITE(6,*) "Subtracting cube           ",TRIM(input%cube_name_subtract)
      CASE("-xyz")
        input%do_xyz=.TRUE.
        CALL getarg(I+1,input%xyz_name)
        WRITE(6,*) "Reading from xyz           ",TRIM(input%xyz_name)
      CASE("-center")
        input%do_center=.TRUE.
        CALL getarg(I+1,nextarg)
        SELECT CASE(nextarg)
        CASE("geo")
            input%do_center_geo=.TRUE.
            WRITE(6,*) "Centering around the geometric center "
        CASE("point")
            input%do_center_point=.TRUE.
            WRITE(6,*) "Centering around a given point center "
        CASE DEFAULT
            input%do_center_atom=.TRUE.
            READ(nextarg,*) input%atom_center
            WRITE(6,*) "Centering around atom      ",input%atom_center
        END SELECT
      CASE("-point")
        CALL getarg(I+1,nextarg)
        READ(nextarg,*) input%center_point(1)
        CALL getarg(I+2,nextarg)
        READ(nextarg,*) input%center_point(2)
        CALL getarg(I+3,nextarg)
        READ(nextarg,*) input%center_point(3)
      CASE("-slice")
        input%do_slice=.TRUE.
        CALL getarg(I+1,nextarg)
        READ(nextarg,*) input%slice(1)
        CALL getarg(I+2,nextarg)
        READ(nextarg,*) input%slice(2)
        CALL getarg(I+3,nextarg)
        READ(nextarg,*) input%slice(3)
      CASE("-iso")
        input%do_iso=.TRUE.
        CALL getarg(I+1,nextarg)
        READ(nextarg,*) input%iso_level(1)
        CALL getarg(I+2,nextarg)
        READ(nextarg,*) input%iso_level(2)
        CALL getarg(I+3,nextarg)
        READ(nextarg,*) input%iso_level(3)
        WRITE(6,'(A,3(1PE12.4))') "Calculate an iso-surface from cube, at value ",input%iso_level
      CASE("-isorange")
        CALL getarg(I+1,nextarg)
        READ(nextarg,*) input%iso_hmin
        CALL getarg(I+2,nextarg)
        READ(nextarg,*) input%iso_hmax
      CASE("-isocurrent")
        input%iso_current = .TRUE.
        CALL getarg(I+1,nextarg)
        READ(nextarg,*) input%bwidth
        WRITE(6,'(A,f10.5)') "Correct the iso-density surface to a better approx. of iso-current, W= ",input%bwidth
      CASE("-espot")
        input%espot_over_iso=.TRUE.
        WRITE(6,'(A)') "Assign the ES-pot values to the iso-surface"
      CASE("-espot_cube")
        input%have_espot_cube=.TRUE.
        CALL getarg(I+1,input%cube_name_espot)
        WRITE(6,*) "Load electrostatic potential cube file ",TRIM(input%cube_name_espot)
      CASE("-iso_delta_grid")
        CALL getarg(I+1,nextarg)
        READ(nextarg,*) input%iso_delta_grid
        WRITE(6,*) " Integration parameter to generate the 2D grid", input%iso_delta_grid
      CASE("-iso_from_top")
        input%from_top=.TRUE.
        WRITE(6,*) " Search for the iso-current surface starting from top of the box "
      CASE("-iso_stop_at_atom")
        input%check_atom_from_down=.TRUE.
        WRITE(6,*) " Does not allow the tip to go through atoms "
      CASE("-1d_profile")
        input%do_1d_profile=.TRUE.
        CALL getarg(I+1,nextarg)
        READ(nextarg,*) input%index_profile
        CALL getarg(I+2,nextarg)
        READ(nextarg,*) input%delta_profile
      CASE("-stride")
        CALL getarg(I+1,nextarg)
        READ(nextarg,*) input%stride
        WRITE(6,*) "write cube with stride     ",input%stride
      CASE("-multiply")
        CALL getarg(I+1,nextarg)
        READ(nextarg,*) input%mult_fact
        WRITE(6,*) "multiply cube values by factor ",input%mult_fact
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

MODULE cube_post_process

!  USE cubecruncher, ONLY: cube_type
   USE cubecruncher

  IMPLICIT NONE
!  REAL(KIND=dp), PARAMETER :: a2au = 1.0/0.5291772083D0 ! multiply angstrom by a2au to get bohr
  REAL(KIND=dp), PARAMETER :: au2a = 0.5291772083D0 ! multiply bohr by au2a to get angstrom 
  REAL(KIND=dp), PARAMETER :: au2eV = 27.211384D0

CONTAINS

  SUBROUTINE compute_profile(cube,index,delta)

   TYPE(cube_type), POINTER :: cube
   INTEGER, INTENT(IN) :: index
   REAL(KIND=dp), INTENT(IN) :: delta

   CHARACTER(LEN=80) :: filename
   INTEGER :: i, ip_low, ip_up,j,k, ncount, np, npp
   REAL(KIND=dp),DIMENSION(:), ALLOCATABLE :: pos, av_val
   REAL(KIND=dp) :: length, p, val, dh1, dh2

   ALLOCATE(pos(cube%npoints(index)))
   ALLOCATE(av_val(cube%npoints(index)))
   DO i=1,cube%npoints(index)
     pos(i) = (i-1)*cube%dh(index,index)*au2a+cube%origin(index)*au2a
   END DO
   IF(index==3) THEN

     IF(cube%dh(1,3)/=0.0d0 .OR. cube%dh(2,3)/=0.0D0) THEN
      write(*,*) ' Calculation of profile along a non orthogonal direction not implemented'
      STOP
     END IF
     ncount = cube%npoints(1)*cube%npoints(2)
!dbg
!     dh1=sqrt(cube%dh(1,1)*cube%dh(1,1)*au2a*au2a+cube%dh(2,1)*cube%dh(2,1)*au2a*au2a)
!     dh2=sqrt(cube%dh(1,2)*cube%dh(1,2)*au2a*au2a+cube%dh(2,2)*cube%dh(2,2)*au2a*au2a)
!     dh1=cube%dh(1,1)*au2a
!     dh2=cube%dh(2,2)*au2a
!  write(*,*) ' dh ', dh1, dh2, dh1*dh2,ncount,cube%dh(1,1)*cube%dh(2,2)*au2a*au2a
!dbg
     DO i=1,cube%npoints(3)
       val = 0.0D0
       DO j=1,cube%npoints(2)
         DO k=1,cube%npoints(1)
           val =val + cube%grid(k,j,i)! *dh1*dh2
         END DO
       END DO
       av_val(i) = val /REAL(ncount,KIND=dp)
     END DO

   ELSE IF(index==1) THEN

     IF(cube%dh(2,1)/=0.0d0 .OR. cube%dh(3,1)/=0.0D0) THEN
      write(*,*) ' Calculation of profile along a non orthogonal direction not implemented'
      STOP
     END IF
     ncount = cube%npoints(2)*cube%npoints(3)
     DO i=1,cube%npoints(1)
       val = 0.0D0
       DO j=1,cube%npoints(3)
         DO k=1,cube%npoints(2)
           val =val + cube%grid(i,k,j)
         END DO
       END DO
       av_val(i) = val/REAL(ncount,KIND=dp)
     END DO

   ELSE IF(index==2) THEN

    IF(cube%dh(2,1)/=0.0d0 .OR. cube%dh(2,3)/=0.0D0) THEN
      write(*,*) ' Calculation of profile along a non orthogonal direction not implemented'
      STOP
     END IF
     ncount = cube%npoints(1)*cube%npoints(3)
     DO i=1,cube%npoints(2)
       val = 0.0D0
       DO j=1,cube%npoints(3)
         DO k=1,cube%npoints(1)
           val =val + cube%grid(k,i,j)
         END DO
       END DO
       av_val(i) = val/REAL(ncount,KIND=dp)
     END DO


   END IF

   WRITE(filename,'(A8,I1,A4)')  'profile_',index,'.dat'
   OPEN(101,FILE=TRIM(filename)) 
   DO i =1,cube%npoints(index)
     WRITE(101,*) pos(i), av_val(i)
   END DO
   CLOSE(101)

   WRITE(filename,'(A12,I1,A4)')  'profile_int_',index,'.dat'
   OPEN(101,FILE=TRIM(filename)) 
   write(*,*) 'delta ', delta
   length = pos(cube%npoints(index))-pos(1)
   write(*,*) 'length ', length
   np = INT(length/delta) + 1
   npp = INT(delta/(cube%dh(index,index)*au2a)) +1
   write(*,*) 'integration intervals and points per interval ', np, npp

   DO j=1,np
     ip_low = (j-1)*npp +1
     ip_up  = MIN(j*npp,cube%npoints(index)) 
     p = pos(ip_low)+delta/2.0D0
     val = 0.0D0
     DO i =ip_low,ip_up
       val = val + av_val(i)
     END DO
     val = val/REAL(ip_up-ip_low+1,KIND=dp)
     WRITE(101,*) p, val
     IF(ip_up==cube%npoints(index)) EXIT
   END DO
   CLOSE(101)
   

   DEALLOCATE(pos, av_val)
  END SUBROUTINE compute_profile

  SUBROUTINE  compute_slice(cube,dir,height)

   TYPE(cube_type), POINTER :: cube
   INTEGER, INTENT(IN) ::dir 
   REAL(KIND=dp), INTENT(IN) :: height

   CHARACTER(LEN=80) :: filename
   INTEGER :: index_plane, d1,d2
   INTEGER :: i, j,k
   REAL(KIND=dp) :: r, r1, r2 , val


   d1=dir+1; if(d1>3) d1=1
   d2=dir-1; if(d2<1) d2=3
   IF(d1>d2) THEN
     d2=d2+2; d1=d1-2
   END IF

   DO i = 1,cube%npoints(dir)
      r=sqrt((real(i-1,KIND=dp)*cube%dh(1,dir)-cube%origin(1))**2&
               +(real(i-1,KIND=dp)*cube%dh(2,dir)-cube%origin(2))**2+&
                (real(i-1,KIND=dp)*cube%dh(3,dir)-cube%origin(3))**2)
      IF(r*au2a>=height) THEN
        index_plane = i
        EXIT
      END IF
   END DO
   write(*,*) ' selected slice near the grid plane labelled ', Index_plane
   

   WRITE(filename,'(A6,I1,A2,I3,A4)')  'slice_',dir,'_h',index_plane,'.dat'
   OPEN(101,FILE=TRIM(filename))

   DO i = 1,cube%npoints(d1)
     DO j = 1,cube%npoints(d2)

       IF(dir==3)THEN
         val = (cube%grid(i,j,index_plane-1)+cube%grid(i,j,index_plane))*0.5D0
         r1 = real(i-1,KIND=dp)*cube%dh(1,1)+real(j-1,KIND=dp)*cube%dh(1,2)+real(index_plane-1,KIND=dp)*cube%dh(1,3)
         r2 = real(i-1,KIND=dp)*cube%dh(2,1)+real(j-1,KIND=dp)*cube%dh(2,2)+real(index_plane-1,KIND=dp)*cube%dh(2,3)
       ELSEIF(dir==1) THEN
         val = (cube%grid(index_plane-1,i,j)+cube%grid(index_plane,i,j))*0.5D0
         r1 = real(index_plane-1,KIND=dp)*cube%dh(2,1)+real(i-1,KIND=dp)*cube%dh(2,2)+real(j-1,KIND=dp)*cube%dh(2,3)
         r2 = real(index_plane-1,KIND=dp)*cube%dh(3,1)+real(i-1,KIND=dp)*cube%dh(3,2)+real(j-1,KIND=dp)*cube%dh(3,3)
       ELSEIF(dir==2) THEN
         val = (cube%grid(i,index_plane-1,j)+cube%grid(i,index_plane,j))*0.5D0
         r1 = real(i-1,KIND=dp)*cube%dh(1,1)+real(index_plane-1,KIND=dp)*cube%dh(1,2)+real(j-1,KIND=dp)*cube%dh(1,3)
         r2 = real(i-1,KIND=dp)*cube%dh(3,1)+real(index_plane-1,KIND=dp)*cube%dh(3,2)+real(j-1,KIND=dp)*cube%dh(3,3)
       END IF

       WRITE(101,*) r1*au2a, r2*au2a, val 
     END DO
   END DO

   CLOSE(101)
  END SUBROUTINE compute_slice

  SUBROUTINE compute_iso_surf(cube,dir,iso_val_inp,h_min, h_max, grid_delta, &
                   use_espot, iso_current, from_top, check_atom_from_down, bwidth, cube_espot)

   TYPE(cube_type), POINTER :: cube
   INTEGER, INTENT(IN) ::dir
   REAL(KIND=dp), INTENT(IN) :: iso_val_inp, h_min, h_max, grid_delta
   LOGICAL :: use_espot, iso_current
   LOGICAL :: check_atom_from_down,  from_top
   REAL(KIND=dp), INTENT(IN), OPTIONAL ::  bwidth 
   TYPE(cube_type), OPTIONAL, POINTER :: cube_espot

   CHARACTER(LEN=80) :: filename
   LOGICAL :: found_val
   INTEGER ::  d1,d2, iold, jold, kold
   INTEGER :: i, iat, ipoint, j,k, kval, index_min, index_max
   INTEGER :: i1,i2,incr
   REAL(KIND=dp) :: d12, iso_val
   REAL(KIND=dp) :: alpha, fac1, fac2, r, r1, r2 , val1, val2, val1_e, val2_e
   REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: ah
   REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:) :: data_grid
   REAL(KIND=dp) :: delta, det, height, L1, L2, r1_shift, r2_shift, tmp
   INTEGER :: ii, jj, im, ip, jm, jp, np1, np2 


   found_val = .false.
   d1=dir+1; if(d1>3) d1=1
   d2=dir-1; if(d2<1) d2=3
   IF(d1>d2)THEN
      d2=d2+2; d1=d1-2
   END IF

   allocate(data_grid(3,cube%npoints(d1),cube%npoints(d2)))
   data_grid = 0.0D0
   iold = 1
   jold = 1
   kold = 1

   !Find starting plane
   IF(h_min<-10000.0D0) THEN
     index_min = 1
   ELSE
     DO i = 1,cube%npoints(dir)
!       r=sqrt((real(i-1,KIND=dp)*cube%dh(1,dir)+cube%origin(1))**2&
!               +(real(i-1,KIND=dp)*cube%dh(2,dir)+cube%origin(2))**2+&
!                (real(i-1,KIND=dp)*cube%dh(3,dir)+cube%origin(3))**2)
       r  = cube%origin(dir) + real(i-1,KIND=dp)*cube%dh(dir,dir)
       IF(r*au2a>=h_min)THEN
         index_min = i-1
         EXIT
       END IF
     END DO
   END IF
   !Find starting plane
   IF(h_max<-10000.0D0) THEN
     index_max = cube%npoints(dir) - 1
   ELSE
     DO i = index_min,cube%npoints(dir)
!       r=sqrt((real(i-1,KIND=dp)*cube%dh(1,dir)+cube%origin(1))**2&
!               +(real(i-1,KIND=dp)*cube%dh(2,dir)+cube%origin(2))**2+&
!                (real(i-1,KIND=dp)*cube%dh(3,dir)+cube%origin(3))**2)
       r  = cube%origin(dir) + real(i-1,KIND=dp)*cube%dh(dir,dir)
!dbg
!    write(*,*) index_min, i, r, r*au2a, h_max
       IF(r*au2a>=h_max)THEN
         index_max = i
         EXIT
       END IF
     END DO
   END IF
   IF(index_max==cube%npoints(dir)) index_max = index_max-1

  !dbg
!    write(*,*) index_min, index_max,  h_min, h_max
!   stop
  !dbg
   IF(from_top) THEN
      i1=index_max
      i2=index_min
      incr=-1
   ELSE
      i1=index_min
      i2=index_max
      incr=1
   END IF
   IF(iso_current) THEN
     alpha = 2.D0*0.5123D0*bwidth*sqrt(au2ev)
     i = cube%npoints(d1)/2
     j = cube%npoints(d2)/2
!dbg
!  write(*,*) 'center grid ', i, j, ' k range ', i1,i2, incr
!dbg
     DO k = i1,i2,incr
       IF(dir==3)THEN
         val1 = cube%grid(i,j,k)
         val2 = cube%grid(i,j,k+1)
!dbg
!  write(*,*) 'cube ', i,j,k, val1, val2
!dbg
       ELSEIF(dir==2)THEN
         val1 = cube%grid(i,k,j)
         val2 = cube%grid(i,k+1,j)
       ELSEIF(dir==1)THEN
         val1 = cube%grid(k,i,j)
         val2 = cube%grid(k+1,i,j)
       END IF
       IF(val1<val2) THEN
           IF(iso_val_inp<val2 .AND. iso_val_inp>val1) THEN
             found_val = .TRUE.
             fac1 = (iso_val_inp-val1)/(val2-val1)
             fac2 = (val2-iso_val_inp)/(val2-val1)
           END IF
       ELSE
           IF(iso_val_inp<val1 .AND. iso_val_inp>val2) THEN
             found_val = .TRUE.
             fac2 = (iso_val_inp-val2)/(val1-val2)
             fac1 = (val1-iso_val_inp)/(val1-val2)
           END IF
       END IF
       IF(found_val) THEN
!dbg
! write(*,*) 'found val ', val1, val2, ' at k ', k, ' weights ', fac1, fac2  
!dbg
         IF(dir==3)THEN
           val1_e = cube_espot%grid(i,j,k) 
           val2_e = cube_espot%grid(i,j,k+1)
         ELSEIF(dir==2) THEN
           val1_e = cube_espot%grid(i,k,j) 
           val2_e = cube_espot%grid(i,k+1,j)
         ELSEIF(dir==1) THEN
           val1_e = cube_espot%grid(k,i,j) 
           val2_e = cube_espot%grid(k+1,i,j)
         END IF 
         iso_val = fac1*(val1*exp(alpha*sqrt(val1_e))) + &
                   fac2*(val2*exp(alpha*sqrt(val2_e)))
!dbg
!  write(*,*) 'iso val ', iso_val, ' with espot ', val1_e, val2_e
!  stop
!dbg
         EXIT
       END IF
     END DO
     IF(.NOT. found_val) THEN
!     iso_val = 
       stop ' reference value for correction not found '
     END IF
   ELSE
     iso_val = iso_val_inp
   END IF


   write(*,*) ' Start search for iso surface at value ', iso_val, ' in the height range between ', h_min, ' and ', h_max
   write(*,*) ' Number of considered planes ', index_max-index_min + 1

   IF(check_atom_from_down) THEN
     allocate(ah(cube%npoints(d1),cube%npoints(d2)))
     ah=h_max/au2a
!dbg
         DO iat = 1,cube%natom
           if(cube%coords(dir,iat)*au2a<h_max) THEN
    write(*,*) iat, cube%coords(dir,iat)*au2a, h_max
           endif
       enddo
! stop
!dbg
     DO i = 1,cube%npoints(d1)
       DO j = 1,cube%npoints(d2)
         r1 = real(i-1,KIND=dp)*cube%dh(1,d1)+cube%origin(1)+real(j-1,KIND=dp)*cube%dh(1,d2)
         r2 = real(i-1,KIND=dp)*cube%dh(2,d1)+real(j-1,KIND=dp)*cube%dh(2,d2)+cube%origin(2)
         DO iat = 1,cube%natom
           if(cube%coords(dir,iat)*au2a<h_max) THEN
             d12=sqrt((cube%coords(d1,iat)-r1)**2+(cube%coords(d2,iat)-r2)**2)
! write(*,*) i,j,iat,r1,r2,cube%coords(d1,iat),cube%coords(d2,iat),d12
             if(d12<2.D0 .AND. cube%coords(dir,iat)<ah(i,j)) THEN
               ah(i,j) = cube%coords(dir,iat)
!dbg
!   write(*,*) ' For grid point ', i, j, r1*au2a,r2*au2a, ' found an atom below hmax, atom ', iat, &
!              cube%coords(d1,iat)*au2a, cube%coords(d2,iat)*au2a, cube%coords(dir,iat)*au2a,&
!                 ', isosurf limit ', ah(i,j)*au2a 
!dbg
             end if
           end if
         END DO
       END DO
     END DO
   END IF

   jold = 1
   kold = index_min
   DO i = 1,cube%npoints(d1)
!   DO i = cube%npoints(d1),1,-1 !cube%npoints(d1)
     DO j = 1,cube%npoints(d2)
!     DO j = cube%npoints(d2),1,-1 !cube%npoints(d2)
       found_val = .FALSE. 
       DO k = i1,i2,incr
         IF(dir==3)THEN
           val1 = cube%grid(i,j,k)
           val2 = cube%grid(i,j,k+1)
         ELSEIF(dir==2)THEN
           val1 = cube%grid(i,k,j)
           val2 = cube%grid(i,k+1,j)
         ELSEIF(dir==1)THEN
           val1 = cube%grid(k,i,j)
           val2 = cube%grid(k+1,i,j)
         END IF
         IF(iso_current) THEN
           IF(dir==3)THEN
             val1_e = cube_espot%grid(i,j,k) 
             val2_e = cube_espot%grid(i,j,k+1)
           ELSEIF(dir==2) THEN
             val1_e = cube_espot%grid(i,k,j) 
             val2_e = cube_espot%grid(i,k+1,j)
           ELSEIF(dir==1) THEN
             val1_e = cube_espot%grid(k,i,j) 
             val2_e = cube_espot%grid(k+1,i,j)
           END IF 
           val1 = (val1*exp(alpha*sqrt(val1_e)))
           val2 = (val2*exp(alpha*sqrt(val2_e)))
         END IF

         IF(val1<val2) THEN 
             IF(iso_val<val2 .AND. iso_val>val1) THEN
               found_val = .TRUE.
               fac1 = (iso_val-val1)/(val2-val1)
               fac2 = (val2-iso_val)/(val2-val1)
             END IF
         ELSE
             IF(iso_val<val1 .AND. iso_val>val2) THEN
               found_val = .TRUE.
               fac2 = (iso_val-val2)/(val1-val2)
               fac1 = (val1-iso_val)/(val1-val2)
             END IF
         END IF
         IF(found_val) THEN
             r  = cube%origin(dir) + real(i-1,KIND=dp)*cube%dh(dir,d1)+real(j-1,KIND=dp)*cube%dh(dir,d2)&
               +real(k-1,KIND=dp)*cube%dh(dir,dir)*fac1+real(k,KIND=dp)*cube%dh(dir,dir)*fac2
             IF(check_atom_from_down) THEN
                 if(r>ah(i,j)) THEN
                   write(*,'(A,2F10.5,A,F16.8,A, F16.8)') ' At ', r1*au2a, r2*au2a, 'tip goes behind atoms, position moved from ', &
                                                        r*au2a, ' to ', (ah(i,j)-0.5D0)*au2a
!                   r=ah(i,j)-0.1D0
                   r=data_grid(3,iold,jold)
                   kval = NINT((ah(i,j)-2.2D0-cube%origin(dir))/cube%dh(dir,dir))+1
                   kval = kold
                 else
                   kval=k
                 endif
             ELSE
               kval = k
             END IF
             r1 = cube%origin(d1) + real(i-1,KIND=dp)*cube%dh(d1,d1) +real(j-1,KIND=dp)*cube%dh(d1,d2) &
                +real(k-1,KIND=dp)*cube%dh(d1,dir) *fac1+real(k,KIND=dp)*cube%dh(d1,dir)*fac2
             r2 = cube%origin(d2) + real(i-1,KIND=dp)*cube%dh(d2,d1) +real(j-1,KIND=dp)*cube%dh(d2,d2) &
                +real(kval-1,KIND=dp)*cube%dh(d2,dir) *fac1+real(kval,KIND=dp)*cube%dh(d2,dir)*fac2
             data_grid(1,i,j) = r1
             data_grid(2,i,j) = r2
             IF(use_espot) THEN
               IF(dir==3)THEN
                 data_grid(3,i,j) = (cube_espot%grid(i,j,kval-20)*fac1+cube_espot%grid(i,j,kval+1-20)*fac2)!/2.0D0
               ELSEIF(dir==2)THEN
                 data_grid(3,i,j) = (cube_espot%grid(i,kval,j)+cube_espot%grid(i,kval+1,j))/2.0D0
               ELSEIF(dir==1) THEN
                 data_grid(3,i,j) = (cube_espot%grid(kval,i,j)+cube_espot%grid(kval+1,i,j))/2.0D0
               END IF
             ELSE
               data_grid(3,i,j) = r
             END IF
            
             jold = j
             iold = i
             kold = kval
             EXIT
         END IF 
       END DO ! k
       IF(.NOT. found_val) THEN
          r1 = cube%origin(d1) + real(i-1,KIND=dp)*cube%dh(d1,d1) +real(j-1,KIND=dp)*cube%dh(d1,d2) &
               +real(kold-1,KIND=dp)*cube%dh(d1,dir)
          r2 = cube%origin(d2) + real(i-1,KIND=dp)*cube%dh(d2,d1) +real(j-1,KIND=dp)*cube%dh(d2,d2) &
               +real(kold-1,KIND=dp)*cube%dh(d2,dir)
          r = data_grid(3,iold,jold)
          data_grid(1,i,j) = r1
          data_grid(2,i,j) = r2
          data_grid(3,i,j) = r
          jold = j
          iold = i
         write(*,'(A,2F10.5,A,F16.8,A, F16.8)') ' At grid point ', r1*au2a, r2*au2a, ' no value close to ', iso_val, &
                    ' has been found; assigned : ',r
       END IF
     END DO ! j
   END DO  ! i


   OPEN(102,FILE="gird_data_index.dat")
     
   DO i = 1,cube%npoints(d1)
     DO j = 1,cube%npoints(d2)
       IF(use_espot) THEN
         write(102,*) i, j, data_grid(3,i,j)*au2eV
       ELSE
         write(102,*) i, j, data_grid(3,i,j)*au2a 
       END IF
     END DO
     write(102,*) " "
   END DO
   CLOSE(102)

   IF (use_espot .AND. iso_current) THEN
     WRITE(filename,'(A12,I1,A2,F10.8,A8)')  'es_on_iso_I_',dir,'_v',iso_val,'_gnu.dat'
   ELSE IF(use_espot) THEN
     WRITE(filename,'(A14,I1,A2,F10.8,A8)')  'es_on_iso_rho_',dir,'_v',iso_val,'_gnu.dat'
   ELSEIF(iso_current) THEN
     WRITE(filename,'(A8,I1,A2,F10.8,A8)')  'h_iso_I_',dir,'_v',iso_val,'_gnu.dat'
   ELSE
     WRITE(filename,'(A10,I1,A2,F10.8,A8)')  'h_iso_rho_',dir,'_v',iso_val,'_gnu.dat'
   END IF
   OPEN(103,FILE=TRIM(filename))

   L1 = cube%dh(d1,d1) * (cube%npoints(d1)-1) + cube%dh(d1,d2) * (cube%npoints(d2)-1)
   L2 = cube%dh(d2,d1) * (cube%npoints(d1)-1) + cube%dh(d2,d2) * (cube%npoints(d2)-1)

   delta = grid_delta

   np1 = INT(L1/delta) + 1
   np2 = INT(2.0D0*L2/delta) + 1

   det = (cube%dh(d1,d1)*cube%dh(d2,d2) - cube%dh(d1,d2)*cube%dh(d2,d1))
   det = 1.0D0/det
   DO ii = 1,np1
     r1 = delta*real(ii-1,KIND=dp)
     DO jj = 1,np2
       r2 = delta*real(jj-1,KIND=dp)
       tmp = (cube%dh(d2,d2)*r1-cube%dh(d1,d2)*r2)*det
       i = int(tmp)+1
       tmp = (cube%dh(d1,d1)*r2-cube%dh(d2,d1)*r1)*det
       j = int(tmp)+1
       IF(i<1) THEN
         r1_shift = r1+cube%npoints(d1)*cube%dh(d1,d1)
         r2_shift = r2+cube%npoints(d1)*cube%dh(d2,d1)
         tmp = (cube%dh(d2,d2)*r1_shift-cube%dh(d1,d2)*r2_shift)*det
         i = int(tmp)+1
!        tmp = (cube%dh(d1,d1)*r2_shift-cube%dh(d2,d1)*r1_shift)*det
!        j = int(tmp)+1
       ELSEIF(i>cube%npoints(d1)) THEN
         r1_shift = r1-cube%npoints(d1)*cube%dh(d1,d1)
         r2_shift = r2-cube%npoints(d1)*cube%dh(d2,d1)
         tmp = (cube%dh(d2,d2)*r1_shift-cube%dh(d1,d2)*r2_shift)*det
         i = int(tmp)+1
!        tmp = (cube%dh(d1,d1)*r2_shift-cube%dh(d2,d1)*r1_shift)*det
!        j = int(tmp)+1
       END IF

       IF(j<1) THEN
         r1_shift = r1+cube%npoints(d2)*cube%dh(d1,d2)
         r2_shift = r2+cube%npoints(d2)*cube%dh(d2,d2)
         tmp = (cube%dh(d1,d1)*r2_shift-cube%dh(d2,d1)*r1_shift)*det
         j = int(tmp)+1
       ELSEIF(j>cube%npoints(d2)) THEN
         r1_shift = r1-cube%npoints(d2)*cube%dh(d1,d2)
         r2_shift = r2-cube%npoints(d2)*cube%dh(d2,d2)
         tmp = (cube%dh(d1,d1)*r2_shift-cube%dh(d2,d1)*r1_shift)*det
         j = int(tmp)+1
       END IF

       im = i 
       ip = i
       jm = j 
       jp = j
       IF(i>1) THEN
           im = i-1
       ELSEIF(i==1) THEN
           im = cube%npoints(d1)
       END IF

       IF(j>1) THEN
           jm = j-1
       ELSE IF (j==1) THEN
           jm = cube%npoints(d2)
       END IF

       IF(i<cube%npoints(d1)) THEN
               ip = i+1
       ELSEIF(i==cube%npoints(d1)) THEN
             ip = 1
       END IF

       IF(j<cube%npoints(d2)) THEN
           jp = j+1
       ELSEIF(j==cube%npoints(d2))THEN
           jp = 1
       END IF
       height = (data_grid(3,i,j)+data_grid(3,im,j)+data_grid(3,ip,j)+&
                 data_grid(3,i,jm)+data_grid(3,i,jp))/5.D0
       IF(use_espot) THEN
         write(103,'(3f16.8)') (r1+cube%origin(d1))*au2a, (r2+cube%origin(d2))*au2a,height*au2eV
       ELSE
         write(103,'(3f16.8)') (r1+cube%origin(d1))*au2a, (r2+cube%origin(d2))*au2a,height*au2a 
       END IF
     END DO
     write(103,*) " "
   END DO
   CLOSE(103)

   DEALLOCATE(data_grid)

  END SUBROUTINE compute_iso_surf


END MODULE cube_post_process


PROGRAM main
   USE cubecruncher
   USE periodic_table
   USE command_line_tools
   USE cube_post_process
   IMPLICIT NONE

   TYPE(cube_type), POINTER :: cube,cube_subtract, cube_espot

   INTEGER :: iunit_cube_in,iunit_cube_out,iunit_xyz,narg
   INTEGER :: dir
   LOGICAL :: did_something
   REAL(KIND=dp) :: height, iso_val
   REAL(KIND=dp), DIMENSION(3) :: center
   TYPE(input_type) :: input
  
   CALL init_input(input)
   CALL parse_command_line(input,narg)
   did_something=.FALSE.
   !
   ! do simple checking of input and output
   !
   IF (input%write_version) THEN
      write(6,*) "This is cubecruncher version 1.1"
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
   
     iunit_cube_in=117
     iunit_cube_out=118
     iunit_xyz=119

     ! read cube
     write(6,FMT='(A)',ADVANCE="No") "Reading cube ... "
     OPEN(iunit_cube_in,FILE=TRIM(input%cube_name_in)) 
     CALL read_cube(cube,iunit_cube_in)
     CLOSE(iunit_cube_in)
     write(6,*) "Done"

     IF (input%do_subtract) THEN
        ! read cube
        ALLOCATE(cube_subtract)
        CALL init_cube(cube_subtract)
        write(6,FMT='(A)',ADVANCE="No") "Reading cube to subtract ... "
        OPEN(iunit_cube_in,FILE=TRIM(input%cube_name_subtract)) 
        CALL read_cube(cube_subtract,iunit_cube_in)
        CLOSE(iunit_cube_in)
        write(6,*) "Done"
        ! subtract from input
        CALL subtract_cubes(cube,cube_subtract)
        CALL deallocate_cube(cube_subtract)
     ENDIF

     IF (input%espot_over_iso .OR. input%iso_current) THEN
       IF(.NOT. input%have_espot_cube) THEN
         STOP " Name of ES-pot cube file needed "
       END IF
               ! read cube
        ALLOCATE(cube_espot)
        CALL init_cube(cube_espot)
        write(6,FMT='(A)',ADVANCE="No") "Reading cube with espot ... "
        OPEN(iunit_cube_in,FILE=TRIM(input%cube_name_espot))
        CALL read_cube(cube_espot,iunit_cube_in)
        CLOSE(iunit_cube_in)
        write(6,*) "Done"
     ELSE
        NULLIFY(cube_espot)
     END IF

     IF (input%do_xyz) THEN  
        ! read xyz
        write(6,FMT='(A)',ADVANCE="No") "Reading xyz ... "
        OPEN(iunit_xyz,FILE=TRIM(input%xyz_name))
        CALL read_xyz(cube,iunit_xyz)
        CLOSE(iunit_xyz)
        write(6,*) "Done"
     ENDIF

     CALL boxmat(cube)
  
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
        IF (input%do_center_point) THEN
           center=input%center_point
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

     IF (input%do_1d_profile) THEN
        write(6,FMT='(A, I4)') "Calculation of the profile along direction ",&
                 input%index_profile
        CALL  compute_profile(cube,input%index_profile,input%delta_profile)
        write(6,*) "Done"
     ENDIF

     IF (input%do_slice) THEN
        IF(input%slice(1)/=0.0D0) THEN
           dir = 1
           height = input%slice(1)
        ELSE IF(input%slice(2)/=0.0D0) THEN
           dir = 2
           height = input%slice(2)
        ELSE IF(input%slice(3)/=0.0D0) THEN
           dir = 3
           height = input%slice(3)
        END IF
        write(6,FMT='(A, I4, A, F10.4)') "Calculation of a cube slice perpendicular to direction",&
                dir, ' at height ', height 
        CALL  compute_slice(cube,dir,height)
        write(6,*) "Done"
     ENDIF

     IF(input%do_iso) THEN
        IF(input%iso_level(1)/=0.0D0) THEN
           dir = 1
           iso_val = input%iso_level(1)
        ELSE IF(input%iso_level(2)/=0.0D0) THEN
           dir = 2
           iso_val = input%iso_level(2)
        ELSE IF(input%iso_level(3)/=0.0D0) THEN
           dir = 3
           iso_val = input%iso_level(3)
        END IF
        write(6,FMT='(A, I4, A, 1PE20.4)') "Calculation of a iso surface, recording the height along direction",&
                dir, ' at level ', iso_val 
        IF(input%espot_over_iso .OR. input%iso_current ) THEN
            CALL compute_iso_surf(cube, dir, iso_val, input%iso_hmin ,input%iso_hmax, input%iso_delta_grid,&
                  input%espot_over_iso,  input%iso_current, input%from_top, input%check_atom_from_down, input%bwidth,  cube_espot)
        ELSE
            CALL compute_iso_surf(cube, dir, iso_val, input%iso_hmin ,input%iso_hmax, input%iso_delta_grid, &
                use_espot=.FALSE., iso_current=.FALSE., from_top=input%from_top, check_atom_from_down=input%check_atom_from_down)
        END IF
        write(6,*) "Done"
     END IF


     IF (ASSOCIATED(cube_espot)) THEN
        CALL deallocate_cube(cube_espot)
     END IF

     ! final out
     write(6,FMT='(A)',ADVANCE="No") "Writing cube ... "
     OPEN(iunit_cube_out,FILE=TRIM(input%cube_name_out)) 
     CALL write_cube(cube,iunit_cube_out,input%stride,input%mult_fact)
     CLOSE(iunit_cube_out)
     write(6,*) "Done"

     ! deallocate stuff
     CALL deallocate_cube(cube)

  ENDIF
END PROGRAM
