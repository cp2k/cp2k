!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_arpack [1.0] *
!!
!!   NAME
!!     qs_arpack
!!
!!   FUNCTION
!!     provide an interface to (p)arpack
!!
!!   AUTHOR
!!     Joost VandeVondele (05.2002)
!!
!!   NOTES
!!     - This module needs the MPI version of PARPACK to linked 
!!     - Notice that the routines PXLAMCH conflict with the scalapack routines, 
!!       and have to be renamed througout the (P)ARPACK code (see the 
!!       tools directory ....). The same is true for the version of SECOND used
!!       hence any downloaded binary of (P)ARPACK is not likely to work
!!     - notice that the code will only work if __ARPACK is defined
!!       otherwise produces a more or less empty body in this way, people
!!       can compile more easily the code for other parts of cp2k without
!!       going through the effort of compiling a library for themselves
!!
!!   MODIFICATION HISTORY
!!     None
!!
!!   SEE ALSO
!!     http://www.netlib.org/arpack/
!!     http://www.caam.rice.edu/software/ARPACK/index.html
!!
!!****

module qs_arpack

  use qs_arpack_diag_sy,   only : arpack_diag_sy
  use qs_arpack_utils,     only : arpack_confirm
  
  implicit none

  private
  public  :: arpack_diag_sy, arpack_confirm

end module qs_arpack
