!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 -2006 CP2K developers group                            !
!-----------------------------------------------------------------------------!
!!****** cp2k/header [1.0] *
!!
!!   NAME
!!     header
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     APSI & CJM & JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE header

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: fist_header, tbmd_header, qs_header, tddfpt_header
  PUBLIC :: kg_header, qmmm_header

  !!*****
  !******************************************************************************

CONTAINS

  !******************************************************************************
  SUBROUTINE qmmm_header (iw, ionode)
    INTEGER, INTENT(in)                      :: iw
    LOGICAL, INTENT(in)                      :: ionode

    IF (.NOT.ionode) RETURN
    WRITE ( iw, '( / )' )
    WRITE ( iw, '( 21(3x,a,/) )' ) &
         "  ___                    _                     __  __           _           ",&
         " / _ \ _   _  __ _ _ __ | |_ _   _ _ __ ___   |  \/  | ___  ___| |__   ___  ",&
         "| | | | | | |/ _` | '_ \| __| | | | '_ ` _ \  | |\/| |/ _ \/ __| '_ \ / __| ",&
         "| |_| | |_| | (_| | | | | |_| |_| | | | | | | | |  | |  __/ (__| | | |\__ \ ",&
         " \__\_\\__,_|\__,_|_| |_|\__|\__,_|_| |_| |_| |_|  |_|\___|\___|_| |_||___/o",&
         "                                                                            ",&
         " __  __       _                 _              __  __           _           ",&
         "|  \/  | ___ | | ___  ___ _   _| | __ _ _ __  |  \/  | ___  ___| |__   ___  ",&
         "| |\/| |/ _ \| |/ _ \/ __| | | | |/ _` | '__| | |\/| |/ _ \/ __| '_ \ / __| ",&
         "| |  | | (_) | |  __/ (__| |_| | | (_| | |    | |  | |  __/ (__| | | |\__ \ ",&
         "|_|  |_|\___/|_|\___|\___|\__,_|_|\__,_|_|    |_|  |_|\___|\___|_| |_||___/o",&
         "                                                                            ",&
         "----------------------------------------------------------------------------",&
         "| Please cite the following reference(s) if you publish or present results |",&
         "| using the QM/MM module of CP2K.                                          |",&
         "|--------------------------------------------------------------------------|",&
         "|                                                                          |",&
         "|   T.Laino, F.Mohamed, A.Laio, M. Parrinello                              |",&
         "|   An Efficient Real Space Multigrid QM/MM Electrostatic Coupling         |",&
         "|   J. Chem. Th. Comp.  1(6) (2005) pp.1176-1184                           |",&
         "----------------------------------------------------------------------------"
         

  END SUBROUTINE qmmm_header
  
  SUBROUTINE tddfpt_header(iw, ionode)
    INTEGER, INTENT(in)                      :: iw
    LOGICAL, INTENT(in)                      :: ionode

    IF (.NOT.ionode) RETURN       
    WRITE (unit=iw, fmt="(7(2X,A,/))")&
         " _____ _                  ____                            _            _   ",&
         "|_   _(_)_ __ ___   ___  |  _ \  ___ _ __   ___ _ __   __| | ___ _ __ | |_ ",&
         "  | | | | '_ ` _ \ / _ \ | | | |/ _ \ '_ \ / _ \ '_ \ / _` |/ _ \ '_ \| __|",&
         "  | | | | | | | | |  __/ | |_| |  __/ |_) |  __/ | | | (_| |  __/ | | | |_ ",&
         "  |_| |_|_| |_| |_|\___| |____/ \___| .__/ \___|_| |_|\__,_|\___|_| |_|\__|",&
         "                                    |_|                                    ",&
         "                                                       CALCULATION STARTED!"
  END SUBROUTINE tddfpt_header

  SUBROUTINE fist_header ( iw, ionode )
    INTEGER, INTENT(in)                      :: iw
    LOGICAL, INTENT(in)                      :: ionode

    IF (.NOT.ionode) RETURN
    WRITE ( iw, '( / )' )
    WRITE ( iw, '( 14(20x,a,/) )' ) &
         '  ************* ******** ***********   ', &
         ' ************* **********************  ', &
         ' ****     **** ****      ***********   ', &
         ' ******** ****  *******     ****       ', &
         ' *******  ****   *******   ****        ', &
         ' ****     ****       ****  ****        ', &
         ' ****     **** **********  ****        ', &
         '  ****     **** ********    ****       ', &
         ' FRONTIERS IN SIMULATION TECHNOLOGY    ', &
         '                                       ', &
         ' C.J. Mundy, S. Balasubramanian,       ', &
         ' Ken Bagchi, J. Hutter, Ari Seitsonen  ', &
         '             Version 0.3               ', &
         '                                                       '

  END SUBROUTINE fist_header

  SUBROUTINE tbmd_header ( iw, ionode )
    INTEGER, INTENT(in)                      :: iw
    LOGICAL, INTENT(in)                      :: ionode

    IF (.NOT.ionode) RETURN

    WRITE ( iw, '( / )' )
    WRITE ( iw, '( 14(12x,a,/) )' ) &
         '  **********************   *************   ********    ', &
         ' ************************ *************** **********   ', &
         '  ***************     *** *************** ***********  ', &
         '     ****    ****    **** ****       **** ***     **** ', &
         '    ****     ***********  ****       **** ***     **** ', &
         '    ****     ****     *** ****       **** ***********  ', &
         '    ****     ***********  ****       **** **********   ', &
         '     ****     *********    ****      ****  *******     ', &
         '                                                       ', &
         '                 University of Zurich                  ', &
         '                        2000                           ', &
         '                                                       ', &
         '                     Version 0.0                       ', &
         '                                                       '

  END SUBROUTINE tbmd_header

  SUBROUTINE qs_header(output_unit, ionode)

    !   Purpose: Write the Quickstep program banner with the copyright and the
    !            program version information to the output unit.
    !   JVdV it is also a great place to provide references to QS related papers

    !   History: - Creation (15.06.2000, Matthias Krack)

    !   ***************************************************************************

    INTEGER                                  :: output_unit
    LOGICAL                                  :: ionode

    IF (.NOT.ionode) RETURN
    WRITE (UNIT=output_unit,FMT="(/,(T2,A79))")&
         "*******************************************************************************",&
         "*******************************************************************************",&
         "**                                                                           **",&
         "**     #####                         ##              ##                      **",&
         "**    ##   ##            ##          ##              ##                      **",&
         "**   ##     ##                       ##            ######                    **",&
         "**   ##     ##  ##   ##  ##   #####  ##  ##   ####   ##    #####    #####    **",&
         "**   ##     ##  ##   ##  ##  ##      ## ##   ##      ##   ##   ##  ##   ##   **",&
         "**   ##  ## ##  ##   ##  ##  ##      ####     ###    ##   ######   ######    **",&
         "**    ##  ###   ##   ##  ##  ##      ## ##      ##   ##   ##       ##        **",&
         "**     #######   #####   ##   #####  ##  ##  ####    ##    #####   ##        **",&
         "**           ##                                                    ##        **",&
         "**                                                                           **",&
         "**                                                ... make the atoms dance   **",&
         "**                                                                           **",&
         "**            Copyright (C) by CP2K Developers Group (2000 - 2006)           **",&
         "**                                                                           **"
         !
         ! add QS related method papers here. Also update the webpages.
         !
    WRITE (UNIT=output_unit,FMT="(/,(T2,A79))")&
         "*******************************************************************************",&
         "** G. Lippert, J. Hutter, and M. Parrinello,                                 **",&
         "** A hybrid Gaussian and plane wave density functional scheme,               **",&
         "** Molecular Physics, Vol. 92, Issue 3, 477-487 (1997)                       **",&
         "** http://dx.doi.org/10.1080/002689797170220                                 **",&
         "*******************************************************************************",&
         "** G. Lippert, J. Hutter, and M. Parrinello                                  **",&
         "** The Gaussian and augmented-plane-wave density functional method for       **",&
         "** ab initio molecular dynamics simulations                                  **",&
         "** Theoretical Chemistry Accounts, Vol. 103, Issue 2, 124-140 (1999)         **",&
         "** http://dx.doi.org/10.1007/s002140050523                                   **",&
         "*******************************************************************************",&
         "** M. Krack and M. Parrinello                                                **",&
         "** All-electron ab-initio molecular dynamics                                 **",&
         "** Physical Chemistry Chemical Physics, Vol. 2, Issue 10, 2105-2112 (2000)   **",&
         "** http://dx.doi.org/10.1039/b001167n                                        **",&
         "*******************************************************************************",&
         "** J. VandeVondele and J. Hutter                                             **",&
         "** An efficient orbital transformation method for                            **",&
         "** electronic structure calculations.                                        **",&
         "** J. Chem. Phys. 118, 4365 (2003)                                           **",&
         "** http://dx.doi.org/10.1063/1.1543154                                       **",&
         "*******************************************************************************",&
         "** J. VandeVondele, M. Krack, F. Mohamed,                                    **",&
         "** M. Parrinello, T. Chassaing and J. Hutter                                 **",&
         "** Quickstep: fast and accurate density functional calculations using        **",&
         "** a mixed Gaussian and plane waves approach.                                **",&
         "** Comp. Phys. Comm. 167, 103 (2005).                                        **",&
         "** http://dx.doi.org/10.1016/j.cpc.2004.12.014                               **",&
         "*******************************************************************************"

  END SUBROUTINE qs_header


  SUBROUTINE kg_header ( iw, ionode )
    INTEGER                                  :: iw
    LOGICAL, INTENT(in)                      :: ionode

    IF (.NOT.ionode) RETURN

    WRITE ( iw, '( / )' )
    WRITE ( iw, '( 14(12x,a,/) )' ) &
         '  ****  ****   *******                                 ', &
         ' ****  ****   **********                               ', &
         ' **** ****   *****                                     ', &
         ' ********    ****  ****   ***   ***  ****   **** **    ', &
         ' ********    ****  ***** ***** ** ** ** ** **    **    ', &
         ' **** ****   ****   **** ** ** ** ** ** ** ***** **    ', &
         ' ****  ****   *********  ** ** ** ** ** ** **    **    ', &
         '  ****  ****   *******   ** **  ***  ****   **** ***** ', &
         '                                                       ', &
         '                     Version 0.01                      ', &
         '                                                       '

  END SUBROUTINE kg_header


END MODULE header

