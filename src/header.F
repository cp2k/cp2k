!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      none
!> \author APSI & CJM & JGH
! *****************************************************************************
MODULE header
  USE f77_blas

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: fist_header, se_header, dftb_header, qs_header, tddfpt_header
  PUBLIC :: kg_header, qmmm_header, atom_header, atom_footer

CONTAINS

! *****************************************************************************
  SUBROUTINE qmmm_header (iw, ionode)
    INTEGER, INTENT(in)                      :: iw
    LOGICAL, INTENT(in)                      :: ionode

    IF (.NOT.ionode) RETURN
    WRITE ( iw, '( / )' )
    WRITE ( iw, '( 13(3x,a,/) )' ) &
         "----------------------------------------------------------------------------",&
         "|                                                                          |",&
         "| ,ad8888ba,    88b           d88       88b           d88  88b           d88",&
         " d8""'    `""8b   888b         d888       888b         d888  888b         d888",&
         "d8'        `8b  88`8b       d8'88       88`8b       d8'88  88`8b       d8'88",& 
         "88          88  88 `8b     d8' 88       88 `8b     d8' 88  88 `8b     d8' 88",& 
         "88          88  88  `8b   d8'  88       88  `8b   d8'  88  88  `8b   d8'  88",& 
         "Y8,    ""88,,8P  88   `8b d8'   88       88   `8b d8'   88  88   `8b d8'   88",&
         " Y8a.    Y88P   88    `888'    88  888  88    `888'    88  88    `888'    88",& 
         "| `""Y8888Y""Y8a  88     `8'     88  888  88     `8'     88  88     `8'     88",&
         "|                                                                          |",&
         "|                                                     CALCULATION STARTED..|",&
         "----------------------------------------------------------------------------"
         
  END SUBROUTINE qmmm_header
  
! *****************************************************************************
  SUBROUTINE tddfpt_header(iw, ionode)
    INTEGER, INTENT(in)                      :: iw
    LOGICAL, INTENT(in)                      :: ionode

    IF (.NOT.ionode) RETURN       
    WRITE (unit=iw, fmt="(13(2X,A,/))")&
         "----------------------------------------------------------------------------",&
         "                                                                            ",&
         "888888888888  88888888ba,          88888888ba,    88888888888  888888888888",&
         "     88       88      `""8b         88      `""8b   88                88     ",&
         "     88       88        `8b        88        `8b  88                88     ",&
         "     88       88         88        88         88  88aaaaa           88     ",&
         "     88       88         88        88         88  88""""""""""           88     ",&
         "     88       88         8P        88         8P  88                88     ",&
         "     88       88      .a8P         88      .a8P   88                88     ",&
         "     88       88888888Y""'          88888888Y""'    88                88     ",& 
         "                                                                            ",&
         "                                                       CALCULATION STARTED!",&
         "----------------------------------------------------------------------------"

  END SUBROUTINE tddfpt_header

! *****************************************************************************
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
         ' IFW Kuo, T. Laino, J. VandeVondele    ', &
         '             Version 0.31415           ', &
         '                                       '

  END SUBROUTINE fist_header

! *****************************************************************************
  SUBROUTINE dftb_header ( iw, ionode )
    INTEGER, INTENT(in)                      :: iw
    LOGICAL, INTENT(in)                      :: ionode

    IF (.NOT.ionode) RETURN

    WRITE ( iw, '( / )' )
    WRITE ( iw, '( 11(12x,a,/) )' ) &
         '  #####   #####        # ######  ####### ####### ######   ',&
         ' #     # #     #      #  #     # #          #    #     #  ',&
         ' #     # #           #   #     # #          #    #     #  ',&
         ' #     #  #####     #    #     # #####      #    ######   ',&
         ' #   # #       #   #     #     # #          #    #     #  ',&
         ' #    #  #     #  #      #     # #          #    #     #  ',&
         '  #### #  #####  #       ######  #          #    ######   ',&
         '                                                       ', &
         '                     Version 0.1                       ', &
         '                J. Hutter and T. Heine                 ', &
         '                                                       '
  END SUBROUTINE dftb_header

! *****************************************************************************
  SUBROUTINE se_header ( iw, ionode )
    INTEGER, INTENT(in)                      :: iw
    LOGICAL, INTENT(in)                      :: ionode

    IF (.NOT.ionode) RETURN

    WRITE ( iw, '( / )' )
    WRITE ( iw, '( 14(12x,a,/) )' ) &
         '  #####   #####        # #     # ######  ######   #####  ', &
         ' #     # #     #      #  ##    # #     # #     # #     # ', &
         ' #     # #           #   # #   # #     # #     # #     # ', &
         ' #     #  #####     #    #  #  # #     # #     # #     # ', &
         ' #   # #       #   #     #   # # #     # #     # #     # ', &
         ' #    #  #     #  #      #    ## #     # #     # #     # ', &
         '  #### #  #####  #       #     # ######  ######   #####  ', &
         '                                                         ', &
         '                   University of Zurich                  ', &
         '                       2005 -- 2009                      ', &
         '                                                         ', &
         '                       Version 1.0                       ', &
         '                                                         '

  END SUBROUTINE se_header

! *****************************************************************************
  SUBROUTINE qs_header(output_unit, ionode)
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
         "**            Copyright (C) by CP2K Developers Group (2000 - 2009)           **",&
         "**                                                                           **",&
         "*******************************************************************************"
  END SUBROUTINE qs_header

! *****************************************************************************
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

  SUBROUTINE atom_header ( iw )
    INTEGER                                  :: iw

    WRITE ( iw, '( / )' )
    WRITE ( iw, '( 14(25x,a,/) )' ) &
         '  ****  ******  ****   ****   ', &
         ' **  ** ****** **  ** ******  ', &
         ' ******   **   **  ** **  **  ', &
         ' **  **   **    ****  **  **  ', &
         '                              ', &
         '    University of Zurich      ', &
         '           2009               ', &
         '                              ', &
         '        Version 0.0           ', &
         '                                          '
  END SUBROUTINE atom_header
  
  SUBROUTINE atom_footer ( iw )
    INTEGER                                  :: iw

    WRITE ( iw, '( / )' )
    WRITE ( iw, '( 14(25x,a,/) )' ) &
         '    NORMAL TERMINATION OF     ', &
         '                              ', &
         '  ****  ******  ****   ****   ', &
         ' **  ** ****** **  ** ******  ', &
         ' ******   **   **  ** **  **  ', &
         ' **  **   **    ****  **  **  '

  END SUBROUTINE atom_footer

END MODULE header

