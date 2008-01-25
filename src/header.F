!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
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
  PUBLIC :: fist_header, se_header, dftb_header, qs_header, tddfpt_header
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
         '                       2005 -- 2008                      ', &
         '                                                         ', &
         '                       Version 0.2                       ', &
         '                                                         '

  END SUBROUTINE se_header

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
         "**            Copyright (C) by CP2K Developers Group (2000 - 2008)           **",&
         "**                                                                           **",&
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

