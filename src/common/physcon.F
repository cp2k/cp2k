!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Definition of physical constants:
!>
!>      a_bohr      : Bohr radius [m]
!>      a_fine      : Fine-structure constant
!>      a_mass      : Atomic mass unit [kg]; conversion factor [u] -> [kg]
!>      angstrom    : Conversion factor [Bohr] -> [Angstrom]
!>      bar         : Conversion factor [a.u.] -> [bar]
!>      bohr        : Conversion factor [Angstrom] -> [Bohr]
!>      boltzmann   : Boltzmann constant [J/K]
!>      c_light     : Speed of light in vacuum [m/s]
!>      e_charge    : Elementary charge [C]
!>      e_mass      : Electron mass [kg]
!>      e_gfactor   : Electron g factor [ ]
!>      esu         : Conversion factors [a.u.] -> [esu]
!>      evolt       : Conversion factor [a.u.] -> [eV]
!>      femtoseconds: Conversion factor [a.u.] -> [fs]
!>      h_bar       : Planck constant [J*s]
!>      h_planck    : Planck constant [J*s]
!>      hertz       : Conversion factor [a.u.] -> [Hz]
!>      joule       : Conversion factor [a.u.] -> [J]
!>      kcalmol     : Conversion factor [a.u.] -> [kcal/mol]
!>      kelvin      : Conversion factor [a.u.] -> [K]
!>      kjmol       : Conversion factor [a.u.] -> [kJ/mol]
!>      massunit    : Conversion factor [u] -> [a.u.]
!>      mu_perm     : Magnetic constant or permeability of vacuum [N/A**2]
!>      n_avogadro  : Avogadro constant [1/mol]
!>      newton      : Conversion factor [a.u.] -> [N]
!>      pascal      : Conversion factor [a.u.] -> [Pa]
!>      permittivity: Electric constant or permittivity of vacuum [F/m]
!>      picoseconds : Conversion factor [a.u.] -> [ps]
!>      rydberg     : Rydberg constant [1/m]
!>      seconds     : Conversion factor [a.u.] -> [s]
!>      vibfac      : Conversion factor [a.u./Bohr**2] -> [1/cm]
!>      wavenumbers : Conversion factor [a.u.] -> [1/cm]
!>      debye       : Conversion factor [a.u.] -> Debye
!> \note
!>      Fundamental physical constants (SI units)
!>      Literature: - P. J. Mohr and B. N. Taylor,
!>                    "CODATA recommended values of the fundamental physical
!>                     constants: 1998 Rev. Mod. Phys. 72, 351-495 (2000)
!>                  - P. J. Mohr and B. N. Taylor,
!>                    "CODATA recommended values of the fundamental physical
!>                     constants: 2002", Rev. Mod. Phys. 77, 1 (2005).
!>                  - P. J. Mohr, B. N. Taylor, and D. B. Newell,
!>                    "CODATA recommended values of the fundamental physical
!>                     constants: 2006 Rev. Mod. Phys. 80, 633 (2008)
!>                  - P. J. Mohr, B. N. Taylor, and D. B. Newell,
!>                    "CODATA recommended values of the fundamental physical
!>                     constants: 2010", Rev. Mod. Phys. 84, 1527-1605 (2012)
!> \par History
!>      - Adapted for use in CP2K (JGH)
!>      - Updated to CODATA 1998 and cleaned (05.09.2003,MK)
!>      - Updated to CODATA 2006. (26.03.2008,AK)
!>      - Updated to CODATA 2010. (10.12.2012,MK)
!> \author Matthias Krack (MK)
! *****************************************************************************
MODULE physcon

  USE f77_blas
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE periodic_table,                  ONLY: init_periodic_table

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: maxmom = 3

  REAL(KIND = dp) :: a_bohr,a_fine,a_mass,angstrom,atm,bar,bohr,boltzmann,&
                     c_light,debye,e_charge,e_gfactor,e_mass,evolt,&
                     femtoseconds,h_bar,h_planck,hertz,joule,kcalmol,kelvin,&
                     kjmol,massunit,mu_perm,n_avogadro,newton,p_mass,pascal,&
                     permittivity,picoseconds,rydberg,seconds,vibfac,wavenumbers

  REAL(KIND = dp), DIMENSION(maxmom) :: esu

  PUBLIC :: a_bohr, a_fine, a_mass, angstrom, atm, bar, bohr, boltzmann,&
            debye, e_charge, e_gfactor, e_mass, evolt, femtoseconds, h_bar,&
            hertz, joule, kcalmol, kelvin, kjmol, massunit, mu_perm, n_avogadro, newton,&
            p_mass, pascal, picoseconds, seconds, vibfac, wavenumbers

  PUBLIC :: init_physcon,&
            write_physcon

CONTAINS

! *****************************************************************************
!> \brief  Initialize all basic physical constants used by CP2K using
!>         the recommended fundamental constants of physics and chemistry
!>         based on the CODATA adjustments
!> \author  JGH
!> \date    14.11.2000
!> \par History
!>      - Updated to CODATA 1998 and cleaned (05.09.2003,MK)
!>      - Updated to CODATA 2006. (26.03.2008,AK)
!>      - Updated to CODATA 2010. (10.12.2012,MK)
!> \version 1.1
! *****************************************************************************
  SUBROUTINE init_physcon()

    INTEGER                                  :: i

! Exact constants
! Speed of light in vacuum [m/s]

    c_light = 299792458.0_dp

    ! Magnetic constant or permeability of vacuum [N/A**2]
    mu_perm = 4.0_dp*pi*1.0E-7_dp

    ! Electric constant or permittivity of vacuum [F/m]
    permittivity = 1.0_dp/(mu_perm*c_light**2)

#define CODATA_DEFAULT_VERSION 2006

#if !defined CODATA_VERSION
#define CODATA_VERSION CODATA_DEFAULT_VERSION
#endif

#if CODATA_VERSION != 1998 && CODATA_VERSION != 2006 && CODATA_VERSION != 2010
#define CODATA_VERSION CODATA_DEFAULT_VERSION
#endif

#if CODATA_VERSION == 1998
    ! Recommended fundamental constants of physics
    ! and chemistry based on the 1998 adjustment

    ! Planck constant [J*s]
    h_planck = 6.62606876E-34_dp
    h_bar = h_planck/(2.0_dp*pi)

    ! Elementary charge [C]
    e_charge = 1.602176462E-19_dp

    ! Electron mass [kg]
    e_mass = 9.10938188E-31_dp

    ! Proton mass [kg]
    p_mass = 1.67262158Eâ27_dp

    ! Electron g factor [ ]
    e_gfactor = -2.0023193043718_dp

    ! Fine-structure constant
!MK a_fine = 0.5_dp*mu_perm*c_light*e_charge**2/h_planck
    a_fine = 7.297352533E-3_dp

    ! Rydberg constant [1/m]
!MK rydberg = 0.5_dp*e_mass*c_light*a_fine**2/h_planck
    rydberg = 10973731.568549_dp

    ! Avogadro constant [1/mol]
    n_avogadro = 6.02214199E+23_dp

    ! Boltzmann constant [J/K]
    boltzmann = 1.3806503E-23_dp

    ! Atomic mass unit [kg]; conversion factor [u] -> [kg]
    a_mass = 1.66053873E-27_dp

    ! Bohr radius [m]
!MK a_bohr = a_fine/(4.0_dp*pi*rydberg)
    a_bohr = 0.5291772083E-10_dp
#endif

#if CODATA_VERSION == 2006
    ! Recommended fundamental constants of physics
    ! and chemistry based on the 2006 adjustment

    ! Planck constant [J*s]
    h_planck = 6.62606896E-34_dp
    h_bar = h_planck/(2.0_dp*pi)

    ! Elementary charge [C]
    e_charge = 1.602176487E-19_dp

    ! Electron mass [kg]
    e_mass = 9.10938215E-31_dp

    ! Proton mass [kg]
    p_mass = 1.672621637E-27_dp

    ! Electron g factor [ ]
    e_gfactor = -2.0023193043622_dp

    ! Fine-structure constant
!MK a_fine = 0.5_dp*mu_perm*c_light*e_charge**2/h_planck
    a_fine = 7.2973525376E-3_dp

    ! Rydberg constant [1/m]
!MK rydberg = 0.5_dp*e_mass*c_light*a_fine**2/h_planck
    rydberg = 10973731.568527_dp

    ! Avogadro constant [1/mol]
    n_avogadro = 6.02214179E+23_dp

    ! Boltzmann constant [J/K]
    boltzmann = 1.3806504E-23_dp

    ! Atomic mass unit [kg]; conversion factor [u] -> [kg]
    a_mass = 1.660538782E-27_dp

    ! Bohr radius [m]
!MK a_bohr = a_fine/(4.0_dp*pi*rydberg)
    a_bohr = 0.52917720859E-10_dp
#endif

#if CODATA_VERSION == 2010
    ! Recommended fundamental constants of physics
    ! and chemistry based on the 2010 adjustment

    ! Planck constant [J*s]
    h_planck = 6.62606957E-34_dp
    h_bar = h_planck/(2.0_dp*pi)

    ! Elementary charge [C]
    e_charge = 1.602176565E-19_dp

    ! Electron mass [kg]
    e_mass = 9.10938291E-31_dp

    ! Proton mass [kg]
    p_mass = 1.672621777E-27_dp

    ! Electron g factor [ ]
    e_gfactor = -2.00231930436153_dp

    ! Fine-structure constant
!MK a_fine = 0.5_dp*mu_perm*c_light*e_charge**2/h_planck
    a_fine = 7.2973525698E-3_dp

    ! Rydberg constant [1/m]
!MK rydberg = 0.5_dp*e_mass*c_light*a_fine**2/h_planck
    rydberg = 10973731.568539_dp

    ! Avogadro constant [1/mol]
    n_avogadro = 6.02214129E+23_dp

    ! Boltzmann constant [J/K]
    boltzmann = 1.3806488E-23_dp

    ! Atomic mass unit [kg]; conversion factor [u] -> [kg]
    a_mass = 1.660538921E-27_dp

    ! Bohr radius [m]
!MK a_bohr = a_fine/(4.0_dp*pi*rydberg)
    a_bohr = 0.52917721092E-10_dp
#endif

    ! Conversion factors

    ! [u] -> [a.u.]
    massunit = a_mass/e_mass

    ! [Bohr] -> [Angstrom]
    angstrom = 1.0E+10_dp*a_bohr

    ! [Angstrom] -> [Bohr]
    bohr = 1.0_dp/angstrom

    ! [a.u.] -> [s]
    seconds = 1.0_dp/(4.0_dp*pi*rydberg*c_light)

    ! [a.u.] -> [fs]
    femtoseconds = 1.0E+15_dp*seconds

    ! [a.u.] -> [ps]
    picoseconds = 1.0E+12_dp*seconds

    ! [a.u.] -> [J]
    joule = 2.0_dp*rydberg*h_planck*c_light

    ! [a.u.] -> [N]
    newton = joule/a_bohr

    ! [a.u.] -> [K]
    kelvin = joule/boltzmann

    ! [a.u.] -> [kJ/mol]
    kjmol = 0.001_dp*joule*n_avogadro

    ! [a.u.] -> [kcal/mol]
    kcalmol = kjmol/4.184_dp

    ! [a.u.] -> [Pa]
    pascal = joule/a_bohr**3

    ! [a.u.] -> [bar]
    bar = pascal/1.0E+5_dp

    ! [a.u.] -> [atm]
    atm = pascal/1.013250E+5_dp

    ! [a.u.] -> [eV]
    evolt = joule/e_charge

    ! [a.u.] -> [Hz]
    hertz = joule/h_planck

    ! [a.u./Bohr**2] -> [1/cm] (wave numbers)
    vibfac = 5.0_dp*SQRT(kjmol)/(pi*a_bohr*c_light)

    ! [a.u.] -> [1/cm] (wave numbers)
    wavenumbers = 0.02_dp*rydberg

    ! [a.u.] -> [esu] (electrostatic units)
    esu(1) = 1.0E+21_dp*a_bohr*c_light*e_charge
    DO i=2,maxmom
      esu(i) = esu(i-1)/bohr
    END DO

    ! [a.u.] -> [debye] (electrostatic units)
    Debye = esu(1)

    ! Initialize the periodic table
    CALL init_periodic_table()

  END SUBROUTINE init_physcon

! *****************************************************************************
!> \brief  Write all basic physical constants used by CP2K to a logical
!>           output unit.
!> \author  JGH
!> \date    14.11.2000
!> \par History
!>       - Updated to CODATA 1998 and cleaned (05.09.2003,MK)
!>       - Updated to CODATA 2006. (26.03.2008,AK)
!> \version 1.1
! *****************************************************************************
  SUBROUTINE write_physcon(output_unit)

    INTEGER, INTENT(IN) :: output_unit

#if CODATA_VERSION == 1998
    WRITE (UNIT=output_unit,FMT="(T2,/,T2,A,/,/,(T2,A))")&
      "*** Fundamental physical constants (SI units) ***",&
      "*** Literature: B. J. Mohr and B. N. Taylor,",&
      "***             CODATA recommended values of the fundamental physical",&
      "***             constants: 1998, Rev. Mod. Phys. 72(2), 351 (2000)"
#endif

#if CODATA_VERSION == 2006
    WRITE (UNIT=output_unit,FMT="(T2,/,T2,A,/,/,(T2,A))")&
      "*** Fundamental physical constants (SI units) ***",&
      "*** Literature: B. J. Mohr and B. N. Taylor,",&
      "***             CODATA recommended values of the fundamental physical",&
      "***             constants: 2006, Web Version 5.1",&
      "***             http://physics.nist.gov/constants"
#endif

#if CODATA_VERSION == 2010
    WRITE (UNIT=output_unit,FMT="(T2,/,T2,A,/,/,(T2,A))")&
      "*** Fundamental physical constants (SI units) ***",&
      "*** Literature: P. J. Mohr, B. N. Taylor, and D. B. Newell,",&
      "***             CODATA recommended values of the fundamental physical",&
      "***             constants: 2010, Rev. Mod. Phys. 84, 1527-1605 (2012)"
#endif

    WRITE (UNIT=output_unit,FMT="(/,T2,A,T61,ES20.14)")&
      "Speed of light in vacuum [m/s]",c_light
    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "Magnetic constant or permeability of vacuum [N/A**2]",mu_perm
    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "Electric constant or permittivity of vacuum [F/m]",permittivity
    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "Planck constant (h) [J*s]",h_planck
    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "Planck constant (h-bar) [J*s]",h_bar
    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "Elementary charge [C]",e_charge
    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "Electron mass [kg]",e_mass
    WRITE (UNIT=output_unit,FMT="(T2,A,T60,ES21.14)")&
      "Electron g factor [ ]",e_gfactor
    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "Proton mass [kg]",p_mass
    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "Fine-structure constant",a_fine
    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "Rydberg constant [1/m]",rydberg
    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "Avogadro constant [1/mol]",n_avogadro
    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "Boltzmann constant [J/K]",boltzmann
    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "Atomic mass unit [kg]",a_mass
    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "Bohr radius [m]",a_bohr

    ! Conversion factors

    WRITE (UNIT=output_unit,FMT="(/,T2,A,/)")&
      "*** Conversion factors ***"

    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "[u] -> [a.u.]",massunit,&
      "[Angstrom] -> [Bohr] = [a.u.]",bohr,&
      "[a.u.] = [Bohr] -> [Angstrom]",angstrom,&
      "[a.u.] -> [s]",seconds,&
      "[a.u.] -> [fs]",femtoseconds,&
      "[a.u.] -> [J]",joule,&
      "[a.u.] -> [N]",newton,&
      "[a.u.] -> [K]",kelvin,&
      "[a.u.] -> [kJ/mol]",kjmol,&
      "[a.u.] -> [kcal/mol]",kcalmol,&
      "[a.u.] -> [Pa]",pascal,&
      "[a.u.] -> [bar]",bar,&
      "[a.u.] -> [atm]",atm,&
      "[a.u.] -> [eV]",evolt,&
      "[a.u.] -> [Hz]",hertz,&
      "[a.u.] -> [1/cm] (wave numbers)",wavenumbers,&
      "[a.u./Bohr**2] -> [1/cm]",vibfac
    WRITE (UNIT=output_unit,FMT="(T2,A)") ""

  END SUBROUTINE write_physcon

END MODULE physcon
