!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001 - 2003  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/physcon [1.0] *
!!
!!   NAME
!!     physcon
!!
!!   FUNCTION
!!     Definition of physical constants:
!!
!!     a_bohr      : Bohr radius [m]
!!     a_fine      : Fine-structure constant
!!     a_mass      : Atomic mass unit [kg]; conversion factor [u] -> [kg]
!!     angstrom    : Conversion factor [Bohr] -> [Angstrom]
!!     bohr        : Conversion factor [Angstrom] -> [Bohr]
!!     boltzmann   : Boltzmann constant [J/K]
!!     c_light     : Speed of light in vacuum [m/s]
!!     e_charge    : Elementary charge [C]
!!     e_mass      : Electron mass [kg]
!!     esu         : Conversion factors [a.u.] -> [esu]
!!     evolt       : Conversion factor [a.u.] -> [eV]
!!     femtoseconds: Conversion factor [a.u.] -> [fs]
!!     h_bar       : Planck constant [J*s]
!!     h_planck    : Planck constant [J*s]
!!     hz          : Conversion factor [a.u.] -> [Hz]
!!     joule       : Conversion factor [a.u.] -> [J]
!!     kcalmol     : Conversion factor [a.u.] -> [kcal/mol]
!!     kelvin      : Conversion factor [a.u.] -> [K]
!!     kjmol       : Conversion factor [a.u.] -> [kJ/mol]
!!     massunit    : Conversion factor [u] -> [a.u.]
!!     mhz         : Conversion factor [a.u.] -> [MHz]
!!     mu_perm     : Magnetic constant or permeability of vacuum [N/A**2]
!!     n_avogadro  : Avogrado constant [1/mol]
!!     permittivity: Electric constant or permittivity of vacuum [F/m]
!!     rydberg     : Rydberg constant [1/m]
!!     seconds     : Conversion factor [a.u.] -> [s]
!!     vibfac      : Conversion factor [a.u./Bohr**2] -> [1/cm]
!!     wavenum     : Conversion factor [a.u.] -> [1/cm] (wave numbers)
!!     debye       : Conversion factor [a.u.] -> Debye
!!
!!   AUTHOR
!!     Matthias Krack
!!
!!   MODIFICATION HISTORY
!!     - Adapted for use in CP2K (JGH)
!!     - Updated to CODATA 1998 and cleaned (05.09.2003,MK)
!!
!!   NOTES
!!     Fundamental physical constants (SI units)
!!     Literature: - B. J. Mohr and B. N. Taylor,
!!                   "CODATA recommended values of the fundamental physical
!!                    constants: 1998", Rev. Mod. Phys. 72(2), 351 (2000).
!!
!!   SOURCE
!******************************************************************************

MODULE physcon

  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE periodic_table,                  ONLY: init_periodic_table

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  INTEGER, PARAMETER :: maxmom = 3

  REAL(KIND = dp) :: a_bohr,a_fine,a_mass,angstrom,bohr,boltzmann,c_light,e_charge,&
              e_mass,evolt,femtoseconds,h_bar,h_planck,hz,joule,kcalmol,kelvin,&
              kjmol,massunit,mu_perm,mhz,n_avogadro,permittivity,rydberg,&
              seconds,vibfac,wavenum, debye

  REAL(KIND = dp), DIMENSION(maxmom) :: esu

! *** Public variables ***

  PUBLIC :: a_bohr,a_fine,a_mass,angstrom,bohr,boltzmann,c_light,e_charge,&
            e_mass,evolt,femtoseconds,h_bar,h_planck,hz,joule,kcalmol,kelvin,&
            kjmol,massunit,mu_perm,mhz,n_avogadro,permittivity,rydberg,&
            seconds,vibfac,wavenum, debye

! *** Public subroutines ***

  PUBLIC :: init_physcon,&
            write_physcon

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE init_physcon()

!   Purpose: Initialize all basic physical constants used by CP2K using
!            the recommended fundamental constants of physics and chemistry
!            based on the 1998 adjustment (CODATA 1998)
                                                        
!   History: - Creation (14.11.2000,JGH)
!            - Updated to CODATA 1998 and cleaned (05.09.2003,MK)

!   ***************************************************************************


    INTEGER                                  :: i

!   ---------------------------------------------------------------------------
!   *** Exact constants ***
!   *** Speed of light in vacuum [m/s] ***

    c_light = 299792458.0_dp

!   *** Magnetic constant or permeability of vacuum [N/A**2] ***

    mu_perm = 4.0_dp*pi*1.0E-7_dp

!   *** Electric constant or permittivity of vacuum [F/m] ***

    permittivity = 1.0_dp/(mu_perm*c_light**2)

!   ---------------------------------------------------------------------------

!   *** Recommended fundamental constants of physics ***
!   *** and chemistry based on the 1998 adjustment   ***

!   *** Planck constant [J*s] ***

    h_planck = 6.62606876E-34_dp
    h_bar = h_planck/(2.0_dp*pi)

!   *** Elementary charge [C] ***

    e_charge = 1.602176462E-19_dp

!   *** Electron mass [kg] ***

    e_mass = 9.10938188E-31_dp

!   *** Fine-structure constant ***

!MK  a_fine = 0.5_dp*mu_perm*c_light*e_charge**2/h_planck
     a_fine = 7.297352533E-3_dp

!   *** Rydberg constant [1/m] ***

!MK rydberg = 0.5_dp*e_mass*c_light*a_fine**2/h_planck
    rydberg = 10973731.568549_dp

!   *** Avogradro constant [1/mol] ***

    n_avogadro = 6.02214199E+23_dp

!   *** Boltzmann constant [J/K] ***

    boltzmann = 1.3806503E-23_dp

!   *** Atomic mass unit [kg]; conversion factor [u] -> [kg] ***

    a_mass = 1.66053873E-27_dp

!   *** Bohr radius [m];  ***

!MK a_bohr = a_fine/(4.0_dp*pi*rydberg)
    a_bohr = 0.5291772083E-10_dp

!   ---------------------------------------------------------------------------

!   *** Conversion factors ***

!   *** [u] -> [a.u.] ***

    massunit = a_mass/e_mass

!   *** [Bohr] -> [Angstrom] ***

    angstrom = 1.0E+10_dp*a_bohr

!   *** [Angstrom] -> [Bohr] ***

    bohr = 1.0_dp/angstrom

!   *** [a.u.] -> [s] ***

    seconds = 1.0_dp/(4.0_dp*pi*rydberg*c_light)

!   *** [a.u.] -> [fs] ***

    femtoseconds = 1.0E+15_dp*seconds

!   *** [a.u.] -> [J] ***

    joule = 2.0_dp*rydberg*h_planck*c_light

!   *** [a.u.] -> [K] ***

    kelvin = joule/boltzmann

!   *** [a.u.] -> [kJ/mol] ***

    kjmol = 0.001_dp*joule*n_avogadro

!   *** [a.u.] -> [kcal/mol] ***

    kcalmol = kjmol/4.184_dp

!   *** [a.u.] -> [eV] ***

    evolt = joule/e_charge

!   *** [a.u.] -> [Hz] ***

    hz = joule/h_planck

!   *** [a.u.] -> [MHz] ***

    mhz = 1.0E-6_dp*hz

!   *** [a.u./Bohr**2] -> [1/cm] (wave numbers) ***

    vibfac = 5.0_dp*SQRT(kjmol)/(pi*a_bohr*c_light)

!   *** [a.u.] -> [1/cm] (wave numbers) ***

    wavenum = 0.02_dp*rydberg

!   *** [a.u.] -> [esu] (electrostatic units) ***

    esu(1) = 1.0E+21_dp*a_bohr*c_light*e_charge

    DO i=2,maxmom
      esu(i) = esu(i-1)/bohr
    END DO

!   *** [a.u.] -> [debye] (electrostatic units) ***
    Debye = esu(1)

!   *** Initialize the periodic table ***

    CALL init_periodic_table()

  END SUBROUTINE init_physcon

! *****************************************************************************

  SUBROUTINE write_physcon(output_unit)

!   Purpose: Write all basic physical constants used by CP2K to a logical
!            output unit.

!   History: - Creation (14.11.2000,JGH)
!            - Updated to CODATA 1998 and cleaned (05.09.2003,MK)

!   ***************************************************************************

    INTEGER, INTENT(IN)                      :: output_unit

!   ---------------------------------------------------------------------------

    WRITE (UNIT=output_unit,FMT="(T2,A,/,/,(T2,A))")&
      "*** Fundamental physical constants (SI units) ***",&
      "*** Literature: B. J. Mohr and B. N. Taylor,",&
      "***             CODATA recommended values of the fundamental physical",&
      "***             constants: 1998, Rev. Mod. Phys. 72(2), 351 (2000)."

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
    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "Fine-structure constant",a_fine
    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "Rydberg constant [1/m]",rydberg
    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "Avogrado constant [1/mol]",n_avogadro
    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "Boltzmann constant [J/K]",boltzmann
    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "Atomic mass unit [kg]",a_mass
    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "Bohr radius [m]",a_bohr

!   *** Conversion factors ***

    WRITE (UNIT=output_unit,FMT="(/,T2,A,/)")&
      "*** Conversion factors ***"

    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "[u] -> [a.u.]",massunit,&
      "[Angstrom] -> [Bohr] = [a.u.]",bohr,&
      "[a.u.] = [Bohr] -> [Angstrom]",angstrom,&
      "[a.u.] -> [s]",seconds,&
      "[a.u.] -> [fs]",femtoseconds,&
      "[a.u.] -> [J]",joule,&
      "[a.u.] -> [K]",kelvin,&
      "[a.u.] -> [kJ/mol]",kjmol,&
      "[a.u.] -> [kcal/mol]",kcalmol,&
      "[a.u.] -> [eV]",evolt,&
      "[a.u.] -> [Hz]",hz,&
      "[a.u.] -> [MHz]",mhz,&
      "[a.u.] -> [1/cm] (wave numbers)",wavenum,&
      "[a.u./Bohr**2] -> [1/cm]",vibfac


  END SUBROUTINE write_physcon

! *****************************************************************************

END MODULE physcon
