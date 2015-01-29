!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief collects all references to literature in CP2K
!>      as new algorithms / method are included from literature sources
!>      we add a citation, to give credit to the authors
!>      (often ourselves, but this need not to be the case), and as a form
!>      of documentation.
!>
!>      reference_manager.F provides an easy way to cite these references from the code
!>      (e.g. whenever the relevant part of the code is executed)
!>      and to add citations to the input
!>      at the end of the run a bibliography is printed that can be used e.g. in papers
!> \note
!>      references need to be input using the ISI citation format, because it is
!>      uniform, easy to parse, and can be exported for example from web of science.
!>      (mark the reference, and click 'export to reference software', open the resulting file)
!>      Furthermore, it can be easily converted to and from using the bibutils tools
!>      a collection of easy to use conversion programs that can be found at
!>      http://www.scripps.edu/~cdputnam/software/bibutils/ by Chris Putnam
!> \par History
!>      08.2007 [Joost VandeVondele] Created
!>      01.2014 [Lianheng Tong] Added Kantorovich2008, Kantorovich2008a
!> \author Joost VandeVondele
! *****************************************************************************
MODULE bibliography
  
  USE reference_manager,               ONLY: add_reference,&
                                             cite_reference
  USE string_utilities,                ONLY: s2a
#include "cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cite_reference     ! for convenience
  PUBLIC :: add_all_references ! f77_interface only

  ! all reference keys
  INTEGER, PUBLIC, SAVE ::  VandeVondele2005a, VandeVondele2003,       &
       Laino2005, Laino2006, Goedecker1996, Hartwigsen1998, Krack2005, &
       Lippert1997, Lippert1999, Krack2000, Krack2002, Iannuzzi2005,   &
       Iannuzzi2006, Iannuzzi2007, Toukmaji1996, Martyna1999,          &
       VandeVondele2005b, Perdew1981, Avezac2005, Porezag1995,         &
       Seifert1996, Elstner1998, Zhechkov2005, Henkelman1999,          &
       Frigo2005, Nose1984a, Nose1984b, Brelaz1979, Bussi2007,         &
       Genovese2006, Genovese2007, Evans1983, Minary2003, Byrd1995,    &
       VandeVondele2007, Ortiz1994, Becke1988, Perdew1996, Zhang1998,  &
       Perdew2008, Lee1988, Heyd2006, Vydrov2006, Heyd2003, Heyd2004,  &
       Vosko1980, Aguado2003, Essmann1995, Ewald1921, Darden1993,      &
       Siepmann1995, Tersoff1988, Tosi1964a, Tosi1964b, Yamada2000,    &
       Dudarev1997, Dudarev1998, Dewar1977, Dewar1985, Rocha2006,      &
       Stewart1989, Thiel1992, Repasky2002, Stewart2007, Weber2008,    &
       Hunt2003, Guidon2008, Elber1987, Jonsson1998, Jonsson2000_1,    &
       Jonsson2000_2, Wales2004, Stewart1982, E2002, Laino2008,        &
       Blochl1995, Tao2003, VandeVondele2006, Hu2007, Grimme2006,      &
       Branduardi2007, Schenter2008, BeckeRoussel1989, Proynov2007,    &
       VandenCic2006, Becke1997, Foiles1986, Ricci2003, Kolafa2004,    &
       Kuhne2007, Kunert2003, Ceriotti2009, ceriotti2009b, Guidon2009, &
       BarducBus2008, Grimme2010, Guidon2010, Marques2012, Jones2011,  &
       Bernstein2009, Bernstein2012, Dick1958, Mitchell1993,           &
       Devynck2012, VandeVondele2012, Niklasson2003, Shao2003,         &
       VandeVondele2002, Dion2004, Romanperez2009, DelBen2012,         &
       Sabatini2013, DelBen2013, Kikuchi2009, Putrino2000,             &
       Putrino2002, Sebastiani2001, Weber2009, Tran2013, Golze2013,    &
       Tuckerman1992,  Zhao1994, Tozer1996, Goedecker2004,             &
       Khaliullin2013, Hutter2014, Bengtsson1999, Kantorovich2008,     &
       Kantorovich2008a, Wellendorff2012, Niklasson2014, Borstnik2014, &
       Grimme2011, Fattebert2002, Andreussi2012,                       &
       Khaliullin2007, Khaliullin2008, Merlot2014,Lin2009,Lin2013,     &
       QUIP_ref,DelBen2015

CONTAINS

! *****************************************************************************
!> \brief adds references that can later be cited / printed using the key
!> \par History
!>      08.2007 created [Joost VandeVondele]
!> \note
!>     - note that the ISI record is fixed format (line length limited and the
!>       first 3 characters can indicate record type)
!>     - please add DOI whenever available, this will result in a clickable
!>       link in the input reference manual.
! *****************************************************************************
  SUBROUTINE add_all_references()

    CALL add_reference(key=Wellendorff2012, ISI_record=s2a(&
         "AU Wellendorff, J",&
         "   Lundgaard, K",&
         "   Mogelhoj, A",&
         "   Petzold, V",&
         "   Landis, D",&
         "   Norskov, J",&
         "   Bligaard, T",&
         "   Jacobsen, K",&
         "TI Density functionals for surface science: Exchange-correlation model development with Bayesian error estimation",&
         "SO PHYSICAL REVIEW B",&
         "PY 2012",&
         "VL 85",&
         "IS 23",&
         "BP 235149",&
         "ER"),&
          DOI="10.1103/PhysRevB.85.235149")

    CALL add_reference(key=Brelaz1979, ISI_record=s2a(&
         "PT J",&
         "AU Brelaz, D",&
         "TI New methods to color the vertices of a graph",&
         "SO COMMUNICATIONS OF THE ACM",&
         "SN 0001-0782",&
         "PY 1979",&
         "VL 22",&
         "IS 4",&
         "BP 251",&
         "EP 256",&
         "UT WOS:A1979GR73200004",&
         "ER"),&
         DOI="10.1145/359094.359101")

    CALL add_reference(key=Bengtsson1999, ISI_record=s2a(&
         "AU BENGTSSON, L",&
         "TI DIPOLE CORRECTION FOR SURFACE SUPERCELL CALCULATIONE",&
         "SO PHYSICAL REVIEW B",&
         "SN 0163-1829",&
         "PY 1999",&
         "VL 59",&
         "IS 19",&
         "BP 12301",&
         "EP 12304",&
         "ER"),&
          DOI="10.1103/PhysRevB.59.12301")


    CALL add_reference(key=Foiles1986, ISI_record=s2a(&
         "AU FOILES, SM",&
         "   BASKES, MI",&
         "   DAW, MS",&
         "TI EMBEDDED-ATOM-METHOD FUNCTIONS FOR THE FCC METALS CU, AG, AU, NI, PD,",&
         "   PT, AND THEIR ALLOYS",&
         "SO PHYSICAL REVIEW B",&
         "SN 0163-1829",&
         "PD JUN 15",&
         "PY 1986",&
         "VL 33",&
         "IS 12",&
         "PN Part 1",&
         "BP 7983",&
         "EP 7991",&
         "UT ISI:A1986C899400011",&
         "ER"),&
         DOI="10.1103/PhysRevB.33.7983")

    CALL add_reference(key=QUIP_ref, ISI_record=s2a(&
         "AU QUIP", &
         "TI libAtoms/QUIP libraries from http://www.libatoms.org, please cite",&
         "   web site and references for appropriate potential invoked",&
         "SO web site", &
         "PY 2014", &
         "ER"),&
         DOI="")

    CALL add_reference(key=VandenCic2006, ISI_record=s2a(&
         "AU Vanden-Eijnden, E",&
         "   Ciccotti, G",&
         "AF Vanden-Eijnden, Eric",&
         "   Ciccotti, Giovanni",&
         "TI Second-order integrators for Langevin equations with holonomic",&
         "   constraints",&
         "SO CHEMICAL PHYSICS LETTERS",&
         "SN 0009-2614",&
         "PD SEP 29",&
         "PY 2006",&
         "VL 429",&
         "IS 1-3",&
         "BP 310",&
         "EP 316",&
         "UT ISI:000241018800060",&
         "ER"),&
         DOI="10.1016/j.cplett.2006.07.086")

    CALL add_reference(key=Hu2007, ISI_record=s2a(&
         "AU Hu, H",&
         "   Lu, ZY",&
         "   Elstner, M",&
         "   Hermans, J",&
         "   Yang, WT",&
         "TI Simulating water with the self-consistent-charge density functional",&
         "   tight binding method: From molecular clusters to the liquid state",&
         "SO JOURNAL OF PHYSICAL CHEMISTRY A",&
         "CT 232nd National Meeting of the American-Chemical-Society",&
         "CY SEP 10-14, 2006",&
         "CL San Francisco, CA",&
         "SN 1089-5639",&
         "PD JUL 5",&
         "PY 2007",&
         "VL 111",&
         "IS 26",&
         "BP 5685",&
         "EP 5691",&
         "UT ISI:000247573600013",&
         "ER"),&
         DOI="10.1021/jp070308d")

    CALL add_reference(key=Zhao1994, ISI_record=s2a(&
         "AU ZHAO, QS",&
         "   MORRISON, RC",&
         "   PARR, RG",&
         "TI FOR ELECTRON DENSITIES TO KOHN-SHAM KINETIC ENERGIES, ORBITAL ENERGIES.",&
         "   EXCHANGE-CORRELATION POTENTIALS, AND EXCHANGE-CORRELATION ENERGIES",&
         "SO PHYSICAL REVIEW A",&
         "SN 1050-2947",&
         "PD SEP 1",&
         "PY 1994",&
         "VL 50",&
         "IS 3",&
         "BP 2138",&
         "EP 2142",&
         "UT WOS:A1994PG80900023",&
         "ER"),&
         DOI="10.1103/PhysRevA.50.2138")

    CALL add_reference(key=Tozer1996, ISI_record=s2a(&
         "AU TOZER, DJ",&
         "   INGAMELLS, VE",&
         "   HANDY, NC",&
         "TI EXCHANGE-CORRELATION POTENTIALS",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD NOV 22",&
         "PY 1996",&
         "VL 105",&
         "IS 20",&
         "BP 9200",&
         "EP 9213",&
         "UT WOS:A1996VV04900024",&
         "ER"),&
         DOI="10.1063/1.472753")


    CALL add_reference(key=Blochl1995, ISI_record=s2a(&
         "AU BLOCHL, PE",&
         "TI ELECTROSTATIC DECOUPLING OF PERIODIC IMAGES OF PLANE-WAVE-EXPANDED",&
         "   DENSITIES AND DERIVED ATOMIC POINT CHARGES",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD NOV 1",&
         "PY 1995",&
         "VL 103",&
         "IS 17",&
         "BP 7422",&
         "EP 7428",&
         "UT ISI:A1995TB63300019",&
         "ER"),&
         DOI="10.1063/1.470314")

    CALL add_reference(key=Laino2008, ISI_record=s2a(&
         "AU Laino, T",&
         "   Hutter, J",&
         'TI Notes on "Ewald summation of electrostatic multipole interactions up to',&
         '   quadrupolar level" [J. Chem. Phys. 119, 7471 (2003)]',&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD AUG 21",&
         "PY 2008",&
         "VL 129",&
         "IS 7",&
         "AR 074102",&
         "UT ISI:000258643300004",&
         "ER"),&
         DOI="10.1063/1.2970887")

    CALL add_reference(key=E2002, ISI_record=s2a(&
         "AU E, WN",&
         "   Ren, WQ",&
         "   Vanden-Eijnden, E",&
         "TI String method for the study of rare events",&
         "SO PHYSICAL REVIEW B",&
         "SN 1098-0121",&
         "PD AUG 1",&
         "PY 2002",&
         "VL 66",&
         "IS 5",&
         "AR 052301",&
         "UT ISI:000177873000010",&
         "ER"),&
         DOI="10.1103/PhysRevB.66.052301")

    CALL add_reference(key=Wales2004, ISI_record=s2a(&
         "AU Trygubenko, SA",&
         "   Wales, DJ",&
         "TI A doubly nudged elastic band method for finding transition states",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD FEB 1",&
         "PY 2004",&
         "VL 120",&
         "IS 5",&
         "BP 2082",&
         "EP 2094",&
         "DI 10.1063/1.1636455",&
         "UT ISI:000188498400003",&
         "ER"),&
         DOI="10.1063/1.1636455")

    CALL add_reference(key=Jonsson2000_2, ISI_record=s2a(&
         "AU Henkelman, G",&
         "   Uberuaga, BP",&
         "   Jonsson, H",&
         "TI A climbing image nudged elastic band method for finding saddle points",&
         "   and minimum energy paths",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD DEC 8",&
         "PY 2000",&
         "VL 113",&
         "IS 22",&
         "BP 9901",&
         "EP 9904",&
         "UT ISI:000165584900005",&
         "ER"),&
         DOI="10.1063/1.1329672")

    CALL add_reference(key=Jonsson2000_1, ISI_record=s2a(&
         "AU Henkelman, G",&
         "   Jonsson, H",&
         "TI Improved tangent estimate in the nudged elastic band method for finding",&
         "   minimum energy paths and saddle points",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD DEC 8",&
         "PY 2000",&
         "VL 113",&
         "IS 22",&
         "BP 9978",&
         "EP 9985",&
         "UT ISI:000165584900014",&
         "ER"),&
         DOI="10.1063/1.1323224")

    CALL add_reference(key=Jonsson1998, ISI_record=s2a(&
         "AU JONSSON, H",&
         "   MILLS, G",&
         "   JACOBSEN, K W",&
         "TI Nudged elastic band method for finding minimum energy paths ",&
         "   of transitions",&
         "SO Classical and Quantum Dynamics in Condensed Phase Simulations",&
         "PY 1998",&
         "BP 385",&
         "EP 404",&
         "ER"),&
         DOI="")

    CALL add_reference(key=Elber1987, ISI_record=s2a(&
         "AU ELBER, R",&
         "   KARPLUS, M",&
         "TI A METHOD FOR DETERMINING REACTION PATHS IN LARGE MOLECULES -",&
         "   APPLICATION TO MYOGLOBIN",&
         "SO CHEMICAL PHYSICS LETTERS",&
         "SN 0009-2614",&
         "PD SEP 4",&
         "PY 1987",&
         "VL 139",&
         "IS 5",&
         "BP 375",&
         "EP 380",&
         "UT ISI:A1987K321300001",&
         "ER"),&
         DOI="10.1016/0009-2614(87)80576-6")

    CALL add_reference(key=Weber2008, ISI_record=s2a(&
         "AU Weber, V",&
         "   VandeVondele, J",&
         "   Hutter, J",&
         "   Niklasson, AMN",&
         "AF Weber, Valery",&
         "   VandeVondele, Joost",&
         "   Hutter, Juerg",&
         "   Niklasson, Anders M. N.",&
         "TI Direct energy functional minimization under orthogonality constraints",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD FEB 28",&
         "PY 2008",&
         "VL 128",&
         "IS 8",&
         "AR 084113",&
         "DI ARTN 084113",&
         "UT ISI:000254047200013",&
         "ER"),&
          DOI="10.1063/1.2841077")

    CALL add_reference(key=Stewart2007, ISI_record=s2a(&
         "AU Stewart, JJP",&
         "AF Stewart, James J. P.",&
         "TI Optimization of parameters for semiempirical methods V: Modification of",&
         "   NDDO approximations and application to 70 elements",&
         "SO JOURNAL OF MOLECULAR MODELING",&
         "SN 1610-2940",&
         "PD DEC",&
         "PY 2007",&
         "VL 13",&
         "IS 12",&
         "BP 1173",&
         "EP 1213",&
         "UT ISI:000250308000001",&
         "ER"),&
         DOI="10.1007/s00894-007-0233-4")

    CALL add_reference(key=Repasky2002, ISI_record=s2a(&
         "AU Repasky, MP",&
         "   Chandrasekhar, J",&
         "   Jorgensen, WL",&
         "TI PDDG/PM3 and PDDG/MNDO: Improved semiempirical methods",&
         "SO JOURNAL OF COMPUTATIONAL CHEMISTRY",&
         "SN 0192-8651",&
         "PD DEC",&
         "PY 2002",&
         "VL 23",&
         "IS 16",&
         "BP 1601",&
         "EP 1622",&
         "UT ISI:000179022700011",&
         "ER"),&
         DOI="10.1002/jcc.10162")

    CALL add_reference(key=Thiel1992, ISI_record=s2a(&
         "AU THIEL, W",&
         "   VOITYUK, AA",&
         "TI EXTENSION OF THE MNDO FORMALISM TO D-ORBITALS - INTEGRAL APPROXIMATIONS",&
         "   AND PRELIMINARY NUMERICAL RESULTS",&
         "SO THEORETICA CHIMICA ACTA",&
         "SN 0040-5744",&
         "PD FEB",&
         "PY 1992",&
         "VL 81",&
         "IS 6",&
         "BP 391",&
         "EP 404",&
         "UT ISI:A1992HF04700004",&
         "ER"),&
         DOI="10.1007/BF01134863")

    CALL add_reference(key=Stewart1989, ISI_record=s2a(&
         "AU STEWART, JJP",&
         "TI OPTIMIZATION OF PARAMETERS FOR SEMIEMPIRICAL METHODS .1. METHOD",&
         "SO JOURNAL OF COMPUTATIONAL CHEMISTRY",&
         "SN 0192-8651",&
         "PD MAR",&
         "PY 1989",&
         "VL 10",&
         "IS 2",&
         "BP 209",&
         "EP 220",&
         "UT ISI:A1989T690100007",&
         "ER"),&
         DOI="10.1002/jcc.540100208")

    CALL add_reference(key=Rocha2006, ISI_record=s2a(&
         "AU Rocha, GB",&
         "   Freire, RO",&
         "   Simas, AM",&
         "   Stewart, JJP",&
         "TI RMI: A reparameterization of AM1 for H, C, N, O, P, S, F, Cl, Br, and I",&
         "SO JOURNAL OF COMPUTATIONAL CHEMISTRY",&
         "SN 0192-8651",&
         "PD JUL 30",&
         "PY 2006",&
         "VL 27",&
         "IS 10",&
         "BP 1101",&
         "EP 1111",&
         "UT ISI:000238158900007",&
         "ER"),&
         DOI="10.1002/jcc.20425")

    CALL add_reference(key=Dewar1985, ISI_record=s2a(&
         "AU DEWAR, MJS",&
         "   ZOEBISCH, EG",&
         "   HEALY, EF",&
         "   STEWART, JJP",&
         "TI THE DEVELOPMENT AND USE OF QUANTUM-MECHANICAL MOLECULAR-MODELS .76. AM1",&
         "   - A NEW GENERAL-PURPOSE QUANTUM-MECHANICAL MOLECULAR-MODEL",&
         "SO JOURNAL OF THE AMERICAN CHEMICAL SOCIETY",&
         "SN 0002-7863",&
         "PY 1985",&
         "VL 107",&
         "IS 13",&
         "BP 3902",&
         "EP 3909",&
         "UT ISI:A1985ALC3500024",&
         "ER"),&
         DOI="10.1021/ja00299a024")

    CALL add_reference(key=Dewar1977, ISI_record=s2a(&
         "AU DEWAR, MJS",&
         "   THIEL, W",&
         "TI GROUND-STATES OF MOLECULES .38. MNDO METHOD - APPROXIMATIONS AND",&
         "   PARAMETERS",&
         "SO JOURNAL OF THE AMERICAN CHEMICAL SOCIETY",&
         "SN 0002-7863",&
         "PY 1977",&
         "VL 99",&
         "IS 15",&
         "BP 4899",&
         "EP 4907",&
         "UT ISI:A1977DN51700004",&
         "ER"),&
         DOI="10.1021/ja00457a004")

    CALL add_reference(key=Henkelman1999, ISI_record=s2a(&
         "PT J",&
         "AU Henkelman, G",&
         "   Jonsson, H",&
         "TI A dimer method for finding saddle points on high dimensional potential",&
         "   surfaces using only first derivatives",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD OCT 15",&
         "PY 1999",&
         "VL 111",&
         "IS 15",&
         "BP 7010",&
         "EP 7022",&
         "UT ISI:000083051600040",&
         "ER"),&
         DOI="10.1063/1.480097")

    CALL add_reference(key=Aguado2003, ISI_record=s2a(&
         "PT J",&
         "AU Aguado, A",&
         "   Madden, PA",&
         "TI Ewald summation of electrostatic multipole interactions up to the",&
         "   quadrupolar level",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD OCT 8",&
         "PY 2003",&
         "VL 119",&
         "IS 14",&
         "BP 7471",&
         "EP 7483",&
         "UT ISI:000185575700055",&
         "ER"),&
         DOI="10.1063/1.1605941")

    CALL add_reference(key=Yamada2000, ISI_record=s2a(&
         "PT J",&
         "AU Yamada, K",&
         "   Kurosaki, K",&
         "   Uno, M",&
         "   Yamanaka, S",&
         "TI Evaluation of thermal properties of uranium dioxide by molecular",&
         "   dynamics",&
         "SO JOURNAL OF ALLOYS AND COMPOUNDS",&
         "SN 0925-8388",&
         "PD JUL 14",&
         "PY 2000",&
         "VL 307",&
         "BP 10",&
         "EP 16",&
         "UT ISI:000088590900003",&
         "ER"),&
         DOI="10.1016/S0925-8388(00)00806-9")

    CALL add_reference(key=Tosi1964a, ISI_record=s2a(&
         "PT J",&
         "AU FUMI, FG",&
         "   TOSI, MP",&
         "TI IONIC SIZES + BORN REPULSIVE PARAMETERS IN NACL-TYPE ALKALI HALIDES .I.",&
         "   HUGGINS-MAYER + PAULING FORMS",&
         "SO JOURNAL OF PHYSICS AND CHEMISTRY OF SOLIDS",&
         "SN 0022-3697",&
         "PY 1964",&
         "VL 25",&
         "IS 1",&
         "BP 31",&
         "EP 43",&
         "UT ISI:A19646597B00004",&
         "ER"),&
         DOI="10.1016/0022-3697(64)90159-3")

    CALL add_reference(key=Tosi1964b, ISI_record=s2a(&
         "PT J",&
         "AU TOSI, MP",&
         "   FUMI, FG",&
         "TI IONIC SIZES + BORN REPULSIVE PARAMETERS IN NACL-TYPE ALKALI HALIDES .2.",&
         "   GENERALIZED",&
         "SO JOURNAL OF PHYSICS AND CHEMISTRY OF SOLIDS",&
         "SN 0022-3697",&
         "PY 1964",&
         "VL 25",&
         "IS 1",&
         "BP 45",&
         "EP 52",&
         "UT ISI:A19646597B00017",&
         "ER"),&
         DOI="10.1016/0022-3697(64)90160-X")

    CALL add_reference(key=Tersoff1988, ISI_record=s2a(&
         "PT J",&
         "AU TERSOFF, J",&
         "TI EMPIRICAL INTERATOMIC POTENTIAL FOR SILICON WITH IMPROVED ELASTIC",&
         "   PROPERTIES",&
         "SO PHYSICAL REVIEW B",&
         "SN 0163-1829",&
         "PD NOV 15",&
         "PY 1988",&
         "VL 38",&
         "IS 14",&
         "BP 9902",&
         "EP 9905",&
         "UT ISI:A1988Q924400066",&
         "ER"),&
         DOI="10.1103/PhysRevB.38.9902")

    CALL add_reference(key=Siepmann1995, ISI_record=s2a(&
         "PT J",&
         "AU SIEPMANN, JI",&
         "   SPRIK, M   ",&
         "TI INFLUENCE OF SURFACE TOPOLOGY AND ELECTROSTATIC POTENTIAL ON",&
         "   WATER/ELECTRODE SYSTEMS",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD JAN 1",&
         "PY 1995",&
         "VL 102",&
         "IS 1",&
         "BP 511",&
         "EP 524",&
         "UT ISI:A1995PZ71400055",&
         "ER"),&
         DOI="10.1063/1.469429")

    CALL add_reference(key=Bussi2007,ISI_record=s2a(&
         "PT J",&
         "AU Bussi, G",&
         "   Donadio, D",&
         "   Parrinello, M",&
         "AF Bussi, Giovanni",&
         "   Donadio, Davide",&
         "   Parrinello, Michele",&
         "TI Canonical sampling through velocity rescaling",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD JAN 7",&
         "PY 2007",&
         "VL 126",&
         "IS 1",&
         "AR 014101",&
         "DI ARTN 014101",&
         "UT ISI:000243380000005",&
         "ER"),&
         DOI="10.1063/1.2408420")

    CALL add_reference(key=Nose1984a,ISI_record=s2a(&
         "PT J",&
         "AU NOSE, S",&
         "TI A UNIFIED FORMULATION OF THE CONSTANT TEMPERATURE MOLECULAR-DYNAMICS",&
         "   METHODS",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PY 1984",&
         "VL 81",&
         "IS 1",&
         "BP 511",&
         "EP 519",&
         "UT ISI:A1984TA66100062",&
         "ER"),&
         DOI="10.1063/1.447334")

    CALL add_reference(key=Nose1984b,ISI_record=s2a(&
         "PT J",&
         "AU NOSE, S",&
         "TI A MOLECULAR-DYNAMICS METHOD FOR SIMULATIONS IN THE CANONICAL ENSEMBLE",&
         "SO MOLECULAR PHYSICS",&
         "LA English",&
         "DT Article",&
         "NR 17",&
         "TC 1593",&
         "PU TAYLOR & FRANCIS LTD",&
         "PI LONDON",&
         "PA ONE GUNDPOWDER SQUARE, LONDON, ENGLAND EC4A 3DE",&
         "SN 0026-8976",&
         "J9 MOL PHYS",&
         "JI Mol. Phys.",&
         "PY 1984",&
         "VL 52",&
         "IS 2",&
         "BP 255",&
         "EP 268",&
         "PG 14",&
         "SC Physics, Atomic, Molecular & Chemical",&
         "GA SV647",&
         "UT ISI:A1984SV64700001",&
         "ER"),&
         DOI="10.1080/00268978400101201")

    CALL add_reference(key=VandeVondele2005a,ISI_record=s2a(&
         "PT J",&
         "AU VandeVondele, J",&
         "   Krack, M",&
         "   Mohamed, F",&
         "   Parrinello, M",&
         "   Chassaing, T",&
         "   Hutter, J",&
         "TI QUICKSTEP: Fast and accurate density functional calculations using a",&
         "   mixed Gaussian and plane waves approach",&
         "SO COMPUTER PHYSICS COMMUNICATIONS",&
         "SN 0010-4655",&
         "PD APR 15",&
         "PY 2005",&
         "VL 167",&
         "IS 2",&
         "BP 103",&
         "EP 128",&
         "UT ISI:000228421500005",&
         "ER"),&
         DOI="10.1016/j.cpc.2004.12.014")

    CALL add_reference(key=VandeVondele2003,ISI_record=s2a(&
         "PT J",&
         "AU VandeVondele, J",&
         "   Hutter, J",&
         "TI An efficient orbital transformation method for electronic structure",&
         "   calculations",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD MAR 8",&
         "PY 2003",&
         "VL 118",&
         "IS 10",&
         "BP 4365",&
         "EP 4369",&
         "UT ISI:000181065800006",&
         "ER"),&
         DOI="10.1063/1.1543154")

    CALL add_reference(key=Laino2005,ISI_record=s2a(&
         "PT J",&
         "AU Laino, T",&
         "   Mohamed, F",&
         "   Laio, A",&
         "   Parrinello, M",&
         "TI An efficient real space multigrid OM/MM electrostatic coupling",&
         "SO JOURNAL OF CHEMICAL THEORY AND COMPUTATION",&
         "SN 1549-9618",&
         "PD NOV-DEC",&
         "PY 2005",&
         "VL 1",&
         "IS 6",&
         "BP 1176",&
         "EP 1184",&
         "UT ISI:000233260300012",&
         "ER"),&
         DOI="10.1021/ct050123f")

    CALL add_reference(key=Laino2006,ISI_record=s2a(&
         "PT J",&
         "AU Laino, T",&
         "   Mohamed, F",&
         "   Laio, A",&
         "   Parrinello, M",&
         "AF Laino, Teodoro",&
         "   Mohamed, Fawzi",&
         "   Laio, Alessandro",&
         "   Parrinello, Michele",&
         "TI An efficient linear-scaling electrostatic coupling for treating",&
         "   periodic boundary conditions in QM/MM simulations",&
         "SO JOURNAL OF CHEMICAL THEORY AND COMPUTATION",&
         "SN 1549-9618",&
         "PD SEP 12",&
         "PY 2006",&
         "VL 2",&
         "IS 5",&
         "BP 1370",&
         "EP 1378",&
         "UT ISI:000240437000015",&
         "ER"),&
         DOI="10.1021/ct6001169")

    CALL add_reference(key=Goedecker1996,ISI_record=s2a(&
         "PT J",&
         "AU Goedecker, S",&
         "   Teter, M",&
         "   Hutter, J",&
         "TI Separable dual-space Gaussian pseudopotentials",&
         "SO PHYSICAL REVIEW B",&
         "SN 1098-0121",&
         "PD JUL 15",&
         "PY 1996",&
         "VL 54",&
         "IS 3",&
         "BP 1703",&
         "EP 1710",&
         "UT ISI:A1996UZ86100053",&
         "ER"),&
         DOI="10.1103/PhysRevB.54.1703")

    CALL add_reference(key=Hartwigsen1998,ISI_record=s2a(&
         "PT J",&
         "AU Hartwigsen, C",&
         "   Goedecker, S",&
         "   Hutter, J",&
         "TI Relativistic separable dual-space Gaussian pseudopotentials from H to Rn",&
         "SO PHYSICAL REVIEW B",&
         "SN 0163-1829",&
         "PD AUG 15",&
         "PY 1998",&
         "VL 58",&
         "IS 7",&
         "BP 3641",&
         "EP 3662",&
         "UT ISI:000075616800043",&
         "ER"),&
         DOI="10.1103/PhysRevB.58.3641")

    CALL add_reference(key=Krack2005,ISI_record=s2a(&
         "PT J",&
         "AU Krack, M",&
         "TI Pseudopotentials for H to Kr optimized for gradient-corrected",&
         "   exchange-correlation functionals",&
         "SO THEORETICAL CHEMISTRY ACCOUNTS",&
         "SN 1432-881X",&
         "PD SEP",&
         "PY 2005",&
         "VL 114",&
         "IS 1-3",&
         "BP 145",&
         "EP 152",&
         "UT ISI:000232053800020",&
         "ER"),&
         DOI="10.1007/s00214-005-0655-y")

    CALL add_reference(key=Lippert1997,ISI_record=s2a(&
         "PT J",&
         "AU Lippert, G",&
         "   Hutter, J",&
         "   Parrinello, M",&
         "TI A hybrid Gaussian and plane wave density functional scheme",&
         "SO MOLECULAR PHYSICS",&
         "SN 0026-8976",&
         "PD OCT 20",&
         "PY 1997",&
         "VL 92",&
         "IS 3",&
         "BP 477",&
         "EP 487",&
         "UT ISI:A1997YC60700017",&
         "ER"),&
         DOI="10.1080/002689797170220")

    CALL add_reference(key=Lippert1999,ISI_record=s2a(&
         "PT J",&
         "AU Lippert, G",&
         "   Hutter, J",&
         "   Parrinello, M",&
         "TI The Gaussian and augmented-plane-wave density functional method for ab",&
         "   initio molecular dynamics simulations",&
         "SO THEORETICAL CHEMISTRY ACCOUNTS",&
         "SN 1432-881X",&
         "PD DEC",&
         "PY 1999",&
         "VL 103",&
         "IS 2",&
         "BP 124",&
         "EP 140",&
         "UT ISI:000084308100005",&
         "ER"),&
         DOI="10.1007/s002140050523")

    CALL add_reference(key=Krack2002,ISI_record=s2a(&
         "PT J",&
         "AU Krack, M",&
         "   Gambirasio, A",&
         "   Parrinello, M",&
         "TI Ab initio x-ray scattering of liquid water",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD NOV 22",&
         "PY 2002",&
         "VL 117",&
         "IS 20",&
         "BP 9409",&
         "EP 9412",&
         "UT ISI:000179042300036",&
         "ER"),&
         DOI="10.1063/1.1517040")

    CALL add_reference(key=Krack2000,ISI_record=s2a(&
         "PT J",&
         "AU Krack, M",&
         "   Parrinello, M",&
         "TI All-electron ab-initio molecular dynamics",&
         "SO PHYSICAL CHEMISTRY CHEMICAL PHYSICS",&
         "SN 1463-9076",&
         "PY 2000",&
         "VL 2",&
         "IS 10",&
         "BP 2105",&
         "EP 2112",&
         "UT ISI:000086928000007",&
         "ER"),&
         DOI="10.1039/b001167n")

    CALL add_reference(key=Iannuzzi2007,ISI_record=s2a(&
         "PT J",&
         "AU Iannuzzi, M",&
         "   Hutter, J",&
         "AF Iannuzzi, Marcella",&
         "   Hutter, Juerg",&
         "TI Inner-shell spectroscopy by the Gaussian and augmented plane wave method",&
         "SO PHYSICAL CHEMISTRY CHEMICAL PHYSICS",&
         "SN 1463-9076",&
         "PY 2007",&
         "VL 9",&
         "IS 13",&
         "BP 1599",&
         "EP 1610",&
         "UT ISI:000245633400008",&
         "ER"),&
         DOI="10.1039/b615522g")

    CALL add_reference(key=Iannuzzi2006,ISI_record=s2a(&
         "PT J",&
         "AU Iannuzzi, M",&
         "   Kirchner, B",&
         "   Hutter, J",&
         "TI Density functional embedding for molecular systems",&
         "SO CHEMICAL PHYSICS LETTERS",&
         "SN 0009-2614",&
         "PD APR 3",&
         "PY 2006",&
         "VL 421",&
         "IS 1-3",&
         "BP 16",&
         "EP 20",&
         "UT ISI:000236662600004",&
         "ER"),&
         DOI="10.1016/j.cplett.2005.08.155")

    CALL add_reference(key=Iannuzzi2005,ISI_record=s2a(&
         "PT J",&
         "AU Iannuzzi, M",&
         "   Chassaing, T",&
         "   Wallman, T",&
         "   Hutter, J",&
         "TI Ground and excited state density functional calculations with the",&
         "   Gaussian and augmented-plane-wave method",&
         "SO CHIMIA",&
         "SN 0009-4293",&
         "PY 2005",&
         "VL 59",&
         "IS 7-8",&
         "BP 499",&
         "EP 503",&
         "UT ISI:000231615300005",&
         "ER"),&
         DOI="10.2533/000942905777676164")

    CALL add_reference(key=Toukmaji1996,ISI_record=s2a(&
         "PT J",&
         "AU Toukmaji, AY",&
         "   Board, JA",&
         "TI Ewald summation techniques in perspective: A survey",&
         "SO COMPUTER PHYSICS COMMUNICATIONS",&
         "SN 0010-4655",&
         "PD JUN",&
         "PY 1996",&
         "VL 95",&
         "IS 2-3",&
         "BP 73",&
         "EP 92",&
         "UT ISI:A1996UQ97300001",&
         "ER"),&
         DOI="10.1016/0010-4655(96)00016-1")

    CALL add_reference(key=Martyna1999,ISI_record=s2a(&
         "PT J",&
         "AU Martyna, GJ",&
         "   Tuckerman, ME",&
         "TI A reciprocal space based method for treating long range interactions in",&
         "   ab initio and force-field-based calculations in clusters",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD FEB 8",&
         "PY 1999",&
         "VL 110",&
         "IS 6",&
         "BP 2810",&
         "EP 2821",&
         "UT ISI:000078379800010",&
         "ER"),&
         DOI="10.1063/1.477923")

    CALL add_reference(key=VandeVondele2005b,ISI_record=s2a(&
         "PT J",&
         "AU VandeVondele, J",&
         "   Sprik, M",&
         "TI A molecular dynamics study of the hydroxyl radical in solution applying",&
         "   self-interaction-corrected density functional methods",&
         "SO PHYSICAL CHEMISTRY CHEMICAL PHYSICS",&
         "SN 1463-9076",&
         "PY 2005",&
         "VL 7",&
         "IS 7",&
         "BP 1363",&
         "EP 1367",&
         "UT ISI:000227707200005",&
         "ER"),&
         DOI="10.1039/b501603g")

    CALL add_reference(key=Perdew1981,ISI_record=s2a(&
         "PT J",&
         "AU PERDEW, JP",&
         "   ZUNGER, A",&
         "TI SELF-INTERACTION CORRECTION TO DENSITY-FUNCTIONAL APPROXIMATIONS FOR",&
         "   MANY-ELECTRON SYSTEMS",&
         "SO PHYSICAL REVIEW B",&
         "SN 0163-1829",&
         "PY 1981",&
         "VL 23",&
         "IS 10",&
         "BP 5048",&
         "EP 5079",&
         "UT ISI:A1981LR24700023",&
         "ER"),&
         DOI="10.1103/PhysRevB.23.5048")

    CALL add_reference(key=Avezac2005,ISI_record=s2a(&
         "PT J",&
         'AU d''Avezac, M',&
         "   Calandra, M",&
         "   Mauri, F",&
         "TI Density functional theory description of hole-trapping in SiO2: A",&
         "   self-interaction-corrected approach",&
         "SO PHYSICAL REVIEW B",&
         "SN 1098-0121",&
         "PD MAY",&
         "PY 2005",&
         "VL 71",&
         "IS 20",&
         "AR 205210",&
         "DI ARTN 205210",&
         "UT ISI:000230244300043",&
         "ER"),&
         DOI="10.1103/PhysRevB.71.205210")

    CALL add_reference(key=Zhechkov2005,ISI_record=s2a(&
         "PT J",&
         "AU Zhechkov, L",&
         "   Heine, T",&
         "   Patchkovskii, S",&
         "   Seifert, G",&
         "   Duarte, HA",&
         "TI An efficient a Posteriori treatment for dispersion interaction in",&
         "   density-functional-based tight binding",&
         "SO JOURNAL OF CHEMICAL THEORY AND COMPUTATION",&
         "SN 1549-9618",&
         "PD SEP-OCT",&
         "PY 2005",&
         "VL 1",&
         "IS 5",&
         "BP 841",&
         "EP 847",&
         "UT ISI:000231871500009",&
         "ER"),&
         DOI="10.1021/ct050065y")

    CALL add_reference(key=Elstner1998,ISI_record=s2a(&
         "PT J",&
         "AU Elstner, M",&
         "   Porezag, D",&
         "   Jungnickel, G",&
         "   Elsner, J",&
         "   Haugk, M",&
         "   Frauenheim, T",&
         "   Suhai, S",&
         "   Seifert, G",&
         "TI Self-consistent-charge density-functional tight-binding method for",&
         "   simulations of complex materials properties",&
         "SO PHYSICAL REVIEW B",&
         "SN 1098-0121",&
         "PD SEP 15",&
         "PY 1998",&
         "VL 58",&
         "IS 11",&
         "BP 7260",&
         "EP 7268",&
         "UT ISI:000076058800081",&
         "ER"),&
         DOI="10.1103/PhysRevB.58.7260")

    CALL add_reference(key=Seifert1996,ISI_record=s2a(&
         "PT J",&
         "AU Seifert, G",&
         "   Porezag, D",&
         "   Frauenheim, T",&
         "TI Calculations of molecules, clusters, and solids with a simplified",&
         "   LCAO-DFT-LDA scheme",&
         "SO INTERNATIONAL JOURNAL OF QUANTUM CHEMISTRY",&
         "SN 0020-7608",&
         "PD APR 15",&
         "PY 1996",&
         "VL 58",&
         "IS 2",&
         "BP 185",&
         "EP 192",&
         "UT ISI:A1996UB73600008",&
         "ER"),&
         DOI="10.1002/(SICI)1097-461X(1996)58:2<185::AID-QUA7>3.0.CO;2-U")

    CALL add_reference(key=Porezag1995,ISI_record=s2a(&
         "PT J",&
         "AU POREZAG, D",&
         "   FRAUENHEIM, T",&
         "   KOHLER, T",&
         "   SEIFERT, G",&
         "   KASCHNER, R",&
         "TI CONSTRUCTION OF TIGHT-BINDING-LIKE POTENTIALS ON THE BASIS OF",&
         "   DENSITY-FUNCTIONAL THEORY - APPLICATION TO CARBON",&
         "SO PHYSICAL REVIEW B",&
         "SN 0163-1829",&
         "PD MAY 15",&
         "PY 1995",&
         "VL 51",&
         "IS 19",&
         "BP 12947",&
         "EP 12957",&
         "UT ISI:A1995QZ16500002",&
         "ER"),&
         DOI="10.1103/PhysRevB.51.12947")

    CALL add_reference(key=Frigo2005,ISI_record=s2a(&
         "PT J",&
         "AU Frigo, M",&
         "   Johnson, SG",&
         "TI The design and implementation of FFTW3",&
         "SO PROCEEDINGS OF THE IEEE",&
         "SN 0018-9219",&
         "PD FEB",&
         "PY 2005",&
         "VL 93",&
         "IS 2",&
         "BP 216",&
         "EP 231",&
         "UT ISI:000226542300002",&
         "ER"),&
         DOI="10.1109/JPROC.2004.840301")

    CALL add_reference(key=Genovese2006,ISI_record=s2a(&
         "PT J",&
         "AU Genovese, L",&
         "   Deutsch, T",&
         "   Neelov, A",&
         "   Goedecker, S",&
         "   Beylkin, G",&
         "AF Genovese, Luigi",&
         "   Deutsch, Thierry",&
         "   Neelov, Alexey",&
         "   Goedecker, Stefan",&
         "   Beylkin, Gregory",&
         "TI Efficient solution of Poisson's equation with free boundary conditions",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD AUG 21",&
         "PY 2006",&
         "VL 125",&
         "IS 7",&
         "AR 074105",&
         "DI ARTN 074105",&
         "UT ISI:000239914800005",&
         "ER"),&
         DOI="10.1063/1.2335442")

    CALL add_reference(key=Genovese2007,ISI_record=s2a(&
         "PT J",&
         "AU Genovese, L",&
         "   Deutsch, T",&
         "   Goedecker, S",&
         "AF Genovese, Luigi",&
         "   Deutsch, Thierry",&
         "   Goedecker, Stefan",&
         "TI Efficient and accurate three-dimensional Poisson solver for surface",&
         "   problems",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD AUG 7",&
         "PY 2007",&
         "VL 127",&
         "IS 5",&
         "AR 054704",&
         "DI ARTN 054704",&
         "UT ISI:000248594200029",&
         "ER"),&
          DOI="10.1063/1.2754685")

    CALL add_reference(key=Minary2003,ISI_record=s2a(&
         "PT J",&
         "AU Minary, P",&
         "   Martyna, GJ",&
         "   Tuckerman, ME",&
         "TI Algorithms and novel applications based on the isokinetic ensemble. I.",&
         "   Biophysical and path integral molecular dynamics",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD FEB 8",&
         "PY 2003",&
         "VL 118",&
         "IS 6",&
         "BP 2510",&
         "EP 2526",&
         "UT ISI:000180564800009",&
         "ER"),&
         DOI="10.1063/1.1534582")

    CALL add_reference(key=Evans1983,ISI_record=s2a(&
         "PT J",&
         "AU EVANS, DJ",&
         "   HOOVER, WG",&
         "   FAILOR, BH",&
         "   MORAN, B",&
         "   LADD, AJC",&
         "TI NON-EQUILIBRIUM MOLECULAR-DYNAMICS VIA GAUSS PRINCIPLE OF LEAST",&
         "   CONSTRAINT",&
         "SO PHYSICAL REVIEW A",&
         "SN 1050-2947",&
         "PY 1983",&
         "VL 28",&
         "IS 2",&
         "BP 1016",&
         "EP 1021",&
         "UT ISI:A1983RB52800057",&
         "ER"),&
         DOI="10.1103/PhysRevA.28.1016")

    CALL add_reference(key=Byrd1995,ISI_record=s2a(&
         "PT J",&
         "TI A LIMITED MEMORY ALGORITHM FOR BOUND CONSTRAINED OPTIMIZATION",&
         "AU BYRD, RH",&
         "   LU, PH",&
         "   NOCEDAL, J",&
         "   ZHU, CY",&
         "SO SIAM JOURNAL ON SCIENTIFIC COMPUTING",&
         "VL 16",&
         "IS 5",&
         "BP 1190",&
         "EP 1208",&
         "PY 1995",&
         "TC 214",&
         "UT WOS:A1995RR54100011",&
         "ER"),&
         DOI="10.1137/0916069")

    CALL add_reference(key=VandeVondele2007,ISI_record=s2a(&
         "PT J",&
         "AU VandeVondele, J",&
         "   Hutter, J",&
         "AF VandeVondele, Joost",&
         "   Hutter, Juerg",&
         "TI Gaussian basis sets for accurate calculations on molecular systems in",&
         "   gas and condensed phases",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD SEP 21",&
         "PY 2007",&
         "VL 127",&
         "IS 11",&
         "AR 114105",&
         "DI ARTN 114105",&
         "UT ISI:000249667400011",&
         "ER"),&
         DOI="10.1063/1.2770708")

    CALL add_reference(key=Ortiz1994,ISI_record=s2a(&
         "PT J",&
         "AU ORTIZ, G",&
         "   BALLONE, P",&
         "TI CORRELATION-ENERGY, STRUCTURE FACTOR, RADIAL-DISTRIBUTION FUNCTION, AND",&
         "   MOMENTUM DISTRIBUTION OF THE SPIN-POLARIZED UNIFORM ELECTRON-GAS",&
         "SO PHYSICAL REVIEW B",&
         "SN 0163-1829",&
         "PD JUL 15",&
         "PY 1994",&
         "VL 50",&
         "IS 3",&
         "BP 1391",&
         "EP 1405",&
         "UT ISI:A1994NZ15300007",&
         "ER"),&
         DOI="10.1103/PhysRevB.50.1391")

    CALL add_reference(key=Becke1988,ISI_record=s2a(&
         "PT J",&
         "AU BECKE, AD",&
         "TI DENSITY-FUNCTIONAL EXCHANGE-ENERGY APPROXIMATION WITH CORRECT",&
         "   ASYMPTOTIC-BEHAVIOR",&
         "SO PHYSICAL REVIEW A",&
         "SN 1050-2947",&
         "PD SEP 15",&
         "PY 1988",&
         "VL 38",&
         "IS 6",&
         "BP 3098",&
         "EP 3100",&
         "UT ISI:A1988Q146900044",&
         "ER"),&
         DOI="10.1103/PhysRevA.38.3098")

    CALL add_reference(key=Perdew1996,ISI_record=s2a(&
         "PT J",&
         "AU Perdew, JP",&
         "   Burke, K",&
         "   Ernzerhof, M",&
         "TI Generalized gradient approximation made simple",&
         "SO PHYSICAL REVIEW LETTERS",&
         "SN 0031-9007",&
         "PD OCT 28",&
         "PY 1996",&
         "VL 77",&
         "IS 18",&
         "BP 3865",&
         "EP 3868",&
         "UT ISI:A1996VP22500044",&
         "ER"),&
         DOI="10.1103/PhysRevLett.77.3865")

    CALL add_reference(key=Zhang1998,ISI_record=s2a(&
         "PT J",&
         "AU Zhang, YK",&
         "   Yang, WT",&
         "TI Comment on Generalized gradient approximation made simple",&
         "SO PHYSICAL REVIEW LETTERS",&
         "SN 0031-9007",&
         "PD JAN 26",&
         "PY 1998",&
         "VL 80",&
         "IS 4",&
         "BP 890",&
         "EP 890",&
         "UT ISI:000071717100066",&
         "ER"),&
         DOI="10.1103/PhysRevLett.80.890")

    CALL add_reference(key=Perdew2008,ISI_record=s2a(&
         "PT J",&
         "AU Perdew, JP",&
         "   Ruzsinszky, A",&
         "   Csonka, GI",&
         "   Vydrov, OA",&
         "   Scuseria, GE",&
         "   Constantin, LA",&
         "   Zhou, X",&
         "   Burke, K",&
         "TI Restoring the Density-Gradient Expansion for Exchange in Solids and Surfaces",&
         "SO PHYSICAL REVIEW LETTERS",&
         "SN ",&
         "PD APR 04",&
         "PY 2008",&
         "VL 100",&
         "IS 13",&
         "BP 136406",&
         "EP 136409",&
         "UT ",&
         "ER"),&
         DOI="10.1103/PhysRevLett.100.136406")

    CALL add_reference(key=Lee1988,ISI_record=s2a(&
         "PT J",&
         "AU LEE, CT",&
         "   YANG, WT",&
         "   PARR, RG",&
         "TI DEVELOPMENT OF THE COLLE-SALVETTI CORRELATION-ENERGY FORMULA INTO A",&
         "   FUNCTIONAL OF THE ELECTRON-DENSITY",&
         "SO PHYSICAL REVIEW B",&
         "SN 0163-1829",&
         "PD JAN 15",&
         "PY 1988",&
         "VL 37",&
         "IS 2",&
         "BP 785",&
         "EP 789",&
         "UT ISI:A1988L976200011",&
         "ER"),&
         DOI="10.1103/PhysRevB.37.785")

    CALL add_reference(key=Heyd2004,ISI_record=s2a(&
         "PT J",&
         "AU Heyd, J",&
         "   Scuseria, GE",&
         "TI Assessment and validation of a screened Coulomb hybrid density",&
         "   functional",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD APR 22",&
         "PY 2004",&
         "VL 120",&
         "IS 16",&
         "BP 7274",&
         "EP 7280",&
         "UT ISI:000220676000005",&
         "ER"),&
         DOI="10.1063/1.1668634")

    CALL add_reference(key=Heyd2003,ISI_record=s2a(&
         "PT J",&
         "AU Heyd, J",&
         "   Scuseria, GE",&
         "   Ernzerhof, M",&
         "TI Hybrid functionals based on a screened Coulomb potential",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD MAY 8",&
         "PY 2003",&
         "VL 118",&
         "IS 18",&
         "BP 8207",&
         "EP 8215",&
         "UT ISI:000182454100010",&
         "ER"),&
         DOI="10.1063/1.1564060")

    CALL add_reference(key=Heyd2006,ISI_record=s2a(&
         "PT J",&
         "AU Heyd, J",&
         "   Scuseria, GE",&
         "   Ernzerhof, M",&
         "TI Hybrid functionals based on a screened Coulomb potential (vol 118, pg",&
         "   8207, 2003)",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD JUN 7",&
         "PY 2006",&
         "VL 124",&
         "IS 21",&
         "AR 219906",&
         "DI ARTN 219906",&
         "UT ISI:000238758700058",&
         "ER"),&
         DOI="10.1063/1.2204597")

    CALL add_reference(key=Vydrov2006,ISI_record=s2a(&
         "PT J",&
         "AU Vydrov, OA",&
         "   Heyd, J",&
         "   Krukau, AV",&
         "   Scuseria, GE",&
         "AF Vydrov, Oleg A.",&
         "   Heyd, Jochen",&
         "   Krukau, Aliaksandr V.",&
         "   Scuseria, Gustavo E.",&
         "TI Importance of short-range versus long-range Hartree-Fock exchange for",&
         "   the performance of hybrid density functionals",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD AUG 21",&
         "PY 2006",&
         "VL 125",&
         "IS 7",&
         "AR 074106",&
         "DI ARTN 074106",&
         "UT ISI:000239914800006",&
         "ER"),&
         DOI="10.1063/1.2244560")

    CALL add_reference(key=Vosko1980,ISI_record=s2a(&
         "PT J",&
         "AU VOSKO, SH",&
         "   WILK, L",&
         "   NUSAIR, M",&
         "TI ACCURATE SPIN-DEPENDENT ELECTRON LIQUID CORRELATION ENERGIES FOR LOCAL",&
         "   SPIN-DENSITY CALCULATIONS - A CRITICAL ANALYSIS",&
         "SO CANADIAN JOURNAL OF PHYSICS",&
         "SN 0008-4204",&
         "PY 1980",&
         "VL 58",&
         "IS 8",&
         "BP 1200",&
         "EP 1211",&
         "UT ISI:A1980KE76300015",&
         "ER"),&
         DOI="10.1139/p80-159") ! this doi should be correct later this week, but does not yet resolve

    CALL add_reference(key=Essmann1995,ISI_record=s2a(&
         "PT J",&
         "AU ESSMANN, U",&
         "   PERERA, L",&
         "   BERKOWITZ, ML",&
         "   DARDEN, T",&
         "   LEE, H",&
         "   PEDERSEN, LG",&
         "TI A SMOOTH PARTICLE MESH EWALD METHOD",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD NOV 15",&
         "PY 1995",&
         "VL 103",&
         "IS 19",&
         "BP 8577",&
         "EP 8593",&
         "UT ISI:A1995TE36400026",&
         "ER"),&
         DOI="10.1063/1.470117")

    CALL add_reference(key=Ewald1921,ISI_record=s2a(&
         "PT J",&
         "AU Ewald, PP",&
         "TI Die Berechnung optischer und elektrostatischer Gitterpotentiale",&
         "SO ANNALEN DER PHYSIK",&
         "SN 0003-3804",&
         "PD FEB",&
         "PY 1921",&
         "VL 369",&
         "IS 3",&
         "BP 253",&
         "EP 287",&
         "UT ISI:000201909400003",&
         "ER"),&
         DOI="10.1002/andp.19213690304")

    CALL add_reference(key=Darden1993,ISI_record=s2a(&
         "PT J",&
         "AU DARDEN, T",&
         "   YORK, D",&
         "   PEDERSEN, L",&
         "TI PARTICLE MESH EWALD - AN N.LOG(N) METHOD FOR EWALD SUMS IN LARGE SYSTEMS",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD JUN 15",&
         "PY 1993",&
         "VL 98",&
         "IS 12",&
         "BP 10089",&
         "EP 10092",&
         "UT ISI:A1993LG10100091",&
         "ER"),&
         DOI="10.1063/1.464397")

    CALL add_reference(key=Dudarev1997,ISI_record=s2a(&
         "PT J",&
         "AU Dudarev, SL",&
         "   Manh, DN",&
         "   Sutton, AP",&
         "TI Effect of Mott-Hubbard correlations on the electronic",&
         "   structure and structural stability of uranium dioxide",&
         "SO PHILOSOPHICAL MAGAZINE B",&
         "SN 0141-8637",&
         "PD MAY",&
         "PY 1997",&
         "VL 75",&
         "IS 5",&
         "BP 613",&
         "EP 628",&
         "UT ISI:A1997WX94300001",&
         "ER"),&
         DOI="10.1080/13642819708202343")

    CALL add_reference(key=Dudarev1998,ISI_record=s2a(&
         "PT J",&
         "AU Dudarev, SL",&
         "   Botton, GA",&
         "   Savrasov, SY",&
         "   Humphreys, CJ",&
         "   Sutton, AP",&
         "TI Electron-energy-loss spectra and the structural stability of",&
         "   nickel oxide: An LSDA+U study",&
         "SO PHYSICAL REVIEW B",&
         "SN 0163-1829",&
         "PD JAN 15",&
         "PY 1998",&
         "VL 57",&
         "IS 3",&
         "BP 1505",&
         "EP 1509",&
         "UT ISI:000071716800040",&
         "ER"),&
         DOI="10.1103/PhysRevB.57.1505")

    CALL add_reference(key=Hunt2003,ISI_record=s2a(&
         "PT J",&
         "AU Hunt, P",&
         "   Sprik, M",&
         "   Vuilleumier, R",&
         "TI Thermal versus electronic broadening in the density of states of liquid",&
         "   water",&
         "SO CHEMICAL PHYSICS LETTERS",&
         "SN 0009-2614",&
         "PD JUL 17",&
         "PY 2003",&
         "VL 376",&
         "IS 1-2",&
         "BP 68",&
         "EP 74",&
         "UT ISI:000184315600011",&
         "ER"),&
         DOI="10.1016/S0009-2614(03)00954-0")

    CALL add_reference(key=Guidon2008,ISI_record=s2a(&
         "PT J",&
         "AU Guidon, M",&
         "   Schiffmann, F",&
         "   Hutter, J",&
         "   VandeVondele, J",&
         "AF Guidon, Manuel",&
         "   Schiffmann, Florian",&
         "   Hutter, Juerg",&
         "   VandeVondele, Joost",&
         "TI Ab initio molecular dynamics using hybrid density functionals",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD JUN 7",&
         "PY 2008",&
         "VL 128",&
         "IS 21",&
         "AR 214104",&
         "DI ARTN 214104",&
         "UT ISI:000256527500005",&
         "ER"),&
         DOI="10.1063/1.2931945")

    CALL add_reference(key=Stewart1982,ISI_record=s2a(&
         "PT J",&
         "AU Stewart, JJP",&
         "   Csaszar, P",&
         "   Pulay, P",&
         "TI Fast semi-empirical calculations",&
         "SO JOURNAL OF COMPUTATIONAL CHEMISTRY",&
         "SN 0192-8651",&
         "PY 1982",&
         "VL 3",&
         "IS 2",&
         "BP 227",&
         "EP 228",&
         "UT ISI:A1982NN22800013",&
         "ER"),&
         DOI="10.1002/jcc.540030214")

    CALL add_reference(key=Tao2003,ISI_record=s2a(&
         "PT J",&
         "AU Tao, JM",&
         "   Perdew, JP",&
         "   Staroverov, VN",&
         "   Scuseria, GE",&
         "TI Climbing the density functional ladder: Nonempirical meta-generalized",&
         "   gradient approximation designed for molecules and solids",&
         "SO PHYSICAL REVIEW LETTERS",&
         "PD OCT 3",&
         "PY 2003",&
         "VL 91",&
         "IS 14",&
         "AR 146401",&
         "DI 10.1103/PhysRevLett.91.146401",&
         "UT ISI:000185719500032",&
         "ER"),&
         DOI="10.1103/PhysRevLett.91.146401")

    CALL add_reference(key=VandeVondele2006,ISI_record=s2a(&
         "PT C",&
         "AU VandeVondele, J",&
         "   Iannuzzi, M",&
         "   Hutter, J",&
         "ED Ferrario, M; Ciccotti, G; Binder, K",&
         "TI Large scale condensed matter calculations using the gaussian and",&
         "   augmented plane waves method",&
         "SO Computer Simulations in Condensed Matter Systems: From Materials to",&
         "   Chemical Biology, Vol 1",&
         "SE LECTURE NOTES IN PHYSICS",&
         "LA English",&
         "DT Proceedings Paper",&
         "CT Conference on Computer Simulations in Condensed Matter Systems",&
         "CY JUL, 2005",&
         "CL Erice, ITALY",&
         "SP Ctr Sci Culture",&
         "HO Ettore Majorana Fdn",&
         "NR 35",&
         "TC 2",&
         "PU SPRINGER",&
         "PI NEW YORK",&
         "PA 233 SPRING STREET, NEW YORK, NY 10013, UNITED STATES",&
         "SN 0075-8450",&
         "BN 978-3-540-35270-9",&
         "J9 LECT NOTES PHYS",&
         "PY 2006",&
         "VL 703",&
         "BP 287",&
         "EP 314",&
         "DI 10.1007/3-540-35273-2_8",&
         "PG 28",&
         "SC Physics, Multidisciplinary",&
         "GA BFQ13",&
         "UT ISI:000243725000009",&
         "ER"),&
         DOI="10.1007/3-540-35273-2_8")

    CALL add_reference(key=Grimme2006,ISI_record=s2a(&
         "PT J",&
         "AU Grimme, S",&
         "AF Grimme, Stefan",&
         "TI Semiempirical GGA-type density functional constructed with a long-range",&
         "   dispersion correction",&
         "SO JOURNAL OF COMPUTATIONAL CHEMISTRY",&
         "SN 0192-8651",&
         "PD NOV 30",&
         "PY 2006",&
         "VL 27",&
         "IS 15",&
         "BP 1787",&
         "EP 1799",&
         "DI 10.1002/jcc.20495",&
         "UT ISI:000241477200003",&
         "ER"),&
         DOI="10.1002/jcc.20495")

    CALL add_reference(key=Grimme2010,ISI_record=s2a(&
         "PT J",&
         "AU Grimme, S",&
         "   Antony, J",&
         "   Ehrlich, S",&
         "   Krieg, H",&
         "AF Grimme, Stefan",&
         "   Antony, Jens",&
         "   Ehrlich, Stephan",&
         "   Krieg, Helge",&
         "TI A consistent and accurate ab initio parametrization of density",&
         "   functional dispersion correction (DFT-D) for the 94 elements H-Pu",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD APR 21",&
         "PY 2010",&
         "VL 132",&
         "IS 15",&
         "AR 154104",&
         "DI 10.1063/1.3382344",&
         "UT ISI:000276971500005",&
         "ER"),&
         DOI="10.1063/1.3382344")

    CALL add_reference(key=Grimme2011,ISI_record=s2a(&
         "PT J",&
         "AU Grimme, S",&
         "   Ehrlich, S",&
         "   Goerigk, L",&
         "AF Grimme, Stefan",&
         "   Ehrlich, Stephan",&
         "   Goerigk, Lars",&
         "TI Effect of the damping function in dispersion corrected density functional theory",&
         "SO JOURNAL OF COMPUTATIONAL CHEMISTRY",&
         "PY 2011",&
         "VL 32",&
         "BP 1456",&
         "DI 10.1002/jcc.21759",&
         "ER"),&
         DOI="10.1002/jcc.21759")

    CALL add_reference(key=Branduardi2007,ISI_record=s2a(&
         "PT J",&
         "AU Branduardi, D",&
         "   Gervasio, FL",&
         "   Parrinello, M",&
         "AF Branduardi, Davide",&
         "   Gervasio, Francesco Luigi",&
         "   Parrinello, Michele",&
         "TI From A to B in free energy space",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD FEB 7",&
         "PY 2007",&
         "VL 126",&
         "IS 5",&
         "AR 054103",&
         "DI 10.1063/1.2432340",&
         "UT ISI:000244044600009",&
         "ER"),&
         DOI="10.1063/1.2432340")

    CALL add_reference(key=Schenter2008,ISI_record=s2a(&
         "PT J",&
         "AU Chang, DT",&
         "   Schenter, GK",&
         "   Garrett, BC ",&
         "TI Self-consistent polarization neglect of diatomic differential overlap: Applications to water clusters",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD APR 28",&
         "PY 2008",&
         "VL 128",&
         "IS 16",&
         "AR 164111",&
         "DI 10.1063/1.2905230",&
         "ER"),&
         DOI="10.1063/1.2905230")

    CALL add_reference(key=Proynov2007, ISI_record=s2a(&
         "AU Proynov, E",&
         "   Gan, Z",&
         "   Kong, J",&
         "AF Proynov, Emil",&
         "   Gan, Zhenting",&
         "   Kong, Jing",&
         "TI Analytical representation of the Becke-Roussel exchange functional",&
         "SO CHEMICAL PHYSICS LETTERS",&
         "SN 0009-2614",&
         "PD MAR 31",&
         "PY 2008",&
         "VL 455",&
         "IS 1-3",&
         "BP 103",&
         "EP 109",&
         "DI 10.1016/j.cplett.2008.02.039",&
         "UT ISI:000254425400020",&
         "ER"),&
         DOI="10.1016/j.cplett.2008.02.039")

    CALL add_reference(key=BeckeRoussel1989, ISI_record=s2a(&
         "AU BECKE, AD",&
         "   ROUSSEL, MR",&
         "TI EXCHANGE HOLES IN INHOMOGENEOUS SYSTEMS - A COORDINATE-SPACE MODEL",&
         "SO PHYSICAL REVIEW A",&
         "DT Article",&
         "RP BECKE, AD, QUEENS UNIV,DEPT CHEM,KINGSTON K7L 3N6,ONTARIO,CANADA.",&
         "NR 17",&
         "TC 114",&
         "PU AMERICAN PHYSICAL SOC",&
         "PI COLLEGE PK",&
         "PA ONE PHYSICS ELLIPSE, COLLEGE PK, MD 20740-3844 USA",&
         "SN 1050-2947",&
         "J9 PHYS REV A",&
         "JI Phys. Rev. A",&
         "PD APR 15",&
         "PY 1989",&
         "VL 39",&
         "IS 8",&
         "BP 3761",&
         "EP 3767",&
         "PG 7",&
         "SC Optics; Physics, Atomic, Molecular & Chemical",&
         "GA U3175",&
         "UT ISI:A1989U317500001",&
         "ER"),&
         DOI="10.1103/PhysRevA.39.3761")

    CALL add_reference(key=Becke1997, ISI_record=s2a(&
         "AU Becke, AD", &
         "PY 1997", &
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606", &
         "TI Density-functional thermochemistry .", &
         "   5. Systematic optimization of exchange-correlation functionals", &
         "VL 107", &
         "IS 20", &
         "SP 8554", &
         "EP 8560"),&
         DOI="10.1063/1.475007")

    CALL add_reference(key=Ricci2003, ISI_record=s2a(&
         "AU Ricci, A",&
         "   Ciccotti, G",&
         "TI Algorithms for Brownian dynamics",&
         "SO MOLECULAR PHYSICS",&
         "SN 0026-8976",&
         "J9 MOL PHYS",&
         "JI Mol. Phys.",&
         "PD JUN 20",&
         "PY 2003",&
         "VL 101",&
         "IS 12",&
         "BP 1927",&
         "EP 1931",&
         "UT ISI:000184158800014",&
         "ER"),&
         DOI="10.1080/0026897031000108113")

    CALL add_reference(key=Kolafa2004, ISI_record=s2a(&
         "AU Kolafa, J",&
         "TI Time-reversible always stable predictor-corrector method for molecular dynamics of polarizable molecules",&
         "SO JOURNAL OF COMPUTATIONAL CHEMISTRY",&
         "SN 0192-8651",&
         "J9 J COMPUT CHEM",&
         "JI J. Comput. Chem.",&
         "PD FEB",&
         "PY 2004",&
         "VL 25",&
         "IS 3",&
         "BP 335",&
         "EP 342",&
         "UT ISI:000188245900004",&
         "ER"),&
         DOI="10.1002/jcc.10385")

    CALL add_reference(key=Kuhne2007, ISI_record=s2a(&
         "AU Kuhne, TD",&
         "   Krack, M",&
         "   Mohamed, FR",&
         "   Parrinello, M",&
         "TI Efficient and accurate Car-Parrinello-like approach to Born-Oppenheimer molecular dynamics",&
         "SO PHYSICAL REVIEW LETTERS",&
         "SN 0031-9007",&
         "J9 PHYS REV LETT",&
         "JI Phys. Rev. Lett.",&
         "PD FEB 9",&
         "PY 2007",&
         "VL 98",&
         "IS 6",&
         "AR 066401",&
         "PG 4",&
         "UT ISI:000244117400046",&
         "ER"),&
         DOI="10.1103/PhysRevLett.98.066401")

    CALL add_reference(key=Kunert2003, ISI_record=s2a(&
         "AU Kunert, T",&
         "   Schmidt, R",&
         "TI Non-adiabatic quantum molecular dynamics: General formalism and case",&
         "   study H-2(+) in strong laser fields",&
         "SO EUROPEAN PHYSICAL JOURNAL D",&
         "SN 1434-6060",&
         "J9 EUR PHYS J D",&
         "JI Eur. Phys. J. D",&
         "PD JUL",&
         "PY 2003",&
         "VL 25",&
         "IS 1",&
         "PG 10",&
         "UT ISI:000184771500002",&
         "ER"),&
         DOI="10.1140/epjd/e2003-00086-8")

    CALL add_reference(key=Ceriotti2009, ISI_record=s2a(&
         "AU Ceriotti, M",&
         "   Bussi, G",&
         "   Parrinello, M",&
         "TI Langevin equation with colored noise for constant-temperature molecular dynamics simulations",&
         "SO PHYSICAL REVIEW LETTERS",&
         "SN 0031-9007",&
         "J9 PHYS REV LETT",&
         "JI Phys. Rev. Lett.",&
         "PD JAN 16",&
         "PY 2009",&
         "VL 102",&
         "IS 2",&
         "AR 020601",&
         "PG 4",&
         "UT ISI:000244117400046",&
         "ER"),&
         DOI="10.1103/PhysRevLett.102.020601")

   CALL add_reference(key=Ceriotti2009b, ISI_record=s2a(&
          "AU Ceriotti, M",&
          "   Bussi, G",&
          "   Parrinello, M",&
          "AF Ceriotti, Michele",&
          "   Bussi, Giovanni",&
          "   Parrinello, Michele",&
          "TI Nuclear Quantum Effects in Solids Using a Colored-Noise Thermostat",&
          "SO PHYSICAL REVIEW LETTERS",&
          "LA English",&
          "J9 PHYS REV LETT",&
          "JI Phys. Rev. Lett.",&
          "PD JUL 17",&
          "PY 2009",&
          "VL 103",&
          "IS 3",&
          "AR 030603",&
          "PG 4",&
          "ER"),&
          DOI="10.1103/PhysRevLett.103.030603")

    CALL add_reference(key=Guidon2009, ISI_record=s2a(&
         "AU Guidon, M",&
         "   Hutter, J",&
         "   VandeVondele, J",&
         "AF Guidon, Manuel",&
         "   Hutter, Juerg",&
         "   VandeVondele, Joost",&
         "TI Robust Periodic Hartree-Fock Exchange for Large-Scale Simulations Using",&
         "   Gaussian Basis Sets",&
         "SO JOURNAL OF CHEMICAL THEORY AND COMPUTATION",&
         "SN 1549-9618",&
         "PD NOV",&
         "PY 2009",&
         "VL 5",&
         "IS 11",&
         "BP 3010",&
         "EP 3021",&
         "DI 10.1021/ct900494g",&
         "UT ISI:000271522500011",&
         "ER"),&
          DOI="10.1021/ct900494g")

    CALL add_reference(key=BarducBus2008, ISI_record=s2a(&
         "AU Barducci, A",&
         "   Bussi, G",&
         "   Parrinello, M",&
         "TI Well-Tempered Metadynamics: A Smoothly Converging and Tunable Free-Energy Method",&
         "SO PHYSICAL REVIEW LETTERS",&
         "SN 0031-9007",&
         "J9 PHYS REV LETT",&
         "JI Phys. Rev. Lett.",&
         "PD JAN 18",&
         "PY 2008",&
         "VL 100",&
         "IS 2",&
         "AR 020603",&
         "PG 4",&
         "UT ISI:000252471200008",&
         "ER"),&
         DOI="10.1103/PhysRevLett.100.020603")

    CALL add_reference(key=Guidon2010, ISI_record=s2a(&
         "AU Guidon, M",&
         "   Hutter, J",&
         "   VandeVondele, J",&
         "AF Guidon, Manuel",&
         "   Hutter, Jurg",&
         "   VandeVondele, Joost",&
         "TI Auxiliary Density Matrix Methods for Hartree-Fock Exchange Calculations",&
         "SO JOURNAL OF CHEMICAL THEORY AND COMPUTATION",&
         "LA English",&
         "DT Article",&
         "NR 65",&
         "TC 0",&
         "PU AMER CHEMICAL SOC",&
         "PI WASHINGTON",&
         "PA 1155 16TH ST, NW, WASHINGTON, DC 20036 USA",&
         "SN 1549-9618",&
         "J9 J CHEM THEORY COMPUT",&
         "JI J. Chem. Theory Comput.",&
         "PD AUG",&
         "PY 2010",&
         "VL 6",&
         "IS 8",&
         "BP 2348",&
         "EP 2364",&
         "DI 10.1021/ct1002225",&
         "PG 17",&
         "SC Chemistry, Physical; Physics, Atomic, Molecular & Chemical",&
         "GA 635NE",&
         "UT ISI:000280661600009",&
         "ER"),&
         DOI="10.1021/ct1002225")

    CALL add_reference(key=Marques2012, ISI_record=s2a(&
         "AU Marques, MAL",&
         "   Oliveira, MJT",&
         "   Burnus, T",&
         "AF Marques, Miguel A. L.",&
         "   Oliveira, Micael J. T.",&
         "   Burnus, Tobias",&
         "TI LIBXC: A library of exchange and correlation functionals for density",&
         "   functional theory",&
         "SO COMPUTER PHYSICS COMMUNICATIONS",&
         "SN 0010-4655",&
         "PD OCT",&
         "PY 2012",&
         "VL 183",&
         "IS 10",&
         "BP 2272",&
         "EP 2281",&
         "DI 10.1016/j.cpc.2012.05.007",&
         "UT WOS:000306771900024",&
         "ER"),&
         DOI="10.1016/j.cpc.2012.05.007")

    CALL add_reference(key=Jones2011, ISI_record=s2a(&
         "AU Jones, Andrew",&
         "   Leimkuhler, Ben",&
         "TI Adaptive stochastic methods for sampling driven molecular systems",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "VL 135",&
         "IS 8",&
         "AR 084125",&
         "DI 10.1063/1.3626941",&
         "PD AUG 28",&
         "PY 2011",&
         "TC 0",&
         "SN 0021-9606",&
         "UT WOS:000294484700031",&
         "ER"),&
         DOI="10.1063/1.3626941")

     CALL add_reference(key=Bernstein2012, ISI_record=s2a(&
         "AU Bernstein, Noam",&
         "   Varnai, Csilla",&
         "   Solt, Ivan",&
         "   Winfield, Steven A",&
         "   Payne, Mike C",&
         "   Simon, Istvan",&
         "   Fuxreiter, Monika",&
         "   Csanyi, Gabor",&
         "TI QM/MM simulation of liquid water with an adaptive quantum region",&
         "SO PHYSICAL CHEMISTRY CHEMICAL PHYSICS",&
         "VL 14",&
         "IS 2",&
         "BP 646",&
         "EP 656",&
         "DI 10.1039/c1cp22600b",&
         "PY 2012",&
         "TC 0",&
         "ZB 0",&
         "Z8 0",&
         "Z9 0",&
         "SN 1463-9076",&
         "UT WOS:000298552800027",&
         "ER"),&
         DOI="10.1039/c1cp22600b")

     CALL add_reference(key=Bernstein2009, ISI_record=s2a(&
         "AU Bernstein, N",&
         "   Kermode, J R",&
         "   Csanyi, G",&
         "TI Hybrid atomistic simulation methods for materials systems",&
         "SO REPORTS ON PROGRESS IN PHYSICS",&
         "VL 72",&
         "IS 2",&
         "AR 026501",&
         "DI 10.1088/0034-4885/72/2/026501",&
         "PD FEB 1",&
         "PY 2009",&
         "TC 25",&
         "ZB 0",&
         "Z8 0",&
         "Z9 25",&
         "SN 0034-4885",&
         "UT WOS:000262785900001",&
         "ER"),&
         DOI="10.1088/0034-4885/72/2/026501")

    CALL add_reference(key=Dick1958, ISI_record=s2a(&
         "PT J",&
         "AU Dick, BG",&
         "   Overhauser, AW",&
         "TI Theory of the Dielectric Constants of Alkali Halide Crystals",&
         "SO Phys. Rev.",&
         "VL 112",&
         "IS 1",&
         "SP 90",&
         "EP 103",&
         "PD OCT 1",&
         "PY 1958",&
         "UR http://link.aps.org/doi/10.1103/PhysRev.112.90",&
         "ER"),&
         DOI="10.1103/PhysRev.112.90")

    CALL add_reference(key=Mitchell1993, ISI_record=s2a(&
         "PT J",&
         "AU Mitchell, PJ",&
         "   Fincham D",&
         "TI Shell model simulations by adiabatic dynamics",&
         "SO J. Phys.: Condens. Matter",&
         "VL 5",&
         "IS 8",&
         "SP 1031",&
         "EP 1038",&
         "PD FEB 22",&
         "PY 1993",&
         "UR http://stacks.iop.org/0953-8984/5/i=8/a=006",&
         "ER"),&
         DOI="10.1088/0953-8984/5/8/006")

    CALL add_reference(key=Devynck2012, ISI_record=s2a(&
         "PT J",&
         "AU Devynck, F",&
         "   Iannuzzi, M",&
         "   Krack, M",&
         "AF Devynck, Fabien",&
         "   Iannuzzi, Marcella",&
         "   Krack, Matthias",&
         "TI Frenkel pair recombinations in UO2: Importance of explicit description of",&
         "   polarizability in core-shell molecular dynamics simulations",&
         "SO Phys. Rev. B",&
         "VL 85",&
         "IS 18",&
         "SP 184103",&
         "EP",&
         "PD MAY 15",&
         "PY 2012",&
         "UR http://link.aps.org/doi/10.1103/PhysRevB.85.184103",&
         "ER"),&
         DOI="10.1103/PhysRevB.85.184103")

     CALL add_reference(key=VandeVondele2012, ISI_record=s2a(&
          "PT J",&
          "AU VandeVondele, J",&
          "   Borstnik, U",&
          "   Hutter, J",&
          "AF VandeVondele, Joost",&
          "   Borstnik, Urban",&
          "   Hutter, Juerg",&
          "TI Linear Scaling Self-Consistent Field Calculations with Millions of Atoms",&
          "   in the Condensed Phase",&
          "SO JOURNAL OF CHEMICAL THEORY AND COMPUTATION",&
          "SN 1549-9618",&
          "PD OCT",&
          "PY 2012",&
          "VL 8",&
          "IS 10",&
          "BP 3565",&
          "EP 3573",&
          "DI 10.1021/ct200897x",&
          "UT WOS:000309505600014",&
          "ER"),&
          DOI="10.1021/ct200897x")

      CALL add_reference(key=Niklasson2003, ISI_record=s2a(&
           "PT J",&
           "AU Niklasson, AMN",&
           "   Tymczak, CJ",&
           "   Challacombe, M",&
           "AF Niklasson, AMN",&
           "   Tymczak, CJ",&
           "   Challacombe, M",&
           "TI Trace resetting density matrix purification in O(N)",&
           "   self-consistent-field theory",&
           "SO JOURNAL OF CHEMICAL PHYSICS",&
           "SN 0021-9606",&
           "PD MAY 15",&
           "PY 2003",&
           "VL 118",&
           "IS 19",&
           "BP 8611",&
           "EP 8620",&
           "DI 10.1063/1.1559913",&
           "UT WOS:000182603000010",&
           "ER"),&
           DOI="10.1063/1.1559913")

      CALL add_reference(key=Shao2003, ISI_record=s2a(&
            "PT J",&
            "AU Shao, Y",&
            "   Saravanan, C",&
            "   Head-Gordon, M",&
            "   White, CA",&
            "AF Shao, Y",&
            "   Saravanan, C",&
            "   Head-Gordon, M",&
            "   White, CA",&
            "TI Curvy steps for density matrix-based energy minimization: Application to",&
            "   large-scale self-consistent-field calculations",&
            "SO JOURNAL OF CHEMICAL PHYSICS",&
            "RI White, Christopher/B-7390-2009",&
            "SN 0021-9606",&
            "PD APR 8",&
            "PY 2003",&
            "VL 118",&
            "IS 14",&
            "BP 6144",&
            "EP 6151",&
            "DI 10.1063/1.1558476",&
            "UT WOS:000181834900004",&
            "ER"),&
          DOI="10.1063/1.1558476")

      CALL add_reference(key=VandeVondele2002, ISI_record=s2a(&
            "PT J",&
            "AU VandeVondele, J",&
            "   Rothlisberger, U",&
            "AF VandeVondele, J",&
            "   Rothlisberger, U",&
            "TI Canonical adiabatic free energy sampling (CAFES): A novel method for the",&
            "   exploration of free energy surfaces",&
            "SO JOURNAL OF PHYSICAL CHEMISTRY B",&
            "SN 1520-6106",&
            "PD JAN 10",&
            "PY 2002",&
            "VL 106",&
            "IS 1",&
            "BP 203",&
            "EP 208",&
            "DI 10.1021/jp013346k",&
            "UT WOS:000173256700030",&
            "ER"),&
            DOI="10.1021/jp013346k")

    CALL add_reference(key=Dion2004, ISI_record=s2a(&
         "PT J",&
         "AU Dion, M",&
         "   Rydberg, H",&
         "   Schroder, E",&
         "   Langreth, DC",&
         "   Lundqvist, BI",&
         "TI Van der Waals density functional for general geometries",&
         "SO Phys. Rev. Lett.",&
         "PD JUN 18",&
         "PY 2004",&
         "VL 92",&
         "IS 24",&
         "AR 246401",&
         "UR http://prl.aps.org/abstract/PRL/v92/i24/e246401",&
         "ER"),&
         DOI="10.1103/PhysRevLett.92.246401")

    CALL add_reference(key=RomanPerez2009, ISI_record=s2a(&
         "PT J",&
         "AU Roman-Perez, G",&
         "   Soler, JM",&
         "TI Efficient Implementation of a van der Waals Density Functional:",&
         "   Application to Double-Wall Carbon Nanotubes",&
         "SO Phys. Rev. Lett.",&
         "PD AUG 28",&
         "PY 2009",&
         "VL 103",&
         "IS 9",&
         "AR 096102",&
         "UR http://prl.aps.org/abstract/PRL/v103/i9/e096102",&
         "ER"),&
         DOI="10.1103/PhysRevLett.103.096102")

     CALL add_reference(key=DelBen2012, ISI_record=s2a(&
         "PT J",&
         "AU Del Ben, M",&
         "   Hutter, J",&
         "   VandeVondele, J",&
         "AF Del Ben, Mauro",&
         "   Hutter, Juerg",&
         "   VandeVondele, Joost",&
         "TI Second-Order Moller-Plesset Perturbation Theory in the Condensed Phase:",&
         "   An Efficient and Massively Parallel Gaussian and Plane Waves Approach",&
         "SO JOURNAL OF CHEMICAL THEORY AND COMPUTATION",&
         "SN 1549-9618",&
         "PD NOV",&
         "PY 2012",&
         "VL 8",&
         "IS 11",&
         "BP 4177",&
         "EP 4188",&
         "DI 10.1021/ct300531w",&
         "UT WOS:000311191900019",&
         "ER"),&
         DOI="10.1021/ct300531w")

    CALL add_reference(key=Sabatini2013, ISI_record=s2a(&
         "PT J",&
         "AU Sabatini, R",&
         "   Gorni, T",&
         "   de Gironcoli, S",&
         "AF Sabatini, Riccardo",&
         "   Gorni, Tommaso",&
         "   de Gironcoli, Stefano",&
         "TI Nonlocal van der Waals density functional made simple and efficient",&
         "SO Phys. Rev. B",&
         "PD JAN 15",&
         "PY 2013",&
         "VL 87",&
         "DI 10.1103/PhysRevB.87.041108",&
         "IS 4",&
         "AR 041108(R)",&
         "ER"),&
         DOI="10.1103/PhysRevB.87.041108")

    CALL add_reference(key=Delben2013, ISI_record=s2a(&
         "PT J",&
         "AU Del Ben, M",&
         "   Hutter, J",&
         "   VandeVondele, J",&
         "AF Del Ben, Mauro",&
         "   Hutter, Juerg",&
         "   VandeVondele, Joost",&
         "TI Electron Correlation in the Condensed Phase from a Resolution of",&
         "   Identity Approach Based on the Gaussian and Plane Waves Scheme",&
         "SO JOURNAL OF CHEMICAL THEORY AND COMPUTATION",&
         "SN 1549-9618",&
         "PD JUN",&
         "PY 2013",&
         "VL 9",&
         "IS 6",&
         "BP 2654",&
         "EP 2671",&
         "DI 10.1021/ct4002202",&
         "UT WOS:000320484500014",&
         "ER"),&
         DOI="10.1021/ct4002202")

    CALL add_reference(key=Kikuchi2009, ISI_record=s2a(&
         "PT J",&
         "AU Kikuchi, Y",&
         "   Imamura, Y",&
         "   Nakai, H",&
         "AF Kikuchi, Yasuaki",&
         "   Imamura, Yutaka",&
         "   Nakai, Hiromi",&
         "TI One-Body Energy Decomposition Schemes Revisited: Assessment of",&
         "   Mulliken-, Grid-, and Conventional Energy Density Analyses",&
         "SO INTERNATIONAL JOURNAL OF QUANTUM CHEMISTRY",&
         "TC 3",&
         "Z9 3",&
         "SN 0020-7608",&
         "PD SEP",&
         "PY 2009",&
         "VL 109",&
         "IS 11",&
         "BP 2464",&
         "EP 2473",&
         "DI 10.1002/qua.22017",&
         "ER"),&
         DOI="10.1002/qua.22017")

    CALL add_reference(key=Putrino2000, ISI_record=s2a(&
         "PT J",&
         "AU Putrino, A",&
         "   Sebastiani, D",&
         "   Parrinello, M",&
         "AF Putrino, Anna",&
         "   Sebastiani, Daniel",&
         "   Parrinello, Michele",&
         "TI Generalized Variational Density Functional Perturbation Theory",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD NOV",&
         "PY 2000",&
         "VL 113",&
         "IS 17",&
         "BP 7102",&
         "EP 7109",&
         "DI 10.1063/1.1312830",&
         "ER"),&
         DOI="10.1063/1.1312830")

    CALL add_reference(key=Tran2013, ISI_record=s2a(&
         "PT J",&
         "AU Tran, F",&
         "   Hutter, J",&
         "AF Tran, Fabien",&
         "   Hutter, Juerg",&
         "TI Nonlocal van der Waals functionals: The case of rare-gas dimers and",&
         "   solids",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "TC 2",&
         "Z9 2",&
         "SN 0021-9606",&
         "PD MAY 28",&
         "PY 2013",&
         "VL 138",&
         "IS 20",&
         "AR 204103",&
         "DI 10.1063/1.4807332",&
         "UT WOS:000320131100005",&
         "ER"),&
         DOI="10.1063/1.4807332")

    CALL add_reference(key=Putrino2002, ISI_record=s2a(&
         "PT J",&
         "AU Putrino, A",&
         "   Parrinello, M",&
         "AF Putrino, Anna",&
         "   Parrinello, Michele",&
         "TI Anharmonic Raman Spectra in High-Pressure Ice from Ab Initio Simulations",&
         "SO PHYSICAL REVIEW LETTERS",&
         "SN 0031-9007",&
         "PD APR 29",&
         "PY 2002",&
         "VL 88",&
         "IS 17",&
         "AR 176401",&
         "DI 10.1103/PhysRevLett.88.176401",&
         "ER"),&
         DOI="10.1103/PhysRevLett.88.176401")

    CALL add_reference(key=Sebastiani2001, ISI_record=s2a(&
         "PT J",&
         "AU Sebastiani, D",&
         "   Parrinello, M",&
         "AF Sebastiani, Daniel",&
         "   Parrinello, Michele",&
         "TI A New Ab Initio Approach for NMR Chemical Shifts in Periodic Systems",&
         "SO THE JOURNAL OF PHYSICAL CHEMISTRY A",&
         "PD MAR",&
         "PY 2001",&
         "VL 105",&
         "IS 10",&
         "BP 1951",&
         "EP 1958",&
         "DI 10.1021/jp002807j",&
         "ER"),&
         DOI="10.1021/jp002807j")

    CALL add_reference(key=Weber2009, ISI_record=s2a(&
         "PT J",&
         "AU Weber, V",&
         "   Iannuzzi, M",&
         "   Giani, S",&
         "   Hutter, S",&
         "   Declerck, R",&
         "   Waroduier, M",&
         "AF Weber, Valery",&
         "   Iannuzzi, Marcella",&
         "   Giani, Samuele",&
         "   Hutter, Juerg",&
         "   Declerck, Reinout",&
         "   Waroduier, Michel",&
         "TI Magnetic Linear Response Properties Calculations with the ",&
         "   Gaussian and Augmanted-Plane-Wave Method",&
         "SO THE JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD JUL 7",&
         "PY 2009",&
         "VL 131",&
         "IS 1",&
         "AR 014106",&
         "DI 10.1063/1.3156803",&
         "ER"),&
         DOI="10.1063/1.3156803")

    CALL add_reference(key=Golze2013, ISI_record=s2a(&
         "PT J",&
         "AU Golze, D",&
         "   Iannuzzi, M",&
         "   Nguyen, M-T",&
         "   Passerone, D",&
         "   Hutter, J",&
         "AF Golze, Dorothea",&
         "   Iannuzzi, Marcella",&
         "   Nguyen, Manh-Thuong",&
         "   Passerone, Daniele",&
         "   Hutter, Juerg",&
         "TI Simulation of Adsorption Processes at Metallic Interfaces:",&
         "   An Image Charge Augmented QM/MM Approach",&
         "SO Journal of Chemical Theory and Computation",&
         "SN 1549-9618^",&
         "PD NOV 12",&
         "PY 2013",&
         "VL 9",&
         "IS 11",&
         "BP 5086",&
         "EP 5097",&
         "DI 10.1021/ct400698y",&
         "ER"),&
         DOI="10.1021/ct400698y")

    CALL add_reference(key=Fattebert2002, ISI_record=s2a(&
         "PT J",&
         "AU Fattebert, JL",&
         "   Gygi, F",&
         "AF Fattebert, Jean-Luc",&
         "   Gygi, Francois",&
         "TI Density functional theory for efficient ab initio molecular dynamics",&
         "   simulations in solution",&
         "SO J. Comput. Chem.",&
         "VL 23",&
         "IS 6",&
         "SP 662",&
         "EP 666",&
         "PD MAR 18",&
         "PY 2002",&
         "UR http://dx.doi.org/10.1002/jcc.10069",&
         "ER"),&
         DOI="10.1002/jcc.10069")

    CALL add_reference(key=Andreussi2012, ISI_record=s2a(&
         "PT J",&
         "AU Andreussi, O",&
         "   Dabo, I",&
         "   Marzari, N",&
         "AF Andreussi, Oliviero",&
         "   Dabo, Ismaila",&
         "   Marzari, Nicola",&
         "TI Revised self-consistent continuum solvation in electronic-structure calculations",&
         "SO J. Chem. Phys.",&
         "VL 136",&
         "IS 6",&
         "AR 064102",&
         "PD FEB 8",&
         "PY 2012",&
         "UR http://dx.doi.org/10.1063/1.3676407",&
         "ER"),&
         DOI="10.1063/1.3676407")

    CALL add_reference(key=Tuckerman1992, ISI_record=s2a(&
         "PT J",&
         "AU TUCKERMAN, M",&
         "   BERNE, BJ",&
         "   MARTYNA, GJ",&
         "AF TUCKERMAN, M",&
         "   BERNE, BJ",&
         "   MARTYNA, GJ",&
         "TI REVERSIBLE MULTIPLE TIME SCALE MOLECULAR-DYNAMICS",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD AUG 1",&
         "PY 1992",&
         "VL 97",&
         "IS 3",&
         "BP 1990",&
         "EP 2001",&
         "DI 10.1063/1.463137",&
         "UT WOS:A1992JE89100044",&
         "ER"),&
         DOI="10.1063/1.463137")

    CALL add_reference(key=Goedecker2004, ISI_record=s2a(&
         "PT J",&
         "AU Goedecker, S",&
         "AF Goedecker, Stefan",&
         "TI Minima hopping: An efficient search method for the global minimum",&
         "   of the potential energy surface of complex molecular systems",&
         "SO Journal of Chemical Physics",&
         "PY 2004",&
         "VL 120",&
         "IS 21",&
         "BP 9911",&
         "EP 9917",&
         "DI 10.1063/1.1724816",&
         "ER"),&
         DOI="10.1063/1.1724816")

    CALL add_reference(key=Khaliullin2007, ISI_record=s2a(&
         "PT J",&
         "AU Khaliullin, RZ",&
         "   Cobar, EA",&
         "   Lochan, RC",&
         "   Bell, AT",&
         "   Head-Gordon, M",&
         "AF Khaliullin, Rustam Z.",&
         "   Cobar, Erika A.",&
         "   Lochan, Rohini C.",&
         "   Bell, Alexis T.",&
         "   Head-Gordon, Martin",&
         "TI Unravelling the origin of intermolecular interactions using absolutely localized molecular orbitals",&
         "SO Journal of Physical Chemistry A",&
         "PY 2007",&
         "VL 111",&
         "IS 36",&
         "BP 8753",&
         "EP 8765",&
         "DI 10.1021/jp073685z",&
         "ER"),&
         DOI="10.1021/jp073685z")

    CALL add_reference(key=Khaliullin2008, ISI_record=s2a(&
         "PT J",&
         "AU Khaliullin, RZ",&
         "   Bell, AT",&
         "   Head-Gordon, M",&
         "AF Khaliullin, Rustam Z.",&
         "   Bell, Alexis T.",&
         "   Head-Gordon, Martin",&
         "TI Analysis of charge transfer effects in molecular complexes based on absolutely localized molecular orbitals",&
         "SO Journal of Chemical Physics",&
         "SN 0021-9606",&
         "PY 2008",&
         "VL 128",&
         "IS 18",&
         "AR 184112",&
         "DI 10.1063/1.2912041",&
         "ER"),&
         DOI="10.1063/1.2912041")

    CALL add_reference(key=Khaliullin2013, ISI_record=s2a(&
         "PT J",&
         "AU Khaliullin, RZ",&
         "   VandeVondele, J",&
         "   Hutter, J",&
         "AF Khaliullin, Rustam Z.",&
         "   VandeVondele, Joost",&
         "   Hutter, Juerg",&
         "TI Efficient Linear-Scaling Density Functional Theory for Molecular Systems",&
         "SO JOURNAL OF CHEMICAL THEORY AND COMPUTATION",&
         "SN 1549-9618",&
         "PD OCT",&
         "PY 2013",&
         "VL 9",&
         "IS 10",&
         "BP 4421",&
         "EP 4427",&
         "DI 10.1021/ct400595k",&
         "UT WOS:000326355100014",&
         "ER"),&
         DOI="10.1021/ct400595k")

    CALL add_reference(key=Hutter2014, ISI_record=s2a(&
         "PT J",&
         "AU Hutter, J",&
         "   Iannuzzi, M",&
         "   Schiffmann, F",&
         "   VandeVondele, J",&
         "AF Hutter, Juerg",&
         "   Iannuzzi, Marcella",&
         "   Schiffmann, Florian",&
         "   VandeVondele, Joost",&
         "TI CP2K: atomistic simulations of condensed matter systems",&
         "SO WILEY INTERDISCIPLINARY REVIEWS-COMPUTATIONAL MOLECULAR SCIENCE",&
         "SN 1759-0876",&
         "PD JAN",&
         "PY 2014",&
         "VL 4",&
         "IS 1",&
         "BP 15",&
         "EP 25",&
         "DI 10.1002/wcms.1159",&
         "UT WOS:000328646100002",&
         "ER"),&
         DOI="10.1002/wcms.1159")

    CALL add_reference(key=Kantorovich2008, ISI_record=s2a(&
         "PT J",&
         "AU Kantorovich, L",&
         "AF Kantorovich, Lev",&
         "TI Generalized Langevin equation for solids. I. ",&
         "   Rigorous derivation and main properties",&
         "SO PHYSICAL REVIEW B",&
         "SN 1098-0121",&
         "PD SEP",&
         "PY 2008",&
         "VL 78",&
         "IS 9",&
         "AR 094304",&
         "DI 10.1103/PhysRevB.78.094304",&
         "ER"),&
         DOI="10.1103/PhysRevB.78.094304")

    CALL add_reference(key=Kantorovich2008a, ISI_record=s2a(&
         "PT J",&
         "AU Kantorovich, L",&
         "   Rompotis, N",&
         "AF Kantorovich, Lev",&
         "   Rompotis, N",&
         "TI Generalized Langevin equation for solids. II. ",&
         "   Stochastic boundary conditions for nonequilibrium molecular dynamics simulations",&
         "SO PHYSICAL REVIEW B",&
         "SN 1098-0121",&
         "PD SEP",&
         "PY 2008",&
         "VL 78",&
         "IS 9",&
         "AR 094305",&
         "DI 10.1103/PhysRevB.78.094305",&
         "ER"),&
         DOI="10.1103/PhysRevB.78.094305")

    CALL add_reference(key=Niklasson2014, ISI_record=s2a(&
         "PT J",&
         "AU Rubensson, E.",&
         "AU Niklasson, A.",&
         "TI Interior Eigenvalues from Density Matrix Expansions in Quantum",&
         "    Mechanical Molecular Dynamics",&
         "SO SIAM Journal on Scientific Computing",&
         "PD MAR",&
         "PY 2014",&
         "VL 36",&
         "IS 2",&
         "BP B147",&
         "EP B170",&
         "DI 10.1137/130911585",&
         "ER" ),&
         DOI="10.1137/130911585")

    CALL add_reference(key=Borstnik2014, ISI_record=s2a(&
         "PT J",&
         "AU Borstnik, U",&
         "   VandeVondele, J",&
         "   Weber, V",&
         "   Hutter, J",&
         "AF Borstnik, Urban",&
         "   VandeVondele, Joost",&
         "   Weber, Valery",&
         "   Hutter, Juerg",&
         "TI Sparse matrix multiplication: The distributed block-compressed sparse",&
         "   row library",&
         "SO PARALLEL COMPUTING",&
         "SN 0167-8191",&
         "EI 1872-7336",&
         "PD MAY",&
         "PY 2014",&
         "VL 40",&
         "IS 5-6",&
         "BP 47",&
         "EP 58",&
         "DI 10.1016/j.parco.2014.03.012",&
         "UT WOS:000338614300001",&
         "ER"),&
         DOI="10.1016/j.parco.2014.03.012")

    CALL add_reference(key=Merlot2014, ISI_record=s2a(&
         "PT J",&
         "AU Merlot, P",&
         "   Izsak, R",&
         "   Borgo, A",&
         "   Kjaergaard, T",&
         "   Helgaker, T",&
         "   Reine, S",&
         "AF Merlot, Patrick",&
         "   Izsak, Robert",&
         "   Borgo, Alex",&
         "   Kjaergaard, Thomas",&
         "   Helgaker, Trygve",&
         "   Reine, Simen",&
         "TI Charge-constrained auxiliary-density-matrix methods for the",&
         "   Hartree-Fock exchange contribution",&
         "SO JOURNAL OF CHEMICAL PHYSICS",&
         "SN 0021-9606",&
         "PD SEP",&
         "PY 2014",&
         "VL 141",&
         "IS 9" ,&
         "AR 094104",&
         "DI /10.1063/1.4894267",&
         "ER"),&
         DOI="10.1063/1.4894267")

    CALL add_reference(key=Lin2009, ISI_record=s2a(&
          "PT J",&
          "AU Lin, L",&
          "   Lu, JF",&
          "   Ying, LX",&
          "   Car, R",&
          "   E, WN",&
          "AF Lin, Lin",&
          "   Lu, Jianfeng",&
          "   Ying, Lexing",&
          "   Car, Roberto",&
          "   E, Weinan",&
          "TI FAST ALGORITHM FOR EXTRACTING THE DIAGONAL OF THE INVERSE MATRIX WITH",&
          "   APPLICATION TO THE ELECTRONIC STRUCTURE ANALYSIS OF METALLIC SYSTEMS",&
          "SO COMMUNICATIONS IN MATHEMATICAL SCIENCES",&
          "VL 7",&
          "SN 1539-6746",&
          "PD SEP",&
          "PY 2009",&
          "RI Lin, Lin/I-2726-2012",&
          "OI Lin, Lin/0000-0001-7738-5947",&
          "TC 20",&
          "ZB 0",&
          "Z8 3",&
          "ZS 0",&
          "Z9 21",&
          "VL 7",&
          "IS 3",&
          "BP 755",&
          "EP 777",&
          "UT WOS:000271568500012",&
          "ER"),&
          DOI="10.4310/CMS.2009.v7.n3.a12")

     CALL add_reference(key=Lin2013, ISI_record=s2a(&
          "PT J",&
          "AU Lin, Lin",&
          "   Chen, Mohan",&
          "   Yang, Chao",&
          "   He, Lixin",&
          "TI Accelerating atomic orbital-based electronic structure calculation via",&
          "   pole expansion and selected inversion",&
          "SO JOURNAL OF PHYSICS-CONDENSED MATTER",&
          "VL 25",&
          "IS 29",&
          "AR 295501",&
          "PD JUL 24 2013",&
          "PY 2013",&
          "TC 3",&
          "ZB 1",&
          "Z8 0",&
          "ZS 0",&
          "Z9 3",&
          "SN 0953-8984",&
          "UT WOS:000321219200006",&
          "PM 23803312",&
          "ER"),&
          DOI="10.1088/0953-8984/25/29/295501")

     CALL add_reference(key=DelBen2015, ISI_record=s2a(&
          "PT J",&
          "AU Del Ben, M",&
          "   Schuett, O",&
          "   Wentz, T",&
          "   Messmer, P",&
          "   Hutter, J",&
          "   VandeVondele, J",&
          "AF Del Ben, Mauro",&
          "   Schuett, Ole",&
          "   Wentz, Tim",&
          "   Messmer, Peter",&
          "   Hutter, Juerg",&
          "   VandeVondele, Joost",&
          "TI Enabling simulation at the fifth rung of DFT: Large scale RPA",&
          "   calculations with excellent time to solution",&
          "SO COMPUTER PHYSICS COMMUNICATIONS",&
          "SN 0010-4655",&
          "EI 1879-2944",&
          "PD FEB",&
          "PY 2015",&
          "VL 187",&
          "BP 120",&
          "EP 129",&
          "DI 10.1016/j.cpc.2014.10.021",&
          "UT WOS:000346954200014",&
          "ER"),&
          DOI="10.1016/j.cpc.2014.10.021")

  END SUBROUTINE add_all_references

END MODULE bibliography
