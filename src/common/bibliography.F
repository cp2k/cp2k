!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief collects all references to literature in CP2K
!>      as new algorithms / method are included from literature sources
!>      we add a citation, to give credit to the authors
!>      (often ourselves, but this need not to be the case), and as a form
!>      of documentation.
!>
!>      reference_manager.F provides an easy way to cite these references from the code
!>      (e.g., whenever the relevant part of the code is executed)
!>      and to add citations to the input
!>      at the end of the run a bibliography is printed that can be used, e.g., in papers
!> \note
!>      references need to be input using the ISI citation format, because it is
!>      uniform, easy to parse, and can be exported for example from web of science.
!>      (mark the reference, and click 'export to reference software', open the resulting file)
!>      Furthermore, it can be easily converted to and from using the bibutils tools
!>      a collection of easy to use conversion programs that can be found at
!>      https://ctan.org/pkg/bibutils by Chris Putnam
!> \author Joost VandeVondele
! **************************************************************************************************
MODULE bibliography

   USE reference_manager,               ONLY: add_reference,&
                                              cite_reference
   USE string_utilities,                ONLY: s2a
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: cite_reference ! for convenience
   PUBLIC :: add_all_references ! f77_interface only

   ! all reference keys
   INTEGER, PUBLIC, SAVE :: Richters2018, VandeVondele2005a, VandeVondele2003, &
                            Laino2005, Laino2006, Goedecker1996, Hartwigsen1998, Krack2005, &
                            Lippert1997, Lippert1999, Krack2000, Krack2002, Iannuzzi2005, &
                            Iannuzzi2006, Iannuzzi2007, Toukmaji1996, Martyna1999, &
                            VandeVondele2005b, Perdew1981, Avezac2005, Porezag1995, &
                            Seifert1996, Elstner1998, Zhechkov2005, Henkelman1999, Henkelman2014, &
                            Frigo2005, Nose1984a, Nose1984b, Brelaz1979, Bussi2007, &
                            Genovese2006, Genovese2007, Evans1983, Minary2003, Byrd1995, &
                            VandeVondele2007, Ortiz1994, Becke1988, Perdew1996, Zhang1998, &
                            Perdew2008, Lee1988, Heyd2006, Vydrov2006, Heyd2003, Heyd2004, &
                            Vosko1980, Aguado2003, Essmann1995, Ewald1921, Darden1993, &
                            Siepmann1995, Tersoff1988, Tosi1964a, Tosi1964b, Yamada2000, &
                            Dudarev1997, Dudarev1998, Dewar1977, Dewar1985, Rocha2006, &
                            Stewart1989, Thiel1992, Repasky2002, Stewart2007, Weber2008, &
                            Hunt2003, Guidon2008, Elber1987, Jonsson1998, Jonsson2000_1, &
                            Jonsson2000_2, Wales2004, Stewart1982, E2002, Laino2008, &
                            Blochl1995, Tao2003, VandeVondele2006, VanVoorhis2015, Hu2007, Grimme2006, &
                            Branduardi2007, Schenter2008, BeckeRoussel1989, Proynov2007, &
                            VandenCic2006, Becke1997, Foiles1986, Ricci2003, Kolafa2004, &
                            Kuhne2007, Kunert2003, Ceriotti2009, ceriotti2009b, Guidon2009, &
                            BarducBus2008, Grimme2010, Guidon2010, Marques2012, Jones2011, &
                            Bernstein2009, Bernstein2012, Dick1958, Mitchell1993, &
                            Devynck2012, VandeVondele2012, Niklasson2003, Shao2003, &
                            VandeVondele2002, Dion2004, Romanperez2009, DelBen2012, &
                            Sabatini2013, DelBen2013, Kikuchi2009, Putrino2000, &
                            Putrino2002, Sebastiani2001, Weber2009, Tran2013, Golze2013, &
                            Golze2015, Golze2017a, Golze2017b, Tuckerman1992, Zhao1994, &
                            Tozer1996, Goedecker2004, Rengaraj2020, &
                            Khaliullin2013, Kruse2012, Hutter2014, Bengtsson1999, Kantorovich2008, &
                            Kantorovich2008a, Wellendorff2012, Niklasson2014, Borstnik2014, &
                            Rayson2009, Grimme2011, Fattebert2002, Andreussi2012, &
                            Khaliullin2007, Khaliullin2008, Merlot2014, Lin2009, Lin2013, &
                            QUIP_ref, Batzner2022, DelBen2015, Souza2002, Umari2002, Stengel2009, &
                            Luber2014, Berghold2011, DelBen2015b, Campana2009, &
                            Schiffmann2015, Bruck2014, Rappe1992, Ceriotti2012, &
                            Ceriotti2010, Walewski2014, Monkhorst1976, MacDonald1978, &
                            Gilbert2008, Schonherr2014, Ceriotti2014, BaniHashemian2016, &
                            Kapil2016, Heinzmann1976, Ehrhardt1985, Rybkin2016, West2006, &
                            Bates2013, Andermatt2016, Zhu2016, Schuett2016, Lu2004, &
                            Becke1988b, Migliore2009, Mavros2015, Holmberg2017, Marek2014, &
                            Stoychev2016, Futera2017, Bailey2006, Papior2017, Lehtola2018, &
                            Brieuc2016, Barca2018, Scheiber2018, Huang2011, Heaton_Burgess2007, &
                            Schuett2018, Holmberg2018, Togo2018, Staub2019, Grimme2013, Grimme2016, &
                            Grimme2017, Kondov2007, Clabaut2020, Clabaut2021, &
                            Ren2011, Ren2013, Cohen2000, Rogers2002, Filippetti2000, &
                            Limpanuparb2011, Martin2003, Yin2017, Goerigk2017, &
                            Wilhelm2016a, Wilhelm2016b, Wilhelm2017, Wilhelm2018, Wilhelm2021, Lass2018, &
                            cp2kqs2020, Behler2007, Behler2011, Schran2020a, Schran2020b, &
                            Rycroft2009, Thomas2015, Brehm2018, Brehm2020, Shigeta2001, Heinecke2016, &
                            Brehm2021, Bussy2021a, Bussy2021b, Ditler2021, Ditler2022, Mattiat2019, &
                            Mattiat2022, Belleflamme2023, Knizia2013, Musaelian2023, Eriksen2020, &
                            Bussy2023, Bussy2024, Wang2018, Zeng2023, Graml2024, Solca2024, &
                            Caldeweyher2017, Caldeweyher2019, Caldeweyher2020, Freeman1977, Gruneis2009, &
                            Stein2022, Stein2024, &
                            Blase2018, Blase2020, Bruneval2015, Golze2019, Gui2018, Jacquemin2017, Liu2020, &
                            Sander2015, Schreiber2008, Hueser2013, vanSetten2015, Setyawan2010

CONTAINS

! **************************************************************************************************
!> \brief adds references that can later be cited / printed using the key
!> \par History
!>      08.2007 created [Joost VandeVondele]
!> \note
!>     - note that the ISI record is fixed format (line length limited and the
!>       first 3 characters can indicate record type)
!>     - please add DOI whenever available, this will result in a clickable
!>       link in the input reference manual.
! **************************************************************************************************
   SUBROUTINE add_all_references()
      CALL add_reference(key=Ceriotti2012, &
                         authors=s2a("Ceriotti, M", "Manolopoulos, D"), &
                         title="Efficient First-Principles Calculation "// &
                         "of the Quantum Kinetic Energy and Momentum Distribution of Nuclei", &
                         source="PHYSICAL REVIEW LETTERS", volume="109", issue="10", pages="100604", &
                         year=2012, doi="10.1103/PhysRevLett.109.100604")

      CALL add_reference(key=Ceriotti2010, &
                         authors=s2a("Ceriotti, M", "Parrinello, M", "Markland, T", "Manolopoulos, D"), &
                         title="Efficient stochastic thermostatting of path integral molecular dynamics", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="133", issue="12", pages="124104", &
                         year=2010, doi="10.1063/1.3489925")

      CALL add_reference(key=Wellendorff2012, &
                         authors=s2a("Wellendorff, J", "Lundgaard, K", "Mogelhoj, A", "Petzold, V", &
                                     "Landis, D", "Norskov, J", "Bligaard, T", "Jacobsen, K"), &
                         title="Density functionals for surface science: "// &
                         "Exchange-correlation model development with Bayesian error estimation", &
                         source="PHYSICAL REVIEW B", volume="85", issue="23", pages="235149", &
                         year=2012, doi="10.1103/PhysRevB.85.235149")

      CALL add_reference(key=Brelaz1979, &
                         authors=s2a("Brelaz, D"), &
                         title="New methods to color the vertices of a graph", &
                         source="COMMUNICATIONS OF THE ACM", volume="22", issue="4", pages="251-256", &
                         year=1979, doi="10.1145/359094.359101")

      CALL add_reference(key=Bengtsson1999, &
                         authors=s2a("BENGTSSON, L"), &
                         title="DIPOLE CORRECTION FOR SURFACE SUPERCELL CALCULATIONE", &
                         source="PHYSICAL REVIEW B", volume="59", issue="19", pages="12301-12304", &
                         year=1999, doi="10.1103/PhysRevB.59.12301")

      CALL add_reference(key=Foiles1986, &
                         authors=s2a("FOILES, SM", "BASKES, MI", "DAW, MS"), &
                         title="EMBEDDED-ATOM-METHOD FUNCTIONS FOR THE FCC METALS CU, AG, AU, NI, PD, PT, AND THEIR ALLOYS", &
                         source="PHYSICAL REVIEW B", volume="33", issue="12", pages="7983-7991", &
                         year=1986, month=6, day=15, doi="10.1103/PhysRevB.33.7983")

      CALL add_reference(key=QUIP_ref, &
                         authors=s2a("QUIP"), &
                         title="libAtoms/QUIP libraries from http://www.libatoms.org, "// &
                         "please cite web site and references for appropriate potential invoked", &
                         source="web site", &
                         year=2014)

      CALL add_reference(key=Batzner2022, &
                         authors=s2a("Batzner, S", "Musaelian, A", "Sun, L", "Geiger, M", "Mailoa, JP", &
                                     "Kornbluth, M", "Molinari, N", "Smidt, TE", "Kozinsky, B"), &
                         title="E(3)-equivariant graph neural networks for data-efficient and accurate interatomic potentials", &
                         source="Nature Communications", volume="13", pages="2453", &
                         year=2022, doi="10.1038/s41467-022-29939-5")

      CALL add_reference(key=VandenCic2006, &
                         authors=s2a("Vanden-Eijnden, E", "Ciccotti, G"), &
                         title="Second-order integrators for Langevin equations with holonomic constraints", &
                         source="CHEMICAL PHYSICS LETTERS", volume="429", issue="1-3", pages="310-316", &
                         year=2006, month=9, day=29, doi="10.1016/j.cplett.2006.07.086")

      CALL add_reference(key=Hu2007, &
                         authors=s2a("Hu, H", "Lu, ZY", "Elstner, M", "Hermans, J", "Yang, WT"), &
                         title="Simulating water with the self-consistent-charge "// &
                         "density functional tight binding method: From molecular clusters to the liquid state", &
                         source="JOURNAL OF PHYSICAL CHEMISTRY A", volume="111", issue="26", pages="5685-5691", &
                         year=2007, month=7, day=5, doi="10.1021/jp070308d")

      CALL add_reference(key=Zhao1994, &
                         authors=s2a("ZHAO, QS", "MORRISON, RC", "PARR, RG"), &
                         title="From electron densities to Kohn-Sham kinetic energies, "// &
                         "orbital energies, exchange-correlation potentials, and exchange-correlation energies", &
                         source="PHYSICAL REVIEW A", volume="50", issue="3", pages="2138-2142", &
                         year=1994, month=9, day=1, doi="10.1103/PhysRevA.50.2138")

      CALL add_reference(key=Tozer1996, &
                         authors=s2a("TOZER, DJ", "INGAMELLS, VE", "HANDY, NC"), &
                         title="EXCHANGE-CORRELATION POTENTIALS", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="105", issue="20", pages="9200-9213", &
                         year=1996, month=11, day=22, doi="10.1063/1.472753")

      CALL add_reference(key=Blochl1995, &
                         authors=s2a("BLOCHL, PE"), &
                         title="ELECTROSTATIC DECOUPLING OF PERIODIC IMAGES OF "// &
                         "PLANE-WAVE-EXPANDED DENSITIES AND DERIVED ATOMIC POINT CHARGES", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="103", issue="17", pages="7422-7428", &
                         year=1995, month=11, day=1, doi="10.1063/1.470314")

      CALL add_reference(key=Laino2008, &
                         authors=s2a("Laino, T", "Hutter, J"), &
                         title='Notes on "Ewald summation of electrostatic multipole interactions '// &
                         'up to quadrupolar level [J. Chem. Phys. 119, 7471 (2003)]', &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="129", issue="7", pages="074102", &
                         year=2008, month=8, day=21, doi="10.1063/1.2970887")

      CALL add_reference(key=E2002, &
                         authors=s2a("E, WN", "Ren, WQ", "Vanden-Eijnden, E"), &
                         title="String method for the study of rare events", &
                         source="PHYSICAL REVIEW B", volume="66", issue="5", pages="052301", &
                         year=2002, month=8, day=1, doi="10.1103/PhysRevB.66.052301")

      CALL add_reference(key=Wales2004, &
                         authors=s2a("Trygubenko, SA", "Wales, DJ"), &
                         title="A doubly nudged elastic band method for finding transition states", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="120", issue="5", pages="2082-2094", &
                         year=2004, month=2, day=1, doi="10.1063/1.1636455")

      CALL add_reference(key=Jonsson2000_2, &
                         authors=s2a("Henkelman, G", "Uberuaga, BP", "Jonsson, H"), &
                         title="A climbing image nudged elastic band method for finding "// &
                         "saddle points and minimum energy paths", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="113", issue="22", pages="9901-9904", &
                         year=2000, month=12, day=8, doi="10.1063/1.1329672")

      CALL add_reference(key=Jonsson2000_1, &
                         authors=s2a("Henkelman, G", "Jonsson, H"), &
                         title="Improved tangent estimate in the nudged elastic band method for "// &
                         "finding minimum energy paths and saddle points", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="113", issue="22", pages="9978-9985", &
                         year=2000, month=12, day=8, doi="10.1063/1.1323224")

      CALL add_reference(key=Jonsson1998, &
                         authors=s2a("JONSSON, H", "MILLS, G", "JACOBSEN, K W"), &
                         title="Nudged elastic band method for finding minimum energy paths of transitions", &
                         source="Classical and Quantum Dynamics in Condensed Phase Simulations", pages="385-404", &
                         year=1998)

      CALL add_reference(key=Elber1987, &
                         authors=s2a("ELBER, R", "KARPLUS, M"), &
                         title="A METHOD FOR DETERMINING REACTION PATHS IN LARGE MOLECULES - APPLICATION TO MYOGLOBIN", &
                         source="CHEMICAL PHYSICS LETTERS", volume="139", issue="5", pages="375-380", &
                         year=1987, month=9, day=4, doi="10.1016/0009-2614(87)80576-6")

      CALL add_reference(key=Weber2008, &
                         authors=s2a("Weber, V", "VandeVondele, J", "Hutter, J", "Niklasson, AMN"), &
                         title="Direct energy functional minimization under orthogonality constraints", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="128", issue="8", pages="084113", &
                         year=2008, month=2, day=28, doi="10.1063/1.2841077")

      CALL add_reference(key=Stewart2007, &
                         authors=s2a("Stewart, JJP"), &
                         title="Optimization of parameters for semiempirical methods V: "// &
                         "Modification of NDDO approximations and application to 70 elements", &
                         source="JOURNAL OF MOLECULAR MODELING", volume="13", issue="12", pages="1173-1213", &
                         year=2007, month=12, doi="10.1007/s00894-007-0233-4")

      CALL add_reference(key=Repasky2002, &
                         authors=s2a("Repasky, MP", "Chandrasekhar, J", "Jorgensen, WL"), &
                         title="PDDG/PM3 and PDDG/MNDO: Improved semiempirical methods", &
                         source="JOURNAL OF COMPUTATIONAL CHEMISTRY", volume="23", issue="16", pages="1601-1622", &
                         year=2002, month=12, doi="10.1002/jcc.10162")

      CALL add_reference(key=Thiel1992, &
                         authors=s2a("THIEL, W", "VOITYUK, AA"), &
                         title="EXTENSION OF THE MNDO FORMALISM TO D-ORBITALS - "// &
                         "INTEGRAL APPROXIMATIONS AND PRELIMINARY NUMERICAL RESULTS", &
                         source="THEORETICA CHIMICA ACTA", volume="81", issue="6", pages="391-404", &
                         year=1992, month=2, doi="10.1007/BF01134863")

      CALL add_reference(key=Stewart1989, &
                         authors=s2a("STEWART, JJP"), &
                         title="OPTIMIZATION OF PARAMETERS FOR SEMIEMPIRICAL METHODS .1. METHOD", &
                         source="JOURNAL OF COMPUTATIONAL CHEMISTRY", volume="10", issue="2", pages="209-220", &
                         year=1989, month=3, doi="10.1002/jcc.540100208")

      CALL add_reference(key=Rocha2006, &
                         authors=s2a("Rocha, GB", "Freire, RO", "Simas, AM", "Stewart, JJP"), &
                         title="RMI: A reparameterization of AM1 for H, C, N, O, P, S, F, Cl, Br, and I", &
                         source="JOURNAL OF COMPUTATIONAL CHEMISTRY", volume="27", issue="10", pages="1101-1111", &
                         year=2006, month=7, day=30, doi="10.1002/jcc.20425")

      CALL add_reference(key=Dewar1985, &
                         authors=s2a("DEWAR, MJS", "ZOEBISCH, EG", "HEALY, EF", "STEWART, JJP"), &
                         title="THE DEVELOPMENT AND USE OF QUANTUM-MECHANICAL MOLECULAR-MODELS .76. AM1 - "// &
                         "A NEW GENERAL-PURPOSE QUANTUM-MECHANICAL MOLECULAR-MODEL", &
                         source="JOURNAL OF THE AMERICAN CHEMICAL SOCIETY", volume="107", issue="13", pages="3902-3909", &
                         year=1985, doi="10.1021/ja00299a024")

      CALL add_reference(key=Dewar1977, &
                         authors=s2a("DEWAR, MJS", "THIEL, W"), &
                         title="GROUND-STATES OF MOLECULES .38. MNDO METHOD - APPROXIMATIONS AND PARAMETERS", &
                         source="JOURNAL OF THE AMERICAN CHEMICAL SOCIETY", volume="99", issue="15", pages="4899-4907", &
                         year=1977, doi="10.1021/ja00457a004")

      CALL add_reference(key=Henkelman1999, &
                         authors=s2a("Henkelman, G", "Jonsson, H"), &
                         title="A dimer method for finding saddle points on high dimensional "// &
                         "potential surfaces using only first derivatives", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="111", issue="15", pages="7010-7022", &
                         year=1999, month=10, day=15, doi="10.1063/1.480097")

      CALL add_reference(key=Henkelman2014, &
                         authors=s2a("Xiao, P, Wu, Q, Henkelman, G"), &
                         title="Basin constrained k-dimer method for saddle point finding", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="141", issue="16", pages="164111", &
                         year=2014, month=10, day=15, doi="10.1063/1.4898664")

      CALL add_reference(key=Aguado2003, &
                         authors=s2a("Aguado, A", "Madden, PA"), &
                         title="Ewald summation of electrostatic multipole interactions up to the quadrupolar level", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="119", issue="14", pages="7471-7483", &
                         year=2003, month=10, day=8, doi="10.1063/1.1605941")

      CALL add_reference(key=Yamada2000, &
                         authors=s2a("Yamada, K", "Kurosaki, K", "Uno, M", "Yamanaka, S"), &
                         title="Evaluation of thermal properties of uranium dioxide by molecular dynamics", &
                         source="JOURNAL OF ALLOYS AND COMPOUNDS", volume="307", pages="10-16", &
                         year=2000, month=7, day=14, doi="10.1016/S0925-8388(00)00806-9")

      CALL add_reference(key=Tosi1964a, &
                         authors=s2a("FUMI, FG", "TOSI, MP"), &
                         title="IONIC SIZES + BORN REPULSIVE PARAMETERS IN NACL-TYPE ALKALI HALIDES "// &
                         ".I. HUGGINS-MAYER + PAULING FORMS", &
                         source="JOURNAL OF PHYSICS AND CHEMISTRY OF SOLIDS", volume="25", issue="1", pages="31-43", &
                         year=1964, doi="10.1016/0022-3697(64)90159-3")

      CALL add_reference(key=Tosi1964b, &
                         authors=s2a("TOSI, MP", "FUMI, FG"), &
                         title="IONIC SIZES + BORN REPULSIVE PARAMETERS IN NACL-TYPE ALKALI HALIDES .2. GENERALIZED", &
                         source="JOURNAL OF PHYSICS AND CHEMISTRY OF SOLIDS", volume="25", issue="1", pages="45-52", &
                         year=1964, doi="10.1016/0022-3697(64)90160-X")

      CALL add_reference(key=Tersoff1988, &
                         authors=s2a("TERSOFF, J"), &
                         title="EMPIRICAL INTERATOMIC POTENTIAL FOR SILICON WITH IMPROVED ELASTIC PROPERTIES", &
                         source="PHYSICAL REVIEW B", volume="38", issue="14", pages="9902-9905", &
                         year=1988, month=11, day=15, doi="10.1103/PhysRevB.38.9902")

      CALL add_reference(key=Siepmann1995, &
                         authors=s2a("SIEPMANN, JI", "SPRIK, M"), &
                         title="INFLUENCE OF SURFACE TOPOLOGY AND ELECTROSTATIC POTENTIAL ON WATER/ELECTRODE SYSTEMS", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="102", issue="1", pages="511-524", &
                         year=1995, month=1, day=1, doi="10.1063/1.469429")

      CALL add_reference(key=Bussi2007, &
                         authors=s2a("Bussi, G", "Donadio, D", "Parrinello, M"), &
                         title="Canonical sampling through velocity rescaling", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="126", issue="1", pages="014101", &
                         year=2007, month=1, day=7, doi="10.1063/1.2408420")

      CALL add_reference(key=Nose1984a, &
                         authors=s2a("NOSE, S"), &
                         title="A UNIFIED FORMULATION OF THE CONSTANT TEMPERATURE MOLECULAR-DYNAMICS METHODS", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="81", issue="1", pages="511-519", &
                         year=1984, doi="10.1063/1.447334")

      CALL add_reference(key=Nose1984b, &
                         authors=s2a("NOSE, S"), &
                         title="A MOLECULAR-DYNAMICS METHOD FOR SIMULATIONS IN THE CANONICAL ENSEMBLE", &
                         source="MOLECULAR PHYSICS", volume="52", issue="2", pages="255-268", &
                         year=1984, doi="10.1080/00268978400101201")

      CALL add_reference(key=VandeVondele2005a, &
                         authors=s2a("VandeVondele, J", "Krack, M", "Mohamed, F", "Parrinello, M", &
                                     "Chassaing, T", "Hutter, J"), &
                         title="QUICKSTEP: Fast and accurate density functional calculations "// &
                         "using a mixed Gaussian and plane waves approach", &
                         source="COMPUTER PHYSICS COMMUNICATIONS", volume="167", issue="2", pages="103-128", &
                         year=2005, month=4, day=15, doi="10.1016/j.cpc.2004.12.014")

      CALL add_reference(key=VandeVondele2003, &
                         authors=s2a("VandeVondele, J", "Hutter, J"), &
                         title="An efficient orbital transformation method for electronic structure calculations", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="118", issue="10", pages="4365-4369", &
                         year=2003, month=3, day=8, doi="10.1063/1.1543154")

      CALL add_reference(key=Laino2005, &
                         authors=s2a("Laino, T", "Mohamed, F", "Laio, A", "Parrinello, M"), &
                         title="An efficient real space multigrid QM/MM electrostatic coupling", &
                         source="JOURNAL OF CHEMICAL THEORY AND COMPUTATION", volume="1", issue="6", pages="1176-1184", &
                         year=2005, month=11, doi="10.1021/ct050123f")

      CALL add_reference(key=Laino2006, &
                         authors=s2a("Laino, T", "Mohamed, F", "Laio, A", "Parrinello, M"), &
                         title="An efficient linear-scaling electrostatic coupling for treating "// &
                         "periodic boundary conditions in QM/MM simulations", &
                         source="JOURNAL OF CHEMICAL THEORY AND COMPUTATION", volume="2", issue="5", pages="1370-1378", &
                         year=2006, month=9, day=12, doi="10.1021/ct6001169")

      CALL add_reference(key=Goedecker1996, &
                         authors=s2a("Goedecker, S", "Teter, M", "Hutter, J"), &
                         title="Separable dual-space Gaussian pseudopotentials", &
                         source="PHYSICAL REVIEW B", volume="54", issue="3", pages="1703-1710", &
                         year=1996, month=7, day=15, doi="10.1103/PhysRevB.54.1703")

      CALL add_reference(key=Hartwigsen1998, &
                         authors=s2a("Hartwigsen, C", "Goedecker, S", "Hutter, J"), &
                         title="Relativistic separable dual-space Gaussian pseudopotentials from H to Rn", &
                         source="PHYSICAL REVIEW B", volume="58", issue="7", pages="3641-3662", &
                         year=1998, month=8, day=15, doi="10.1103/PhysRevB.58.3641")

      CALL add_reference(key=Krack2005, &
                         authors=s2a("Krack, M"), &
                         title="Pseudopotentials for H to Kr optimized for "// &
                         "gradient-corrected exchange-correlation functionals", &
                         source="THEORETICAL CHEMISTRY ACCOUNTS", volume="114", issue="1-3", pages="145-152", &
                         year=2005, month=9, doi="10.1007/s00214-005-0655-y")

      CALL add_reference(key=Lippert1997, &
                         authors=s2a("Lippert, G", "Hutter, J", "Parrinello, M"), &
                         title="A hybrid Gaussian and plane wave density functional scheme", &
                         source="MOLECULAR PHYSICS", volume="92", issue="3", pages="477-487", &
                         year=1997, month=10, day=20, doi="10.1080/002689797170220")

      CALL add_reference(key=Lippert1999, &
                         authors=s2a("Lippert, G", "Hutter, J", "Parrinello, M"), &
                         title="The Gaussian and augmented-plane-wave density functional method for "// &
                         "ab initio molecular dynamics simulations", &
                         source="THEORETICAL CHEMISTRY ACCOUNTS", volume="103", issue="2", pages="124-140", &
                         year=1999, month=12, doi="10.1007/s002140050523")

      CALL add_reference(key=Krack2002, &
                         authors=s2a("Krack, M", "Gambirasio, A", "Parrinello, M"), &
                         title="Ab initio x-ray scattering of liquid water", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="117", issue="20", pages="9409-9412", &
                         year=2002, month=11, day=22, doi="10.1063/1.1517040")

      CALL add_reference(key=Krack2000, &
                         authors=s2a("Krack, M", "Parrinello, M"), &
                         title="All-electron ab-initio molecular dynamics", &
                         source="PHYSICAL CHEMISTRY CHEMICAL PHYSICS", volume="2", issue="10", pages="2105-2112", &
                         year=2000, doi="10.1039/b001167n")

      CALL add_reference(key=Iannuzzi2007, &
                         authors=s2a("Iannuzzi, M", "Hutter, J"), &
                         title="Inner-shell spectroscopy by the Gaussian and augmented plane wave method", &
                         source="PHYSICAL CHEMISTRY CHEMICAL PHYSICS", volume="9", issue="13", pages="1599-1610", &
                         year=2007, doi="10.1039/b615522g")

      CALL add_reference(key=Iannuzzi2006, &
                         authors=s2a("Iannuzzi, M", "Kirchner, B", "Hutter, J"), &
                         title="Density functional embedding for molecular systems", &
                         source="CHEMICAL PHYSICS LETTERS", volume="421", issue="1-3", pages="16-20", &
                         year=2006, month=4, day=3, doi="10.1016/j.cplett.2005.08.155")

      CALL add_reference(key=Iannuzzi2005, &
                         authors=s2a("Iannuzzi, M", "Chassaing, T", "Wallman, T", "Hutter, J"), &
                         title="Ground and excited state density functional calculations with the "// &
                         "Gaussian and augmented-plane-wave method", &
                         source="CHIMIA", volume="59", issue="7-8", pages="499-503", &
                         year=2005, doi="10.2533/000942905777676164")

      CALL add_reference(key=Toukmaji1996, &
                         authors=s2a("Toukmaji, AY", "Board, JA"), &
                         title="Ewald summation techniques in perspective: A survey", &
                         source="COMPUTER PHYSICS COMMUNICATIONS", volume="95", issue="2-3", pages="73-92", &
                         year=1996, month=6, doi="10.1016/0010-4655(96)00016-1")

      CALL add_reference(key=Martyna1999, &
                         authors=s2a("Martyna, GJ", "Tuckerman, ME"), &
                         title="A reciprocal space based method for treating long range interactions in "// &
                         "ab initio and force-field-based calculations in clusters", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="110", issue="6", pages="2810-2821", &
                         year=1999, month=2, day=8, doi="10.1063/1.477923")

      CALL add_reference(key=VandeVondele2005b, &
                         authors=s2a("VandeVondele, J", "Sprik, M"), &
                         title="A molecular dynamics study of the hydroxyl radical in solution "// &
                         "applying self-interaction-corrected density functional methods", &
                         source="PHYSICAL CHEMISTRY CHEMICAL PHYSICS", volume="7", issue="7", pages="1363-1367", &
                         year=2005, doi="10.1039/b501603g")

      CALL add_reference(key=Perdew1981, &
                         authors=s2a("PERDEW, JP", "ZUNGER, A"), &
                         title="SELF-INTERACTION CORRECTION TO DENSITY-FUNCTIONAL APPROXIMATIONS "// &
                         "FOR MANY-ELECTRON SYSTEMS", &
                         source="PHYSICAL REVIEW B", volume="23", issue="10", pages="5048-5079", &
                         year=1981, doi="10.1103/PhysRevB.23.5048")

      CALL add_reference(key=Avezac2005, &
                         authors=s2a("d'Avezac, M", "Calandra, M", "Mauri, F"), &
                         title="Density functional theory description of hole-trapping in SiO2: "// &
                         "A self-interaction-corrected approach", &
                         source="PHYSICAL REVIEW B", volume="71", issue="20", pages="205210", &
                         year=2005, month=5, doi="10.1103/PhysRevB.71.205210")

      CALL add_reference(key=Zhechkov2005, &
                         authors=s2a("Zhechkov, L", "Heine, T", "Patchkovskii, S", "Seifert, G", "Duarte, HA"), &
                         title="An efficient a Posteriori treatment for dispersion interaction in "// &
                         "density-functional-based tight binding", &
                         source="JOURNAL OF CHEMICAL THEORY AND COMPUTATION", volume="1", issue="5", pages="841-847", &
                         year=2005, month=9, doi="10.1021/ct050065y")

      CALL add_reference(key=Elstner1998, &
                         authors=s2a("Elstner, M", "Porezag, D", "Jungnickel, G", "Elsner, J", &
                                     "Haugk, M", "Frauenheim, T", "Suhai, S", "Seifert, G"), &
                         title="Self-consistent-charge density-functional tight-binding method for "// &
                         "simulations of complex materials properties", &
                         source="PHYSICAL REVIEW B", volume="58", issue="11", pages="7260-7268", &
                         year=1998, month=9, day=15, doi="10.1103/PhysRevB.58.7260")

      CALL add_reference(key=Seifert1996, &
                         authors=s2a("Seifert, G", "Porezag, D", "Frauenheim, T"), &
                         title="Calculations of molecules, clusters, and solids with a simplified LCAO-DFT-LDA scheme", &
                         source="INTERNATIONAL JOURNAL OF QUANTUM CHEMISTRY", volume="58", issue="2", pages="185-192", &
                         year=1996, month=4, day=15, doi="10.1002/(SICI)1097-461X(1996)58:2<185::AID-QUA7>3.0.CO;2-U")

      CALL add_reference(key=Porezag1995, &
                         authors=s2a("POREZAG, D", "FRAUENHEIM, T", "KOHLER, T", "SEIFERT, G", "KASCHNER, R"), &
                         title="CONSTRUCTION OF TIGHT-BINDING-LIKE POTENTIALS ON THE BASIS OF "// &
                         "DENSITY-FUNCTIONAL THEORY - APPLICATION TO CARBON", &
                         source="PHYSICAL REVIEW B", volume="51", issue="19", pages="12947-12957", &
                         year=1995, month=5, day=15, doi="10.1103/PhysRevB.51.12947")

      CALL add_reference(key=Frigo2005, &
                         authors=s2a("Frigo, M", "Johnson, SG"), &
                         title="The design and implementation of FFTW3", &
                         source="PROCEEDINGS OF THE IEEE", volume="93", issue="2", pages="216-231", &
                         year=2005, month=2, doi="10.1109/JPROC.2004.840301")

      CALL add_reference(key=Genovese2006, &
                         authors=s2a("Genovese, L", "Deutsch, T", "Neelov, A", "Goedecker, S", "Beylkin, G"), &
                         title="Efficient solution of Poisson's equation with free boundary conditions", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="125", issue="7", pages="074105", &
                         year=2006, month=8, day=21, doi="10.1063/1.2335442")

      CALL add_reference(key=Genovese2007, &
                         authors=s2a("Genovese, L", "Deutsch, T", "Goedecker, S"), &
                         title="Efficient and accurate three-dimensional Poisson solver for surface problems", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="127", issue="5", pages="054704", &
                         year=2007, month=8, day=7, doi="10.1063/1.2754685")

      CALL add_reference(key=Minary2003, &
                         authors=s2a("Minary, P", "Martyna, GJ", "Tuckerman, ME"), &
                         title="Algorithms and novel applications based on the isokinetic ensemble. "// &
                         "I. Biophysical and path integral molecular dynamics", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="118", issue="6", pages="2510-2526", &
                         year=2003, month=2, day=8, doi="10.1063/1.1534582")

      CALL add_reference(key=Evans1983, &
                         authors=s2a("EVANS, DJ", "HOOVER, WG", "FAILOR, BH", "MORAN, B", "LADD, AJC"), &
                         title="NON-EQUILIBRIUM MOLECULAR-DYNAMICS VIA GAUSS PRINCIPLE OF LEAST CONSTRAINT", &
                         source="PHYSICAL REVIEW A", volume="28", issue="2", pages="1016-1021", &
                         year=1983, doi="10.1103/PhysRevA.28.1016")

      CALL add_reference(key=Byrd1995, &
                         authors=s2a("BYRD, RH", "LU, PH", "NOCEDAL, J", "ZHU, CY"), &
                         title="A LIMITED MEMORY ALGORITHM FOR BOUND CONSTRAINED OPTIMIZATION", &
                         source="SIAM JOURNAL ON SCIENTIFIC COMPUTING", volume="16", issue="5", pages="1190-1208", &
                         year=1995, doi="10.1137/0916069")

      CALL add_reference(key=VandeVondele2007, &
                         authors=s2a("VandeVondele, J", "Hutter, J"), &
                         title="Gaussian basis sets for accurate calculations on molecular systems "// &
                         "in gas and condensed phases", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="127", issue="11", pages="114105", &
                         year=2007, month=9, day=21, doi="10.1063/1.2770708")

      CALL add_reference(key=Ortiz1994, &
                         authors=s2a("ORTIZ, G", "BALLONE, P"), &
                         title="CORRELATION-ENERGY, STRUCTURE FACTOR, RADIAL-DISTRIBUTION FUNCTION, "// &
                         "AND MOMENTUM DISTRIBUTION OF THE SPIN-POLARIZED UNIFORM ELECTRON-GAS", &
                         source="PHYSICAL REVIEW B", volume="50", issue="3", pages="1391-1405", &
                         year=1994, month=7, day=15, doi="10.1103/PhysRevB.50.1391")

      CALL add_reference(key=Becke1988, &
                         authors=s2a("BECKE, AD"), &
                         title="DENSITY-FUNCTIONAL EXCHANGE-ENERGY APPROXIMATION WITH "// &
                         "CORRECT ASYMPTOTIC-BEHAVIOR", &
                         source="PHYSICAL REVIEW A", volume="38", issue="6", pages="3098-3100", &
                         year=1988, month=9, day=15, doi="10.1103/PhysRevA.38.3098")

      CALL add_reference(key=Perdew1996, &
                         authors=s2a("Perdew, JP", "Burke, K", "Ernzerhof, M"), &
                         title="Generalized gradient approximation made simple", &
                         source="PHYSICAL REVIEW LETTERS", volume="77", issue="18", pages="3865-3868", &
                         year=1996, month=10, day=28, doi="10.1103/PhysRevLett.77.3865")

      CALL add_reference(key=Zhang1998, &
                         authors=s2a("Zhang, YK", "Yang, WT"), &
                         title="Comment on Generalized gradient approximation made simple", &
                         source="PHYSICAL REVIEW LETTERS", volume="80", issue="4", pages="890-890", &
                         year=1998, month=1, day=26, doi="10.1103/PhysRevLett.80.890")

      CALL add_reference(key=Perdew2008, &
                         authors=s2a("Perdew, JP", "Ruzsinszky, A", "Csonka, GI", "Vydrov, OA", &
                                     "Scuseria, GE", "Constantin, LA", "Zhou, X", "Burke, K"), &
                         title="Restoring the Density-Gradient Expansion for Exchange in Solids and Surfaces", &
                         source="PHYSICAL REVIEW LETTERS", volume="100", issue="13", pages="136406-136409", &
                         year=2008, month=4, day=4, doi="10.1103/PhysRevLett.100.136406")

      CALL add_reference(key=Lee1988, &
                         authors=s2a("LEE, CT", "YANG, WT", "PARR, RG"), &
                         title="DEVELOPMENT OF THE COLLE-SALVETTI CORRELATION-ENERGY FORMULA INTO A "// &
                         "FUNCTIONAL OF THE ELECTRON-DENSITY", &
                         source="PHYSICAL REVIEW B", volume="37", issue="2", pages="785-789", &
                         year=1988, month=1, day=15, doi="10.1103/PhysRevB.37.785")

      CALL add_reference(key=Heyd2004, &
                         authors=s2a("Heyd, J", "Scuseria, GE"), &
                         title="Assessment and validation of a screened Coulomb hybrid density functional", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="120", issue="16", pages="7274-7280", &
                         year=2004, month=4, day=22, doi="10.1063/1.1668634")

      CALL add_reference(key=Heyd2003, &
                         authors=s2a("Heyd, J", "Scuseria, GE", "Ernzerhof, M"), &
                         title="Hybrid functionals based on a screened Coulomb potential", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="118", issue="18", pages="8207-8215", &
                         year=2003, month=5, day=8, doi="10.1063/1.1564060")

      CALL add_reference(key=Heyd2006, &
                         authors=s2a("Heyd, J", "Scuseria, GE", "Ernzerhof, M"), &
                         title="Hybrid functionals based on a screened Coulomb potential (vol 118, pg 8207, 2003)", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="124", issue="21", pages="219906", &
                         year=2006, month=6, day=7, doi="10.1063/1.2204597")

      CALL add_reference(key=Vydrov2006, &
                         authors=s2a("Vydrov, OA", "Heyd, J", "Krukau, AV", "Scuseria, GE"), &
                         title="Importance of short-range versus long-range Hartree-Fock exchange "// &
                         "for the performance of hybrid density functionals", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="125", issue="7", pages="074106", &
                         year=2006, month=8, day=21, doi="10.1063/1.2244560")

      CALL add_reference(key=Vosko1980, &
                         authors=s2a("VOSKO, SH", "WILK, L", "NUSAIR, M"), &
                         title="ACCURATE SPIN-DEPENDENT ELECTRON LIQUID CORRELATION ENERGIES FOR "// &
                         "LOCAL SPIN-DENSITY CALCULATIONS - A CRITICAL ANALYSIS", &
                         source="CANADIAN JOURNAL OF PHYSICS", volume="58", issue="8", pages="1200-1211", &
                         year=1980, doi="10.1139/p80-159")

      CALL add_reference(key=Essmann1995, &
                         authors=s2a("ESSMANN, U", "PERERA, L", "BERKOWITZ, ML", "DARDEN, T", &
                                     "LEE, H", "PEDERSEN, LG"), &
                         title="A SMOOTH PARTICLE MESH EWALD METHOD", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="103", issue="19", pages="8577-8593", &
                         year=1995, month=11, day=15, doi="10.1063/1.470117")

      CALL add_reference(key=Ewald1921, &
                         authors=s2a("Ewald, PP"), &
                         title="Die Berechnung optischer und elektrostatischer Gitterpotentiale", &
                         source="ANNALEN DER PHYSIK", volume="369", issue="3", pages="253-287", &
                         year=1921, month=2, doi="10.1002/andp.19213690304")

      CALL add_reference(key=Darden1993, &
                         authors=s2a("DARDEN, T", "YORK, D", "PEDERSEN, L"), &
                         title="PARTICLE MESH EWALD - AN N.LOG(N) METHOD FOR EWALD SUMS IN LARGE SYSTEMS", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="98", issue="12", pages="10089-10092", &
                         year=1993, month=6, day=15, doi="10.1063/1.464397")

      CALL add_reference(key=Dudarev1997, &
                         authors=s2a("Dudarev, SL", "Manh, DN", "Sutton, AP"), &
                         title="Effect of Mott-Hubbard correlations on the electronic structure "// &
                         "and structural stability of uranium dioxide", &
                         source="PHILOSOPHICAL MAGAZINE B", volume="75", issue="5", pages="613-628", &
                         year=1997, month=5, doi="10.1080/13642819708202343")

      CALL add_reference(key=Dudarev1998, &
                         authors=s2a("Dudarev, SL", "Botton, GA", "Savrasov, SY", "Humphreys, CJ", &
                                     "Sutton, AP"), &
                         title="Electron-energy-loss spectra and the structural stability of "// &
                         "nickel oxide: An LSDA+U study", &
                         source="PHYSICAL REVIEW B", volume="57", issue="3", pages="1505-1509", &
                         year=1998, month=1, day=15, doi="10.1103/PhysRevB.57.1505")

      CALL add_reference(key=Hunt2003, &
                         authors=s2a("Hunt, P", "Sprik, M", "Vuilleumier, R"), &
                         title="Thermal versus electronic broadening in the density of states of liquid water", &
                         source="CHEMICAL PHYSICS LETTERS", volume="376", issue="1-2", pages="68-74", &
                         year=2003, month=7, day=17, doi="10.1016/S0009-2614(03)00954-0")

      CALL add_reference(key=Guidon2008, &
                         authors=s2a("Guidon, M", "Schiffmann, F", "Hutter, J", "VandeVondele, J"), &
                         title="Ab initio molecular dynamics using hybrid density functionals", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="128", issue="21", pages="214104", &
                         year=2008, month=6, day=7, doi="10.1063/1.2931945")

      CALL add_reference(key=Stewart1982, &
                         authors=s2a("Stewart, JJP", "Csaszar, P", "Pulay, P"), &
                         title="Fast semi-empirical calculations", &
                         source="JOURNAL OF COMPUTATIONAL CHEMISTRY", volume="3", issue="2", pages="227-228", &
                         year=1982, doi="10.1002/jcc.540030214")

      CALL add_reference(key=Tao2003, &
                         authors=s2a("Tao, JM", "Perdew, JP", "Staroverov, VN", "Scuseria, GE"), &
                         title="Climbing the density functional ladder: Nonempirical meta-generalized "// &
                         "gradient approximation designed for molecules and solids", &
                         source="PHYSICAL REVIEW LETTERS", volume="91", issue="14", pages="146401", &
                         year=2003, month=10, day=3, doi="10.1103/PhysRevLett.91.146401")

      CALL add_reference(key=VandeVondele2006, &
                         authors=s2a("VandeVondele, J", "Iannuzzi, M", "Hutter, J"), &
                         title="Large scale condensed matter calculations using the gaussian and "// &
                         "augmented plane waves method", &
                         source="Computer Simulations in Condensed Matter Systems: "// &
                         "From Materials to Chemical Biology, Vol 1", volume="703", pages="287-314", &
                         year=2006, doi="10.1007/3-540-35273-2_8")

      CALL add_reference(key=Grimme2006, &
                         authors=s2a("Grimme, S"), &
                         title="Semiempirical GGA-type density functional constructed with "// &
                         "a long-range dispersion correction", &
                         source="JOURNAL OF COMPUTATIONAL CHEMISTRY", volume="27", issue="15", pages="1787-1799", &
                         year=2006, month=11, day=30, doi="10.1002/jcc.20495")

      CALL add_reference(key=Grimme2010, &
                         authors=s2a("Grimme, S", "Antony, J", "Ehrlich, S", "Krieg, H"), &
                         title="A consistent and accurate ab initio parametrization of density "// &
                         "functional dispersion correction (DFT-D) for the 94 elements H-Pu", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="132", issue="15", pages="154104", &
                         year=2010, month=4, day=21, doi="10.1063/1.3382344")

      CALL add_reference(key=Grimme2011, &
                         authors=s2a("Grimme, S", "Ehrlich, S", "Goerigk, L"), &
                         title="Effect of the damping function in dispersion corrected density functional theory", &
                         source="JOURNAL OF COMPUTATIONAL CHEMISTRY", volume="32", pages="1456", &
                         year=2011, doi="10.1002/jcc.21759")

      CALL add_reference(key=Grimme2013, &
                         authors=s2a("Grimme, S"), &
                         title="A simplified Tamm-Dancoff density functional approach for the "// &
                         "electronic excitation spectra of very large molecules", &
                         source="The Journal of Chemical Physics", volume="138", pages="244104", &
                         year=2013, doi="10.1063/1.4811331")

      CALL add_reference(key=Grimme2016, &
                         authors=s2a("Grimme, S", "Bannwarth, C"), &
                         title="Ultra-fast computation of electronic spectra for large systems by "// &
                         "tight-binding based simplified Tamm-Dancoff approximation (sTDA-xTB)", &
                         source="The Journal of Chemical Physics", volume="145", pages="054103", &
                         year=2016, doi="10.1063/1.4959605")

      CALL add_reference(key=Grimme2017, &
                         authors=s2a("Grimme, S", "Bannwarth, C", "Shushkov, P"), &
                         title="A Robust and Accurate Tight-Binding Quantum Chemical Method for "// &
                         "Structures, Vibrational Frequencies, and Noncovalent Interactions of "// &
                         "Large Molecular Systems Parametrized for All spd-Block Elements (Z = 1-86)", &
                         source="Journal of Chemical Theory and Computation", volume="13", pages="1989", &
                         year=2017, doi="10.1021/acs.jctc.7b00118")

      CALL add_reference(key=Branduardi2007, &
                         authors=s2a("Branduardi, D", "Gervasio, FL", "Parrinello, M"), &
                         title="From A to B in free energy space", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="126", issue="5", pages="054103", &
                         year=2007, month=2, day=7, doi="10.1063/1.2432340")

      CALL add_reference(key=Schenter2008, &
                         authors=s2a("Chang, DT", "Schenter, GK", "Garrett, BC"), &
                         title="Self-consistent polarization neglect of diatomic differential overlap: "// &
                         "Applications to water clusters", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="128", issue="16", pages="164111", &
                         year=2008, month=4, day=28, doi="10.1063/1.2905230")

      CALL add_reference(key=Proynov2007, &
                         authors=s2a("Proynov, E", "Gan, Z", "Kong, J"), &
                         title="Analytical representation of the Becke-Roussel exchange functional", &
                         source="CHEMICAL PHYSICS LETTERS", volume="455", issue="1-3", pages="103-109", &
                         year=2008, month=3, day=31, doi="10.1016/j.cplett.2008.02.039")

      CALL add_reference(key=BeckeRoussel1989, &
                         authors=s2a("BECKE, AD", "ROUSSEL, MR"), &
                         title="EXCHANGE HOLES IN INHOMOGENEOUS SYSTEMS - A COORDINATE-SPACE MODEL", &
                         source="PHYSICAL REVIEW A", volume="39", issue="8", pages="3761-3767", &
                         year=1989, month=4, day=15, doi="10.1103/PhysRevA.39.3761")

      CALL add_reference(key=Becke1997, &
                         authors=s2a("Becke, AD"), &
                         title="Density-functional thermochemistry . 5. "// &
                         "Systematic optimization of exchange-correlation functionals", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="107", issue="20", &
                         year=1997, doi="10.1063/1.475007")

      CALL add_reference(key=Ricci2003, &
                         authors=s2a("Ricci, A", "Ciccotti, G"), &
                         title="Algorithms for Brownian dynamics", &
                         source="MOLECULAR PHYSICS", volume="101", issue="12", pages="1927-1931", &
                         year=2003, month=6, day=20, doi="10.1080/0026897031000108113")

      CALL add_reference(key=Kolafa2004, &
                         authors=s2a("Kolafa, J"), &
                         title="Time-reversible always stable predictor-corrector method for "// &
                         "molecular dynamics of polarizable molecules", &
                         source="JOURNAL OF COMPUTATIONAL CHEMISTRY", volume="25", issue="3", pages="335-342", &
                         year=2004, month=2, doi="10.1002/jcc.10385")

      CALL add_reference(key=Kuhne2007, &
                         authors=s2a("Kuhne, TD", "Krack, M", "Mohamed, FR", "Parrinello, M"), &
                         title="Efficient and accurate Car-Parrinello-like approach to "// &
                         "Born-Oppenheimer molecular dynamics", &
                         source="PHYSICAL REVIEW LETTERS", volume="98", issue="6", pages="066401", &
                         year=2007, month=2, day=9, doi="10.1103/PhysRevLett.98.066401")

      CALL add_reference(key=Rengaraj2020, &
                         authors=s2a("Rengaraj, V", "Lass, M", "Plessl, C", "Kuhne, TD"), &
                         title="Accurate Sampling with Noisy Forces from Approximate Computing", &
                         source="COMPUTATION", volume="8", issue="2", pages="39", &
                         year=2020, month=4, doi="10.3390/computation8020039")

      CALL add_reference(key=Kunert2003, &
                         authors=s2a("Kunert, T", "Schmidt, R"), &
                         title="Non-adiabatic quantum molecular dynamics: "// &
                         "General formalism and case study H-2(+) in strong laser fields", &
                         source="EUROPEAN PHYSICAL JOURNAL D", volume="25", issue="1", &
                         year=2003, month=7, doi="10.1140/epjd/e2003-00086-8")

      CALL add_reference(key=Ceriotti2009, &
                         authors=s2a("Ceriotti, M", "Bussi, G", "Parrinello, M"), &
                         title="Langevin equation with colored noise for constant-temperature "// &
                         "molecular dynamics simulations", &
                         source="PHYSICAL REVIEW LETTERS", volume="102", issue="2", pages="020601", &
                         year=2009, month=1, day=16, doi="10.1103/PhysRevLett.102.020601")

      CALL add_reference(key=Ceriotti2009b, &
                         authors=s2a("Ceriotti, M", "Bussi, G", "Parrinello, M"), &
                         title="Nuclear Quantum Effects in Solids Using a Colored-Noise Thermostat", &
                         source="PHYSICAL REVIEW LETTERS", volume="103", issue="3", pages="030603", &
                         year=2009, month=7, day=17, doi="10.1103/PhysRevLett.103.030603")

      CALL add_reference(key=Guidon2009, &
                         authors=s2a("Guidon, M", "Hutter, J", "VandeVondele, J"), &
                         title="Robust Periodic Hartree-Fock Exchange for Large-Scale Simulations "// &
                         "Using Gaussian Basis Sets", &
                         source="JOURNAL OF CHEMICAL THEORY AND COMPUTATION", volume="5", issue="11", pages="3010-3021", &
                         year=2009, month=11, doi="10.1021/ct900494g")

      CALL add_reference(key=BarducBus2008, &
                         authors=s2a("Barducci, A", "Bussi, G", "Parrinello, M"), &
                         title="Well-Tempered Metadynamics: A Smoothly Converging and Tunable Free-Energy Method", &
                         source="PHYSICAL REVIEW LETTERS", volume="100", issue="2", pages="020603", &
                         year=2008, month=1, day=18, doi="10.1103/PhysRevLett.100.020603")

      CALL add_reference(key=Guidon2010, &
                         authors=s2a("Guidon, M", "Hutter, J", "VandeVondele, J"), &
                         title="Auxiliary Density Matrix Methods for Hartree-Fock Exchange Calculations", &
                         source="JOURNAL OF CHEMICAL THEORY AND COMPUTATION", volume="6", issue="8", pages="2348-2364", &
                         year=2010, month=8, doi="10.1021/ct1002225")

      CALL add_reference(key=Marques2012, &
                         authors=s2a("Marques, MAL", "Oliveira, MJT", "Burnus, T"), &
                         title="LIBXC: A library of exchange and correlation functionals for density functional theory", &
                         source="COMPUTER PHYSICS COMMUNICATIONS", volume="183", issue="10", pages="2272-2281", &
                         year=2012, month=10, doi="10.1016/j.cpc.2012.05.007")

      CALL add_reference(key=Lehtola2018, &
                         authors=s2a("Lehtola, S", "Steigemann, C", "Oliveira, MJT", "Marques, MAL"), &
                         title="Recent developments in libxc - A comprehensive library of "// &
                         "functionals for density functional theory", &
                         source="SoftwareX", volume="7", pages="1-5", &
                         year=2018, month=1, doi="10.1016/j.softx.2017.11.002")

      CALL add_reference(key=Jones2011, &
                         authors=s2a("Jones, Andrew", "Leimkuhler, Ben"), &
                         title="Adaptive stochastic methods for sampling driven molecular systems", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="135", issue="8", pages="084125", &
                         year=2011, month=8, day=28, doi="10.1063/1.3626941")

      CALL add_reference(key=Bernstein2012, &
                         authors=s2a("Bernstein, Noam", "Varnai, Csilla", "Solt, Ivan", "Winfield, Steven A", &
                                     "Payne, Mike C", "Simon, Istvan", "Fuxreiter, Monika", "Csanyi, Gabor"), &
                         title="QM/MM simulation of liquid water with an adaptive quantum region", &
                         source="PHYSICAL CHEMISTRY CHEMICAL PHYSICS", volume="14", issue="2", pages="646-656", &
                         year=2012, doi="10.1039/c1cp22600b")

      CALL add_reference(key=Bernstein2009, &
                         authors=s2a("Bernstein, N", "Kermode, J R", "Csanyi, G"), &
                         title="Hybrid atomistic simulation methods for materials systems", &
                         source="REPORTS ON PROGRESS IN PHYSICS", volume="72", issue="2", pages="026501", &
                         year=2009, month=2, day=1, doi="10.1088/0034-4885/72/2/026501")

      CALL add_reference(key=Dick1958, &
                         authors=s2a("Dick, BG", "Overhauser, AW"), &
                         title="Theory of the Dielectric Constants of Alkali Halide Crystals", &
                         source="Phys. Rev.", volume="112", issue="1", &
                         year=1958, month=10, day=1, doi="10.1103/PhysRev.112.90")

      CALL add_reference(key=Mitchell1993, &
                         authors=s2a("Mitchell, PJ", "Fincham D"), &
                         title="Shell model simulations by adiabatic dynamics", &
                         source="J. Phys.: Condens. Matter", volume="5", issue="8", &
                         year=1993, month=2, day=22, doi="10.1088/0953-8984/5/8/006")

      CALL add_reference(key=Devynck2012, &
                         authors=s2a("Devynck, F", "Iannuzzi, M", "Krack, M"), &
                         title="Frenkel pair recombinations in UO2: Importance of explicit "// &
                         "description of polarizability in core-shell molecular dynamics simulations", &
                         source="Phys. Rev. B", volume="85", issue="18", &
                         year=2012, month=5, day=15, doi="10.1103/PhysRevB.85.184103")

      CALL add_reference(key=VandeVondele2012, &
                         authors=s2a("VandeVondele, J", "Borstnik, U", "Hutter, J"), &
                         title="Linear Scaling Self-Consistent Field Calculations with "// &
                         "Millions of Atoms in the Condensed Phase", &
                         source="JOURNAL OF CHEMICAL THEORY AND COMPUTATION", volume="8", issue="10", pages="3565-3573", &
                         year=2012, month=10, doi="10.1021/ct200897x")

      CALL add_reference(key=Niklasson2003, &
                         authors=s2a("Niklasson, AMN", "Tymczak, CJ", "Challacombe, M"), &
                         title="Trace resetting density matrix purification in O(N) self-consistent-field theory", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="118", issue="19", pages="8611-8620", &
                         year=2003, month=5, day=15, doi="10.1063/1.1559913")

      CALL add_reference(key=Shao2003, &
                         authors=s2a("Shao, Y", "Saravanan, C", "Head-Gordon, M", "White, CA"), &
                         title="Curvy steps for density matrix-based energy minimization: "// &
                         "Application to large-scale self-consistent-field calculations", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="118", issue="14", pages="6144-6151", &
                         year=2003, month=4, day=8, doi="10.1063/1.1558476")

      CALL add_reference(key=VandeVondele2002, &
                         authors=s2a("VandeVondele, J", "Rothlisberger, U"), &
                         title="Canonical adiabatic free energy sampling (CAFES): "// &
                         "A novel method for the exploration of free energy surfaces", &
                         source="JOURNAL OF PHYSICAL CHEMISTRY B", volume="106", issue="1", pages="203-208", &
                         year=2002, month=1, day=10, doi="10.1021/jp013346k")

      CALL add_reference(key=Dion2004, &
                         authors=s2a("Dion, M", "Rydberg, H", "Schroder, E", "Langreth, DC", "Lundqvist, BI"), &
                         title="Van der Waals density functional for general geometries", &
                         source="Phys. Rev. Lett.", volume="92", issue="24", pages="246401", &
                         year=2004, month=6, day=18, doi="10.1103/PhysRevLett.92.246401")

      CALL add_reference(key=RomanPerez2009, &
                         authors=s2a("Roman-Perez, G", "Soler, JM"), &
                         title="Efficient Implementation of a van der Waals Density Functional: "// &
                         "Application to Double-Wall Carbon Nanotubes", &
                         source="Phys. Rev. Lett.", volume="103", issue="9", pages="096102", &
                         year=2009, month=8, day=28, doi="10.1103/PhysRevLett.103.096102")

      CALL add_reference(key=DelBen2012, &
                         authors=s2a("Del Ben, M", "Hutter, J", "VandeVondele, J"), &
                         title="Second-Order Moller-Plesset Perturbation Theory in the Condensed Phase: "// &
                         "An Efficient and Massively Parallel Gaussian and Plane Waves Approach", &
                         source="JOURNAL OF CHEMICAL THEORY AND COMPUTATION", volume="8", issue="11", pages="4177-4188", &
                         year=2012, month=11, doi="10.1021/ct300531w")

      CALL add_reference(key=Sabatini2013, &
                         authors=s2a("Sabatini, R", "Gorni, T", "de Gironcoli, S"), &
                         title="Nonlocal van der Waals density functional made simple and efficient", &
                         source="Phys. Rev. B", volume="87", issue="4", pages="041108(R)", &
                         year=2013, month=1, day=15, doi="10.1103/PhysRevB.87.041108")

      CALL add_reference(key=Walewski2014, &
                         authors=s2a("Walewski, L", "Forbert, H", "Marx, D"), &
                         title="Reactive path integral quantum simulations of molecules solvated in superfluid helium", &
                         source="COMPUTER PHYSICS COMMUNICATIONS", volume="185", issue="3", pages="884-899", &
                         year=2014, month=3, doi="10.1016/j.cpc.2013.12.011")

      CALL add_reference(key=Delben2013, &
                         authors=s2a("Del Ben, M", "Hutter, J", "VandeVondele, J"), &
                         title="Electron Correlation in the Condensed Phase from a Resolution of "// &
                         "Identity Approach Based on the Gaussian and Plane Waves Scheme", &
                         source="JOURNAL OF CHEMICAL THEORY AND COMPUTATION", volume="9", issue="6", pages="2654-2671", &
                         year=2013, month=6, doi="10.1021/ct4002202")

      CALL add_reference(key=Kikuchi2009, &
                         authors=s2a("Kikuchi, Y", "Imamura, Y", "Nakai, H"), &
                         title="One-Body Energy Decomposition Schemes Revisited: Assessment of "// &
                         "Mulliken-, Grid-, and Conventional Energy Density Analyses", &
                         source="INTERNATIONAL JOURNAL OF QUANTUM CHEMISTRY", volume="109", issue="11", pages="2464-2473", &
                         year=2009, month=9, doi="10.1002/qua.22017")

      CALL add_reference(key=Putrino2000, &
                         authors=s2a("Putrino, A", "Sebastiani, D", "Parrinello, M"), &
                         title="Generalized Variational Density Functional Perturbation Theory", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="113", issue="17", pages="7102-7109", &
                         year=2000, month=11, doi="10.1063/1.1312830")

      CALL add_reference(key=Tran2013, &
                         authors=s2a("Tran, F", "Hutter, J"), &
                         title="Nonlocal van der Waals functionals: The case of rare-gas dimers and solids", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="138", issue="20", pages="204103", &
                         year=2013, month=5, day=28, doi="10.1063/1.4807332")

      CALL add_reference(key=Putrino2002, &
                         authors=s2a("Putrino, A", "Parrinello, M"), &
                         title="Anharmonic Raman Spectra in High-Pressure Ice from Ab Initio Simulations", &
                         source="PHYSICAL REVIEW LETTERS", volume="88", issue="17", pages="176401", &
                         year=2002, month=4, day=29, doi="10.1103/PhysRevLett.88.176401")

      CALL add_reference(key=Sebastiani2001, &
                         authors=s2a("Sebastiani, D", "Parrinello, M"), &
                         title="A New Ab Initio Approach for NMR Chemical Shifts in Periodic Systems", &
                         source="THE JOURNAL OF PHYSICAL CHEMISTRY A", volume="105", issue="10", pages="1951-1958", &
                         year=2001, month=3, doi="10.1021/jp002807j")

      CALL add_reference(key=Weber2009, &
                         authors=s2a("Weber, V", "Iannuzzi, M", "Giani, S", "Hutter, J", "Declerck, R", "Waroduier, M"), &
                         title="Magnetic Linear Response Properties Calculations with the "// &
                         "Gaussian and Augmanted-Plane-Wave Method", &
                         source="THE JOURNAL OF CHEMICAL PHYSICS", volume="131", issue="1", pages="014106", &
                         year=2009, month=7, day=7, doi="10.1063/1.3156803")

      CALL add_reference(key=Golze2013, &
                         authors=s2a("Golze, D", "Iannuzzi, M", "Nguyen, M-T", "Passerone, D", "Hutter, J"), &
                         title="Simulation of Adsorption Processes at Metallic Interfaces: "// &
                         "An Image Charge Augmented QM/MM Approach", &
                         source="Journal of Chemical Theory and Computation", volume="9", issue="11", pages="5086-5097", &
                         year=2013, month=11, day=12, doi="10.1021/ct400698y")

      CALL add_reference(key=Golze2015, &
                         authors=s2a("Golze, D", "Hutter, J", "Iannuzzi, M"), &
                         title="Wetting of water on hexagonal boron nitride@Rh(111): "// &
                         "a QM/MM model based on atomic charges derived for nano-structured substrates", &
                         source="Physical Chemistry Chemical Physics", volume="17", issue="22", pages="14307-14316", &
                         year=2015, month=6, day=14, doi="10.1039/C4CP04638B")

      CALL add_reference(key=Golze2017a, &
                         authors=s2a("Golze, D", "Benedikter, N", "Iannuzzi, M", "Wilhelm, J", "Hutter, J"), &
                         title="Fast evaluation of solid harmonic Gaussian integrals for local "// &
                         "resolution-of-the-identity methods and range-separated hybrid functionals", &
                         source="The Journal of Chemical Physics", volume="146", issue="3", pages="034105", &
                         year=2017, month=1, day=17, doi="10.1063/1.4973510")

      CALL add_reference(key=Golze2017b, &
                         authors=s2a("Golze, D", "Iannuzzi, M", "Hutter, J"), &
                         title="Local Fitting of the Kohn-Sham Density in a Gaussian and "// &
                         "Plane Waves Scheme for Large-Scale Density Functional Theory Simulations", &
                         source="Journal of Chemical Theory and Computation", volume="13", issue="5", pages="2202-2214", &
                         year=2017, month=5, day=17, doi="10.1021/acs.jctc.7b00148")

      CALL add_reference(key=Fattebert2002, &
                         authors=s2a("Fattebert, JL", "Gygi, F"), &
                         title="Density functional theory for efficient ab initio "// &
                         "molecular dynamics simulations in solution", &
                         source="J. Comput. Chem.", volume="23", issue="6", &
                         year=2002, month=3, day=18, doi="10.1002/jcc.10069")

      CALL add_reference(key=Andreussi2012, &
                         authors=s2a("Andreussi, O", "Dabo, I", "Marzari, N"), &
                         title="Revised self-consistent continuum solvation in electronic-structure calculations", &
                         source="J. Chem. Phys.", volume="136", issue="6", pages="064102", &
                         year=2012, month=2, day=8, doi="10.1063/1.3676407")

      CALL add_reference(key=Tuckerman1992, &
                         authors=s2a("TUCKERMAN, M", "BERNE, BJ", "MARTYNA, GJ"), &
                         title="REVERSIBLE MULTIPLE TIME SCALE MOLECULAR-DYNAMICS", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="97", issue="3", pages="1990-2001", &
                         year=1992, month=8, day=1, doi="10.1063/1.463137")

      CALL add_reference(key=Goedecker2004, &
                         authors=s2a("Goedecker, S"), &
                         title="Minima hopping: An efficient search method for the global minimum of "// &
                         "the potential energy surface of complex molecular systems", &
                         source="Journal of Chemical Physics", volume="120", issue="21", pages="9911-9917", &
                         year=2004, doi="10.1063/1.1724816")

      CALL add_reference(key=Khaliullin2007, &
                         authors=s2a("Khaliullin, RZ", "Cobar, EA", "Lochan, RC", &
                                     "Bell, AT", "Head-Gordon, M"), &
                         title="Unravelling the origin of intermolecular interactions using "// &
                         "absolutely localized molecular orbitals", &
                         source="Journal of Physical Chemistry A", volume="111", issue="36", pages="8753-8765", &
                         year=2007, doi="10.1021/jp073685z")

      CALL add_reference(key=Khaliullin2008, &
                         authors=s2a("Khaliullin, RZ", "Bell, AT", "Head-Gordon, M"), &
                         title="Analysis of charge transfer effects in molecular complexes "// &
                         "based on absolutely localized molecular orbitals", &
                         source="Journal of Chemical Physics", volume="128", issue="18", pages="184112", &
                         year=2008, doi="10.1063/1.2912041")

      CALL add_reference(key=Khaliullin2013, &
                         authors=s2a("Khaliullin, RZ", "VandeVondele, J", "Hutter, J"), &
                         title="Efficient Linear-Scaling Density Functional Theory for Molecular Systems", &
                         source="JOURNAL OF CHEMICAL THEORY AND COMPUTATION", volume="9", issue="10", pages="4421-4427", &
                         year=2013, month=10, doi="10.1021/ct400595k")

      CALL add_reference(key=Hutter2014, &
                         authors=s2a("Hutter, J", "Iannuzzi, M", "Schiffmann, F", "VandeVondele, J"), &
                         title="CP2K: atomistic simulations of condensed matter systems", &
                         source="WIREs Comput Mol Sci.", volume="4", issue="1", pages="15-25", &
                         year=2014, month=1, doi="10.1002/wcms.1159")

      CALL add_reference(key=Kantorovich2008, &
                         authors=s2a("Kantorovich, L"), &
                         title="Generalized Langevin equation for solids. I. Rigorous derivation and main properties", &
                         source="PHYSICAL REVIEW B", volume="78", issue="9", pages="094304", &
                         year=2008, month=9, doi="10.1103/PhysRevB.78.094304")

      CALL add_reference(key=Kantorovich2008a, &
                         authors=s2a("Kantorovich, L", "Rompotis, N"), &
                         title="Generalized Langevin equation for solids. II. Stochastic boundary "// &
                         "conditions for nonequilibrium molecular dynamics simulations", &
                         source="PHYSICAL REVIEW B", volume="78", issue="9", pages="094305", &
                         year=2008, month=9, doi="10.1103/PhysRevB.78.094305")

      CALL add_reference(key=Niklasson2014, &
                         authors=s2a("Rubensson, E.", "Niklasson, A."), &
                         title="Interior Eigenvalues from Density Matrix Expansions in "// &
                         "Quantum Mechanical Molecular Dynamics", &
                         source="SIAM Journal on Scientific Computing", volume="36", issue="2", pages="B147-B170", &
                         year=2014, month=3, doi="10.1137/130911585")

      CALL add_reference(key=Borstnik2014, &
                         authors=s2a("Borstnik, U", "VandeVondele, J", "Weber, V", "Hutter, J"), &
                         title="Sparse matrix multiplication: The distributed block-compressed sparse row library", &
                         source="PARALLEL COMPUTING", volume="40", issue="5-6", pages="47-58", &
                         year=2014, month=5, doi="10.1016/j.parco.2014.03.012")

      CALL add_reference(key=Rayson2009, &
                         authors=s2a("Rayson, M. J.", "Briddon, P. R."), &
                         title="Highly efficient method for Kohn-Sham density functional "// &
                         "calculations of 500-10 000 atom systems", &
                         source="PHYSICAL REVIEW B", volume="80", issue="20", pages="205104", &
                         year=2009, month=11, doi="10.1103/PhysRevB.80.205104")

      CALL add_reference(key=Merlot2014, &
                         authors=s2a("Merlot, P", "Izsak, R", "Borgo, A", &
                                     "Kjaergaard, T", "Helgaker, T", "Reine, S"), &
                         title="Charge-constrained auxiliary-density-matrix methods for the "// &
                         "Hartree-Fock exchange contribution", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="141", issue="9", pages="094104", &
                         year=2014, month=9, doi="10.1063/1.4894267")

      CALL add_reference(key=Lin2009, &
                         authors=s2a("Lin, L", "Lu, JF", "Ying, LX", "Car, R", "E, WN"), &
                         title="FAST ALGORITHM FOR EXTRACTING THE DIAGONAL OF THE INVERSE MATRIX "// &
                         "WITH APPLICATION TO THE ELECTRONIC STRUCTURE ANALYSIS OF METALLIC SYSTEMS", &
                         source="COMMUNICATIONS IN MATHEMATICAL SCIENCES", volume="7", issue="3", pages="755-777", &
                         year=2009, month=9, doi="10.4310/CMS.2009.v7.n3.a12")

      CALL add_reference(key=Lin2013, &
                         authors=s2a("Lin, Lin", "Chen, Mohan", "Yang, Chao", "He, Lixin"), &
                         title="Accelerating atomic orbital-based electronic structure calculation "// &
                         "via pole expansion and selected inversion", &
                         source="JOURNAL OF PHYSICS-CONDENSED MATTER", volume="25", issue="29", pages="295501", &
                         year=2013, month=7, day=24, doi="10.1088/0953-8984/25/29/295501")

      CALL add_reference(key=DelBen2015, &
                         authors=s2a("Del Ben, M", "Schuett, O", "Wentz, T", "Messmer, P", "Hutter, J", "VandeVondele, J"), &
                         title="Enabling simulation at the fifth rung of DFT: "// &
                         "Large scale RPA calculations with excellent time to solution", &
                         source="COMPUTER PHYSICS COMMUNICATIONS", volume="187", pages="120-129", &
                         year=2015, month=2, doi="10.1016/j.cpc.2014.10.021")

      CALL add_reference(key=Souza2002, &
                         authors=s2a("Souza, I", "Iniguez, J", "Vanderbilt, D"), &
                         title="First-principles approach to insulators in finite electric fields", &
                         source="PHYSICAL REVIEW LETTERS", volume="89", issue="11", pages="117602", &
                         year=2002, month=9, doi="10.1103/PhysRevLett.89.117602")

      CALL add_reference(key=Umari2002, &
                         authors=s2a("Umari, P", "Pasquarello, A"), &
                         title="Ab initio molecular dynamics in a finite homogeneous electric field", &
                         source="PHYSICAL REVIEW LETTERS", volume="89", issue="15", pages="157602", &
                         year=2002, month=10, doi="10.1103/PhysRevLett.89.157602")

      CALL add_reference(key=Stengel2009, &
                         authors=s2a("Stengel, Massimiliano", "Spaldin, Nicola A.", "Vanderbilt, David"), &
                         title="Ab initio molecular dynamics in a finite homogeneous electric field "// &
                         "Electric displacement as the fundamental variable in electronic-structure calculations", &
                         source="NATURE PHYSICS", volume="5", issue="4", pages="304-308", &
                         year=2009, month=4, doi="10.1038/NPHYS1185")

      CALL add_reference(key=Luber2014, &
                         authors=s2a("Luber, S", "Iannuzzi, M", "Hutter, J"), &
                         title="Raman spectra from ab initio molecular dynamics and its "// &
                         "application to liquid S-methyloxirane", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="141", issue="9", pages="094503", &
                         year=2014, month=9, day=7, doi="10.1063/1.4894425")

      CALL add_reference(key=Berghold2011, &
                         authors=s2a("Berghold, G", "Parrinello, M", "Hutter, J"), &
                         title="Polarized atomic orbitals for linear scaling methods", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="116", issue="5", pages="1800-1810", &
                         year=2002, month=2, day=1, doi="10.1063/1.1431270")

      CALL add_reference(key=DelBen2015b, &
                         authors=s2a("Del Ben, M", "Hutter, J", "VandeVondele, J"), &
                         title="Forces and stress in second order Moller-Plesset perturbation theory "// &
                         "for condensed phase systems within the resolution-of-identity "// &
                         "Gaussian and plane waves approach", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="143", issue="10", pages="102803", &
                         year=2015, month=9, day=14, doi="10.1063/1.4919238")

      CALL add_reference(key=Campana2009, &
                         authors=s2a("Campana, C", "Mussard, B", "Woo, T K"), &
                         title="Electrostatic Potential Derived Atomic Charges for "// &
                         "Periodic Systems Using a Modified Error Functional", &
                         source="JOURNAL OF CHEMICAL THEORY AND COMPUTATION", volume="5", issue="10", pages="2866", &
                         year=2009, month=10, day=13, doi="10.1021/ct9003405")

      CALL add_reference(key=Schiffmann2015, &
                         authors=s2a("Schiffmann, F", "VandeVondele, J"), &
                         title="Efficient preconditioning of the electronic structure problem in "// &
                         "large scale ab initio molecular dynamics simulations", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="142", issue="24", pages="244117", &
                         year=2015, month=6, day=28, doi="10.1063/1.4922988")

      CALL add_reference(key=Bruck2014, &
                         authors=s2a("Bruck, S", "Calderara, M", "Bani-Hashemian, MH", "VandeVondele, J", "Luisier, M"), &
                         title="Towards ab-initio simulations of nanowire field-effect transistors", &
                         source="2014 INTERNATIONAL WORKSHOP ON COMPUTATIONAL ELECTRONICS (IWCE)", &
                         year=2014, doi="10.1109/IWCE.2014.6865831")

      CALL add_reference(key=Rappe1992, &
                         authors=s2a("Rappe, AK", "Casewit, CJ", "Colwell, KS", "Goddard, WA", "Skiff, WM"), &
                         title="UFF, a full periodic table force field for molecular mechanics "// &
                         "and molecular dynamics simulations", &
                         source="JOURNAL OF THE AMERICAN CHEMICAL SOCIETY", volume="114", issue="25", pages="10024-10035", &
                         year=1992, month=12, day=1, doi="10.1021/ja00051a040")

      CALL add_reference(key=Monkhorst1976, &
                         authors=s2a("Monkhorst, HJ", "Pack, JD"), &
                         title="Special points for Brillouin-zone integrations", &
                         source="PHYSICAL REVIEW B", volume="13", issue="12", pages="5188-5192", &
                         year=1976, doi="10.1103/PhysRevB.13.5188")

      CALL add_reference(key=MacDonald1978, &
                         authors=s2a("MacDonald, AH"), &
                         title="Comment on special points for Brillouin-zone integrations", &
                         source="PHYSICAL REVIEW B", volume="18", issue="10", pages="5897-5899", &
                         year=1978, doi="10.1103/PhysRevB.18.5897")

      CALL add_reference(key=Gilbert2008, &
                         authors=s2a("Gilbert, ATB", "Besley, NA", "Gill, PMW"), &
                         title="Self-consistent field calculations of excited states using the "// &
                         "maximum overlap method (MOM)", &
                         source="THE JOURNAL OF PHYSICAL CHEMISTRY A", volume="112", pages="13164-13171", &
                         year=2008, month=8, day=26, doi="10.1021/jp801738f")

      CALL add_reference(key=Barca2018, &
                         authors=s2a("Barca, GMJ", "Gilbert, ATB", "Gill, PMW"), &
                         title="Simple models for difficult electronic excitations", &
                         source="JOURNAL OF CHEMICAL THEORY AND COMPUTATION", volume="14", pages="1501-1509", &
                         year=2018, doi="10.1021/acs.jctc.7b00994")

      CALL add_reference(key=Schonherr2014, &
                         authors=s2a("Schonherr, M", "Slater, B", "Hutter, J", "VandeVondele, J"), &
                         title="Dielectric Properties of Water Ice, the Ice Ih/XI Phase Transition, "// &
                         "and an Assessment of Density Functional Theory", &
                         source="JOURNAL OF PHYSICAL CHEMISTRY B", volume="118", issue="2", pages="590-596", &
                         year=2014, month=1, day=16, doi="10.1021/jp4103355")

      CALL add_reference(key=Ceriotti2014, &
                         authors=s2a("Ceriotti, M", "More, J", "Manolopoulos, DE"), &
                         title="i-PI: A Python interface for ab initio path integral molecular dynamics simulations", &
                         source="COMPUTER PHYSICS COMMUNICATIONS", volume="185", issue="3", pages="1019-1026", &
                         year=2014, month=3, doi="10.1016/j.cpc.2013.10.027")

      CALL add_reference(key=BaniHashemian2016, &
                         authors=s2a("Bani-Hashemian, MH", "Bruck, S", "Luisier, M", "VandeVondele, J"), &
                         title="A generalized Poisson solver for first-principles device simulations", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="144", issue="4", pages="044113", &
                         year=2016, month=1, day=28, doi="10.1063/1.4940796")

      CALL add_reference(key=Kapil2016, &
                         authors=s2a("Kapil, V", "VandeVondele, J", "Ceriotti, M"), &
                         title="Accurate molecular dynamics and nuclear quantum effects at low cost by multiple steps "// &
                         "in real and imaginary time: Using density functional theory to accelerate wavefunction methods", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="144", issue="5", pages="054111", &
                         year=2016, month=2, day=7, doi="10.1063/1.4941091")

      CALL add_reference(key=Heinzmann1976, &
                         authors=s2a("Heinzmann, R", "Ahlrichs, R"), &
                         title="Population analysis based on occupation numbers of modified atomic orbitals (MAOs)", &
                         source="Theoret. Chim. Acta", volume="42", issue="1", pages="33-45", &
                         year=1976, doi="10.1007/BF00548289")

      CALL add_reference(key=Ehrhardt1985, &
                         authors=s2a("Ehrhardt, C", "Ahlrichs, R"), &
                         title="Population analysis based on occupation numbers II. Relationship between shared "// &
                         "electron numbers and bond energies and characterization of hypervalent contributions", &
                         source="Theoret. Chim. Acta", volume="68", issue="3", pages="231-245", &
                         year=1985, doi="10.1007/BF00526774")

      CALL add_reference(key=Rybkin2016, &
                         authors=s2a("Rybkin, VV", "VandeVondele, J"), &
                         title="Spin-Unrestricted Second-Order Moller-Plesset (MP2) Forces for the Condensed Phase: "// &
                         "From Molecular Radicals to F-Centers in Solids", &
                         source="JOURNAL OF CHEMICAL THEORY AND COMPUTATION", volume="12", issue="5", pages="2214-2223", &
                         year=2016, month=5, doi="10.1021/acs.jctc.6b00015")

      CALL add_reference(key=West2006, &
                         authors=s2a("West, D", "Estreicher, S. K."), &
                         title="First-Principles Calculations of Vibrational Lifetimes and Decay Channels:"// &
                         " Hydrogen-Related Modes in Si", &
                         source="PHYSICAL REVIEW LETTERS", volume="96", issue="11", pages="115504", &
                         year=2006, doi="10.1103/PhysRevLett.96.115504")

      CALL add_reference(key=Bates2013, &
                         authors=s2a("Bates, JE"), &
                         title="Communication: Random phase approximation renormalized many-body perturbation theory", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="139", issue="17", pages="171103", &
                         year=2013, month=11, day=7, doi="10.1063/1.4827254")

      CALL add_reference(key=Andermatt2016, &
                         authors=s2a("Andermatt, S", "Cha, J", "Schiffmann, F", "VandeVondele, J"), &
                         title="Combining Linear-Scaling DFT with Subsystem DFT in Born Oppenheimer "// &
                         "and Ehrenfest Molecular Dynamics Simulations: From Molecules to a Virus in Solution", &
                         source="JOURNAL OF CHEMICAL THEORY AND COMPUTATION", volume="12", issue="7", pages="3214-3227", &
                         year=2016, month=7, doi="10.1021/acs.jctc.6b00398")

      CALL add_reference(key=Zhu2016, &
                         authors=s2a("Zhu, L", "Amsler, M", "Fuhrer, T", "Schaefer, B", "Faraji, S", "Rostami, S", &
                                     "Ghasemi, SA", "Sadeghi, A", "Grauzinyte, M", "Wolverton, C", "Goedecker, S"), &
                         title="A fingerprint based metric for measuring similarities of crystalline structures", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="144", issue="3", pages="034203", &
                         year=2016, month=1, day=21, doi="10.1063/1.4940026")

      CALL add_reference(key=Schuett2016, &
                         authors=s2a("Schuett, Ole", "Messmer, Peter", "Hutter, Juerg", "VandeVondele, Joost"), &
                         title="GPU-Accelerated Sparse Matrix-Matrix Multiplication for Linear Scaling Density Functional Theory", &
                         source="Electronic Structure Calculations on Graphics Processing Units", pages="173-190", &
                         year=2016, doi="10.1002/9781118670712.ch8")

      CALL add_reference(key=Schran2020a, &
                         authors=s2a("Schran, C", "Behler, J", "Marx, D"), &
                         title="Automated Fitting of Neural Network Potentials at Coupled Cluster Accuracy: "// &
                         "Protonated Water Clusters as Testing Ground", &
                         source="Journal of Chemical Theory and Computation", volume="16", issue="1", &
                         year=2020, doi="10.1021/acs.jctc.9b00805")

      CALL add_reference(key=Schran2020b, &
                         authors=s2a("Schran, C", "Brezina, K", "Marsalek, O"), &
                         title="Committee neural network potentials control generalization errors and enable active learning", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="153", issue="10", &
                         year=2020, doi="10.1063/5.0016004")

      CALL add_reference(key=Behler2007, &
                         authors=s2a("Behler, J", "Parrinello, M"), &
                         title="Generalized neural-network representation of high-dimensional potential-energy surfaces", &
                         source="PHYSICAL REVIEW LETTERS", volume="98", issue="14", pages="146401", &
                         year=2007, doi="10.1103/PhysRevLett.98.146401")

      CALL add_reference(key=Behler2011, &
                         authors=s2a("Behler, J"), &
                         title="Atom-centered symmetry functions for constructing high-dimensional neural network potentials", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="134", issue="7", &
                         year=2011, doi="10.1063/1.3553717")

      CALL add_reference(key=Lu2004, &
                         authors=s2a("Lu, WC", "Wang, CZ", "Schmidt, MW", "Bytautas, L", "Ho, KM", "Ruedenberg, K"), &
                         title="Molecule intrinsic minimal basis sets. I. Exact resolution of ab initio "// &
                         "optimized molecular orbitals in terms of deformed atomic minimal-basis orbitals", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="120", issue="6", pages="2629-2637", &
                         year=2004, month=2, day=8, doi="10.1063/1.1638731")

      CALL add_reference(key=Migliore2009, &
                         authors=s2a("Migliore, A"), &
                         title="Full-electron calculation of effective electronic couplings and", &
                         source="The Journal of Chemical Physics", volume="131", issue="11", &
                         year=2009, doi="10.1063/1.3232007")

      CALL add_reference(key=Mavros2015, &
                         authors=s2a("Mavros, MG", " Van Voorhis, T"), &
                         title="Communication: CDFT-CI couplings can be unreliable when there is fractional charge transfer", &
                         source="The Journal of Chemical Physics", volume="143", issue="23", &
                         year=2015, doi="10.1063/1.4938103")

      CALL add_reference(key=Becke1988b, &
                         authors=s2a("Becke, AD"), &
                         title="A multicenter numerical integration scheme for polyatomic molecules", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="88", issue="4", pages="2547-2553", &
                         year=1988, doi="10.1063/1.454033")

      CALL add_reference(key=Holmberg2017, &
                         authors=s2a("Holmberg, N", "Laasonen, K"), &
                         title="Efficient Constrained Density Functional Theory Implementation for Simulation of "// &
                         "Condensed Phase Electron Transfer Reactions", &
                         source="Journal of Chemical Theory and Computation", volume="13", issue="2", &
                         year=2017, doi="10.1021/acs.jctc.6b01085")

      CALL add_reference(key=Marek2014, &
                         authors=s2a("Marek, A", "Blum, V", "Johanni, R", "Havu, V", "Lang, B", &
                                     "Auckenthaler, T", "Heinecke, A", "Bungartz, H", "Lederer, H"), &
                         title="The ELPA library: scalable parallel eigenvalue solutions for "// &
                         "electronic structure theory and computational science", &
                         source="Journal of Physics: Condensed Matter", volume="26", issue="21", &
                         year=2014, doi="10.1088/0953-8984/26/21/213201")

      CALL add_reference(key=VanVoorhis2015, &
                         authors=s2a("VanVoorhis, T", "Welborn, M", "Chen, J", "Wang, L"), &
                         title="Why many semiempirical molecular orbital fail for liquid water and how to fix them", &
                         source="Journal of Computational Chemistry", volume="36", issue="12", &
                         year=2015, doi="10.1002/jcc.23887")

      CALL add_reference(key=Stoychev2016, &
                         authors=s2a("Stoychev, Georgi L.", "Auer, Alexander A.", "Neese, Frank"), &
                         title="Automatic Generation of Auxiliary Basis Sets", &
                         source="Journal of Chemical Theory and Computation", volume="13", issue="2", pages="554-562", &
                         year=2017, doi="10.1021/acs.jctc.6b01041")

      CALL add_reference(key=Kondov2007, &
                         authors=s2a("Kondov, Ivan", "Cizek, Martin", "Benesch, Claudia", "Wang, Haobin", "Thoss, Michael"), &
                         title="Quantum dynamics of photoinduced electron-transfer reactions in dye-semiconductor systems: "// &
                         "First-principles description and application to coumarin 343-TiO2", &
                         source="Journal of Physical Chemistry C", volume="111", issue="32", pages="11970-11981", &
                         year=2007, doi="10.1021/jp072217m")

      CALL add_reference(key=Futera2017, &
                         authors=s2a("Futera, Zdenek", "Blumberger, Jochen"), &
                         title="Electronic Couplings for Charge Transfer across Molecule/Metal and "// &
                         "Molecule/Semiconductor Interfaces: Performance of the Projector Operator-Based Diabatization Approach", &
                         source="Journal Physical Chemistry C", volume="121", issue="36", pages="19677-19689", &
                         year=2017, doi="10.1021/acs.jpcc.7b06566")

      CALL add_reference(key=Bailey2006, &
                         authors=s2a("Rocha, AR", "Garcia-Suarez, VM", "Bailey, S", "Lambert, C", "Ferrer, J", "Sanvito, S"), &
                         title="Spin and molecular electronics in atomically generated orbital landscapes", &
                         source="PHYSICAL REVIEW B", volume="73", pages="085414", &
                         year=2006, doi="10.1103/PhysRevB.73.085414")

      CALL add_reference(key=Papior2017, &
                         authors=s2a("Papior, N", "Lorente, N", "Frederiksen, T", "Garcia, A", "Brandbyge, M"), &
                         title="Improvements on non-equilibrium and transport Green function techniques", &
                         source="COMPUTER PHYSICS COMMUNICATIONS", volume="212", pages="8-24", &
                         year=2017, doi="10.1016/j.cpc.2016.09.022")

      CALL add_reference(key=Brieuc2016, &
                         authors=s2a("Brieuc, F", "Dammak, H", "Hayoun, M"), &
                         title="Quantum thermal Bath for Path Integral Molecular Dynamics Simulation", &
                         source="Journal of Chemical Theory and Computation", volume="12", &
                         year=2016, doi="10.1021/acs.jctc.5b01146")

      CALL add_reference(key=Huang2011, &
                         authors=s2a("Huang, C", "Pavone, M", "Carter, EA"), &
                         title="Quantum mechanical embedding theory based on a unique embedding potential", &
                         source="Journal of Chemical Physics", volume="134", pages="154110", &
                         year=2011, doi="10.1063/1.3577516")

      CALL add_reference(key=Heaton_Burgess2007, &
                         authors=s2a("Heaton-Burgess, T.", "Bulat, FA", "Yang, WT"), &
                         title="Optimized effective potentials in finite basis sets", &
                         source="PHYSICAL REVIEW LETTERS", volume="98", issue="25", pages="256401", &
                         year=2007, doi="10.1103/PhysRevLett.98.256401")

      CALL add_reference(key=Scheiber2018, &
                         authors=s2a("Scheiber, H", "Shi, Y", "Khaliullin, RZ"), &
                         title="Compact orbitals enable low-cost linear-scaling ab initio "// &
                         "molecular dynamics for weakly-interacting systems", &
                         source="The Journal of Chemical Physics", volume="148", pages="231103", &
                         year=2018, month=6, day=21, doi="10.1063/1.5029939")

      CALL add_reference(key=Schuett2018, &
                         authors=s2a("Schuett, O", "VandeVondele, J"), &
                         title="Machine Learning Adaptive Basis Sets for Efficient Large Scale "// &
                         "Density Functional Theory Simulation", &
                         source="Journal of Chemical Theory and Computation", volume="14", &
                         year=2018, doi="10.1021/acs.jctc.8b00378")

      CALL add_reference(key=Holmberg2018, &
                         authors=s2a("Holmberg, Nico", "Laasonen, Kari"), &
                         title="Diabatic model for electrochemical hydrogen evolution based on "// &
                         "constrained DFT configuration interaction", &
                         source="The Journal of Chemical Physics", volume="149", issue="10", pages="104702", &
                         year=2018, doi="10.1063/1.5038959")

      CALL add_reference(key=Togo2018, &
                         authors=s2a("Togo, Atsushi", "Tanaka, Isao"), &
                         title="Spglib : a software library for crystal symmetry search", &
                         source="arXiv", pages="1808.01590", &
                         year=2018)

      CALL add_reference(key=Staub2019, &
                         authors=s2a("Staub, R", "Iannuzzi, M", "Khaliullin, RZ", "Steinmann, SN"), &
                         title="Energy Decomposition Analysis for Metal Surface-Adsorbate Interactions "// &
                         "by Block Localized Wave Functions", &
                         source="Journal of Chemical Theory and Computation", volume="15", &
                         year=2019, doi="10.1021/acs.jctc.8b00957")

      CALL add_reference(key=Clabaut2020, &
                         authors=s2a("Clabaut, P", "Fleurat-Lessard, P", "Michel, C", "Steinmann, SN"), &
                         title="Ten Facets, One Force Field: The GAL19 Force Field for Water-Noble Metal Interfaces", &
                         source="Journal of Chemical Theory and Computation", &
                         year=2020, doi="10.1021/acs.jctc.0c00091")

      CALL add_reference(key=Clabaut2021, &
                         authors=s2a("Clabaut, P"), &
                         title="Solvation and adsorptions at the solid/water interface: Developments and applications", &
                         source="Ph.D. thesis at Université de Lyon", &
                         year=2021)

      CALL add_reference(key=Richters2018, &
                         authors=s2a("Richters, D", "Lass, M", "Walther, A", "Plessl, C", "Kuehne, T D"), &
                         title="A General Algorithm to Calculate the Inverse Principal p-th Root of "// &
                         "Symmetric Positive Definite Matrices", &
                         source="Communications in Computational Physics", volume="25", issue="2", pages="564-585", &
                         year=2018, month=10, doi="10.4208/cicp.OA-2018-0053")

      CALL add_reference(key=Kruse2012, &
                         authors=s2a("Kruse,Holger", "Grimme,Stefan"), &
                         title="A geometrical correction for the inter- and intra-molecular basis set superposition error "// &
                         "in Hartree-Fock and density functional theory calculations for large systems", &
                         source="The Journal of Chemical Physics", volume="136", issue="15", pages="154101", &
                         year=2012, doi="10.1063/1.3700154")

      CALL add_reference(key=Ren2011, &
                         authors=s2a("Ren,Xinguo", "Tkatchenko,Aleksandre", "Rinke,Patrick", "Scheffler,Matthias"), &
                         title="Beyond the Random-Phase Approximation for the Electron Correlation Energy: "// &
                         "The Importance of Single Excitations", &
                         source="PHYSICAL REVIEW LETTERS", volume="106", issue="15", pages="153003", &
                         year=2011, doi="10.1103/PhysRevLett.106.153003")

      CALL add_reference(key=Ren2013, &
                         authors=s2a("Ren,Xinguo", "Rinke,Patrick", "Scuseria,Gustavo", "Scheffler,Matthias"), &
                         title="Renormalized second-order perturbation theory for the electron correlation energy: "// &
                         "Concept, implementation, and benchmarks", &
                         source="PHYSICAL REVIEW B", volume="88", issue="3", pages="035120", &
                         year=2013, doi="10.1103/PhysRevB.88.035120")

      CALL add_reference(key=Martin2003, &
                         authors=s2a("Martin,Richard L."), &
                         title="Natural transition orbitals", &
                         source="The Journal of Chemical Physics", volume="118", issue="11", pages="4775-4777", &
                         year=2003, doi="10.1063/1.1558471")

      CALL add_reference(key=Cohen2000, &
                         authors=s2a("Cohen, Morrel H.", "Frydel, Derek", "Burke, Kieron", "Engel, Eberhard"), &
                         title="Total energy density as an interpretative tool", &
                         source="The Journal of Chemical Physics", volume="113", pages="2990", &
                         year=2000, doi="10.1063/1.1286805")

      CALL add_reference(key=Rogers2002, &
                         authors=s2a("Rogers, Christopher L.", "Rappe, Andrew M."), &
                         title="Geometric formulation of quantum stress fields", &
                         source="PHYSICAL REVIEW B", volume="65", pages="224117", &
                         year=2002, doi="10.1103/PhysRevB.65.224117")

      CALL add_reference(key=Cohen2000, &
                         authors=s2a("Cohen, Morrel H.", "Frydel, Derek", "Burke, Kieron", "Engel, Eberhard"), &
                         title="Total energy density as an interpretative tool", &
                         source="The Journal of Chemical Physics", volume="113", pages="2990", &
                         year=2000, doi="10.1063/1.1286805")

      CALL add_reference(key=Rogers2002, &
                         authors=s2a("Rogers, Christopher L.", "Rappe, Andrew M."), &
                         title="Geometric formulation of quantum stress fields", &
                         source="PHYSICAL REVIEW B", volume="65", pages="224117", &
                         year=2002, doi="10.1103/PhysRevB.65.224117")

      CALL add_reference(key=Cohen2000, &
                         authors=s2a("Cohen, Morrel H.", "Frydel, Derek", "Burke, Kieron", "Engel, Eberhard"), &
                         title="Total energy density as an interpretative tool", &
                         source="The Journal of Chemical Physics", volume="113", pages="2990", &
                         year=2000, doi="10.1063/1.1286805")

      CALL add_reference(key=Rogers2002, &
                         authors=s2a("Rogers, Christopher L.", "Rappe, Andrew M."), &
                         title="Geometric formulation of quantum stress fields", &
                         source="PHYSICAL REVIEW B", volume="65", pages="224117", &
                         year=2002, doi="10.1103/PhysRevB.65.224117")

      CALL add_reference(key=Filippetti2000, &
                         authors=s2a("Filippetti, Alessio", "Fiorentini, Vincenzo"), &
                         title="Theory and applications of the stress density", &
                         source="PHYSICAL REVIEW B", volume="61", pages="8433", &
                         year=2000, doi="10.1103/PhysRevB.61.8433")

      CALL add_reference(key=Limpanuparb2011, &
                         authors=s2a("Limpanuparb, Taweetham", "Gill, Peter M. W."), &
                         title="Resolutions of the Coulomb Operator: V. The Long-Range Ewald Operator", &
                         source="Journal of Chemical Theory and Computation", volume="7", issue="8", &
                         year=2011, doi="10.1021/ct200305n")

      CALL add_reference(key=Yin2017, &
                         authors=s2a("Yin, Wen-Jin", "Krack, Matthias", "Li, Xibo", "Chen, Li-Zhen", "Liu, Li-Min"), &
                         title="Periodic continuum solvation model integrated with "// &
                         "first-principles calculations for solid surfaces", &
                         source="Prog. Nat. Sci.", volume="27", issue="2", pages="283-288", &
                         year=2017, doi="10.1016/j.pnsc.2017.03.003")

      CALL add_reference(key=Goerigk2017, &
                         authors=s2a("Goerigk, Lars", "Hansen, Andreas", "Bauer, Christoph", &
                                     "Ehrlich, Stephan", "Najibi, Asim", "Grimme, Stefan"), &
                         title="A look at the density functional theory zoo with the advanced GMTKN55 database for "// &
                         "general main group thermochemistry, kinetics and noncovalent interactions", &
                         source="Phys. Chem. Chem. Phys.", volume="19", issue="48", &
                         year=2017, doi="10.1039/C7CP04913G")

      CALL add_reference(key=Wilhelm2016a, &
                         authors=s2a("Wilhelm, J", "Del Ben, M", "Hutter, J"), &
                         title="GW in the Gaussian and plane waves scheme with application to linear acenes", &
                         source="JOURNAL OF CHEMICAL THEORY AND COMPUTATION", volume="12", pages="3623-3635", &
                         year=2016, month=2, doi="10.1021/acs.jctc.6b00380")

      CALL add_reference(key=Wilhelm2016b, &
                         authors=s2a("Wilhelm, J", "Seewald, P", "Del Ben, M", "Hutter, J"), &
                         title="Large-Scale Cubic-Scaling Random Phase Approximation Correlation", &
                         source="JOURNAL OF CHEMICAL THEORY AND COMPUTATION", volume="12", pages="5851-5859", &
                         year=2016, month=2, doi="10.1021/acs.jctc.6b00840")

      CALL add_reference(key=Wilhelm2017, &
                         authors=s2a("Wilhelm, J", "Hutter, J"), &
                         title="Periodic GW calculations in the Gaussian and plane-waves scheme", &
                         source="PHYSICAL REVIEW B", volume="95", pages="235123", &
                         year=2017, month=2, doi="10.1103/PhysRevB.95.235123")

      CALL add_reference(key=Wilhelm2018, &
                         authors=s2a("Wilhelm, J", "Golze, D", "Talirz, L", "Hutter, J", "Pignedoli, CA"), &
                         title="Toward GW calculations on thousands of atoms", &
                         source="JOURNAL OF PHYSICAL CHEMISTRY LETTERS", volume="9", pages="306-312", &
                         year=2018, month=2, doi="10.1021/acs.jpclett.7b02740")

      CALL add_reference(key=Wilhelm2021, &
                         authors=s2a("Wilhelm, J", "Seewald, P", "Golze, D"), &
                         title="Low-Scaling GW with Benchmark Accuracy and Application to Phosphorene Nanosheets", &
                         source="JOURNAL OF CHEMICAL THEORY AND COMPUTATION", volume="9", pages="1662-1677", &
                         year=2021, month=2, doi="10.1021/acs.jctc.0c01282")

      CALL add_reference(key=Lass2018, &
                         authors=s2a("Lass, M", "Mohr, S", "Wiebeler, H", "Kuehne, TD", "Plessl, C"), &
                         title="A Massively Parallel Algorithm for the Approximate Calculation of "// &
                         "Inverse P-Th Roots of Large Sparse Matrices", &
                         source="Proceedings of the Platform for Advanced Scientific Computing (PASC) Conference", &
                         year=2018, doi="10.1145/3218176.3218231")

      CALL add_reference(key=cp2kqs2020, &
                         authors=s2a("Kuehne,Thomas D.", "Iannuzzi,Marcella", "Del Ben,Mauro", "Rybkin,Vladimir V.", &
                                     "Seewald,Patrick", "Stein,Frederick", "Laino,Teodoro", "Khaliullin,Rustam Z.", &
                                     "Schuett,Ole", "Schiffmann,Florian", "Golze,Dorothea", "Wilhelm,Jan", &
                                     "Chulkov,Sergey", "Bani-Hashemian,Mohammad Hossein", "Weber,Valery", &
                                     "Borstnik,Urban", "Taillefumier,Mathieu", "Jakobovits,Alice Shoshana", &
                                     "Lazzaro,Alfio", "Pabst,Hans", "Mueller,Tiziano", "Schade,Robert", "Guidon,Manuel", &
                                     "Andermatt,Samuel", "Holmberg,Nico", "Schenter,Gregory K.", "Hehn,Anna", &
                                     "Bussy,Augustin", "Belleflamme,Fabian", "Tabacchi,Gloria", "Gloess,Andreas", &
                                     "Lass,Michael", "Bethune,Iain", "Mundy,Christopher J.", "Plessl,Christian", &
                                     "Watkins,Matt", "VandeVondele,Joost", "Krack,Matthias", "Hutter,Juerg"), &
                         title="CP2K: An electronic structure and molecular dynamics software package - Quickstep: "// &
                         "Efficient and accurate electronic structure calculations", &
                         source="The Journal of Chemical Physics", volume="152", issue="19", &
                         year=2020, doi="10.1063/5.0007045")

      CALL add_reference(key=Rycroft2009, &
                         authors=s2a("Rycroft, Chris H"), &
                         title="VORO++: A three-dimensional Voronoi cell library in C++", &
                         source="Chaos: An Interdisciplinary Journal of Nonlinear Science", volume="19", issue="4", &
                         year=2009, doi="10.1063/1.3215722")

      CALL add_reference(key=Thomas2015, &
                         authors=s2a("Thomas, Martin", "Brehm, Martin", "Kirchner, Barbara"), &
                         title="Voronoi dipole moments for the simulation of bulk phase vibrational spectra", &
                         source="Physical Chemistry Chemical Physics", volume="17", issue="5", &
                         year=2015, doi="10.1039/C4CP05272B")

      CALL add_reference(key=Brehm2018, &
                         authors=s2a("Brehm, Martin", "Thomas, Martin"), &
                         title="An Efficient Lossless Compression Algorithm for Trajectories of "// &
                         "Atom Positions and Volumetric Data", &
                         source="Journal of Chemical Information and Modeling", volume="58", issue="10", &
                         year=2018, doi="10.1021/acs.jcim.8b00501")

      CALL add_reference(key=Brehm2020, &
                         authors=s2a("Brehm, Martin", "Thomas, Martin", "Gehrke, Sascha", "Kirchner, Barbara"), &
                         title="TRAVIS - A free analyzer for trajectories from molecular simulation", &
                         source="The Journal of Chemical Physics", volume="152", issue="16", &
                         year=2020, doi="10.1063/5.0005078")

      CALL add_reference(key=Shigeta2001, &
                         authors=s2a("Y. Shigeta", "A. M. Ferreira", "V. G. Zakrzewski", "J. V. Ortiz"), &
                         title="Electron propagator calculations with Kohn–Sham reference states", &
                         source="International Journal of Quantum Chemistry", volume="85", issue="4-5", &
                         year=2001, doi="10.1002/qua.1543")

      CALL add_reference(key=Bussy2021a, &
                         authors=s2a("Bussy, Augustin", "Hutter, Juerg"), &
                         title="Efficient and low-scaling linear-response time-dependent density "// &
                         "functional theory implementation for core-level spectroscopy of large and periodic systems", &
                         source="Phys. Chem. Chem. Phys.", volume="23", issue="8", pages="4736-4746", &
                         year=2021, doi="10.1039/D0CP06164F")

      CALL add_reference(key=Bussy2021b, &
                         authors=s2a("Bussy, Augustin", "Hutter, Juerg"), &
                         title="First-principles correction scheme for linear-response time-dependent density "// &
                         "functional theory calculations of core electronic states", &
                         source="J. Chem. Phys.", volume="155", issue="3", &
                         year=2021, doi="10.1063/5.0058124")

      CALL add_reference(key=Bussy2023, &
                         authors=s2a("Bussy, Augustin", "Schuett, Ole", "Hutter, Juerg"), &
                         title="Sparse tensor based nuclear gradients for periodic Hartree-Fock and "// &
                         "low-scaling correlated wave function methods in the CP2K software package: "// &
                         "A massively parallel and GPU accelerated implementation.", &
                         source="J. Chem. Phys.", volume="158", issue="16", &
                         year=2023, doi="10.1063/5.0144493")

      CALL add_reference(key=Bussy2024, &
                         authors=s2a("Bussy, Augustin", "Hutter, Juerg"), &
                         title="Efficient periodic resolution-of-the-identity Hartree–Fock exchange "// &
                         "method with k-point sampling and Gaussian basis sets", &
                         source="J. Chem. Phys.", volume="160", issue="6", &
                         year=2024, doi="10.1063/5.0189659")

      CALL add_reference(key=Heinecke2016, &
                         authors=s2a("Heinecke, A", "Henry, G", "Hutchinson, M", "Pabst, H"), &
                         title="LIBXSMM: Accelerating Small Matrix Multiplications by Runtime Code Generation", &
                         source="Proceedings of Intl. Supercomputing Conference", pages="981-991", &
                         year=2016, doi="10.1109/SC.2016.83")

      CALL add_reference(key=Brehm2021, &
                         authors=s2a("Brehm, Martin", "Thomas, Martin"), &
                         title="Optimized Atomic Partial Charges and Radii Defined by "// &
                         "Radical Voronoi Tessellation of Bulk Phase Simulations", &
                         source="Molecules", volume="26", issue="7", &
                         year=2021, doi="10.3390/molecules26071875")

      CALL add_reference(key=Ditler2021, &
                         authors=s2a("Ditler, Edward", "Kumar, Chandan", "Luber, Sandra"), &
                         title="Analytic calculation and analysis of atomic polar tensors for molecules "// &
                         "and materials using the Gaussian and plane waves approach", &
                         source="The Journal of Chemical Physics", volume="154", pages="104121", &
                         year=2021, doi="10.1063/5.0041056")

      CALL add_reference(key=Ditler2022, &
                         authors=s2a("Ditler, Edward", "Zimmermann, Tomas", "Kumar, Chandan", "Luber, Sandra"), &
                         title="Implementation of Nuclear Velocity Perturbation and Magnetic Field "// &
                         "Perturbation Theory in CP2K and Their Application to Vibrational Circular Dichroism", &
                         source="The Journal of Chemical Theory and Computation", volume="18", pages="2448-2461", &
                         year=2022, doi="10.1021/acs.jctc.2c00006")

      CALL add_reference(key=Mattiat2019, &
                         authors=s2a("Mattiat, Johann", "Luber, Sandra"), &
                         title="Vibrational (resonance) Raman optical activity with "// &
                         "real time time dependent density functional theory", &
                         source="The Journal of Chemical Physics", volume="151", pages="234110", &
                         year=2019, doi="10.1063/1.5132294")

      CALL add_reference(key=Mattiat2022, &
                         authors=s2a("Mattiat, Johann", "Luber, Sandra"), &
                         title="Comparison of Length, Velocity, and Symmetric Gauges for the Calculation of "// &
                         "Absorption and Electric Circular Dichroism Spectra with "// &
                         "Real-Time Time-Dependent Density Functional Theory", &
                         source="The Journal of Chemical Theory and Computation", volume="18", pages="5513-5526", &
                         year=2022, doi="10.1021/acs.jctc.2c00644")

      CALL add_reference(key=Belleflamme2023, &
                         authors=s2a("Belleflamme, Fabian", "Hehn, Anna", "Iannuzzi, Marcella", "Hutter, Juerg"), &
                         title="A variational formulation of the Harris functional as a correction to "// &
                         "approximate Kohn–Sham density functional theory", &
                         source="The Journal of Chemical Physics", volume="158", pages="054111", &
                         year=2023, doi="10.1063/5.0122671")

      CALL add_reference(key=Knizia2013, &
                         authors=s2a("Knizia Gerald"), &
                         title="Intrinsic Atomic Orbitals: An Unbiased Bridge between Quantum Theory and Chemical Concepts", &
                         source="Journal of Chemical Theory and Computation", volume="9", pages="4834-4843", &
                         year=2013, doi="10.1021/ct400687b")

      CALL add_reference(key=Musaelian2023, &
                         authors=s2a("Musaelian, A", "Batzner, S", "Johansson, A", &
                                     "Sun, L", "Owen, CJ", "Kornbluth, M", "Kozinsky, B"), &
                         title="Learning local equivariant representations for large-scale atomistic dynamics", &
                         source="Nature Communications", volume="14", pages="579", &
                         year=2023, doi="10.1038/s41467-023-36329-y")

      CALL add_reference(key=Eriksen2020, &
                         authors=s2a("Eriksen, JJ"), &
                         title="Mean-Field density matrix decompositions", &
                         source="The Journal of Chemical Physics", volume="153", pages="214109", &
                         year=2020, doi="10.1063/5.0030764")

      CALL add_reference(key=Graml2024, &
                         authors=s2a("Graml, Maximilian", "Zollner, Klaus", &
                                     "Hernangómez-Pérez, Daniel", "Faria Junior, Paulo Eduardo", "Wilhelm, Jan"), &
                         title="Low-scaling GW algorithm applied to twisted transition-metal dichalcogenide heterobilayers", &
                         source="Journal of Chemical Theory and Computation", volume="20", pages="2202-2208", &
                         year=2024, doi="10.1021/acs.jctc.3c01230")

      CALL add_reference(key=Wang2018, &
                         authors=s2a("Wang, Han", "Zhang, Linfeng", "Han, Jiequn", "E, Weinan"), &
                         title="DeePMD-kit: A Deep Learning Package for Many-body Potential Energy "// &
                         "Representation and Molecular Dynamics", &
                         source="Computer Physics Communications", volume="228", issue="7", &
                         year=2018, doi="10.1016/j.cpc.2018.03.016")

      CALL add_reference(key=Zeng2023, &
                         authors=s2a("Zeng, Jinzhe", "Zhang, Duo", "Lu, Denghui", "Mo, Pinghui", "Li, Zeyu", &
                                     "Chen, Yixiao", "Rynik, Marián", "Huang, Li'ang", "Li, Ziyao", "Shi, Shaochen", &
                                     "Wang, Yingze", "Ye, Haotian", "Tuo, Ping", "Yang, Jiabin", "Ding, Ye", &
                                     "Li, Yifan", "Tisi, Davide", "Zeng, Qiyu", "Bao, Han", "Xia, Yu", &
                                     "Huang, Jiameng", "Muraoka, Koki", "Wang, Yibo", "Chang, Junhan", "Yuan, Fengbo", &
                                     "Bore, Sigbjørn Løland", "Cai, Chun", "Lin, Yinnian", "Wang, Bo", "Xu, Jiayan", &
                                     "Zhu, Jia-Xin", "Luo, Chenxing", "Zhang, Yuzhi", "Goodall, Rhys E. A.", &
                                     "Liang, Wenshuo", "Singh, Anurag Kumar", "Yao, Sikai", "Zhang, Jingchao", &
                                     "Wentzcovitch, Renata", "Han, Jiequn", "Liu, Jie", "Jia, Weile", "York, Darrin M.", &
                                     "E, Weinan", "Car, Roberto", "Zhang, Linfeng", "Wang, Han"), &
                         title="DeePMD-kit v2: A software package for deep potential models", &
                         source="The Journal of Chemical Physics", volume="159", &
                         year=2023, doi="10.1063/5.0155600")

      CALL add_reference(key=Solca2024, &
                         authors=s2a("Solcà, Raffaele", "Simberg, Mikael", "Meli, Rocco", &
                                     "Invernizzi, Alberto", "Reverdell, Auriane", "Biddiscombe, John"), &
                         title="DLA-Future: A Task-Based Linear Algebra Library "// &
                         "Which Provides a GPU-Enabled Distributed Eigensolver", &
                         source="Asynchronous Many-Task Systems and Applications (WAMTA2024)", &
                         year=2024, doi="10.1007/978-3-031-61763-8_13")

      CALL add_reference(key=Caldeweyher2017, &
                         authors=s2a("Caldeweyher, E", "Bannwarth, C", "Grimme, S"), &
                         title="Extension of the D3 dispersion coefficient model", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="147", pages="034112", &
                         year=2017, month=7, day=21, doi="10.1063/1.4993215")

      CALL add_reference(key=Caldeweyher2019, &
                         authors=s2a("Caldeweyher, E", "Ehlert, S", "Hansen, A", "Neugebauer, H", &
                                     "Spicher, S", "Bannwarth, C", "Grimme, S"), &
                         title="A generally applicable atomic-charge dependent London dispersion correction", &
                         source="JOURNAL OF CHEMICAL PHYSICS", volume="150", pages="154122", &
                         year=2019, doi="10.1063/1.5090222")

      CALL add_reference(key=Caldeweyher2020, &
                         authors=s2a("Caldeweyher, E", "Mewes, JM", "Ehlert, S", "Grimme, S"), &
                         title="Extension and evaluation of the D4 London-dispersion model for periodic systems", &
                         source="PHYSICAL CHEMISTRY CHEMICAL PHYSICS", volume="22", pages="8499-8512", &
                         year=2020, doi="10.1039/d0cp00502a")

      CALL add_reference(key=Freeman1977, &
                         authors=s2a("Freeman, David L."), &
                         title="Coupled-cluster expansion applied to the electron gas: Inclusion of ring and exchange effects", &
                         source="Phys. Rev. B", volume="15", issue="12", &
                         year=1977, doi="10.1103/PhysRevB.15.5512")

      CALL add_reference(key=Gruneis2009, &
                         authors=s2a("Grueneis, Andreas", "Marsman, Martijn", "Harl, Judith", &
                                     "Schimka, Laurids", "Kresse, Georg"), &
                         title="Making the random phase approximation to electronic correlation accurate", &
                         source="J. Chem. Phys.", volume="131", issue="15", &
                         year=2009, doi="10.1063/1.3250347")

      CALL add_reference(key=Stein2022, &
                         authors=s2a("Stein, Frederick", "Hutter, Juerg"), &
                         title="Double-hybrid density functionals for the condensed phase: "// &
                         "Gradients, stress tensor, and auxiliary-density matrix method acceleration", &
                         source="J. Chem. Phys.", volume="156", issue="7", &
                         year=2022, doi="10.1063/5.0082327")

      CALL add_reference(key=Stein2024, &
                         authors=s2a("Stein, Frederick", "Hutter, Juerg"), &
                         title="Massively parallel implementation of gradients within the random phase approximation: "// &
                         "Application to the polymorphs of benzene", &
                         source="J. Chem. Phys.", volume="160", issue="2", &
                         year=2024, doi="10.1063/5.0180704")

      CALL add_reference(key=Blase2018, &
                         authors=s2a("Blase, X", "Duchemin, I", "Jacquemin, D"), &
                         title="The Bethe-Salpeter equation in chemistry: relations with TD-DFT, applications and challenges", &
                         source="Chemical Society reviews", volume="47", pages="1022-1043", &
                         year=2018, doi="10.1039/c7cs00049a")

      CALL add_reference(key=Blase2020, &
                         authors=s2a("Blase, X", "Duchemin, I", "Jacquemin, D", "Loos, P"), &
                         title="The Bethe–Salpeter Equation Formalism: From Physics to Chemistry", &
                         source="The Journal of Physical Chemistry Letters", volume="11", pages="7371-7382", &
                         year=2020, doi="10.1021/acs.jpclett.0c01875")

      CALL add_reference(key=Bruneval2015, &
                         authors=s2a("Bruneval, F", "Hamed, SM", "Neaton, JB"), &
                         title="A systematic benchmark of the ab initio Bethe-Salpeter equation approach for "// &
                         "low-lying optical excitations of small organic molecules", &
                         source="The Journal of chemical physics", volume="142", pages="244101", &
                         year=2015, doi="10.1063/1.4922489")

      CALL add_reference(key=Golze2019, &
                         authors=s2a("Golze, D", "Dvorak, M", "Rinke, P"), &
                         title="The GW Compendium: A Practical Guide to Theoretical Photoemission Spectroscopy", &
                         source="Frontiers in chemistry", volume="7", pages="377", &
                         year=2019, doi="10.3389/fchem.2019.00377")

      CALL add_reference(key=Gui2018, &
                         authors=s2a("Gui, X", "Holzer, C", "Klopper, W"), &
                         title="Accuracy Assessment of GW Starting Points for Calculating "// &
                         "Molecular Excitation Energies Using the Bethe-Salpeter Formalism", &
                         source="Journal of chemical theory and computation", volume="14", pages="2127-2136", &
                         year=2018, doi="10.1021/acs.jctc.8b00014")

      CALL add_reference(key=Jacquemin2017, &
                         authors=s2a("Jacquemin, D", "Duchemin, I", "Blase, X"), &
                         title="Is the Bethe-Salpeter Formalism Accurate for Excitation Energies? "// &
                         "Comparisons with TD-DFT, CASPT2, and EOM-CCSD", &
                         source="The Journal of Physical Chemistry Letters", volume="8", pages="1524-1529", &
                         year=2017, doi="10.1021/acs.jpclett.7b00381")

      CALL add_reference(key=Liu2020, &
                         authors=s2a("Liu, C", "Kloppenburg, J", "Yao, Y", "Ren, X", "Appel, H", &
                                     "Kanai, Y", "Blum, V"), &
                         title="All-electron ab initio Bethe-Salpeter equation approach to "// &
                         "neutral excitations in molecules with numeric atom-centered orbitals", &
                         source="The Journal of chemical physics", volume="152", pages="044105", &
                         year=2020, doi="10.1063/1.5123290")

      CALL add_reference(key=Sander2015, &
                         authors=s2a("Sander, T", "Maggio, E", "Kresse, G"), &
                         title="Beyond the Tamm-Dancoff approximation for extended systems using exact diagonalization", &
                         source="Physical Review B", volume="92", &
                         year=2015, doi="10.1103/PhysRevB.92.045209")

      CALL add_reference(key=Schreiber2008, &
                         authors=s2a("Schreiber, M", "Silva-Junior, MR", "Sauer, SPA", "Thiel, W"), &
                         title="Benchmarks for electronically excited states: CASPT2, CC2, CCSD, and CC3", &
                         source="The Journal of chemical physics", volume="128", pages="134110", &
                         year=2008, doi="10.1063/1.2889385")

      CALL add_reference(key=Hueser2013, &
                         authors=s2a("Hueser, F", "Olsen, T", "Thygesen, KS"), &
                         title="Quasiparticle GW calculations for solids, molecules, and two-dimensional materials", &
                         source="Physical Review B", volume="87", &
                         year=2013, doi="10.1103/PhysRevB.87.235132")

      CALL add_reference(key=vanSetten2015, &
                         authors=s2a("van Setten, MJ", "Caruso, F", "Sharifzadeh, S", "Ren, X", "Scheffler, M", &
                                     "Liu, F", "Lischner, J", "Lin, L", "Deslippe, JR", "Louie, SG", "Yang, C", &
                                     "Weigend, F", "Neaton, JB", "Evers, F", "Rinke, P"), &
                         title="GW100: Benchmarking G0W0 for Molecular Systems", &
                         source="Journal of chemical theory and computation", volume="11", pages="5665-5687", &
                         year=2015, doi="10.1021/acs.jctc.5b00453")

      CALL add_reference(key=Setyawan2010, &
                         authors=s2a("Setyawan, W", "Curtarolo, S"), &
                         title="High-throughput electronic band structure calculations: Challenges and tools", &
                         source="Comput. Mater. Sci.", volume="49", pages="299-312", &
                         year=2010, doi="10.1016/j.commatsci.2010.05.010")

   END SUBROUTINE add_all_references

END MODULE bibliography
