!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
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
!>      add references here as an easy way to cite these references from the code
!>      (e.g., whenever the relevant part of the code is executed)
!>      and to add citations to the input
!>      at the end of the run a bibliography is printed that can be used, e.g., in papers
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
                            Grimme2017, Asgeirsson2017, Bannwarth2019, Kondov2007, Clabaut2020, &
                            Clabaut2021, Ren2011, Ren2013, Cohen2000, Rogers2002, Filippetti2000, &
                            Limpanuparb2011, Martin2003, Yin2017, Goerigk2017, &
                            Wilhelm2016a, Wilhelm2016b, Wilhelm2017, Wilhelm2018, Wilhelm2021, Lass2018, &
                            cp2kqs2020, Behler2007, Behler2011, Schran2020a, Schran2020b, &
                            Rycroft2009, Thomas2015, Brehm2018, Brehm2020, Shigeta2001, Heinecke2016, &
                            Brehm2021, Bussy2021a, Bussy2021b, Ditler2021, Ditler2022, Mattiat2019, &
                            Mattiat2022, Belleflamme2023, Knizia2013, Musaelian2023, Eriksen2020, &
                            Bussy2023, Bussy2024, Wang2018, Zeng2023, Graml2024, Solca2024, &
                            Caldeweyher2017, Caldeweyher2019, Caldeweyher2020, Freeman1977, Gruneis2009, &
                            Stein2022, Stein2024, Pracht2019, &
                            Blase2018, Blase2020, Bruneval2015, Golze2019, Gui2018, Jacquemin2017, Liu2020, &
                            Sander2015, Schreiber2008, vanSetten2015, Setyawan2010, Ahart2024, Knysh2024, &
                            Schambeck2024, Mewes2018, Sertcan2024, Drautz2019, Lysogorskiy2021, Bochkarev2024, &
                            VazdaCruz2021, Chen2025

CONTAINS

! **************************************************************************************************
!> \brief adds references that can later be cited / printed using the key
!> \par History
!>      08.2007 created [Joost VandeVondele]
!> \note
!>     - please add DOI whenever available, this will result in a clickable
!>       link in the input reference manual.
!>     - please provide for journal articles:
!>       - all author names, abbreviate the first names
!>       - the title of the article
!>       - the abbreviated journal name (please use the ISO4 standard)
!>       - volume
!>       - pages (in case there is starting and end page: please provide start page-end page;
!>                in case there is article number, e.g. J. Chem. Phys, Phys. Rev., then provide
!>                the article number)
! **************************************************************************************************
   SUBROUTINE add_all_references()
      CALL add_reference(key=Ceriotti2012, &
                         authors=s2a("M. Ceriotti", "D. Manolopoulos"), &
                         title="Efficient First-Principles Calculation "// &
                         "of the Quantum Kinetic Energy and Momentum Distribution of Nuclei", &
                         source="Phys. Rev. Lett.", volume="109", pages="100604", &
                         year=2012, doi="10.1103/PhysRevLett.109.100604")

      CALL add_reference(key=Ceriotti2010, &
                         authors=s2a("M. Ceriotti", "M. Parrinello", "T. Markland", "D. Manolopoulos"), &
                         title="Efficient stochastic thermostatting of path integral molecular dynamics", &
                         source="J. Chem. Phys.", volume="133", pages="124104", &
                         year=2010, doi="10.1063/1.3489925")

      CALL add_reference(key=Wellendorff2012, &
                         authors=s2a("J. Wellendorff", "K. Lundgaard", "A. Mogelhoj", "V. Petzold", &
                                     "D. Landis", "J. Norskov", "T. Bligaard", "K. Jacobsen"), &
                         title="Density functionals for surface science: "// &
                         "Exchange-correlation model development with Bayesian error estimation", &
                         source="Phys. Rev. B", volume="85", pages="235149", &
                         year=2012, doi="10.1103/PhysRevB.85.235149")

      CALL add_reference(key=Brelaz1979, &
                         authors=s2a("D. Brelaz"), &
                         title="New methods to color the vertices of a graph", &
                         source="Commun. ACM", volume="22", pages="251-256", &
                         year=1979, doi="10.1145/359094.359101")

      CALL add_reference(key=Bengtsson1999, &
                         authors=s2a("L. Bengtsson"), &
                         title="Dipole correction for surface supercell calculations", &
                         source="Phys. Rev. B", volume="59", pages="12301-12304", &
                         year=1999, doi="10.1103/PhysRevB.59.12301")

      CALL add_reference(key=Foiles1986, &
                         authors=s2a("S. M. Foiles", "M. I. Baskes", "M. S. Daw"), &
                         title="Embedded-atom-method functions for the fcc metals Cu, Ag, Au, Ni, Pd, Pt, and their alloys", &
                         source="Phys. Rev. B", volume="33", pages="7983-7991", &
                         year=1986, doi="10.1103/PhysRevB.33.7983")

      CALL add_reference(key=QUIP_ref, &
                         authors=s2a("QUIP"), &
                         title="libAtoms/QUIP libraries from http://www.libatoms.org, "// &
                         "please cite web site and references for appropriate potential invoked", &
                         source="web site", &
                         year=2014)

      CALL add_reference(key=Batzner2022, &
                         authors=s2a("S. Batzner", "A. Musaelian", "L. Sun", "M. Geiger", "J. P. Mailoa", &
                                     "M. Kornbluth", "N. Molinari", "T. E. Smidt", "B. Kozinsky"), &
                         title="E(3)-equivariant graph neural networks for data-efficient and accurate interatomic potentials", &
                         source="Nat. Commun.", volume="13", pages="2453", &
                         year=2022, doi="10.1038/s41467-022-29939-5")

      CALL add_reference(key=VandenCic2006, &
                         authors=s2a("E. Vanden-Eijnden", "G. Ciccotti"), &
                         title="Second-order integrators for Langevin equations with holonomic constraints", &
                         source="Chem. Phys. Lett.", volume="429", pages="310-316", &
                         year=2006, doi="10.1016/j.cplett.2006.07.086")

      CALL add_reference(key=Hu2007, &
                         authors=s2a("H. Hu", "Z. Y. Lu", "M. Elstner", "J. Hermans", "W. T. Yang"), &
                         title="Simulating water with the self-consistent-charge "// &
                         "density functional tight binding method: From molecular clusters to the liquid state", &
                         source="J. Phys. Chem. A", volume="111", pages="5685-5691", &
                         year=2007, doi="10.1021/jp070308d")

      CALL add_reference(key=Zhao1994, &
                         authors=s2a("Q. S. Zhao", "R. C. Morrison", "R. G. Parr"), &
                         title="From electron densities to Kohn-Sham kinetic energies, "// &
                         "orbital energies, exchange-correlation potentials, and exchange-correlation energies", &
                         source="Phys. Rev. A", volume="50", pages="2138-2142", &
                         year=1994, doi="10.1103/PhysRevA.50.2138")

      CALL add_reference(key=Tozer1996, &
                         authors=s2a("D. J. Tozer", "V. E. Ingamells", "N. C. Handy"), &
                         title="Exchange-correlation potentials", &
                         source="J. Chem. Phys.", volume="105", pages="9200-9213", &
                         year=1996, doi="10.1063/1.472753")

      CALL add_reference(key=Blochl1995, &
                         authors=s2a("P. E. Blöchl"), &
                         title="Electrostatic decoupling of periodic images of plane‐wave‐expanded"// &
                         " densities and derived atomic point charges", &
                         source="J. Chem. Phys.", volume="103", pages="7422-7428", &
                         year=1995, doi="10.1063/1.470314")

      CALL add_reference(key=Laino2008, &
                         authors=s2a("T. Laino", "J. Hutter"), &
                         title='Notes on "Ewald summation of electrostatic multipole interactions '// &
                         'up to quadrupolar level [J. Chem. Phys. 119, 7471 (2003)]', &
                         source="J. Chem. Phys.", volume="129", pages="074102", &
                         year=2008, doi="10.1063/1.2970887")

      CALL add_reference(key=E2002, &
                         authors=s2a("W. N. E", "W. Q. Ren", "E. Vanden-Eijnden"), &
                         title="String method for the study of rare events", &
                         source="Phys. Rev. B", volume="66", pages="052301", &
                         year=2002, doi="10.1103/PhysRevB.66.052301")

      CALL add_reference(key=Wales2004, &
                         authors=s2a("S. A. Trygubenko", "D. J. Wales"), &
                         title="A doubly nudged elastic band method for finding transition states", &
                         source="J. Chem. Phys.", volume="120", pages="2082-2094", &
                         year=2004, doi="10.1063/1.1636455")

      CALL add_reference(key=Jonsson2000_2, &
                         authors=s2a("G. Henkelman", "B. P. Uberuaga", "H. Jonsson"), &
                         title="A climbing image nudged elastic band method for finding "// &
                         "saddle points and minimum energy paths", &
                         source="J. Chem. Phys.", volume="113", pages="9901-9904", &
                         year=2000, doi="10.1063/1.1329672")

      CALL add_reference(key=Jonsson2000_1, &
                         authors=s2a("G. Henkelman", "H. Jonsson"), &
                         title="Improved tangent estimate in the nudged elastic band method for "// &
                         "finding minimum energy paths and saddle points", &
                         source="J. Chem. Phys.", volume="113", pages="9978-9985", &
                         year=2000, doi="10.1063/1.1323224")

      CALL add_reference(key=Jonsson1998, &
                         authors=s2a("H. Jonsson", "G. Mills", "K. W. Jacobsen"), &
                         title="Nudged elastic band method for finding minimum energy paths of transitions", &
                         source="Classical and Quantum Dynamics in Condensed Phase Simulations", pages="385-404", &
                         year=1998)

      CALL add_reference(key=Elber1987, &
                         authors=s2a("R. Elber", "M. Karplus"), &
                         title="A method for determining reaction paths in large molecules: Application to myoglobin", &
                         source="Chem. Phys. Lett.", volume="139", pages="375-380", &
                         year=1987, doi="10.1016/0009-2614(87)80576-6")

      CALL add_reference(key=Weber2008, &
                         authors=s2a("V. Weber", "J. VandeVondele", "J. Hutter", "A. M. N. Niklasson"), &
                         title="Direct energy functional minimization under orthogonality constraints", &
                         source="J. Chem. Phys.", volume="128", pages="084113", &
                         year=2008, doi="10.1063/1.2841077")

      CALL add_reference(key=Stewart2007, &
                         authors=s2a("J. J. P. Stewart"), &
                         title="Optimization of parameters for semiempirical methods V: "// &
                         "Modification of NDDO approximations and application to 70 elements", &
                         source="J. Mol. Model.", volume="13", pages="1173-1213", &
                         year=2007, doi="10.1007/s00894-007-0233-4")

      CALL add_reference(key=Repasky2002, &
                         authors=s2a("M. P. Repasky", "J. Chandrasekhar", "W. L. Jorgensen"), &
                         title="PDDG/PM3 and PDDG/MNDO: Improved semiempirical methods", &
                         source="J. Comput. Chem.", volume="23", pages="1601-1622", &
                         year=2002, doi="10.1002/jcc.10162")

      CALL add_reference(key=Thiel1992, &
                         authors=s2a("W. Thiel", "A. A. Voityuk"), &
                         title="Extension of the MNDO formalism tod orbitals: Integral approximations "// &
                         "and preliminary numerical results", &
                         source="Theor. Chim. Acta", volume="81", pages="391-404", &
                         year=1992, doi="10.1007/BF01134863")

      CALL add_reference(key=Stewart1989, &
                         authors=s2a("J. J. P. Stewart"), &
                         title="Optimization of parameters for semiempirical methods I. Method", &
                         source="J. Comput. Chem.", volume="10", pages="209-220", &
                         year=1989, doi="10.1002/jcc.540100208")

      CALL add_reference(key=Rocha2006, &
                         authors=s2a("G. B. Rocha", "R. O. Freire", "A. M. Simas", "J. J. P. Stewart"), &
                         title="RMI: A reparameterization of AM1 for H, C, N, O, P, S, F, Cl, Br, and I", &
                         source="J. Comput. Chem.", volume="27", pages="1101-1111", &
                         year=2006, doi="10.1002/jcc.20425")

      CALL add_reference(key=Dewar1985, &
                         authors=s2a("M. J. S. Dewar", "E. G. Zoebisch", "E. F. Healy", "J. J. P. Stewart"), &
                         title="Development and use of quantum mechanical molecular models. 76. AM1: a new "// &
                         "general purpose quantum mechanical molecular model", &
                         source="J. Am. Chem. Soc.", volume="107", pages="3902-3909", &
                         year=1985, doi="10.1021/ja00299a024")

      CALL add_reference(key=Dewar1977, &
                         authors=s2a("M. J. S. Dewar", "W. Thiel"), &
                         title="Ground states of molecules. 38. The MNDO method. Approximations and parameters", &
                         source="J. Am. Chem. Soc.", volume="99", pages="4899-4907", &
                         year=1977, doi="10.1021/ja00457a004")

      CALL add_reference(key=Henkelman1999, &
                         authors=s2a("G. Henkelman", "H. Jonsson"), &
                         title="A dimer method for finding saddle points on high dimensional "// &
                         "potential surfaces using only first derivatives", &
                         source="J. Chem. Phys.", volume="111", pages="7010-7022", &
                         year=1999, doi="10.1063/1.480097")

      CALL add_reference(key=Henkelman2014, &
                         authors=s2a("P. Xiao", "Q. Wu", "G. Henkelman"), &
                         title="Basin constrained k-dimer method for saddle point finding", &
                         source="J. Chem. Phys.", volume="141", pages="164111", &
                         year=2014, doi="10.1063/1.4898664")

      CALL add_reference(key=Aguado2003, &
                         authors=s2a("A. Aguado", "P. A. Madden"), &
                         title="Ewald summation of electrostatic multipole interactions up to the quadrupolar level", &
                         source="J. Chem. Phys.", volume="119", pages="7471-7483", &
                         year=2003, doi="10.1063/1.1605941")

      CALL add_reference(key=Yamada2000, &
                         authors=s2a("K. Yamada", "K. Kurosaki", "M. Uno", "S. Yamanaka"), &
                         title="Evaluation of thermal properties of uranium dioxide by molecular dynamics", &
                         source="J. Alloys Compd.", volume="307", pages="10-16", &
                         year=2000, doi="10.1016/S0925-8388(00)00806-9")

      CALL add_reference(key=Tosi1964a, &
                         authors=s2a("F. G. Fumi", "M. P. Tosi"), &
                         title="Ionic sizes and born repulsive parameters in the NaCl-type alkali halides—I: "// &
                         "The Huggins-Mayer and Pauling forms", &
                         source="J. Phys. Chem. Solids", volume="25", pages="31-43", &
                         year=1964, doi="10.1016/0022-3697(64)90159-3")

      CALL add_reference(key=Tosi1964b, &
                         authors=s2a("M. P. Tosi", "F. G. Fumi"), &
                         title="Ionic sizes and born repulsive parameters in the NaCl-type alkali halides—II: "// &
                         "The generalized Huggins-Mayer form", &
                         source="J. Phys. Chem. Solids", volume="25", pages="45-52", &
                         year=1964, doi="10.1016/0022-3697(64)90160-X")

      CALL add_reference(key=Tersoff1988, &
                         authors=s2a("J. Tersoff"), &
                         title="Empirical interatomic potential for silicon with improved elastic properties", &
                         source="Phys. Rev. B", volume="38", pages="9902-9905", &
                         year=1988, doi="10.1103/PhysRevB.38.9902")

      CALL add_reference(key=Siepmann1995, &
                         authors=s2a("J. I. Siepmann", "M. Sprik"), &
                         title="Influence of surface topology and electrostatic potential on water/electrode systems ", &
                         source="J. Chem. Phys.", volume="102", pages="511-524", &
                         year=1995, doi="10.1063/1.469429")

      CALL add_reference(key=Bussi2007, &
                         authors=s2a("G. Bussi", "D. Donadio", "M. Parrinello"), &
                         title="Canonical sampling through velocity rescaling", &
                         source="J. Chem. Phys.", volume="126", pages="014101", &
                         year=2007, doi="10.1063/1.2408420")

      CALL add_reference(key=Nose1984a, &
                         authors=s2a("S. Nose"), &
                         title="A unified formulation of the constant temperature molecular dynamics methods ", &
                         source="J. Chem. Phys.", volume="81", pages="511-519", &
                         year=1984, doi="10.1063/1.447334")

      CALL add_reference(key=Nose1984b, &
                         authors=s2a("S. Nose"), &
                         title="A molecular dynamics method for simulations in the canonical ensemble", &
                         source="Mol. Phys.", volume="52", pages="255-268", &
                         year=1984, doi="10.1080/00268978400101201")

      CALL add_reference(key=VandeVondele2005a, &
                         authors=s2a("J. VandeVondele", "M. Krack", "F. Mohamed", "M. Parrinello", &
                                     "T. Chassaing", "J. Hutter"), &
                         title="QUICKSTEP: Fast and accurate density functional calculations "// &
                         "using a mixed Gaussian and plane waves approach", &
                         source="Comput. Phys. Commun.", volume="167", pages="103-128", &
                         year=2005, doi="10.1016/j.cpc.2004.12.014")

      CALL add_reference(key=VandeVondele2003, &
                         authors=s2a("J. VandeVondele", "J. Hutter"), &
                         title="An efficient orbital transformation method for electronic structure calculations", &
                         source="J. Chem. Phys.", volume="118", pages="4365-4369", &
                         year=2003, doi="10.1063/1.1543154")

      CALL add_reference(key=Laino2005, &
                         authors=s2a("T. Laino", "F. Mohamed", "A. Laio", "M. Parrinello"), &
                         title="An efficient real space multigrid QM/MM electrostatic coupling", &
                         source="J. Chem. Theory Comput.", volume="1", pages="1176-1184", &
                         year=2005, doi="10.1021/ct050123f")

      CALL add_reference(key=Laino2006, &
                         authors=s2a("T. Laino", "F. Mohamed", "A. Laio", "M. Parrinello"), &
                         title="An efficient linear-scaling electrostatic coupling for treating "// &
                         "periodic boundary conditions in QM/MM simulations", &
                         source="J. Chem. Theory Comput.", volume="2", pages="1370-1378", &
                         year=2006, doi="10.1021/ct6001169")

      CALL add_reference(key=Goedecker1996, &
                         authors=s2a("S. Goedecker", "M. Teter", "J. Hutter"), &
                         title="Separable dual-space Gaussian pseudopotentials", &
                         source="Phys. Rev. B", volume="54", pages="1703-1710", &
                         year=1996, doi="10.1103/PhysRevB.54.1703")

      CALL add_reference(key=Hartwigsen1998, &
                         authors=s2a("C. Hartwigsen", "S. Goedecker", "J. Hutter"), &
                         title="Relativistic separable dual-space Gaussian pseudopotentials from H to Rn", &
                         source="Phys. Rev. B", volume="58", pages="3641-3662", &
                         year=1998, doi="10.1103/PhysRevB.58.3641")

      CALL add_reference(key=Krack2005, &
                         authors=s2a("M. Krack"), &
                         title="Pseudopotentials for H to Kr optimized for "// &
                         "gradient-corrected exchange-correlation functionals", &
                         source="Theor. Chem. Acc.", volume="114", pages="145-152", &
                         year=2005, doi="10.1007/s00214-005-0655-y")

      CALL add_reference(key=Lippert1997, &
                         authors=s2a("G. Lippert", "J. Hutter", "M. Parrinello"), &
                         title="A hybrid Gaussian and plane wave density functional scheme", &
                         source="Mol. Phys.", volume="92", pages="477-487", &
                         year=1997, doi="10.1080/002689797170220")

      CALL add_reference(key=Lippert1999, &
                         authors=s2a("G. Lippert", "J. Hutter", "M. Parrinello"), &
                         title="The Gaussian and augmented-plane-wave density functional method for "// &
                         "ab initio molecular dynamics simulations", &
                         source="Theor. Chem. Acc.", volume="103", pages="124-140", &
                         year=1999, doi="10.1007/s002140050523")

      CALL add_reference(key=Krack2002, &
                         authors=s2a("M. Krack", "A. Gambirasio", "M. Parrinello"), &
                         title="Ab initio x-ray scattering of liquid water", &
                         source="J. Chem. Phys.", volume="117", pages="9409-9412", &
                         year=2002, doi="10.1063/1.1517040")

      CALL add_reference(key=Krack2000, &
                         authors=s2a("M. Krack", "M. Parrinello"), &
                         title="All-electron ab-initio molecular dynamics", &
                         source="Phys. Chem. Chem. Phys.", volume="2", pages="2105-2112", &
                         year=2000, doi="10.1039/b001167n")

      CALL add_reference(key=Iannuzzi2007, &
                         authors=s2a("M. Iannuzzi", "J. Hutter"), &
                         title="Inner-shell spectroscopy by the Gaussian and augmented plane wave method", &
                         source="Phys. Chem. Chem. Phys.", volume="9", pages="1599-1610", &
                         year=2007, doi="10.1039/b615522g")

      CALL add_reference(key=Iannuzzi2006, &
                         authors=s2a("M. Iannuzzi", "B. Kirchner", "J. Hutter"), &
                         title="Density functional embedding for molecular systems", &
                         source="Chem. Phys. Lett.", volume="421", pages="16-20", &
                         year=2006, doi="10.1016/j.cplett.2005.08.155")

      CALL add_reference(key=Iannuzzi2005, &
                         authors=s2a("M. Iannuzzi", "T. Chassaing", "T. Wallman", "J. Hutter"), &
                         title="Ground and excited state density functional calculations with the "// &
                         "Gaussian and augmented-plane-wave method", &
                         source="CHIMIA", volume="59", pages="499-503", &
                         year=2005, doi="10.2533/000942905777676164")

      CALL add_reference(key=Toukmaji1996, &
                         authors=s2a("A. Y. Toukmaji", "J. A. Board"), &
                         title="Ewald summation techniques in perspective: A survey", &
                         source="Comput. Phys. Commun.", volume="95", pages="73-92", &
                         year=1996, doi="10.1016/0010-4655(96)00016-1")

      CALL add_reference(key=Martyna1999, &
                         authors=s2a("G. J. Martyna", "M. E. Tuckerman"), &
                         title="A reciprocal space based method for treating long range interactions in "// &
                         "ab initio and force-field-based calculations in clusters", &
                         source="J. Chem. Phys.", volume="110", pages="2810-2821", &
                         year=1999, doi="10.1063/1.477923")

      CALL add_reference(key=VandeVondele2005b, &
                         authors=s2a("J. VandeVondele", "M. Sprik"), &
                         title="A molecular dynamics study of the hydroxyl radical in solution "// &
                         "applying self-interaction-corrected density functional methods", &
                         source="Phys. Chem. Chem. Phys.", volume="7", pages="1363-1367", &
                         year=2005, doi="10.1039/b501603g")

      CALL add_reference(key=Perdew1981, &
                         authors=s2a("J. P. Perdew", "A. Zunger"), &
                         title="Self-interaction correction to density-functional approximations for many-electron systems", &
                         source="Phys. Rev. B", volume="23", pages="5048-5079", &
                         year=1981, doi="10.1103/PhysRevB.23.5048")

      CALL add_reference(key=Avezac2005, &
                         authors=s2a("M. d'Avezac", "M. Calandra", "F. Mauri"), &
                         title="Density functional theory description of hole-trapping in SiO2: "// &
                         "A self-interaction-corrected approach", &
                         source="Phys. Rev. B", volume="71", pages="205210", &
                         year=2005, doi="10.1103/PhysRevB.71.205210")

      CALL add_reference(key=Zhechkov2005, &
                         authors=s2a("L. Zhechkov", "T. Heine", "S. Patchkovskii", "G. Seifert", "H. A. Duarte"), &
                         title="An efficient a Posteriori treatment for dispersion interaction in "// &
                         "density-functional-based tight binding", &
                         source="J. Chem. Theory Comput.", volume="1", pages="841-847", &
                         year=2005, doi="10.1021/ct050065y")

      CALL add_reference(key=Elstner1998, &
                         authors=s2a("M. Elstner", "D. Porezag", "G. Jungnickel", "J. Elsner", &
                                     "M. Haugk", "T. Frauenheim", "S. Suhai", "G. Seifert"), &
                         title="Self-consistent-charge density-functional tight-binding method for "// &
                         "simulations of complex materials properties", &
                         source="Phys. Rev. B", volume="58", pages="7260-7268", &
                         year=1998, doi="10.1103/PhysRevB.58.7260")

      CALL add_reference(key=Seifert1996, &
                         authors=s2a("G. Seifert", "D. Porezag", "T. Frauenheim"), &
                         title="Calculations of molecules, clusters, and solids with a simplified LCAO-DFT-LDA scheme", &
                         source="Int. J. Quantum Chem.", volume="58", pages="185-192", &
                         year=1996, doi="10.1002/(SICI)1097-461X(1996)58:2<185::AID-QUA7>3.0.CO;2-U")

      CALL add_reference(key=Porezag1995, &
                         authors=s2a("D. Porezag", "T. Frauenheim", "T. Kohler", "G. Seifert", "R. Kaschner"), &
                         title="Construction of tight-binding-like potentials on the basis of density-functional "// &
                         "theory: Application to carbon", &
                         source="Phys. Rev. B", volume="51", pages="12947-12957", &
                         year=1995, doi="10.1103/PhysRevB.51.12947")

      CALL add_reference(key=Frigo2005, &
                         authors=s2a("M. Frigo", "S. G. Johnson"), &
                         title="The design and implementation of FFTW3", &
                         source="Proc. IEEE", volume="93", pages="216-231", &
                         year=2005, doi="10.1109/JPROC.2004.840301")

      CALL add_reference(key=Genovese2006, &
                         authors=s2a("L. Genovese", "T. Deutsch", "A. Neelov", "S. Goedecker", "G. Beylkin"), &
                         title="Efficient solution of Poisson's equation with free boundary conditions", &
                         source="J. Chem. Phys.", volume="125", pages="074105", &
                         year=2006, doi="10.1063/1.2335442")

      CALL add_reference(key=Genovese2007, &
                         authors=s2a("L. Genovese", "T. Deutsch", "S. Goedecker"), &
                         title="Efficient and accurate three-dimensional Poisson solver for surface problems", &
                         source="J. Chem. Phys.", volume="127", pages="054704", &
                         year=2007, doi="10.1063/1.2754685")

      CALL add_reference(key=Minary2003, &
                         authors=s2a("P. Minary", "G. J. Martyna", "M. E. Tuckerman"), &
                         title="Algorithms and novel applications based on the isokinetic ensemble. "// &
                         "I. Biophysical and path integral molecular dynamics", &
                         source="J. Chem. Phys.", volume="118", pages="2510-2526", &
                         year=2003, doi="10.1063/1.1534582")

      CALL add_reference(key=Evans1983, &
                         authors=s2a("D. J. Evans", "W. G. Hoover", "B. H. Failor", "B. Moran", "A. J. C. Ladd"), &
                         title="Nonequilibrium molecular dynamics via Gauss's principle of least constraint", &
                         source="Phys. Rev. A", volume="28", pages="1016-1021", &
                         year=1983, doi="10.1103/PhysRevA.28.1016")

      CALL add_reference(key=Byrd1995, &
                         authors=s2a("R. H. Byrd", "P. H. Lu", "J. Nocedal", "C. Y. Zhu"), &
                         title="A Limited Memory Algorithm for Bound Constrained Optimization", &
                         source="SIAM J. Sci. Comput.", volume="16", pages="1190-1208", &
                         year=1995, doi="10.1137/0916069")

      CALL add_reference(key=VandeVondele2007, &
                         authors=s2a("J. VandeVondele", "J. Hutter"), &
                         title="Gaussian basis sets for accurate calculations on molecular systems "// &
                         "in gas and condensed phases", &
                         source="J. Chem. Phys.", volume="127", pages="114105", &
                         year=2007, doi="10.1063/1.2770708")

      CALL add_reference(key=Ortiz1994, &
                         authors=s2a("G. Ortiz", "P. Ballone"), &
                         title="Correlation energy, structure factor, radial distribution function, "// &
                         "and momentum distribution of the spin-polarized uniform electron gas", &
                         source="Phys. Rev. B", volume="50", pages="1391-1405", &
                         year=1994, doi="10.1103/PhysRevB.50.1391")

      CALL add_reference(key=Becke1988, &
                         authors=s2a("A. D. Becke"), &
                         title="Density-functional exchange-energy approximation with correct asymptotic behavior", &
                         source="Phys. Rev. A", volume="38", pages="3098-3100", &
                         year=1988, doi="10.1103/PhysRevA.38.3098")

      CALL add_reference(key=Perdew1996, &
                         authors=s2a("J. P. Perdew", "K. Burke", "M. Ernzerhof"), &
                         title="Generalized gradient approximation made simple", &
                         source="Phys. Rev. Lett.", volume="77", pages="3865-3868", &
                         year=1996, doi="10.1103/PhysRevLett.77.3865")

      CALL add_reference(key=Zhang1998, &
                         authors=s2a("Y. K. Zhang", "W. T. Yang"), &
                         title="Comment on Generalized gradient approximation made simple", &
                         source="Phys. Rev. Lett.", volume="80", pages="890-890", &
                         year=1998, doi="10.1103/PhysRevLett.80.890")

      CALL add_reference(key=Perdew2008, &
                         authors=s2a("J. P. Perdew", "A. Ruzsinszky", "G. I. Csonka", "O. A. Vydrov", &
                                     "G. E. Scuseria", "L. A. Constantin", "X. Zhou", "K. Burke"), &
                         title="Restoring the Density-Gradient Expansion for Exchange in Solids and Surfaces", &
                         source="Phys. Rev. Lett.", volume="100", pages="136406-136409", &
                         year=2008, doi="10.1103/PhysRevLett.100.136406")

      CALL add_reference(key=Lee1988, &
                         authors=s2a("C. T. Lee", "W. T. Yang", "R. G. Parr"), &
                   title="Development of the Colle-Salvetti correlation-energy formula into a functional of the electron density", &
                         source="Phys. Rev. B", volume="37", pages="785-789", &
                         year=1988, doi="10.1103/PhysRevB.37.785")

      CALL add_reference(key=Heyd2004, &
                         authors=s2a("J. Heyd", "G. E. Scuseria"), &
                         title="Assessment and validation of a screened Coulomb hybrid density functional", &
                         source="J. Chem. Phys.", volume="120", pages="7274-7280", &
                         year=2004, doi="10.1063/1.1668634")

      CALL add_reference(key=Heyd2003, &
                         authors=s2a("J. Heyd", "G. E. Scuseria", "M. Ernzerhof"), &
                         title="Hybrid functionals based on a screened Coulomb potential", &
                         source="J. Chem. Phys.", volume="118", pages="8207-8215", &
                         year=2003, doi="10.1063/1.1564060")

      CALL add_reference(key=Heyd2006, &
                         authors=s2a("J. Heyd", "G. E. Scuseria", "M. Ernzerhof"), &
                       title="Erratum: Hybrid functionals based on a screened Coulomb potential [J. Chem. Phys. 118, 8207, 2003]", &
                         source="J. Chem. Phys.", volume="124", pages="219906", &
                         year=2006, doi="10.1063/1.2204597")

      CALL add_reference(key=Vydrov2006, &
                         authors=s2a("O. A. Vydrov", "J. Heyd", "A. V. Krukau", "G. E. Scuseria"), &
                         title="Importance of short-range versus long-range Hartree-Fock exchange "// &
                         "for the performance of hybrid density functionals", &
                         source="J. Chem. Phys.", volume="125", pages="074106", &
                         year=2006, doi="10.1063/1.2244560")

      CALL add_reference(key=Vosko1980, &
                         authors=s2a("S. H. Vosko", "L. Wilk", "M. Nusair"), &
    title="Accurate spin-dependent electron liquid correlation energies for local spin density calculations: a critical analysis", &
                         source="Can. J. Phys.", volume="58", pages="1200-1211", &
                         year=1980, doi="10.1139/p80-159")

      CALL add_reference(key=Essmann1995, &
                         authors=s2a("U. Essmann", "L. Perera", "M. L. Berkowitz", "T. Darden", &
                                     "H. Lee", "L. G. Pedersen"), &
                         title="A smooth particle mesh Ewald method", &
                         source="J. Chem. Phys.", volume="103", pages="8577-8593", &
                         year=1995, doi="10.1063/1.470117")

      CALL add_reference(key=Ewald1921, &
                         authors=s2a("P. P. Ewald"), &
                         title="Die Berechnung optischer und elektrostatischer Gitterpotentiale", &
                         source="Ann. Phys.", volume="369", pages="253-287", &
                         year=1921, doi="10.1002/andp.19213690304")

      CALL add_reference(key=Darden1993, &
                         authors=s2a("T. Darden", "D. York", "L. Pedersen"), &
                         title="Particle mesh Ewald: An N log(N) method for Ewald sums in large systems ", &
                         source="J. Chem. Phys.", volume="98", pages="10089-10092", &
                         year=1993, doi="10.1063/1.464397")

      CALL add_reference(key=Dudarev1997, &
                         authors=s2a("S. L. Dudarev", "D. N. Manh", "A. P. Sutton"), &
                         title="Effect of Mott-Hubbard correlations on the electronic structure "// &
                         "and structural stability of uranium dioxide", &
                         source="Philos. Mag. B", volume="75", pages="613-628", &
                         year=1997, doi="10.1080/13642819708202343")

      CALL add_reference(key=Dudarev1998, &
                         authors=s2a("S. L. Dudarev", "G. A. Botton", "S. Y. Savrasov", "C. J. Humphreys", &
                                     "A. P. Sutton"), &
                         title="Electron-energy-loss spectra and the structural stability of "// &
                         "nickel oxide: An LSDA+U study", &
                         source="Phys. Rev. B", volume="57", pages="1505-1509", &
                         year=1998, doi="10.1103/PhysRevB.57.1505")

      CALL add_reference(key=Hunt2003, &
                         authors=s2a("P. Hunt", "M. Sprik", "R. Vuilleumier"), &
                         title="Thermal versus electronic broadening in the density of states of liquid water", &
                         source="Chem. Phys. Lett.", volume="376", pages="68-74", &
                         year=2003, doi="10.1016/S0009-2614(03)00954-0")

      CALL add_reference(key=Guidon2008, &
                         authors=s2a("M. Guidon", "F. Schiffmann", "J. Hutter", "J. VandeVondele"), &
                         title="Ab initio molecular dynamics using hybrid density functionals", &
                         source="J. Chem. Phys.", volume="128", pages="214104", &
                         year=2008, doi="10.1063/1.2931945")

      CALL add_reference(key=Stewart1982, &
                         authors=s2a("J. J. P. Stewart", "P. Csaszar", "P. Pulay"), &
                         title="Fast semi-empirical calculations", &
                         source="J. Comput. Chem.", volume="3", pages="227-228", &
                         year=1982, doi="10.1002/jcc.540030214")

      CALL add_reference(key=Tao2003, &
                         authors=s2a("J. M. Tao", "J. P. Perdew", "V. N. Staroverov", "G. E. Scuseria"), &
                         title="Climbing the density functional ladder: Nonempirical meta-generalized "// &
                         "gradient approximation designed for molecules and solids", &
                         source="Phys. Rev. Lett.", volume="91", pages="146401", &
                         year=2003, doi="10.1103/PhysRevLett.91.146401")

      CALL add_reference(key=VandeVondele2006, &
                         authors=s2a("J. VandeVondele", "M. Iannuzzi", "J. Hutter"), &
                         title="Large scale condensed matter calculations using the Gaussian and "// &
                         "augmented plane waves method", &
                         source="Computer Simulations in Condensed Matter Systems: "// &
                         "From Materials to Chemical Biology, Vol 1", volume="703", pages="287-314", &
                         year=2006, doi="10.1007/3-540-35273-2_8")

      CALL add_reference(key=Grimme2006, &
                         authors=s2a("S. Grimme"), &
                         title="Semiempirical GGA-type density functional constructed with "// &
                         "a long-range dispersion correction", &
                         source="J. Comput. Chem.", volume="27", pages="1787-1799", &
                         year=2006, doi="10.1002/jcc.20495")

      CALL add_reference(key=Grimme2010, &
                         authors=s2a("S. Grimme", "J. Antony", "S. Ehrlich", "H. Krieg"), &
                         title="A consistent and accurate ab initio parametrization of density "// &
                         "functional dispersion correction (DFT-D) for the 94 elements H-Pu", &
                         source="J. Chem. Phys.", volume="132", pages="154104", &
                         year=2010, doi="10.1063/1.3382344")

      CALL add_reference(key=Grimme2011, &
                         authors=s2a("S. Grimme", "S. Ehrlich", "L. Goerigk"), &
                         title="Effect of the damping function in dispersion corrected density functional theory", &
                         source="J. Comput. Chem.", volume="32", pages="1456", &
                         year=2011, doi="10.1002/jcc.21759")

      CALL add_reference(key=Grimme2013, &
                         authors=s2a("S. Grimme"), &
                         title="A simplified Tamm-Dancoff density functional approach for the "// &
                         "electronic excitation spectra of very large molecules", &
                         source="J. Chem. Phys.", volume="138", pages="244104", &
                         year=2013, doi="10.1063/1.4811331")

      CALL add_reference(key=Grimme2016, &
                         authors=s2a("S. Grimme", "C. Bannwarth"), &
                         title="Ultra-fast computation of electronic spectra for large systems by "// &
                         "tight-binding based simplified Tamm-Dancoff approximation (sTDA-xTB)", &
                         source="J. Chem. Phys.", volume="145", pages="054103", &
                         year=2016, doi="10.1063/1.4959605")

      CALL add_reference(key=Grimme2017, &
                         authors=s2a("S. Grimme", "C. Bannwarth", "P. Shushkov"), &
                         title="A Robust and Accurate Tight-Binding Quantum Chemical Method for "// &
                         "Structures, Vibrational Frequencies, and Noncovalent Interactions of "// &
                         "Large Molecular Systems Parametrized for All spd-Block Elements (Z = 1-86)", &
                         source="J. Chem. Theory Comput.", volume="13", pages="1989", &
                         year=2017, doi="10.1021/acs.jctc.7b00118")

      CALL add_reference(key=Asgeirsson2017, &
                         authors=s2a("V. Asgeirsson", "C. A. Bauer", "S. Grimme"), &
                         title="Quantum chemical calculation of electron ionization mass spectra "// &
                         "for general organic and inorganic molecules", &
                         source="Chem. Sci.", volume="8", pages="4879-4895", &
                         year=2017, doi="10.1039/C7SC00601B")

      CALL add_reference(key=Bannwarth2019, &
                         authors=s2a("C. Bannwarth", "S. Ehlert", "S. Grimme"), &
                         title="GFN2-xTB--An Accurate and Broadly Parametrized Self-Consistent "// &
                         "Tight-Binding Quantum Chemical Method with Multipole Electrostatics "// &
                         "and Density-Dependent Dispersion Contributions", &
                         source="J. Chem. Theory Comput.", volume="15", pages="1652-1671", &
                         year=2019, doi="10.1021/acs.jctc.8b01176")

      CALL add_reference(key=Branduardi2007, &
                         authors=s2a("D. Branduardi", "F. L. Gervasio", "M. Parrinello"), &
                         title="From A to B in free energy space", &
                         source="J. Chem. Phys.", volume="126", pages="054103", &
                         year=2007, doi="10.1063/1.2432340")

      CALL add_reference(key=Schenter2008, &
                         authors=s2a("D. T. Chang", "G. K. Schenter", "B. C. Garrett"), &
                         title="Self-consistent polarization neglect of diatomic differential overlap: "// &
                         "Applications to water clusters", &
                         source="J. Chem. Phys.", volume="128", pages="164111", &
                         year=2008, doi="10.1063/1.2905230")

      CALL add_reference(key=Proynov2007, &
                         authors=s2a("E. Proynov", "Z. Gan", "J. Kong"), &
                         title="Analytical representation of the Becke-Roussel exchange functional", &
                         source="Chem. Phys. Lett.", volume="455", pages="103-109", &
                         year=2008, doi="10.1016/j.cplett.2008.02.039")

      CALL add_reference(key=BeckeRoussel1989, &
                         authors=s2a("A. D. Becke", "M. R. Roussel"), &
                         title="Exchange holes in inhomogeneous systems: A coordinate-space model", &
                         source="Phys. Rev. A", volume="39", pages="3761-3767", &
                         year=1989, doi="10.1103/PhysRevA.39.3761")

      CALL add_reference(key=Becke1997, &
                         authors=s2a("A. D. Becke"), &
                         title="Density-functional thermochemistry . 5. "// &
                         "Systematic optimization of exchange-correlation functionals", &
                         source="J. Chem. Phys.", volume="107", pages="8554-8560", &
                         year=1997, doi="10.1063/1.475007")

      CALL add_reference(key=Ricci2003, &
                         authors=s2a("A. Ricci", "G. Ciccotti"), &
                         title="Algorithms for Brownian dynamics", &
                         source="Mol. Phys.", volume="101", pages="1927-1931", &
                         year=2003, doi="10.1080/0026897031000108113")

      CALL add_reference(key=Kolafa2004, &
                         authors=s2a("J. Kolafa"), &
                         title="Time-reversible always stable predictor-corrector method for "// &
                         "molecular dynamics of polarizable molecules", &
                         source="J. Comput. Chem.", volume="25", pages="335-342", &
                         year=2004, doi="10.1002/jcc.10385")

      CALL add_reference(key=Kuhne2007, &
                         authors=s2a("T. D. Kühne", "M. Krack", "F. R. Mohamed", "M. Parrinello"), &
                         title="Efficient and accurate Car-Parrinello-like approach to "// &
                         "Born-Oppenheimer molecular dynamics", &
                         source="Phys. Rev. Lett.", volume="98", pages="066401", &
                         year=2007, doi="10.1103/PhysRevLett.98.066401")

      CALL add_reference(key=Rengaraj2020, &
                         authors=s2a("V. Rengaraj", "M. Lass", "C. Plessl", "T. D. Kuhne"), &
                         title="Accurate Sampling with Noisy Forces from Approximate Computing", &
                         source="Computation", volume="8", pages="39", &
                         year=2020, doi="10.3390/computation8020039")

      CALL add_reference(key=Kunert2003, &
                         authors=s2a("T. Kunert", "R. Schmidt"), &
                         title="Non-adiabatic quantum molecular dynamics: "// &
                         "General formalism and case study H-2(+) in strong laser fields", &
                         source="Eur. Phys. J. D", volume="25", pages="15-24", &
                         year=2003, doi="10.1140/epjd/e2003-00086-8")

      CALL add_reference(key=Ceriotti2009, &
                         authors=s2a("M. Ceriotti", "G. Bussi", "M. Parrinello"), &
                         title="Langevin equation with colored noise for constant-temperature "// &
                         "molecular dynamics simulations", &
                         source="Phys. Rev. Lett.", volume="102", pages="020601", &
                         year=2009, doi="10.1103/PhysRevLett.102.020601")

      CALL add_reference(key=Ceriotti2009b, &
                         authors=s2a("M. Ceriotti", "G. Bussi", "M. Parrinello"), &
                         title="Nuclear Quantum Effects in Solids Using a Colored-Noise Thermostat", &
                         source="Phys. Rev. Lett.", volume="103", pages="030603", &
                         year=2009, doi="10.1103/PhysRevLett.103.030603")

      CALL add_reference(key=Guidon2009, &
                         authors=s2a("M. Guidon", "J. Hutter", "J. VandeVondele"), &
                         title="Robust Periodic Hartree-Fock Exchange for Large-Scale Simulations "// &
                         "Using Gaussian Basis Sets", &
                         source="J. Chem. Theory Comput.", volume="5", pages="3010-3021", &
                         year=2009, doi="10.1021/ct900494g")

      CALL add_reference(key=BarducBus2008, &
                         authors=s2a("A. Barducci", "G. Bussi", "M. Parrinello"), &
                         title="Well-Tempered Metadynamics: A Smoothly Converging and Tunable Free-Energy Method", &
                         source="Phys. Rev. Lett.", volume="100", pages="020603", &
                         year=2008, doi="10.1103/PhysRevLett.100.020603")

      CALL add_reference(key=Guidon2010, &
                         authors=s2a("M. Guidon", "J. Hutter", "J. VandeVondele"), &
                         title="Auxiliary Density Matrix Methods for Hartree-Fock Exchange Calculations", &
                         source="J. Chem. Theory Comput.", volume="6", pages="2348-2364", &
                         year=2010, doi="10.1021/ct1002225")

      CALL add_reference(key=Marques2012, &
                         authors=s2a("M. A. L. Marques", "M. J. T. Oliveira", "T. Burnus"), &
                         title="LIBXC: A library of exchange and correlation functionals for density functional theory", &
                         source="Comput. Phys. Commun.", volume="183", pages="2272-2281", &
                         year=2012, doi="10.1016/j.cpc.2012.05.007")

      CALL add_reference(key=Lehtola2018, &
                         authors=s2a("S. Lehtola", "C. Steigemann", "M. J. T. Oliveira", "M. A. L. Marques"), &
                         title="Recent developments in libxc - A comprehensive library of "// &
                         "functionals for density functional theory", &
                         source="SoftwareX", volume="7", pages="1-5", &
                         year=2018, doi="10.1016/j.softx.2017.11.002")

      CALL add_reference(key=Jones2011, &
                         authors=s2a("A. Jones", "B. Leimkuhler"), &
                         title="Adaptive stochastic methods for sampling driven molecular systems", &
                         source="J. Chem. Phys.", volume="135", pages="084125", &
                         year=2011, doi="10.1063/1.3626941")

      CALL add_reference(key=Bernstein2012, &
                         authors=s2a("N. Bernstein", "C. Varnai", "I. Solt", "S. A. Winfield", &
                                     "M. C. Payne", "I. Simon", "M. Fuxreiter", "G. Csanyi"), &
                         title="QM/MM simulation of liquid water with an adaptive quantum region", &
                         source="Phys. Chem. Chem. Phys.", volume="14", pages="646-656", &
                         year=2012, doi="10.1039/c1cp22600b")

      CALL add_reference(key=Bernstein2009, &
                         authors=s2a("N. Bernstein", "J. R. Kermode", "G. Csanyi"), &
                         title="Hybrid atomistic simulation methods for materials systems", &
                         source="Rep. Prog. Phys.", volume="72", pages="026501", &
                         year=2009, doi="10.1088/0034-4885/72/2/026501")

      CALL add_reference(key=Dick1958, &
                         authors=s2a("B. G. Dick", "A. W. Overhauser"), &
                         title="Theory of the Dielectric Constants of Alkali Halide Crystals", &
                         source="Phys. Rev.", volume="112", pages="90", &
                         year=1958, doi="10.1103/PhysRev.112.90")

      CALL add_reference(key=Mitchell1993, &
                         authors=s2a("P. J. Mitchell", "D. Fincham"), &
                         title="Shell model simulations by adiabatic dynamics", &
                         source="J. Phys. Condens. Matter", volume="5", pages="1031", &
                         year=1993, doi="10.1088/0953-8984/5/8/006")

      CALL add_reference(key=Devynck2012, &
                         authors=s2a("F. Devynck", "M. Iannuzzi", "M. Krack"), &
                         title="Frenkel pair recombinations in UO2: Importance of explicit "// &
                         "description of polarizability in core-shell molecular dynamics simulations", &
                         source="Phys. Rev. B", volume="85", pages="184103", &
                         year=2012, doi="10.1103/PhysRevB.85.184103")

      CALL add_reference(key=VandeVondele2012, &
                         authors=s2a("J. VandeVondele", "U. Borstnik", "J. Hutter"), &
                         title="Linear Scaling Self-Consistent Field Calculations with "// &
                         "Millions of Atoms in the Condensed Phase", &
                         source="J. Chem. Theory Comput.", volume="8", pages="3565-3573", &
                         year=2012, doi="10.1021/ct200897x")

      CALL add_reference(key=Niklasson2003, &
                         authors=s2a("A. M. N. Niklasson", "C. J. Tymczak", "M. Challacombe"), &
                         title="Trace resetting density matrix purification in O(N) self-consistent-field theory", &
                         source="J. Chem. Phys.", volume="118", pages="8611-8620", &
                         year=2003, doi="10.1063/1.1559913")

      CALL add_reference(key=Shao2003, &
                         authors=s2a("Y. Shao", "C. Saravanan", "M. Head-Gordon", "C. A. White"), &
                         title="Curvy steps for density matrix-based energy minimization: "// &
                         "Application to large-scale self-consistent-field calculations", &
                         source="J. Chem. Phys.", volume="118", pages="6144-6151", &
                         year=2003, doi="10.1063/1.1558476")

      CALL add_reference(key=VandeVondele2002, &
                         authors=s2a("J. VandeVondele", "U. Rothlisberger"), &
                         title="Canonical adiabatic free energy sampling (CAFES): "// &
                         "A novel method for the exploration of free energy surfaces", &
                         source="J. Phys. Chem. B", volume="106", pages="203-208", &
                         year=2002, doi="10.1021/jp013346k")

      CALL add_reference(key=Dion2004, &
                         authors=s2a("M. Dion", "H. Rydberg", "E. Schroder", "D. C. Langreth", "B. I. Lundqvist"), &
                         title="Van der Waals density functional for general geometries", &
                         source="Phys. Rev. Lett.", volume="92", pages="246401", &
                         year=2004, doi="10.1103/PhysRevLett.92.246401")

      CALL add_reference(key=RomanPerez2009, &
                         authors=s2a("G. Roman-Perez", "J. M. Soler"), &
                         title="Efficient Implementation of a van der Waals Density Functional: "// &
                         "Application to Double-Wall Carbon Nanotubes", &
                         source="Phys. Rev. Lett.", volume="103", pages="096102", &
                         year=2009, doi="10.1103/PhysRevLett.103.096102")

      CALL add_reference(key=DelBen2012, &
                         authors=s2a("M. Del Ben", "J. Hutter", "J. VandeVondele"), &
                         title="Second-Order Moller-Plesset Perturbation Theory in the Condensed Phase: "// &
                         "An Efficient and Massively Parallel Gaussian and Plane Waves Approach", &
                         source="J. Chem. Theory Comput.", volume="8", pages="4177-4188", &
                         year=2012, doi="10.1021/ct300531w")

      CALL add_reference(key=Sabatini2013, &
                         authors=s2a("R. Sabatini", "T. Gorni", "S. de Gironcoli"), &
                         title="Nonlocal van der Waals density functional made simple and efficient", &
                         source="Phys. Rev. B", volume="87", pages="041108(R)", &
                         year=2013, doi="10.1103/PhysRevB.87.041108")

      CALL add_reference(key=Walewski2014, &
                         authors=s2a("L. Walewski", "H. Forbert", "D. Marx"), &
                         title="Reactive path integral quantum simulations of molecules solvated in superfluid helium", &
                         source="Comput. Phys. Commun.", volume="185", pages="884-899", &
                         year=2014, doi="10.1016/j.cpc.2013.12.011")

      CALL add_reference(key=Delben2013, &
                         authors=s2a("M. Del Ben", "J. Hutter", "J. VandeVondele"), &
                         title="Electron Correlation in the Condensed Phase from a Resolution of "// &
                         "Identity Approach Based on the Gaussian and Plane Waves Scheme", &
                         source="J. Chem. Theory Comput.", volume="9", pages="2654-2671", &
                         year=2013, doi="10.1021/ct4002202")

      CALL add_reference(key=Kikuchi2009, &
                         authors=s2a("Y. Kikuchi", "Y. Imamura", "H. Nakai"), &
                         title="One-Body Energy Decomposition Schemes Revisited: Assessment of "// &
                         "Mulliken-, Grid-, and Conventional Energy Density Analyses", &
                         source="Int. J. Quantum Chem.", volume="109", pages="2464-2473", &
                         year=2009, doi="10.1002/qua.22017")

      CALL add_reference(key=Putrino2000, &
                         authors=s2a("A. Putrino", "D. Sebastiani", "M. Parrinello"), &
                         title="Generalized Variational Density Functional Perturbation Theory", &
                         source="J. Chem. Phys.", volume="113", pages="7102-7109", &
                         year=2000, doi="10.1063/1.1312830")

      CALL add_reference(key=Tran2013, &
                         authors=s2a("F. Tran", "J. Hutter"), &
                         title="Nonlocal van der Waals functionals: The case of rare-gas dimers and solids", &
                         source="J. Chem. Phys.", volume="138", pages="204103", &
                         year=2013, doi="10.1063/1.4807332")

      CALL add_reference(key=Putrino2002, &
                         authors=s2a("A. Putrino", "M. Parrinello"), &
                         title="Anharmonic Raman Spectra in High-Pressure Ice from Ab Initio Simulations", &
                         source="Phys. Rev. Lett.", volume="88", pages="176401", &
                         year=2002, doi="10.1103/PhysRevLett.88.176401")

      CALL add_reference(key=Sebastiani2001, &
                         authors=s2a("D. Sebastiani", "M. Parrinello"), &
                         title="A New Ab Initio Approach for NMR Chemical Shifts in Periodic Systems", &
                         source="J. Phys. Chem. A", volume="105", pages="1951-1958", &
                         year=2001, doi="10.1021/jp002807j")

      CALL add_reference(key=Weber2009, &
                         authors=s2a("V. Weber", "M. Iannuzzi", "S. Giani", "J. Hutter", "R. Declerck", "M. Waroduier"), &
                         title="Magnetic Linear Response Properties Calculations with the "// &
                         "Gaussian and Augmanted-Plane-Wave Method", &
                         source="J. Chem. Phys.", volume="131", pages="014106", &
                         year=2009, doi="10.1063/1.3156803")

      CALL add_reference(key=Golze2013, &
                         authors=s2a("D. Golze", "M. Iannuzzi", "M. T. Nguyen", "D. Passerone", "J. Hutter"), &
                         title="Simulation of Adsorption Processes at Metallic Interfaces: "// &
                         "An Image Charge Augmented QM/MM Approach", &
                         source="J. Chem. Theory Comput.", volume="9", pages="5086-5097", &
                         year=2013, doi="10.1021/ct400698y")

      CALL add_reference(key=Golze2015, &
                         authors=s2a("D. Golze", "J. Hutter", "M. Iannuzzi"), &
                         title="Wetting of water on hexagonal boron nitride@Rh(111): "// &
                         "a QM/MM model based on atomic charges derived for nano-structured substrates", &
                         source="Phys. Chem. Chem. Phys.", volume="17", pages="14307-14316", &
                         year=2015, doi="10.1039/C4CP04638B")

      CALL add_reference(key=Golze2017a, &
                         authors=s2a("D. Golze", "N. Benedikter", "M. Iannuzzi", "J. Wilhelm", "J. Hutter"), &
                         title="Fast evaluation of solid harmonic Gaussian integrals for local "// &
                         "resolution-of-the-identity methods and range-separated hybrid functionals", &
                         source="J. Chem. Phys.", volume="146", pages="034105", &
                         year=2017, doi="10.1063/1.4973510")

      CALL add_reference(key=Golze2017b, &
                         authors=s2a("D. Golze", "M. Iannuzzi", "J. Hutter"), &
                         title="Local Fitting of the Kohn-Sham Density in a Gaussian and "// &
                         "Plane Waves Scheme for Large-Scale Density Functional Theory Simulations", &
                         source="J. Chem. Theory Comput.", volume="13", pages="2202-2214", &
                         year=2017, doi="10.1021/acs.jctc.7b00148")

      CALL add_reference(key=Fattebert2002, &
                         authors=s2a("J. L. Fattebert", "F. Gygi"), &
                         title="Density functional theory for efficient ab initio "// &
                         "molecular dynamics simulations in solution", &
                         source="J. Comput. Chem.", volume="23", pages="662-666", &
                         year=2002, doi="10.1002/jcc.10069")

      CALL add_reference(key=Andreussi2012, &
                         authors=s2a("O. Andreussi", "I. Dabo", "N. Marzari"), &
                         title="Revised self-consistent continuum solvation in electronic-structure calculations", &
                         source="J. Chem. Phys.", volume="136", pages="064102", &
                         year=2012, doi="10.1063/1.3676407")

      CALL add_reference(key=Tuckerman1992, &
                         authors=s2a("M. Tuckerman", "B. J. Berne", "G. J. Martyna"), &
                         title="Reversible multiple time scale molecular dynamics", &
                         source="J. Chem. Phys.", volume="97", pages="1990-2001", &
                         year=1992, doi="10.1063/1.463137")

      CALL add_reference(key=Goedecker2004, &
                         authors=s2a("S. Goedecker"), &
                         title="Minima hopping: An efficient search method for the global minimum of "// &
                         "the potential energy surface of complex molecular systems", &
                         source="J. Chem. Phys.", volume="120", pages="9911-9917", &
                         year=2004, doi="10.1063/1.1724816")

      CALL add_reference(key=Khaliullin2007, &
                         authors=s2a("R. Z. Khaliullin", "E. A. Cobar", "R. C. Lochan", &
                                     "A. T. Bell", "M. Head-Gordon"), &
                         title="Unravelling the origin of intermolecular interactions using "// &
                         "absolutely localized molecular orbitals", &
                         source="J. Phys. Chem. A", volume="111", pages="8753-8765", &
                         year=2007, doi="10.1021/jp073685z")

      CALL add_reference(key=Khaliullin2008, &
                         authors=s2a("R. Z. Khaliullin", "A. T. Bell", "M. Head-Gordon"), &
                         title="Analysis of charge transfer effects in molecular complexes "// &
                         "based on absolutely localized molecular orbitals", &
                         source="J. Chem. Phys.", volume="128", pages="184112", &
                         year=2008, doi="10.1063/1.2912041")

      CALL add_reference(key=Khaliullin2013, &
                         authors=s2a("R. Z. Khaliullin", "J. VandeVondele", "J. Hutter"), &
                         title="Efficient Linear-Scaling Density Functional Theory for Molecular Systems", &
                         source="J. Chem. Theory Comput.", volume="9", pages="4421-4427", &
                         year=2013, doi="10.1021/ct400595k")

      CALL add_reference(key=Hutter2014, &
                         authors=s2a("J. Hutter", "M. Iannuzzi", "F. Schiffmann", "J. VandeVondele"), &
                         title="CP2K: atomistic simulations of condensed matter systems", &
                         source="WIREs Comput Mol Sci.", volume="4", pages="15-25", &
                         year=2014, doi="10.1002/wcms.1159")

      CALL add_reference(key=Kantorovich2008, &
                         authors=s2a("L. Kantorovich"), &
                         title="Generalized Langevin equation for solids. I. Rigorous derivation and main properties", &
                         source="Phys. Rev. B", volume="78", pages="094304", &
                         year=2008, doi="10.1103/PhysRevB.78.094304")

      CALL add_reference(key=Kantorovich2008a, &
                         authors=s2a("L. Kantorovich", "N. Rompotis"), &
                         title="Generalized Langevin equation for solids. II. Stochastic boundary "// &
                         "conditions for nonequilibrium molecular dynamics simulations", &
                         source="Phys. Rev. B", volume="78", pages="094305", &
                         year=2008, doi="10.1103/PhysRevB.78.094305")

      CALL add_reference(key=Niklasson2014, &
                         authors=s2a("E. Rubensson", "A. Niklasson"), &
                         title="Interior Eigenvalues from Density Matrix Expansions in "// &
                         "Quantum Mechanical Molecular Dynamics", &
                         source="SIAM J. Sci. Comput.", volume="36", pages="B147-B170", &
                         year=2014, doi="10.1137/130911585")

      CALL add_reference(key=Borstnik2014, &
                         authors=s2a("U. Borstnik", "J. VandeVondele", "V. Weber", "J. Hutter"), &
                         title="Sparse matrix multiplication: The distributed block-compressed sparse row library", &
                         source="Parallel Comput.", volume="40", pages="47-58", &
                         year=2014, doi="10.1016/j.parco.2014.03.012")

      CALL add_reference(key=Rayson2009, &
                         authors=s2a("M. J. Rayson", "P. R. Briddon"), &
                         title="Highly efficient method for Kohn-Sham density functional "// &
                         "calculations of 500-10 000 atom systems", &
                         source="Phys. Rev. B", volume="80", pages="205104", &
                         year=2009, doi="10.1103/PhysRevB.80.205104")

      CALL add_reference(key=Merlot2014, &
                         authors=s2a("P. Merlot", "R. Izsak", "A. Borgo", &
                                     "T. Kjaergaard", "T. Helgaker", "S. Reine"), &
                         title="Charge-constrained auxiliary-density-matrix methods for the "// &
                         "Hartree-Fock exchange contribution", &
                         source="J. Chem. Phys.", volume="141", pages="094104", &
                         year=2014, doi="10.1063/1.4894267")

      CALL add_reference(key=Lin2009, &
                         authors=s2a("L. Lin", "J. F. Lu", "L. X. Ying", "R. Car", "W. N. E"), &
                         title="Fast algorithm for extracting the diagonal of the inverse matrix with "// &
                         "application to the electronic structure analysis of metallic systems", &
                         source="Commun. Math. Sci.", volume="7", pages="755-777", &
                         year=2009, doi="10.4310/CMS.2009.v7.n3.a12")

      CALL add_reference(key=Lin2013, &
                         authors=s2a("L. Lin", "M. Chen", "C. Yang", "L. He"), &
                         title="Accelerating atomic orbital-based electronic structure calculation "// &
                         "via pole expansion and selected inversion", &
                         source="J. Phys. Condens. Matter", volume="25", pages="295501", &
                         year=2013, doi="10.1088/0953-8984/25/29/295501")

      CALL add_reference(key=DelBen2015, &
                         authors=s2a("M. Del Ben", "O. Schuett", "T. Wentz", "P. Messmer", "J. Hutter", "J. VandeVondele"), &
                         title="Enabling simulation at the fifth rung of DFT: "// &
                         "Large scale RPA calculations with excellent time to solution", &
                         source="Comput. Phys. Commun.", volume="187", pages="120-129", &
                         year=2015, doi="10.1016/j.cpc.2014.10.021")

      CALL add_reference(key=Souza2002, &
                         authors=s2a("I. Souza", "J. Iniguez", "D. Vanderbilt"), &
                         title="First-principles approach to insulators in finite electric fields", &
                         source="Phys. Rev. Lett.", volume="89", pages="117602", &
                         year=2002, doi="10.1103/PhysRevLett.89.117602")

      CALL add_reference(key=Umari2002, &
                         authors=s2a("P. Umari", "A. Pasquarello"), &
                         title="Ab initio molecular dynamics in a finite homogeneous electric field", &
                         source="Phys. Rev. Lett.", volume="89", pages="157602", &
                         year=2002, doi="10.1103/PhysRevLett.89.157602")

      CALL add_reference(key=Stengel2009, &
                         authors=s2a("M. Stengel", "N. A. Spaldin", "D. Vanderbilt"), &
                         title="Ab initio molecular dynamics in a finite homogeneous electric field "// &
                         "Electric displacement as the fundamental variable in electronic-structure calculations", &
                         source="Nat. Phys.", volume="5", pages="304-308", &
                         year=2009, doi="10.1038/NPHYS1185")

      CALL add_reference(key=Luber2014, &
                         authors=s2a("S. Luber", "M. Iannuzzi", "J. Hutter"), &
                         title="Raman spectra from ab initio molecular dynamics and its "// &
                         "application to liquid S-methyloxirane", &
                         source="J. Chem. Phys.", volume="141", pages="094503", &
                         year=2014, doi="10.1063/1.4894425")

      CALL add_reference(key=Berghold2011, &
                         authors=s2a("G. Berghold", "M. Parrinello", "J. Hutter"), &
                         title="Polarized atomic orbitals for linear scaling methods", &
                         source="J. Chem. Phys.", volume="116", pages="1800-1810", &
                         year=2002, doi="10.1063/1.1431270")

      CALL add_reference(key=DelBen2015b, &
                         authors=s2a("M. Del Ben", "J. Hutter", "J. VandeVondele"), &
                         title="Forces and stress in second order Moller-Plesset perturbation theory "// &
                         "for condensed phase systems within the resolution-of-identity "// &
                         "Gaussian and plane waves approach", &
                         source="J. Chem. Phys.", volume="143", pages="102803", &
                         year=2015, doi="10.1063/1.4919238")

      CALL add_reference(key=Campana2009, &
                         authors=s2a("C. Campana", "B. Mussard", "T. K. Woo"), &
                         title="Electrostatic Potential Derived Atomic Charges for "// &
                         "Periodic Systems Using a Modified Error Functional", &
                         source="J. Chem. Theory Comput.", volume="5", pages="2866", &
                         year=2009, doi="10.1021/ct9003405")

      CALL add_reference(key=Schiffmann2015, &
                         authors=s2a("F. Schiffmann", "J. VandeVondele"), &
                         title="Efficient preconditioning of the electronic structure problem in "// &
                         "large scale ab initio molecular dynamics simulations", &
                         source="J. Chem. Phys.", volume="142", pages="244117", &
                         year=2015, doi="10.1063/1.4922988")

      CALL add_reference(key=Bruck2014, &
                         authors=s2a("S. Bruck", "M. Calderara", "M. H. Bani-Hashemian", "J. VandeVondele", "M. Luisier"), &
                         title="Towards ab-initio simulations of nanowire field-effect transistors", &
                         source="2014 International Workshop on Computational Electronics (IWCE)", &
                         year=2014, doi="10.1109/IWCE.2014.6865831")

      CALL add_reference(key=Rappe1992, &
                         authors=s2a("A. K. Rappe", "C. J. Casewit", "K. S. Colwell", "W. A. Goddard", "W. M. Skiff"), &
                         title="UFF, a full periodic table force field for molecular mechanics "// &
                         "and molecular dynamics simulations", &
                         source="J. Am. Chem. Soc.", volume="114", pages="10024-10035", &
                         year=1992, doi="10.1021/ja00051a040")

      CALL add_reference(key=Monkhorst1976, &
                         authors=s2a("H. J. Monkhorst", "J. D. Pack"), &
                         title="Special points for Brillouin-zone integrations", &
                         source="Phys. Rev. B", volume="13", pages="5188-5192", &
                         year=1976, doi="10.1103/PhysRevB.13.5188")

      CALL add_reference(key=MacDonald1978, &
                         authors=s2a("A. H. MacDonald"), &
                         title="Comment on special points for Brillouin-zone integrations", &
                         source="Phys. Rev. B", volume="18", pages="5897-5899", &
                         year=1978, doi="10.1103/PhysRevB.18.5897")

      CALL add_reference(key=Gilbert2008, &
                         authors=s2a("A. T. B. Gilbert", "N. A. Besley", "P. M. W. Gill"), &
                         title="Self-consistent field calculations of excited states using the "// &
                         "maximum overlap method (MOM)", &
                         source="J. Phys. Chem. A", volume="112", pages="13164-13171", &
                         year=2008, doi="10.1021/jp801738f")

      CALL add_reference(key=Barca2018, &
                         authors=s2a("G. M. J. Barca", "A. T. B. Gilbert", "P. M. W. Gill"), &
                         title="Simple models for difficult electronic excitations", &
                         source="J. Chem. Theory Comput.", volume="14", pages="1501-1509", &
                         year=2018, doi="10.1021/acs.jctc.7b00994")

      CALL add_reference(key=Schonherr2014, &
                         authors=s2a("M. Schonherr", "B. Slater", "J. Hutter", "J. VandeVondele"), &
                         title="Dielectric Properties of Water Ice, the Ice Ih/XI Phase Transition, "// &
                         "and an Assessment of Density Functional Theory", &
                         source="J. Phys. Chem. B", volume="118", pages="590-596", &
                         year=2014, doi="10.1021/jp4103355")

      CALL add_reference(key=Ceriotti2014, &
                         authors=s2a("M. Ceriotti", "J. More", "D. E. Manolopoulos"), &
                         title="i-PI: A Python interface for ab initio path integral molecular dynamics simulations", &
                         source="Comput. Phys. Commun.", volume="185", pages="1019-1026", &
                         year=2014, doi="10.1016/j.cpc.2013.10.027")

      CALL add_reference(key=BaniHashemian2016, &
                         authors=s2a("M. H. Bani-Hashemian", "S. Bruck", "M. Luisier", "J. VandeVondele"), &
                         title="A generalized Poisson solver for first-principles device simulations", &
                         source="J. Chem. Phys.", volume="144", pages="044113", &
                         year=2016, doi="10.1063/1.4940796")

      CALL add_reference(key=Kapil2016, &
                         authors=s2a("V. Kapil", "J. VandeVondele", "M. Ceriotti"), &
                         title="Accurate molecular dynamics and nuclear quantum effects at low cost by multiple steps "// &
                         "in real and imaginary time: Using density functional theory to accelerate wavefunction methods", &
                         source="J. Chem. Phys.", volume="144", pages="054111", &
                         year=2016, doi="10.1063/1.4941091")

      CALL add_reference(key=Heinzmann1976, &
                         authors=s2a("R. Heinzmann", "R. Ahlrichs"), &
                         title="Population analysis based on occupation numbers of modified atomic orbitals (MAOs)", &
                         source="Theoret. Chim. Acta", volume="42", pages="33-45", &
                         year=1976, doi="10.1007/BF00548289")

      CALL add_reference(key=Ehrhardt1985, &
                         authors=s2a("C. Ehrhardt", "R. Ahlrichs"), &
                         title="Population analysis based on occupation numbers II. Relationship between shared "// &
                         "electron numbers and bond energies and characterization of hypervalent contributions", &
                         source="Theoret. Chim. Acta", volume="68", pages="231-245", &
                         year=1985, doi="10.1007/BF00526774")

      CALL add_reference(key=Rybkin2016, &
                         authors=s2a("V. V. Rybkin", "J. VandeVondele"), &
                         title="Spin-Unrestricted Second-Order Moller-Plesset (MP2) Forces for the Condensed Phase: "// &
                         "From Molecular Radicals to F-Centers in Solids", &
                         source="J. Chem. Theory Comput.", volume="12", pages="2214-2223", &
                         year=2016, doi="10.1021/acs.jctc.6b00015")

      CALL add_reference(key=West2006, &
                         authors=s2a("D. West", "S. K. Estreicher"), &
                         title="First-Principles Calculations of Vibrational Lifetimes and Decay Channels:"// &
                         " Hydrogen-Related Modes in Si", &
                         source="Phys. Rev. Lett.", volume="96", pages="115504", &
                         year=2006, doi="10.1103/PhysRevLett.96.115504")

      CALL add_reference(key=Bates2013, &
                         authors=s2a("J. E. Bates"), &
                         title="Communication: Random phase approximation renormalized many-body perturbation theory", &
                         source="J. Chem. Phys.", volume="139", pages="171103", &
                         year=2013, doi="10.1063/1.4827254")

      CALL add_reference(key=Andermatt2016, &
                         authors=s2a("S. Andermatt", "J. Cha", "F. Schiffmann", "J. VandeVondele"), &
                         title="Combining Linear-Scaling DFT with Subsystem DFT in Born Oppenheimer "// &
                         "and Ehrenfest Molecular Dynamics Simulations: From Molecules to a Virus in Solution", &
                         source="J. Chem. Theory Comput.", volume="12", pages="3214-3227", &
                         year=2016, doi="10.1021/acs.jctc.6b00398")

      CALL add_reference(key=Zhu2016, &
                         authors=s2a("L. Zhu", "M. Amsler", "T. Fuhrer", "B. Schaefer", "S. Faraji", "S. Rostami", &
                                     "S. A. Ghasemi", "A. Sadeghi", "M. Grauzinyte", "C. Wolverton", "S. Goedecker"), &
                         title="A fingerprint based metric for measuring similarities of crystalline structures", &
                         source="J. Chem. Phys.", volume="144", pages="034203", &
                         year=2016, doi="10.1063/1.4940026")

      CALL add_reference(key=Schuett2016, &
                         authors=s2a("O. Schuett", "P. Messmer", "J. Hutter", "J. VandeVondele"), &
                         title="GPU-Accelerated Sparse Matrix-Matrix Multiplication for Linear Scaling Density Functional Theory", &
                         source="Electronic Structure Calculations on Graphics Processing Units", pages="173-190", &
                         year=2016, doi="10.1002/9781118670712.ch8")

      CALL add_reference(key=Schran2020a, &
                         authors=s2a("C. Schran", "J. Behler", "D. Marx"), &
                         title="Automated Fitting of Neural Network Potentials at Coupled Cluster Accuracy: "// &
                         "Protonated Water Clusters as Testing Ground", &
                         source="J. Chem. Theory Comput.", volume="16", pages="88-99", &
                         year=2020, doi="10.1021/acs.jctc.9b00805")

      CALL add_reference(key=Schran2020b, &
                         authors=s2a("C. Schran", "K. Brezina", "O. Marsalek"), &
                         title="Committee neural network potentials control generalization errors and enable active learning", &
                         source="J. Chem. Phys.", volume="153", pages="104105", &
                         year=2020, doi="10.1063/5.0016004")

      CALL add_reference(key=Behler2007, &
                         authors=s2a("J. Behler", "M. Parrinello"), &
                         title="Generalized neural-network representation of high-dimensional potential-energy surfaces", &
                         source="Phys. Rev. Lett.", volume="98", pages="146401", &
                         year=2007, doi="10.1103/PhysRevLett.98.146401")

      CALL add_reference(key=Behler2011, &
                         authors=s2a("J. Behler"), &
                         title="Atom-centered symmetry functions for constructing high-dimensional neural network potentials", &
                         source="J. Chem. Phys.", volume="134", pages="074106", &
                         year=2011, doi="10.1063/1.3553717")

      CALL add_reference(key=Lu2004, &
                         authors=s2a("W. C. Lu", "C. Z. Wang", "M. W. Schmidt", "L. Bytautas", "K. M. Ho", "K. Ruedenberg"), &
                         title="Molecule intrinsic minimal basis sets. I. Exact resolution of ab initio "// &
                         "optimized molecular orbitals in terms of deformed atomic minimal-basis orbitals", &
                         source="J. Chem. Phys.", volume="120", pages="2629-2637", &
                         year=2004, doi="10.1063/1.1638731")

      CALL add_reference(key=Migliore2009, &
                         authors=s2a("A. Migliore"), &
                         title="Full-electron calculation of effective electronic couplings and excitation energies of charge "// &
                         "transfer states: Application to hole transfer in DNA stacks ", &
                         source="J. Chem. Phys.", volume="131", pages="114113", &
                         year=2009, doi="10.1063/1.3232007")

      CALL add_reference(key=Mavros2015, &
                         authors=s2a("M. G. Mavros", "T. Van Voorhis"), &
                         title="Communication: CDFT-CI couplings can be unreliable when there is fractional charge transfer", &
                         source="J. Chem. Phys.", volume="143", pages="231102", &
                         year=2015, doi="10.1063/1.4938103")

      CALL add_reference(key=Becke1988b, &
                         authors=s2a("A. D. Becke"), &
                         title="A multicenter numerical integration scheme for polyatomic molecules", &
                         source="J. Chem. Phys.", volume="88", pages="2547-2553", &
                         year=1988, doi="10.1063/1.454033")

      CALL add_reference(key=Holmberg2017, &
                         authors=s2a("N. Holmberg", "K. Laasonen"), &
                         title="Efficient Constrained Density Functional Theory Implementation for Simulation of "// &
                         "Condensed Phase Electron Transfer Reactions", &
                         source="J. Chem. Theory Comput.", volume="13", pages="587-601", &
                         year=2017, doi="10.1021/acs.jctc.6b01085")

      CALL add_reference(key=Marek2014, &
                         authors=s2a("A. Marek", "V. Blum", "R. Johanni", "V. Havu", "B. Lang", &
                                     "T. Auckenthaler", "A. Heinecke", "H. Bungartz", "H. Lederer"), &
                         title="The ELPA library: scalable parallel eigenvalue solutions for "// &
                         "electronic structure theory and computational science", &
                         source="J. Phys. Condens. Matter", volume="26", pages="213201", &
                         year=2014, doi="10.1088/0953-8984/26/21/213201")

      CALL add_reference(key=VanVoorhis2015, &
                         authors=s2a("T. VanVoorhis", "M. Welborn", "J. Chen", "L. Wang"), &
                         title="Why many semiempirical molecular orbital fail for liquid water and how to fix them", &
                         source="J. Comput. Chem.", volume="36", pages="934-939", &
                         year=2015, doi="10.1002/jcc.23887")

      CALL add_reference(key=Stoychev2016, &
                         authors=s2a("G. L. Stoychev", "A. A. Auer", "F. Neese"), &
                         title="Automatic Generation of Auxiliary Basis Sets", &
                         source="J. Chem. Theory Comput.", volume="13", pages="554-562", &
                         year=2017, doi="10.1021/acs.jctc.6b01041")

      CALL add_reference(key=Kondov2007, &
                         authors=s2a("I. Kondov", "M. Cizek", "C. Benesch", "H. Wang", "M. Thoss"), &
                         title="Quantum dynamics of photoinduced electron-transfer reactions in dye-semiconductor systems: "// &
                         "First-principles description and application to coumarin 343-TiO2", &
                         source="J. Phys. Chem. C", volume="111", pages="11970-11981", &
                         year=2007, doi="10.1021/jp072217m")

      CALL add_reference(key=Futera2017, &
                         authors=s2a("Z. Futera", "J. Blumberger"), &
                         title="Electronic Couplings for Charge Transfer across Molecule/Metal and "// &
                         "Molecule/Semiconductor Interfaces: Performance of the Projector Operator-Based Diabatization Approach", &
                         source="J. Phys. Chem. C", volume="121", pages="19677-19689", &
                         year=2017, doi="10.1021/acs.jpcc.7b06566")

      CALL add_reference(key=Bailey2006, &
                         authors=s2a("A. R. Rocha", "V. M. Garcia-Suarez", "S. Bailey", "C. Lambert", "J. Ferrer", "S. Sanvito"), &
                         title="Spin and molecular electronics in atomically generated orbital landscapes", &
                         source="Phys. Rev. B", volume="73", pages="085414", &
                         year=2006, doi="10.1103/PhysRevB.73.085414")

      CALL add_reference(key=Papior2017, &
                         authors=s2a("N. Papior", "N. Lorente", "T. Frederiksen", "A. Garcia", "M. Brandbyge"), &
                         title="Improvements on non-equilibrium and transport Green function techniques", &
                         source="Comput. Phys. Commun.", volume="212", pages="8-24", &
                         year=2017, doi="10.1016/j.cpc.2016.09.022")

      CALL add_reference(key=Brieuc2016, &
                         authors=s2a("F. Brieuc", "H. Dammak", "M. Hayoun"), &
                         title="Quantum thermal Bath for Path Integral Molecular Dynamics Simulation", &
                         source="J. Chem. Theory Comput.", volume="12", pages="1351-1359", &
                         year=2016, doi="10.1021/acs.jctc.5b01146")

      CALL add_reference(key=Huang2011, &
                         authors=s2a("C. Huang", "M. Pavone", "E. A. Carter"), &
                         title="Quantum mechanical embedding theory based on a unique embedding potential", &
                         source="J. Chem. Phys.", volume="134", pages="154110", &
                         year=2011, doi="10.1063/1.3577516")

      CALL add_reference(key=Heaton_Burgess2007, &
                         authors=s2a("T. Heaton-Burgess", "F. A. Bulat", "W. T. Yang"), &
                         title="Optimized effective potentials in finite basis sets", &
                         source="Phys. Rev. Lett.", volume="98", pages="256401", &
                         year=2007, doi="10.1103/PhysRevLett.98.256401")

      CALL add_reference(key=Scheiber2018, &
                         authors=s2a("H. Scheiber", "Y. Shi", "R. Z. Khaliullin"), &
                         title="Compact orbitals enable low-cost linear-scaling ab initio "// &
                         "molecular dynamics for weakly-interacting systems", &
                         source="J. Chem. Phys.", volume="148", pages="231103", &
                         year=2018, doi="10.1063/1.5029939")

      CALL add_reference(key=Schuett2018, &
                         authors=s2a("O. Schuett", "J. VandeVondele"), &
                         title="Machine Learning Adaptive Basis Sets for Efficient Large Scale "// &
                         "Density Functional Theory Simulation", &
                         source="J. Chem. Theory Comput.", volume="14", pages="4168-4175", &
                         year=2018, doi="10.1021/acs.jctc.8b00378")

      CALL add_reference(key=Holmberg2018, &
                         authors=s2a("N. Holmberg", "K. Laasonen"), &
                         title="Diabatic model for electrochemical hydrogen evolution based on "// &
                         "constrained DFT configuration interaction", &
                         source="J. Chem. Phys.", volume="149", pages="104702", &
                         year=2018, doi="10.1063/1.5038959")

      CALL add_reference(key=Togo2018, &
                         authors=s2a("A. Togo", "I. Tanaka"), &
                         title="Spglib : a software library for crystal symmetry search", &
                         source="arXiv", pages="1808.01590", &
                         year=2018)

      CALL add_reference(key=Staub2019, &
                         authors=s2a("R. Staub", "M. Iannuzzi", "R. Z. Khaliullin", "S. N. Steinmann"), &
                         title="Energy Decomposition Analysis for Metal Surface-Adsorbate Interactions "// &
                         "by Block Localized Wave Functions", &
                         source="J. Chem. Theory Comput.", volume="15", pages="265-275", &
                         year=2019, doi="10.1021/acs.jctc.8b00957")

      CALL add_reference(key=Clabaut2020, &
                         authors=s2a("P. Clabaut", "P. Fleurat-Lessard", "C. Michel", "S. N. Steinmann"), &
                         title="Ten Facets, One Force Field: The GAL19 Force Field for Water-Noble Metal Interfaces", &
                         source="J. Chem. Theory Comput.", volume="16", pages="4565-4578", &
                         year=2020, doi="10.1021/acs.jctc.0c00091")

      CALL add_reference(key=Clabaut2021, &
                         authors=s2a("P. Clabaut"), &
                         title="Solvation and adsorptions at the solid/water interface: Developments and applications", &
                         source="Ph.D. thesis at Université de Lyon", &
                         year=2021)

      CALL add_reference(key=Richters2018, &
                         authors=s2a("D. Richters", "M. Lass", "A. Walther", "C. Plessl", "T. D. Kuehne"), &
                         title="A General Algorithm to Calculate the Inverse Principal p-th Root of "// &
                         "Symmetric Positive Definite Matrices", &
                         source="Commun. Comput. Phys.", volume="25", pages="564-585", &
                         year=2018, doi="10.4208/cicp.OA-2018-0053")

      CALL add_reference(key=Kruse2012, &
                         authors=s2a("H. Kruse", "S. Grimme"), &
                         title="A geometrical correction for the inter- and intra-molecular basis set superposition error "// &
                         "in Hartree-Fock and density functional theory calculations for large systems", &
                         source="J. Chem. Phys.", volume="136", pages="154101", &
                         year=2012, doi="10.1063/1.3700154")

      CALL add_reference(key=Ren2011, &
                         authors=s2a("X. Ren", "A. Tkatchenko", "P. Rinke", "M. Scheffler"), &
                         title="Beyond the Random-Phase Approximation for the Electron Correlation Energy: "// &
                         "The Importance of Single Excitations", &
                         source="Phys. Rev. Lett.", volume="106", pages="153003", &
                         year=2011, doi="10.1103/PhysRevLett.106.153003")

      CALL add_reference(key=Ren2013, &
                         authors=s2a("X. Ren", "P. Rinke", "G. Scuseria", "M. Scheffler"), &
                         title="Renormalized second-order perturbation theory for the electron correlation energy: "// &
                         "Concept, implementation, and benchmarks", &
                         source="Phys. Rev. B", volume="88", pages="035120", &
                         year=2013, doi="10.1103/PhysRevB.88.035120")

      CALL add_reference(key=Martin2003, &
                         authors=s2a("R. L. Martin"), &
                         title="Natural transition orbitals", &
                         source="J. Chem. Phys.", volume="118", pages="4775-4777", &
                         year=2003, doi="10.1063/1.1558471")

      CALL add_reference(key=Cohen2000, &
                         authors=s2a("M. H. Cohen", "D. Frydel", "K. Burke", "E. Engel"), &
                         title="Total energy density as an interpretative tool", &
                         source="J. Chem. Phys.", volume="113", pages="2990", &
                         year=2000, doi="10.1063/1.1286805")

      CALL add_reference(key=Rogers2002, &
                         authors=s2a("C. L. Rogers", "A. M. Rappe"), &
                         title="Geometric formulation of quantum stress fields", &
                         source="Phys. Rev. B", volume="65", pages="224117", &
                         year=2002, doi="10.1103/PhysRevB.65.224117")

      CALL add_reference(key=Filippetti2000, &
                         authors=s2a("A. Filippetti", "V. Fiorentini"), &
                         title="Theory and applications of the stress density", &
                         source="Phys. Rev. B", volume="61", pages="8433", &
                         year=2000, doi="10.1103/PhysRevB.61.8433")

      CALL add_reference(key=Limpanuparb2011, &
                         authors=s2a("T. Limpanuparb", "P. M. W. Gill"), &
                         title="Resolutions of the Coulomb Operator: V. The Long-Range Ewald Operator", &
                         source="J. Chem. Theory Comput.", volume="7", pages="2353-2357", &
                         year=2011, doi="10.1021/ct200305n")

      CALL add_reference(key=Yin2017, &
                         authors=s2a("W. Yin", "M. Krack", "X. Li", "L. Chen", "L. Liu"), &
                         title="Periodic continuum solvation model integrated with "// &
                         "first-principles calculations for solid surfaces", &
                         source="Prog. Nat. Sci.", volume="27", pages="283-288", &
                         year=2017, doi="10.1016/j.pnsc.2017.03.003")

      CALL add_reference(key=Goerigk2017, &
                         authors=s2a("L. Goerigk", "A. Hansen", "C. Bauer", &
                                     "S. Ehrlich", "A. Najibi", "S. Grimme"), &
                         title="A look at the density functional theory zoo with the advanced GMTKN55 database for "// &
                         "general main group thermochemistry, kinetics and noncovalent interactions", &
                         source="Phys. Chem. Chem. Phys.", volume="19", pages="32184-32215", &
                         year=2017, doi="10.1039/C7CP04913G")

      CALL add_reference(key=Wilhelm2016a, &
                         authors=s2a("J. Wilhelm", "M. Del Ben", "J. Hutter"), &
                         title="GW in the Gaussian and plane waves scheme with application to linear acenes", &
                         source="J. Chem. Theory Comput.", volume="12", pages="3623-3635", &
                         year=2016, doi="10.1021/acs.jctc.6b00380")

      CALL add_reference(key=Wilhelm2016b, &
                         authors=s2a("J. Wilhelm", "P. Seewald", "M. Del Ben", "J. Hutter"), &
                         title="Large-Scale Cubic-Scaling Random Phase Approximation Correlation", &
                         source="J. Chem. Theory Comput.", volume="12", pages="5851-5859", &
                         year=2016, doi="10.1021/acs.jctc.6b00840")

      CALL add_reference(key=Wilhelm2017, &
                         authors=s2a("J. Wilhelm", "J. Hutter"), &
                         title="Periodic GW calculations in the Gaussian and plane-waves scheme", &
                         source="Phys. Rev. B", volume="95", pages="235123", &
                         year=2017, doi="10.1103/PhysRevB.95.235123")

      CALL add_reference(key=Wilhelm2018, &
                         authors=s2a("J. Wilhelm", "D. Golze", "L. Talirz", "J. Hutter", "C. A. Pignedoli"), &
                         title="Toward GW calculations on thousands of atoms", &
                         source="J. Phys. Chem. Lett.", volume="9", pages="306-312", &
                         year=2018, doi="10.1021/acs.jpclett.7b02740")

      CALL add_reference(key=Wilhelm2021, &
                         authors=s2a("J. Wilhelm", "P. Seewald", "D. Golze"), &
                         title="Low-Scaling GW with Benchmark Accuracy and Application to Phosphorene Nanosheets", &
                         source="J. Chem. Theory Comput.", volume="9", pages="1662-1677", &
                         year=2021, doi="10.1021/acs.jctc.0c01282")

      CALL add_reference(key=Lass2018, &
                         authors=s2a("M. Lass", "S. Mohr", "H. Wiebeler", "T. D. Kuehne", "C. Plessl"), &
                         title="A Massively Parallel Algorithm for the Approximate Calculation of "// &
                         "Inverse P-Th Roots of Large Sparse Matrices", &
                         source="Proceedings of the Platform for Advanced Scientific Computing (PASC) Conference", &
                         year=2018, doi="10.1145/3218176.3218231")

      CALL add_reference(key=cp2kqs2020, &
                         authors=s2a("T. D. Kuehne", "M. Iannuzzi", "M. Del Ben", "V. V. Rybkin", &
                                     "P. Seewald", "F. Stein", "T. Laino", "R. Z. Khaliullin", &
                                     "O. Schuett", "F. Schiffmann", "D. Golze", "J. Wilhelm", &
                                     "S. Chulkov", "M. H. Bani-Hashemian", "V. Weber", &
                                     "U. Borstnik", "M. Taillefumier", "A. S. Jakobovits", &
                                     "A. Lazzaro", "H. Pabst", "T. Mueller", "R. Schade", "M. Guidon", &
                                     "S. Andermatt", "N. Holmberg", "G. K. Schenter", "A. Hehn", &
                                     "A. Bussy", "F. Belleflamme", "G. Tabacchi", "A. Gloess", &
                                     "M. Lass", "I. Bethune", "C. J. Mundy", "C. Plessl", &
                                     "M. Watkins", "J. VandeVondele", "M. Krack", "J. Hutter"), &
                         title="CP2K: An electronic structure and molecular dynamics software package - Quickstep: "// &
                         "Efficient and accurate electronic structure calculations", &
                         source="J. Chem. Phys.", volume="152", pages="194103", &
                         year=2020, doi="10.1063/5.0007045")

      CALL add_reference(key=Rycroft2009, &
                         authors=s2a("C. H. Rycroft"), &
                         title="VORO++: A three-dimensional Voronoi cell library in C++", &
                         source="Chaos", volume="19", pages="041111", &
                         year=2009, doi="10.1063/1.3215722")

      CALL add_reference(key=Thomas2015, &
                         authors=s2a("M. Thomas", "M. Brehm", "B. Kirchner"), &
                         title="Voronoi dipole moments for the simulation of bulk phase vibrational spectra", &
                         source="Phys. Chem. Chem. Phys.", volume="17", pages="3207-3213", &
                         year=2015, doi="10.1039/C4CP05272B")

      CALL add_reference(key=Brehm2018, &
                         authors=s2a("M. Brehm", "M. Thomas"), &
                         title="An Efficient Lossless Compression Algorithm for Trajectories of "// &
                         "Atom Positions and Volumetric Data", &
                         source="J. Chem. Inf. Model.", volume="58", pages="2092-2107", &
                         year=2018, doi="10.1021/acs.jcim.8b00501")

      CALL add_reference(key=Brehm2020, &
                         authors=s2a("M. Brehm", "M. Thomas", "S. Gehrke", "B. Kirchner"), &
                         title="TRAVIS - A free analyzer for trajectories from molecular simulation", &
                         source="J. Chem. Phys.", volume="152", pages="164105", &
                         year=2020, doi="10.1063/5.0005078")

      CALL add_reference(key=Shigeta2001, &
                         authors=s2a("Y. Shigeta", "A. M. Ferreira", "V. G. Zakrzewski", "J. V. Ortiz"), &
                         title="Electron propagator calculations with Kohn–Sham reference states", &
                         source="Int. J. Quantum Chem.", volume="85", pages="411-420", &
                         year=2001, doi="10.1002/qua.1543")

      CALL add_reference(key=Bussy2021a, &
                         authors=s2a("A. Bussy", "J. Hutter"), &
                         title="Efficient and low-scaling linear-response time-dependent density "// &
                         "functional theory implementation for core-level spectroscopy of large and periodic systems", &
                         source="Phys. Chem. Chem. Phys.", volume="23", pages="4736-4746", &
                         year=2021, doi="10.1039/D0CP06164F")

      CALL add_reference(key=Bussy2021b, &
                         authors=s2a("A. Bussy", "J. Hutter"), &
                         title="First-principles correction scheme for linear-response time-dependent density "// &
                         "functional theory calculations of core electronic states", &
                         source="J. Chem. Phys.", volume="155", pages="034108", &
                         year=2021, doi="10.1063/5.0058124")

      CALL add_reference(key=Bussy2023, &
                         authors=s2a("A. Bussy", "O. Schuett", "J. Hutter"), &
                         title="Sparse tensor based nuclear gradients for periodic Hartree-Fock and "// &
                         "low-scaling correlated wave function methods in the CP2K software package: "// &
                         "A massively parallel and GPU accelerated implementation.", &
                         source="J. Chem. Phys.", volume="158", pages="164109", &
                         year=2023, doi="10.1063/5.0144493")

      CALL add_reference(key=Bussy2024, &
                         authors=s2a("A. Bussy", "J. Hutter"), &
                         title="Efficient periodic resolution-of-the-identity Hartree–Fock exchange "// &
                         "method with k-point sampling and Gaussian basis sets", &
                         source="J. Chem. Phys.", volume="160", pages="064116", &
                         year=2024, doi="10.1063/5.0189659")

      CALL add_reference(key=Heinecke2016, &
                         authors=s2a("A. Heinecke", "G. Henry", "M. Hutchinson", "H. Pabst"), &
                         title="LIBXSMM: Accelerating Small Matrix Multiplications by Runtime Code Generation", &
                         source="Proceedings of Intl. Supercomputing Conference", pages="981-991", &
                         year=2016, doi="10.1109/SC.2016.83")

      CALL add_reference(key=Brehm2021, &
                         authors=s2a("M. Brehm", "M. Thomas"), &
                         title="Optimized Atomic Partial Charges and Radii Defined by "// &
                         "Radical Voronoi Tessellation of Bulk Phase Simulations", &
                         source="Molecules", volume="26", pages="1875", &
                         year=2021, doi="10.3390/molecules26071875")

      CALL add_reference(key=Ditler2021, &
                         authors=s2a("E. Ditler", "C. Kumar", "S. Luber"), &
                         title="Analytic calculation and analysis of atomic polar tensors for molecules "// &
                         "and materials using the Gaussian and plane waves approach", &
                         source="J. Chem. Phys.", volume="154", pages="104121", &
                         year=2021, doi="10.1063/5.0041056")

      CALL add_reference(key=Ditler2022, &
                         authors=s2a("E. Ditler", "T. Zimmermann", "C. Kumar", "S. Luber"), &
                         title="Implementation of Nuclear Velocity Perturbation and Magnetic Field "// &
                         "Perturbation Theory in CP2K and Their Application to Vibrational Circular Dichroism", &
                         source="J. Chem. Theory Comput.", volume="18", pages="2448-2461", &
                         year=2022, doi="10.1021/acs.jctc.2c00006")

      CALL add_reference(key=Mattiat2019, &
                         authors=s2a("J. Mattiat", "S. Luber"), &
                         title="Vibrational (resonance) Raman optical activity with "// &
                         "real time time dependent density functional theory", &
                         source="J. Chem. Phys.", volume="151", pages="234110", &
                         year=2019, doi="10.1063/1.5132294")

      CALL add_reference(key=Mattiat2022, &
                         authors=s2a("J. Mattiat", "S. Luber"), &
                         title="Comparison of Length, Velocity, and Symmetric Gauges for the Calculation of "// &
                         "Absorption and Electric Circular Dichroism Spectra with "// &
                         "Real-Time Time-Dependent Density Functional Theory", &
                         source="J. Chem. Theory Comput.", volume="18", pages="5513-5526", &
                         year=2022, doi="10.1021/acs.jctc.2c00644")

      CALL add_reference(key=Belleflamme2023, &
                         authors=s2a("F. Belleflamme", "A. Hehn", "M. Iannuzzi", "J. Hutter"), &
                         title="A variational formulation of the Harris functional as a correction to "// &
                         "approximate Kohn–Sham density functional theory", &
                         source="J. Chem. Phys.", volume="158", pages="054111", &
                         year=2023, doi="10.1063/5.0122671")

      CALL add_reference(key=Knizia2013, &
                         authors=s2a("G. Knizia"), &
                         title="Intrinsic Atomic Orbitals: An Unbiased Bridge between Quantum Theory and Chemical Concepts", &
                         source="J. Chem. Theory Comput.", volume="9", pages="4834-4843", &
                         year=2013, doi="10.1021/ct400687b")

      CALL add_reference(key=Musaelian2023, &
                         authors=s2a("A. Musaelian", "S. Batzner", "A. Johansson", &
                                     "L. Sun", "C. J. Owen", "M. Kornbluth", "B. Kozinsky"), &
                         title="Learning local equivariant representations for large-scale atomistic dynamics", &
                         source="Nat. Commun.", volume="14", pages="579", &
                         year=2023, doi="10.1038/s41467-023-36329-y")

      CALL add_reference(key=Eriksen2020, &
                         authors=s2a("J. J. Eriksen"), &
                         title="Mean-Field density matrix decompositions", &
                         source="J. Chem. Phys.", volume="153", pages="214109", &
                         year=2020, doi="10.1063/5.0030764")

      CALL add_reference(key=Graml2024, &
                         authors=s2a("M. Graml", "K. Zollner", &
                                     "D. Hernangómez-Pérez", "P. E. Faria Junior", "J. Wilhelm"), &
                         title="Low-scaling GW algorithm applied to twisted transition-metal dichalcogenide heterobilayers", &
                         source="J. Chem. Theory Comput.", volume="20", pages="2202-2208", &
                         year=2024, doi="10.1021/acs.jctc.3c01230")

      CALL add_reference(key=Wang2018, &
                         authors=s2a("H. Wang", "L. Zhang", "J. Han", "W. E"), &
                         title="DeePMD-kit: A Deep Learning Package for Many-body Potential Energy "// &
                         "Representation and Molecular Dynamics", &
                         source="Comput. Phys. Commun.", volume="228", pages=" 178-184", &
                         year=2018, doi="10.1016/j.cpc.2018.03.016")

      CALL add_reference(key=Zeng2023, &
                         authors=s2a("J. Zeng", "D. Zhang", "D. Lu", "P. Mo", "Z. Li", &
                                     "Y. Chen", "M. Rynik", "L. Huang", "Z. Li", "S. Shi", &
                                     "Y. Wang", "H. Ye", "P. Tuo", "J. Yang", "Y. Ding", &
                                     "Y. Li", "D. Tisi", "Q. Zeng", "H. Bao", "Y. Xia", &
                                     "J. Huang", "K. Muraoka", "Y. Wang", "J. Chang", "F. Yuan", &
                                     "S. L. Bore", "C. Cai", "Y. Lin", "B. Wang", "J. Xu", &
                                     "J. Zhu", "C. Luo", "Y. Zhang", "R. E. A. Goodall", &
                                     "W. Liang", "A. K. Singh", "S. Yao", "J. Zhang", &
                                     "R. Wentzcovitch", "J. Han", "J. Liu", "W. Jia", "D. M. York", &
                                     "W. E", "R. Car", "L. Zhang", "H. Wang"), &
                         title="DeePMD-kit v2: A software package for deep potential models", &
                         source="J. Chem. Phys.", volume="159", pages="054801", &
                         year=2023, doi="10.1063/5.0155600")

      CALL add_reference(key=Solca2024, &
                         authors=s2a("R. Solcà", "M. Simberg", "R. Meli", &
                                     "A. Invernizzi", "A. Reverdell", "J. Biddiscombe"), &
                         title="DLA-Future: A Task-Based Linear Algebra Library "// &
                         "Which Provides a GPU-Enabled Distributed Eigensolver", &
                         source="Asynchronous Many-Task Systems and Applications (WAMTA2024)", &
                         year=2024, doi="10.1007/978-3-031-61763-8_13")

      CALL add_reference(key=Caldeweyher2017, &
                         authors=s2a("E. Caldeweyher", "C. Bannwarth", "S. Grimme"), &
                         title="Extension of the D3 dispersion coefficient model", &
                         source="J. Chem. Phys.", volume="147", pages="034112", &
                         year=2017, doi="10.1063/1.4993215")

      CALL add_reference(key=Caldeweyher2019, &
                         authors=s2a("E. Caldeweyher", "S. Ehlert", "A. Hansen", "H. Neugebauer", &
                                     "S. Spicher", "C. Bannwarth", "S. Grimme"), &
                         title="A generally applicable atomic-charge dependent London dispersion correction", &
                         source="J. Chem. Phys.", volume="150", pages="154122", &
                         year=2019, doi="10.1063/1.5090222")

      CALL add_reference(key=Caldeweyher2020, &
                         authors=s2a("E. Caldeweyher", "J. M. Mewes", "S. Ehlert", "S. Grimme"), &
                         title="Extension and evaluation of the D4 London-dispersion model for periodic systems", &
                         source="Phys. Chem. Chem. Phys.", volume="22", pages="8499-8512", &
                         year=2020, doi="10.1039/d0cp00502a")

      CALL add_reference(key=Freeman1977, &
                         authors=s2a("D. L. Freeman"), &
                         title="Coupled-cluster expansion applied to the electron gas: Inclusion of ring and exchange effects", &
                         source="Phys. Rev. B", volume="15", pages="5512", &
                         year=1977, doi="10.1103/PhysRevB.15.5512")

      CALL add_reference(key=Gruneis2009, &
                         authors=s2a("A. Grueneis", "M. Marsman", "J. Harl", &
                                     "L. Schimka", "G. Kresse"), &
                         title="Making the random phase approximation to electronic correlation accurate", &
                         source="J. Chem. Phys.", volume="131", pages="154115", &
                         year=2009, doi="10.1063/1.3250347")

      CALL add_reference(key=Stein2022, &
                         authors=s2a("F. Stein", "J. Hutter"), &
                         title="Double-hybrid density functionals for the condensed phase: "// &
                         "Gradients, stress tensor, and auxiliary-density matrix method acceleration", &
                         source="J. Chem. Phys.", volume="156", pages="024120", &
                         year=2022, doi="10.1063/5.0082327")

      CALL add_reference(key=Stein2024, &
                         authors=s2a("F. Stein", "J. Hutter"), &
                         title="Massively parallel implementation of gradients within the random phase approximation: "// &
                         "Application to the polymorphs of benzene", &
                         source="J. Chem. Phys.", volume="160", pages="024120", &
                         year=2024, doi="10.1063/5.0180704")

      CALL add_reference(key=Blase2018, &
                         authors=s2a("X. Blase", "I. Duchemin", "D. Jacquemin"), &
                         title="The Bethe-Salpeter equation in chemistry: relations with TD-DFT, applications and challenges", &
                         source="Chem. Soc. Rev.", volume="47", pages="1022-1043", &
                         year=2018, doi="10.1039/c7cs00049a")

      CALL add_reference(key=Blase2020, &
                         authors=s2a("X. Blase", "I. Duchemin", "D. Jacquemin", "P. Loos"), &
                         title="The Bethe–Salpeter Equation Formalism: From Physics to Chemistry", &
                         source="J. Phys. Chem. Lett", volume="11", pages="7371-7382", &
                         year=2020, doi="10.1021/acs.jpclett.0c01875")

      CALL add_reference(key=Bruneval2015, &
                         authors=s2a("F. Bruneval", "S. M. Hamed", "J. B. Neaton"), &
                         title="A systematic benchmark of the ab initio Bethe-Salpeter equation approach for "// &
                         "low-lying optical excitations of small organic molecules", &
                         source="J. Chem. Phys.", volume="142", pages="244101", &
                         year=2015, doi="10.1063/1.4922489")

      CALL add_reference(key=Golze2019, &
                         authors=s2a("D. Golze", "M. Dvorak", "P. Rinke"), &
                         title="The GW Compendium: A Practical Guide to Theoretical Photoemission Spectroscopy", &
                         source="Front. Chem.", volume="7", pages="377", &
                         year=2019, doi="10.3389/fchem.2019.00377")

      CALL add_reference(key=Gui2018, &
                         authors=s2a("X. Gui", "C. Holzer", "W. Klopper"), &
                         title="Accuracy Assessment of GW Starting Points for Calculating "// &
                         "Molecular Excitation Energies Using the Bethe-Salpeter Formalism", &
                         source="J. Chem. Theory Comput.", volume="14", pages="2127-2136", &
                         year=2018, doi="10.1021/acs.jctc.8b00014")

      CALL add_reference(key=Jacquemin2017, &
                         authors=s2a("D. Jacquemin", "I. Duchemin", "X. Blase"), &
                         title="Is the Bethe-Salpeter Formalism Accurate for Excitation Energies? "// &
                         "Comparisons with TD-DFT, CASPT2, and EOM-CCSD", &
                         source="J. Phys. Chem. Lett.", volume="8", pages="1524-1529", &
                         year=2017, doi="10.1021/acs.jpclett.7b00381")

      CALL add_reference(key=Liu2020, &
                         authors=s2a("C. Liu", "J. Kloppenburg", "Y. Yao", "X. Ren", "H. Appel", &
                                     "Y. Kanai", "V. Blum"), &
                         title="All-electron ab initio Bethe-Salpeter equation approach to "// &
                         "neutral excitations in molecules with numeric atom-centered orbitals", &
                         source="J. Chem. Phys.", volume="152", pages="044105", &
                         year=2020, doi="10.1063/1.5123290")

      CALL add_reference(key=Sander2015, &
                         authors=s2a("T. Sander", "E. Maggio", "G. Kresse"), &
                         title="Beyond the Tamm-Dancoff approximation for extended systems using exact diagonalization", &
                         source="Phys. Rev. B", volume="92", pages="045209", &
                         year=2015, doi="10.1103/PhysRevB.92.045209")

      CALL add_reference(key=Schreiber2008, &
                         authors=s2a("M. Schreiber", "M. R. Silva-Junior", "S. P. A. Sauer", "W. Thiel"), &
                         title="Benchmarks for electronically excited states: CASPT2, CC2, CCSD, and CC3", &
                         source="J. Chem. Phys.", volume="128", pages="134110", &
                         year=2008, doi="10.1063/1.2889385")

      CALL add_reference(key=vanSetten2015, &
                         authors=s2a("M. J. van Setten", "F. Caruso", "S. Sharifzadeh", "X. Ren", "M. Scheffler", &
                                     "F. Liu", "J. Lischner", "L. Lin", "J. R. Deslippe", "S. G. Louie", "C. Yang", &
                                     "F. Weigend", "J. B. Neaton", "F. Evers", "P. Rinke"), &
                         title="GW100: Benchmarking G0W0 for Molecular Systems", &
                         source="J. Chem. Theory Comput.", volume="11", pages="5665-5687", &
                         year=2015, doi="10.1021/acs.jctc.5b00453")

      CALL add_reference(key=Setyawan2010, &
                         authors=s2a("W. Setyawan", "S. Curtarolo"), &
                         title="High-throughput electronic band structure calculations: Challenges and tools", &
                         source="Comput. Mater. Sci.", volume="49", pages="299-312", &
                         year=2010, doi="10.1016/j.commatsci.2010.05.010")

      CALL add_reference(key=Ahart2024, &
                         authors=s2a("C. Ahart", "S. Chulkov", "C. Cucinotta"), &
                         title="Enabling Ab Initio Molecular Dynamics under Bias: The CP2K+SMEAGOL Interface "// &
                         "for Integrating Density Functional Theory and Non-Equilibrium Green Functions", &
                         source="J. Chem. Theory Comput.", volume="20", pages="6772-6780", &
                         year=2024, doi="10.1021/acs.jctc.4c00371")

      CALL add_reference(key=Knysh2024, &
                         authors=s2a("I. Knysh", "F. Lipparini", "A. Blondel", "I. Duchemin", "X. Blase", &
                                     "P.-F. Loos", "D. Jacquemin"), &
                         title="Reference CC3 Excitation Energies for Organic Chromophores: "// &
                         "Benchmarking TD-DFT, BSE/GW, and Wave Function Methods", &
                         source="J. Chem. Theory Comput.", volume="20", pages="8152-8174", &
                         year=2024, doi="10.1021/acs.jctc.4c00906")

      CALL add_reference(key=Schambeck2024, &
                         authors=s2a("M. Schambeck", "D. Golze", "J. Wilhelm"), &
                         title="Solving multipole challenges in the GW 100 benchmark enables precise low-scaling GW calculations", &
                         source="Phys. Rev. B", volume="110", pages="125146", &
                         year=2024, doi="10.1103/PhysRevB.110.125146")

      CALL add_reference(key=Pracht2019, &
                         authors=s2a("P. Pracht", "E. Caldeweyher", "S. Ehlert", "S. Grimme"), &
                         title="A Robust Non-Self-Consistent Tight-Binding Quantum Chemistry Method for large Molecules", &
                         source="ChemRxiv", volume="", pages="", &
                         year=2019, doi="10.26434/chemrxiv.8326202.v1")

      CALL add_reference(key=Mewes2018, &
                         authors=s2a("S. Mewes", "F. Plasser", "A. Krylov", "A. Dreuw"), &
                         title="Benchmarking Excited-State Calculations Using Exciton Properties", &
                         source="J. Chem. Theory Comput.", volume="14", pages="710-725", &
                         year=2018, doi="10.1021/acs.jctc.7b01145 ")

      CALL add_reference(key=Sertcan2024, &
                         authors=s2a("B. Sertcan-Gokmen", "J. Hutter", "A. Hehn"), &
                         title="Excited-State Forces with the Gaussian and Augmented Plane Wave Method "// &
                         "for the Tamm–Dancoff Approximation of Time-Dependent Density Functional Theory", &
                         source="J. Chem. Theory Comput.", volume="20", pages="8494-8504", &
                         year=2024, doi="10.1021/acs.jctc.4c00614")

      CALL add_reference(key=Drautz2019, &
                         authors=s2a("R. Drautz"), &
                         title="Atomic cluster expansion for accurate and transferable interatomic potentials", &
                         source="Phys. Rev. B", volume="99", pages="014104", &
                         year=2019, doi="10.1103/PhysRevB.99.014104")

      CALL add_reference(key=Lysogorskiy2021, &
                         authors=s2a("Y. Lysogorskiy", "C. van der Ooord", "A. Bochkarev", "S. Menon", &
                                     "M. Rinaldi", "T. Hammerschmidt", "M. Mrovec", "A. Thompson", "G. Csanyi", &
                                     "C. Ortner", "R. Drautz"), &
                   title="Performant implementation of the atomic cluster expansion (PACE) and application to copper and silicon", &
                         source="npj Comput. Mater.", volume="7", pages="97", &
                         year=2021, doi="10.1038/s41524-021-00559-9")

      CALL add_reference(key=Bochkarev2024, &
                         authors=s2a("A. Bochkarev", "Y. Lysogorskiy", "R. Drautz"), &
                         title="Graph Atomic Cluster Expansion for semilocal interactions beyond equivariant message passing", &
                         source="Phys. Rev. X", volume="14", pages="021036", &
                         year=2024, doi="10.1103/PhysRevX.14.021036")

      CALL add_reference(key=VazdaCruz2021, &
                         authors=s2a("V. Vaz da Cruz", "S. Eckert", "A. Fohlisch"), &
                         title="TD-DFT simulations of K-edge resonant inelastic X-ray scattering within "// &
                         "the restricted subspace approximation", &
                         source="Phys. Chem. Chem. Phys.", volume="23", pages="1835-1848", &
                         year=2021, doi="10.1039/d0cp04726k")

      CALL add_reference(key=Chen2025, &
                         authors=s2a("Z. Chen", "Y. Yang"), &
                         title="Periodic Constrained Nuclear-Electronic Orbital Density Functional Theory for "// &
                         "Nuclear Quantum Effects: Method Development and Application to Hydrogen Adsorption on Pt(111)", &
                         source="J. Chem. Theory Comput.", volume="21", pages="7865-7877", &
                         year=2025, doi="10.1021/acs.jctc.5c00837")

   END SUBROUTINE add_all_references

END MODULE bibliography
