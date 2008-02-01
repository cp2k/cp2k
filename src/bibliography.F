!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2007 - 2008 CP2K developers group                                 !
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
!> \author Joost VandeVondele
! *****************************************************************************
MODULE bibliography
  USE reference_manager,               ONLY: add_reference,&
                                             cite_reference
  USE string_utilities,                ONLY: s2a

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cite_reference ! for convenience
  PUBLIC :: add_all_references ! f77_interface only

  ! all reference keys
  INTEGER, PUBLIC, SAVE ::  VandeVondele2005a, VandeVondele2003, Laino2005, Laino2006, &
             Goedecker1996, Hartwigsen1998, Krack2005, Lippert1997, Lippert1999, Krack2000,&
             Krack2002, Iannuzzi2005, Iannuzzi2006, Iannuzzi2007, Toukmaji1996, Martyna1999,&
             VandeVondele2005b, Perdew1981,Avezac2005, Porezag1995, &
             Seifert1996, Elstner1998, Zhechkov2005, Henkelman1999,&
             Frigo2005, Nose1984a, Nose1984b, Bussi2007,Genovese2006,&
             Genovese2007,Evans1983,Minary2003,Byrd1995,&
             VandeVondele2007,Ortiz1994,Becke1988,Perdew1996,Zhang1998,Lee1988, &
             Heyd2006,Vydrov2006,Heyd2003,Heyd2004,Vosko1980, Aguado2003,&
             Essmann1995,Ewald1921,Darden1993, Tersoff1988, Tosi1964a, Tosi1964b,&
             Yamada2000,Dudarev1997,Dudarev1998

CONTAINS

! *****************************************************************************
!> \brief adds references that can later be cited / printed using the key
!> \note
!>     - note that the ISI record is fixed format (line length limited and the
!>       first 3 characters can indicate record type)
!>     - please add DOI whenever available, this will result in a clickable
!>       link in the input reference manual.
!> \par History
!>      08.2007 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE add_all_references()

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
            DOI="10.1139/P80-159")       ! this doi should be correct later this week, but does not yet resolve
           !DOI="10.1139/cjp-58-8-1200") ! this doi is on the publishers web site but does not resolve

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

   END SUBROUTINE add_all_references

END MODULE
