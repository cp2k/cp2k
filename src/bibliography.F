!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2007  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****h* cp2k/bibliography *
!!
!!   NAME
!!     bibliography
!!
!!   FUNCTION
!!     collects all references to literature in CP2K
!!     as new algorithms / method are included from literature sources
!!     we add a citation, to give credit to the authors
!!     (often ourselves, but this need not to be the case), and as a form
!!     of documentation.
!!
!!     reference_manager.F provides an easy way to cite these references from the code
!!     (e.g. whenever the relevant part of the code is executed)
!!     and to add citations to the input
!!     at the end of the run a bibliography is printed that can be used e.g. in papers
!!
!!   NOTES
!!     references need to be input using the ISI citation format, because it is
!!     uniform, easy to parse, and can be exported for example from web of science.
!!     (mark the reference, and click 'export to reference software', open the resulting file)
!!     Furthermore, it can be easily converted to and from using the bibutils tools
!!     a collection of easy to use conversion programs that can be found at
!!     http://www.scripps.edu/~cdputnam/software/bibutils/ by Chris Putnam
!!
!!   AUTHOR
!!     Joost VandeVondele
!!     
!!
!!   MODIFICATION HISTORY
!!     08.2007 [Joost VandeVondele] Created
!!
!!   SOURCE
!****************************************************************************
MODULE bibliography
  USE reference_manager, ONLY: add_reference, cite_reference 
  USE string_utilities, ONLY: s2a
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cite_reference ! for convenience
  PUBLIC :: add_all_references ! f77_interface only

  ! all reference keys
  INTEGER, PUBLIC, SAVE ::  VandeVondele2005a, VandeVondele2003, Laino2005, Laino2006, &
             Goedecker1996, Hartwigsen1998, Krack2005, Lippert1997, Lippert1999, Krack2000,&
             Krack2002, Iannuzzi2005, Iannuzzi2006, Iannuzzi2007, Toukmaji1996, Martyna1999,&
             VandeVondele2005b, Perdew1981,Avezac2005, Porezag1995, Seifert1996, Elstner1998, Zhechkov2005

CONTAINS
!!****f* bibliography/add_all_references *
!!
!!   NAME
!!    add_all_references
!!
!!   FUNCTION
!!    adds references that can later be cited / printed using the key
!!
!!   NOTES
!!    - note that the ISI record is fixed format (line length limited and the
!!      first 3 characters can indicate record type)
!!    - please add DOI whenever available, this will result in a clickable
!!      link in the input reference manual.
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     08.2007 created [Joost VandeVondele]
!!
!!*** **********************************************************************
  SUBROUTINE add_all_references()

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

  END SUBROUTINE add_all_references

END MODULE
