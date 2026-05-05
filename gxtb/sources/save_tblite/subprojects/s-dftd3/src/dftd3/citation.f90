! This file is part of s-dftd3.
! SPDX-Identifier: LGPL-3.0-or-later
!
! s-dftd3 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! s-dftd3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.

!> Module for handling citation data
module dftd3_citation
   implicit none
   private

   public :: citation_type, author_name, new_citation, is_citation_present
   public :: format_bibtex, get_citation, same_citation
   public :: doi_dftd3_0, doi_dftd3_bj, doi_dftd3_m, doi_dftd3_op, doi_dftd3_cso, &
      & doi_gmtkn30_0, doi_gmtkn30_bj, doi_gmtkn55, doi_dsd, doi_dsdpbep86, &
      & doi_drpa, doi_revdsd, doi_pw91_d3, doi_r2scan_d4, doi_scan_d3, &
      & doi_pbeh3c, doi_hse3c, doi_b973c, doi_hf3c, doi_gcp, doi_d3pbc, &
      & doi_r2scan_hyb, doi_r2scan_dhdf, doi_minnesota_d3, doi_b97m_d3, &
      & doi_wb97x_d3, doi_hse06_d3, doi_joss, doi_cf22d, doi_skala

   !> Represents an author to allow creating author lists
   type :: author_type
      !> Name of the author
      character(len=:), allocatable :: name
   end type author_type

   !> Represents a citation
   type :: citation_type
      !> Title of the publication
      character(len=:), allocatable :: title
      !> Authors of the publication
      type(author_type), allocatable :: author(:)
      !> Issue of the publication
      character(len=:), allocatable :: issue
      !> Journal of the publication
      character(len=:), allocatable :: journal
      !> Volume of the publication
      character(len=:), allocatable :: volume
      !> Year of the publication
      character(len=:), allocatable :: year
      !> Page numbers of the publication
      character(len=:), allocatable :: pages
      !> Digital Object Identifier (DOI) of the publication
      character(len=:), allocatable :: doi
   end type citation_type

   character(len=*), parameter :: nl = new_line('a')

   character(len=*), parameter :: &
      & doi_dftd3_0 = "10.1063/1.3382344", &
      & doi_dftd3_bj = "10.1002/jcc.21759", &
      & doi_dftd3_m = "10.1021/acs.jpclett.6b00780", &
      & doi_dftd3_op = "10.1021/acs.jctc.7b00176", &
      & doi_gmtkn30_0 = "10.1021/ct100466k", &
      & doi_gmtkn30_bj = "10.1039/c0cp02984j", &
      & doi_gmtkn55 = "10.1039/c7cp04913g", &
      & doi_dsd = "10.1002/jcc.23391", &
      & doi_dsdpbep86 = "10.1039/c1cp22592h", &
      & doi_drpa = "10.1021/acs.jpca.1c01295", &
      & doi_revdsd = "10.1021/acs.jpca.9b03157", &
      & doi_pw91_d3 = "10.1073/pnas.1516984112", &
      & doi_r2scan_d4 = "10.1063/5.0041008", &
      & doi_scan_d3 = "10.1103/physrevb.94.115144", &
      & doi_pbeh3c = "10.1063/1.4927476", &
      & doi_hse3c = "10.1039/c6cp01697a", &
      & doi_b973c = "10.1063/1.5012601", &
      & doi_hf3c = "10.1002/jcc.23317", &
      & doi_gcp = "10.1063/1.3700154", &
      & doi_d3pbc = "10.1002/cphc.201100521", &
      & doi_r2scan_hyb = "10.1063/5.0086040", &
      & doi_r2scan_dhdf = "10.1063/5.0174988", &
      & doi_minnesota_d3 = "10.1021/acs.jpclett.5b01591", &
      & doi_b97m_d3 = "10.1021/acs.jctc.8b00842", &
      & doi_wb97x_d3 = "10.1021/ct300715s", &
      & doi_hse06_d3 = "10.1021/jp501237c", &
      & doi_joss = "10.21105/joss.07169", &
      & doi_dftd3_cso = "10.1021/acs.jctc.5b00400", &
      & doi_cf22d = "10.1038/s43588-022-00371-5", &
      & doi_skala = "10.48550/arXiv.2506.14665"

contains

!> Create a new citation
pure function new_citation(title, author, journal, issue, volume, year, pages, doi) result(citation)
   !> Title of the publication
   character(len=*), intent(in) :: title
   !> Authors of the publication
   type(author_type), intent(in) :: author(:)
   !> Journal of the publication
   character(len=*), intent(in) :: journal
   !> Issue of the publication
   character(len=*), intent(in), optional :: issue
   !> Volume of the publication
   character(len=*), intent(in) :: volume
   !> Year of the publication
   character(len=*), intent(in) :: year
   !> Page numbers of the publication
   character(len=*), intent(in) :: pages
   !> Digital Object Identifier (DOI) of the publication
   character(len=*), intent(in) :: doi
   !> The new citation
   type(citation_type) :: citation

   citation%title = title
   citation%author = author
   citation%journal = journal
   citation%volume = volume
   citation%year = year
   citation%pages = pages
   citation%doi = doi
   if (present(issue)) citation%issue = issue
end function new_citation

!> Create an author
pure function author_name(name) result(author)
   !> Name of the author
   character(len=*), intent(in) :: name
   !> The new author
   type(author_type) :: author
   author%name = name
end function author_name

!> Check if citation data is present
pure function is_citation_present(citation) result(is_present)
   !> The citation to check
   type(citation_type), intent(in) :: citation
   !> Whether the citation data is present
   logical :: is_present

   is_present = allocated(citation%doi) &
      & .and.allocated(citation%title) &
      & .and.allocated(citation%author) &
      & .and.allocated(citation%journal) &
      & .and.allocated(citation%volume) &
      & .and.allocated(citation%year) &
      & .and.allocated(citation%pages)
end function is_citation_present

!> Check if two citations are the same
pure function same_citation(lhs, rhs) result(same)
   !> The first citation
   type(citation_type), intent(in) :: lhs
   !> The second citation
   type(citation_type), intent(in) :: rhs
   !> Whether the citations are the same
   logical :: same

   if (is_citation_present(lhs) .and. is_citation_present(rhs)) then
      same = lhs%doi == rhs%doi
   else
      same = .false.
   end if
end function same_citation

!> Format a citation as a BibTeX entry
subroutine format_bibtex(string, citation)
   !> The formatted BibTeX entry
   character(len=:), allocatable, intent(out) :: string
   !> The citation to format
   class(citation_type), intent(in) :: citation
   integer :: idx

   if (.not.is_citation_present(citation)) return

   string = &
      & "@article{" // citation%doi // "," // nl // &
      & "  title = {{" // citation%title // "}}," // nl // &
      & "  author = {" // citation%author(1)%name
   do idx = 2, size(citation%author)
      string = string // nl // &
         & "    and " // citation%author(idx)%name
   end do
   string = string // "}," // nl
   if (allocated(citation%issue)) then
      string = string // &
         & "  issue = {" // citation%issue // "}," // nl
   end if
   string = string // &
      & "  volume = {" // citation%volume // "}," // nl // &
      & "  pages = {" // citation%pages // "}," // nl // &
      & "  doi = {" // citation%doi // "}," // nl // &
      & "  url = {https://doi.org/" // citation%doi // "}" // nl // &
      & "}"
end subroutine format_bibtex

!> Get citation data for a given DOI
pure function get_citation(doi) result(citation)
   !> Digital Object Identifier (DOI) of the publication
   character(len=*), intent(in) :: doi
   !> The citation data
   type(citation_type) :: citation

   select case(doi)
   case(doi_dftd3_0)
      citation = new_citation( &
         doi=doi, &
         title="A consistent and accurate ab initio parametrization of density functional &
         & dispersion correction (DFT-D) for the 94 elements H-Pu", &
         author=[ &
         & author_name("Stefan Grimme"), &
         & author_name("Jens Antony"), &
         & author_name("Stephan Ehrlich"), &
         & author_name("Helge Krieg")], &
         journal="J. Chem. Phys.", &
         volume="132", &
         pages="154104", &
         year="2010" &
      )
      
   case(doi_dftd3_bj)
      citation = new_citation( &
         doi=doi, &
         title="Effect of the damping function in dispersion corrected density functional theory", &
         author=[ &
         & author_name("Stefan Grimme"), &
         & author_name("Stephan Ehrlich"), &
         & author_name("Lars Goerigk")], &
         journal="J. Comput. Chem.", &
         issue="7", &
         volume="32", &
         pages="1456-1465", &
         year="2011" &
      )
      
   case(doi_dftd3_m)
      citation = new_citation( &
         doi=doi, &
         title="Revised Damping Parameters for the D3 Dispersion Correction to Density Functional Theory", &
         author=[ &
         & author_name("Daniel G. A. Smith"), &
         & author_name("Lori A. Burns"), &
         & author_name("Konrad Patkowski"), &
         & author_name("C. David Sherrill")], &
         journal="J. Phys. Chem. Lett.", &
         issue="12", &
         volume="7", &
         pages="2197-2203", &
         year="2016" &
      )
      
   case(doi_dftd3_op)
      citation = new_citation( &
         doi=doi, &
         title="Assessing DFT-D3 Damping Functions Across Widely Used Density Functionals: Can We Do Better?", &
         author=[ &
         & author_name("Jonathon Witte"), &
         & author_name("Narbe Mardirossian"), &
         & author_name("Jeffrey B. Neaton"), &
         & author_name("Martin Head-Gordon")], &
         journal="J. Chem. Theory Comput.", &
         issue="5", &
         volume="13", &
         pages="2043-2052", &
         year="2017" &
      )
      
   case(doi_gmtkn30_0)
      citation = new_citation( &
         doi=doi, &
         title="Efficient and Accurate Double-Hybrid-Meta-GGA Density Functionals—Evaluation &
         & with the Extended GMTKN30 Database for General Main Group Thermochemistry, Kinetics, &
         & and Noncovalent Interactions", &
         author=[ &
         & author_name("Lars Goerigk"), &
         & author_name("Stefan Grimme")], &
         journal="J. Chem. Theory Comput.", &
         issue="2", &
         volume="7", &
         pages="291-309", &
         year="2011" &
      )
      
   case(doi_gmtkn30_bj)
      citation = new_citation( &
         doi=doi, &
         title="A thorough benchmark of density functional methods for general main group &
         & thermochemistry, kinetics, and noncovalent interactions", &
         author=[ &
         & author_name("Lars Goerigk"), &
         & author_name("Stefan Grimme")], &
         journal="Phys. Chem. Chem. Phys.", &
         issue="14", &
         volume="13", &
         pages="6670-6688", &
         year="2011" &
      )
      
   case(doi_gmtkn55)
      citation = new_citation( &
         doi=doi, &
         title="A look at the density functional theory zoo with the advanced GMTKN55 &
         & database for general main group thermochemistry, kinetics and noncovalent interactions", &
         author=[ &
         & author_name("Lars Goerigk"), &
         & author_name("Andreas Hansen"), &
         & author_name("Christoph Bauer"), &
         & author_name("Stephan Ehrlich"), &
         & author_name("Asim Najibi"), &
         & author_name("Stefan Grimme")], &
         journal="Phys. Chem. Chem. Phys.", &
         issue="48", &
         volume="19", &
         pages="32184-32215", &
         year="2017" &
      )
      
   case(doi_dsd)
      citation = new_citation( &
         doi=doi, &
         title="Spin-component-scaled double hybrids: An extensive search for the best fifth-rung &
         & functionals blending DFT and perturbation theory", &
         author=[ &
         & author_name("Sebastian Kozuch"), &
         & author_name("Jan M. L. Martin")], &
         journal="J. Comput. Chem.", &
         issue="27", &
         volume="34", &
         pages="2327-2344", &
         year="2013" &
      )
      
   case(doi_dsdpbep86)
      citation = new_citation( &
         doi=doi, &
         title="DSD-PBEP86: in search of the best double-hybrid DFT with spin-component scaled &
         & MP2 and dispersion corrections", &
         author=[ &
         & author_name("Sebastian Kozuch"), &
         & author_name("Jan M. L. Martin")], &
         journal="Phys. Chem. Chem. Phys.", &
         issue="45", &
         volume="13", &
         pages="20104-20107", &
         year="2011" &
      )
      
   case(doi_drpa)
      citation = new_citation( &
         doi=doi, &
         title="Exploring Avenues beyond Revised DSD Functionals: II. Random-Phase Approximation and Scaled MP3 Corrections", &
         author=[ &
         & author_name("Golokesh Santra"), &
         & author_name("Emmanouil Semidalas"), &
         & author_name("Jan M. L. Martin")], &
         journal="J. Phys. Chem. A", &
         issue="21", &
         volume="125", &
         pages="4628-4638", &
         year="2021" &
      )
      
   case(doi_revdsd)
      citation = new_citation( &
         doi=doi, &
         title="Minimally Empirical Double-Hybrid Functionals Trained against the GMTKN55 Database: &
         & revDSD-PBEP86-D4, revDOD-PBE-D4, and DOD-SCAN-D4", &
         author=[ &
         & author_name("Golokesh Santra"), &
         & author_name("Nitai Sylvetsky"), &
         & author_name("Jan M. L. Martin")], &
         journal="J. Phys. Chem. A", &
         issue="24", &
         volume="123", &
         pages="5129-5143", &
         year="2019" &
      )
      
   case(doi_pw91_d3)
      citation = new_citation( &
         doi=doi, &
         title="A priori calculations of the free energy of formation from solution of polymorphic &
         & self-assembled monolayers", &
         author=[ &
         & author_name("Jeffrey R. Reimers"), &
         & author_name("Dwi Panduwinata"), &
         & author_name("Johan Visser"), &
         & author_name("Maxwell J. Crossley")], &
         journal="Proc. Natl. Acad. Sci.", &
         issue="45", &
         volume="112", &
         pages="E6101-E6110", &
         year="2015" &
      )
      
   case(doi_r2scan_d4)
      citation = new_citation( &
         doi=doi, &
         title="r²SCAN-D4: Dispersion corrected meta-generalized gradient approximation for general &
         & chemical applications", &
         author=[ &
         & author_name("Sebastian Ehlert"), &
         & author_name("Uwe Huniar"), &
         & author_name("Jinliang Ning"), &
         & author_name("James W. Furness"), &
         & author_name("Jianwei Sun"), &
         & author_name("Aaron D. Kaplan"), &
         & author_name("John P. Perdew"), &
         & author_name("Jan Gerit Brandenburg")], &
         journal="J. Chem. Phys.", &
         volume="154", &
         pages="061101", &
         year="2021" &
      )
      
   case(doi_scan_d3)
      citation = new_citation( &
         doi=doi, &
         title="Benchmark tests of a strongly constrained semilocal functional with a &
         & long-range dispersion correction", &
         author=[ &
         & author_name("J. G. Brandenburg"), &
         & author_name("J. E. Bates"), &
         & author_name("J. Sun"), &
         & author_name("J. P. Perdew")], &
         journal="Phys. Rev. B", &
         volume="94", &
         pages="115144", &
         year="2016" &
      )
      
   case(doi_pbeh3c)
      citation = new_citation( &
         doi=doi, &
         title="Consistent structures and interactions by density functional theory &
         & with small atomic orbital basis sets", &
         author=[ &
         & author_name("Stefan Grimme"), &
         & author_name("Jan Gerit Brandenburg"), &
         & author_name("Christoph Bannwarth"), &
         & author_name("Andreas Hansen")], &
         journal="J. Chem. Phys.", &
         volume="143", &
         pages="054107", &
         year="2015" &
      )
      
   case(doi_hse3c)
      citation = new_citation( &
         doi=doi, &
         title="Screened exchange hybrid density functional for accurate and efficient &
         & structures and interaction energies", &
         author=[ &
         & author_name("Jan Gerit Brandenburg"), &
         & author_name("Eike Caldeweyher"), &
         & author_name("Stefan Grimme")], &
         journal="Phys. Chem. Chem. Phys.", &
         issue="23", &
         volume="18", &
         pages="15519-15523", &
         year="2016" &
      )
      
   case(doi_b973c)
      citation = new_citation( &
         doi=doi, &
         title="B97-3c: A revised low-cost variant of the B97-D density functional method", &
         author=[ &
         & author_name("Jan Gerit Brandenburg"), &
         & author_name("Christoph Bannwarth"), &
         & author_name("Andreas Hansen"), &
         & author_name("Stefan Grimme")], &
         journal="J. Chem. Phys.", &
         volume=" 148", &
         pages="064104", &
         year="2018" &
      )
      
   case(doi_hf3c)
      citation = new_citation( &
         doi=doi, &
         title="Corrected small basis set Hartree-Fock method for large systems", &
         author=[ &
         & author_name("Rebecca Sure"), &
         & author_name("Stefan Grimme")], &
         journal="J. Comput. Chem.", &
         issue="19", &
         volume="34", &
         pages="1672-1685", &
         year="2013" &
      )
      
   case(doi_gcp)
      citation = new_citation( &
         doi=doi, &
         title="A geometrical correction for the inter- and intra-molecular basis set &
         & superposition error in Hartree-Fock and density functional theory calculations for large systems", &
         author=[ &
         & author_name("Holger Kruse"), &
         & author_name("Stefan Grimme")], &
         journal="J. Chem. Phys.", &
         volume="136", &
         pages="154101", &
         year="2012" &
      )
      
   case(doi_d3pbc)
      citation = new_citation( &
         doi=doi, &
         title="System-Dependent Dispersion Coefficients for the DFT-D3 Treatment of &
         & Adsorption Processes on Ionic Surfaces", &
         author=[ &
         & author_name("Stephan Ehrlich"), &
         & author_name("Jonas Moellmann"), &
         & author_name("Werner Reckien"), &
         & author_name("Thomas Bredow"), &
         & author_name("Stefan Grimme")], &
         journal="ChemPhysChem", &
         issue="17", &
         volume="12", &
         pages="3414", &
         year="2011" &
      )
      
   case(doi_r2scan_hyb)
      citation = new_citation( &
         doi=doi, &
         title="Dispersion corrected r²SCAN based global hybrid functionals: r²SCANh, r²SCAN0, and r²SCAN50", &
         author=[ &
         & author_name("Markus Bursch"), &
         & author_name("Hagen Neugebauer"), &
         & author_name("Sebastian Ehlert"), &
         & author_name("Stefan Grimme")], &
         journal="J. Chem. Phys.", &
         volume="156", &
         pages="134105", &
         year="2022" &
      )
      
   case(doi_r2scan_dhdf)
      citation = new_citation( &
         doi=doi, &
         title="Dispersion-corrected r²SCAN based double-hybrid functionals", &
         author=[ &
         & author_name("Lukas Wittmann"), &
         & author_name("Hagen Neugebauer"), &
         & author_name("Stefan Grimme"), &
         & author_name("Markus Bursch")], &
         journal="J. Chem. Phys.", &
         volume="159", &
         pages="224103", &
         year="2023" &
      )
      
   case(doi_minnesota_d3)
      citation = new_citation( &
         doi=doi, &
         title="Treating London-Dispersion Effects with the Latest Minnesota Density Functionals: &
         & Problems and Possible Solutions", &
         author=[ &
         & author_name("Lars Goerigk")], &
         journal="J. Phys. Chem. Lett.", &
         issue="19", &
         volume="6", &
         pages="3891-3896", &
         year="2015" &
      )
      
   case(doi_b97m_d3)
      citation = new_citation( &
         doi=doi, &
         title="The Nonlocal Kernel in van der Waals Density Functionals as an Additive Correction: &
         & An Extensive Analysis with Special Emphasis on the B97M-V and ωB97M-V Approaches", &
         author=[ &
         & author_name("Asim Najibi"), &
         & author_name("Lars Goerigk")], &
         journal="J. Chem. Theory Comput.", &
         issue="11", &
         volume="14", &
         pages="5725-5738", &
         year="2018" &
      )
      
   case(doi_wb97x_d3)
      citation = new_citation( &
         doi=doi, &
         title="Long-Range Corrected Hybrid Density Functionals with Improved Dispersion Corrections", &
         author=[ &
         & author_name("You-Sheng Lin"), &
         & author_name("Guan-De Li"), &
         & author_name("Shan-Ping Mao"), &
         & author_name("Jeng-Da Chai")], &
         journal="J. Chem. Theory Comput.", &
         issue="1", &
         volume="9", &
         pages="263-272", &
         year="2013" &
      )
      
   case(doi_hse06_d3)
      citation = new_citation( &
         doi=doi, &
         title="DFT-D3 Study of Some Molecular Crystals", &
         author=[ &
         & author_name("Jonas Moellmann"), &
         & author_name("Stefan Grimme")], &
         journal="J. Phys. Chem. C", &
         issue="14", &
         volume="118", &
         pages="7615-7621", &
         year="2014" &
      )
      
   case(doi_joss)
      citation = new_citation( &
         doi=doi, &
         title="Simple DFT-D3: Library first implementation of the D3 dispersion correction", &
         author=[ &
         & author_name("Sebastian Ehlert")], &
         journal="J. Open Source Softw.", &
         issue="103", &
         volume="9", &
         pages="7169", &
         year="2024" &
      )

   case(doi_cf22d)
      citation = new_citation( &
         doi=doi, &
         title="Supervised learning of a chemistry functional with damped dispersion", &
         author=[ &
         & author_name("Yiwei Liu"), &
         & author_name("Cheng Zhang"), &
         & author_name("Zhonghua Liu"), &
         & author_name("Donald G. Truhlar"), &
         & author_name("Ying Wang"), &
         & author_name("Xiao He")], &
         journal="Nat. Comput. Sci.", &
         issue="1", &
         volume="3", &
         pages="48–58", &
         year="2023" &
      )

   case(doi_skala)
      citation = new_citation( &
         doi=doi, &
         title="Accurate and scalable exchange-correlation with deep learning", &
         author=[ &
         & author_name("Giulia Luise"), &
         & author_name("Chin-Wei Huang"), &
         & author_name("Thijs Vogels"), &
         & author_name("Derk P. Kooi"), &
         & author_name("Sebastian Ehlert"), &
         & author_name("Stephanie Lanius"), &
         & author_name("Klaas J. H. Giesbertz"), &
         & author_name("Amir Karton"), &
         & author_name("Deniz Gunceler"), &
         & author_name("Megan Stanley"), &
         & author_name("Wessel P. Bruinsma"), &
         & author_name("Lin Huang"), &
         & author_name("Xinran Wei"), &
         & author_name("Jos{\'e} Garrido Torres"), &
         & author_name("Abylay Katbashev"), &
         & author_name("Rodrigo Chavez Zavaleta"), &
         & author_name("B{\'a}lint M{\'a}t{\'e}"), &
         & author_name("S{\'e}kou-Oumar Kaba"), &
         & author_name("Roberto Sordillo"), &
         & author_name("Yingrong Chen"), &
         & author_name("David B. Williams-Young"), &
         & author_name("Christopher M. Bishop"), &
         & author_name("Jan Hermann"), &
         & author_name("Rianne van den Berg"), &
         & author_name("Paola Gori-Giorgi")], &
         journal="arXiv", &
         volume="2506.14665", &
         pages="2506.14665", &
         year="2025" &
      )

   case(doi_dftd3_cso)
      citation = new_citation( &
         doi=doi, &
         title="Reformulation of the D3(Becke-Johnson) Dispersion Correction "// &
         & "without Resorting to Higher than C6 Dispersion Coefficients", &
         author=[ &
         & author_name('Heiner Schr{\"o}der'), &
         & author_name("Jens Creon"), &
         & author_name("Tobias Schwabe")], &
         journal="J. Chem. Theory Comput.", &
         volume="11", &
         pages="3163--3170", &
         year="2015" &
      )
   end select
end function get_citation

end module dftd3_citation
