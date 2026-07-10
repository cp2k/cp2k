# Foreword and FAQ

It is our great pleasure to present CP2K, an open-source software package for _ab initio_ electronic
structure calculations in atomistic simulations. The code is written in Fortran 2008 and has been
geared towards large-scale, high-performance CPU and GPU computation with multi-threading, MPI, CUDA
and HIP parallelization. For an overview of the capabilities, see
[Features](https://www.cp2k.org/features).

While CP2K started as an implementation of quantum chemical methods (more specifically, the
`QUICKSTEP` module, as presented at [](#VandeVondele2005)) for molecular dynamics simulation,
decades of ensuing development has witnessed a vast team of collaborators with their innumerable
contributions and the ever-growing user base with their valuable feedbacks, to whom we wish to
express sincere gratitude. As CP2K is freely available in various ways and does not requite
registration to use, it is difficult to gather accurate usage stats; but the
[list of publications using CP2K](https://www.cp2k.org/science) speaks for itself.

We would like to ask you, users of CP2K, to acknowledge our work by citing the publications as
listed on the [Bibliography](../bibliography) page and printed as REFERENCES at the end of output
log of the program, in particular the review articles:

- [](#K%C3%BChne2020), on the theoretical background and algorithms;
- [](#Iannuzzi2026), on the practical usage and applications.

We have prepared a list of Q&A for frequently asked things below, which we hope can be helpful for
the experience with the CP2K package and the art of computational chemistry in general.

This program is provided "as-is" without any expressed or implied warranty.

## Firstly, what does the name CP2K stand for?

Simply put, "CP" means Car-Parrinello, the initials of two scientists, and "2K" means year 2000.

Historically there were two formulations developed for _ab initio_ molecular dynamics ({term}`MD`):
the Car-Parrinello Molecular Dynamics ({term}`CPMD`), and the Born-Oppenheimer Molecular Dynamics
({term}`BOMD`). A program named simply also as `CPMD` began its development back in the 1990s
featuring the Car-Parrinello Molecular Dynamics; the sister project, named as CP2K, would become its
spiritual successor in the wake of the New Millenium.

Due to the fact that the original CPMD formulation has not actually been implemented yet, the name
CP2K may sound slightly non-indicative. The BOMD formulation is the major one employed by CP2K and
has seen mainstream applications in a variety of fields in the 21st century.

## Can I try CP2K out somewhere before installation?

Yes. The [CP2K Lab](https://lab.cp2k.com/) is a spin-off commercial platform set up by developers
for building structures, writing input files, executing jobs on the cloud and analyzing the outputs.
Once signed up for free, the free-tier features already allow for experiments with lightweight
computation like the one in [](./first-calculation). For those in need of more resources and
functionalities, the platform offers paid tier, site license and on-premise enterprise support.

## What preliminary knowledge does using CP2K need?

Practically CP2K is built and executed in some Linux-based operating systems, ranging from on
physical high-performance computers for production to in virtual machines for quick small tests.
This implies the need of Linux knowledge including its file system, paths, user privileges and
permissions, environment variables, shells (most commonly Bash and POSIX), utility commands,
stdin/stdout/stderr, piping and redirects, shell scripts, and modules and library files. In
addition, having some experience with Fortran, C, and C++ compilers as well as CMake will be helpful
for configuration and installation. Optionally, learn about upper-level management via job
schedulers and queue systems.

On the science side, introductory courses on chemistry, solid-state physics, and statistical
mechanics are vital prerequisites before carrying out computer simulations just as before performing
experiments in real life. Moreover, it is mandatory to have a clear understanding about the
theoretical methods in simulation; their characteristics and performance, strengths and limitations
should be described in the publications in the original conception and follow-up benchmarks. There
are two polar opposite pitfalls to be avoided: it is easy to overlook the subtleties and adjust
input settings mindlessly hoping that the black box somehow works, but it is also easy to become
absorbed in the maths and spend a lot of time trying to work out the equations that is not the focus
of the actual research project.

```{note}
Worth stressing are two overarching aspects of computer simulation:

- In spite of ever-growing scientific computing power, most of the time it is not affordable to have
  an exact 1:1 computational model of the real-life phenomena of interest. The usual practice is to
  use a much scaled-down model with limited number of atoms and finite length of trajectory for the
  simulation, which should achieve the delicate balance between representativeness and feasibility.
  The discrepancy of space and time scales between simulation and reality can be easily neglected due
  to a lack of awareness of the kinetics, especially for slow, rare events with high energy barrier
  that can only be observed in an extended period of time in real life.
- There is no need to worry if a theoretical method is strictly *ab initio* or not; both styles of
  deriving methods, "starting from physically rigorous and universal first principles" and "taking
  empirical results into account by fitting parameters with extra data", are capable of producing
  useful algorithms and accurate results depending on the case. The real concern is better put on
  the performance of methods on the target system of interest, which should have been benchmarked
  in existing works of the particular subdivision of science; also, the similarity between primitive
  datasets on which empirical parameters of the method (if any) are fitted and the target system
  of interest can be telling.
```

## How do I create the atomistic model for CP2K input?

This is done with external visualization and construction programs. Considering that a task of
geometry and/or cell optimization is usually the very first CP2K job, some general rules are
discussed on the
[relevant documentation page](../methods/optimization/geometry_and_cell_opt.md#starting-structure-and-cell).

If available, *computational* databases and benchmark sets are the most recommended avenue to obtain
structures due to having already been subject to some electronic-structure calculation. Even the
cheap methods and loose thresholds in a high-throughput screening and optimization can make the
structure qualitatively reasonable by chemical and physical intuitions, although further
optimization is still needed.

On the other hand, structures that are from *experimental* characterization are frequently not
"computation-ready", and thus should not be subject to computation without careful validation in
pre-processing. This can be prominent for `cif` and `pdb` structures determined by powder or single-
crystal XRD which can be affected by sample quality and thermal motion.

- Watch out for crystallographic disorder and atoms with low resolution or fractional occupation:
  using the superposition of all atoms as if every occupancy is 1.00 is highly likely to introduce
  contacting or even overlapping atoms.
- Beware of composition: the atomic structure may not match the intended macroscopic, charge-neutral
  chemical formula, owing to missing or duplicated hydrogen atoms, small counter ions, solvent or
  ligand molecules.

Possible resolutions vary from simple manual editing in the modelling stage, to utilization of
supercells and enumeration of special quasirandom structures (common for materials with dopants),
and to more rigorous XRD refinement and application of quantum crystallography methods. It is
believed that further advancements in instrumental analysis and structure resolution techniques
would eventually benefit computational chemistry greatly.

## Do I need PBC for my model?

**Periodic boundary condition** ({term}`PBC`) is a fundamental feature of CP2K, covering the full
range of dimensionalities of translational symmetry from 3D, 2D, 1D to 0D. The key distinction is
how connectivity, neighbor lists and integration grids are generated, how the Poisson solver handles
the electrostatic interaction, and how translational and rotational degrees of freedom of the center
of mass (i.e. collective motion as a whole) are treated.

If the structure involves condensed-phase matter, such as liquid solution, solid crystal, surface
slab and other one- and two-dimensional nano-materials, then generally PBC is used. This is also
applicable to systems with no actual well-defined repeating units like the bulk solutions. A huge
liquid droplet in the gaseous phase, where the diameter is so large that the gas-liquid interface is
almost flat and surface tension is negligible, may just as well be modelled as a combination of a
bulk solution system and an interface between a gaseous/vacuum region and a thin layer of solution,
both of which make use of PBC even though the liquid droplet itself is not periodic. However, it may
be necessary to validate the size of PBC against target properties to confirm that it is
sufficiently large for sampling, sometimes with the minimum image convention in mind.

Isolated molecular clusters in the gaseous phase or vacuum, where external pressure is irrelevant,
can be simulated without PBC. A frequent question is why a molecule optimized in vacuum does not
match its crystal structure; this is because the ordered packing pattern in the crystalline form
creates an environment capable of driving conformational changes. Oftentimes literatures convert a
periodic structure to an isolated model of finite size and apply modifications on the edge in the
form of terminal capping atoms/groups or point charges; these treatments are usually intended to
adapt the structure to quantum chemical softwares with no PBC support, but in CP2K they may not
offer extra advantages over an appropriate PBC for translational symmetry.

In certain cases, the same process can be simulated both with and without PBC. For example, the
reaction between hydroxyl and hydrogen may be modelled as a single $\mathrm{H_2}$ molecule colliding
with a single $\mathrm{OH}$ molecule with different relative orientations, distances and velocities,
which does not need PBC, or modelled as a mixture of numerous $\mathrm{H_2}$ and $\mathrm{OH}$
molecules, which needs PBC. Their behavior regarding responses to external conditions including
temperature, pressure, and any form of energy input may be different, but they provide insights from
distinct perspectives.

## Does CP2K support k-points?

As an essential element for solid-state electronic structure, k-point sampling is supported for some
features in the `QUICKSTEP` module of CP2K, as elaborated on [](../methods/dft/k-points).

## Where can I meet the CP2K community?

Several discussion venues are available:

- The [User Forum](https://groups.google.com/group/cp2k) hosted on Google Groups, with a read-only
  [mirror](https://lists.cp2k.org/listinfo/cp2k-user) and a downloadable
  [archives](https://lists.cp2k.org/archives/cp2k-user/). To use the forum, sign in with a Google
  account, apply to join and then wait for approval.
- The [issues](https://github.com/cp2k/cp2k/issues/) and
  [discussions](https://github.com/cp2k/cp2k/discussions) of the official github repository.
- The
  [Matter Modeling Stack Exchange](https://mattermodeling.stackexchange.com/questions/tagged/cp2k)
  has, among other topics, a tag for CP2K.
- For Chinese users, there is also a CP2K category in the First-principles subforum of the
  [Computational Chemistry Commune](http://bbs.keinsci.com/forum-105-1.html?typeid=42).

Please note that the github issues and discussions are only intended for topics relevant to the
program development and code implementation, such as reproducible bug reports, well-defined feature
requests and revisions to the documentation or manual. For more general help on the usage, as well
as unexpected behaviors that may or may not be bugs, check the other venues first; experts can
handle the questions and determine if they are eligible to be brought to github issues.

## What is the best practice to ask questions?

The general etiquette for requesting tech support online has been summarized nicely by Eric S.
Raymond's [How To Ask Questions The Smart Way](http://www.catb.org/~esr/faqs/smart-questions.html).
(**Disclaimer**: this link does not imply any connection between the original author and the CP2K
developers, nor does it suggest that the original author may be contacted for assistance.)

In the very first place, please refer to the [](./troubleshooting) page for a catalog of well-known
warning and error messages with readily available explanations and suggestions. Search through the
venues mentioned above for similar questions from others, and see if there are recent answers and
advice applicable to the scenario. For the time being, it is not recommended to seek for unofficial
CP2K-specific suggestions from generic large language models (LLM); even if they have been trained
on a refined and verified corpus of CP2K materials one day, they can still hallucinate and generate
superficially convincing but factually incorrect responses. (Unless willing to take responsibility
for the correctness of any content produced by artificial intelligence as with human authors of a
formal academic publication, do not bother mentioning anything from AI in the discussion at all.)

Before submitting a question, please compose it with sufficient details, accuracy, and clarity.
Approach the process in the same way as making a presentation to general audience, or even writing
the "Methods" section in a formal academic publication; this includes giving explanations to
uncommon acronyms (say, the abbreviated name of a specific class of materials, or anything that is
not on the [Acronyms](../acronyms) page) and traceable citations (with publication title, date, and
DOI link, instead of merely showing a screenshot or a paragraph of copy-pasted text).

The release date or git version of CP2K, and custom revisions if any, has to be mentioned at the
beginning. It is encouraged to try out the latest development version from the master branch of the
github repository whenever situation permits, as this is likely containing the resolution patches
already, and if not, works on which will benefit the next release version. Be aware that there are
distinctive sets of manuals, with [](../CP2K_INPUT) for the latest development version and
[](../versions) for the past releases; check twice if a page matches the program actually used prior
to reading it.

For problems related to installation and/or performance, the hardware specification and the
configuration for linked libraries should be explained. The distribution source and means of
preparation of dependencies, like with package managers, environment-controlling modules, or just a
build from source, need clarifying. Faulty libraries are unfortunately very common that problems may
be localized to a machine X or with a dependency Y, or even in a period of time Z with certain
external concurrent processes or other users intervening; try ruling out these factors first.

For error terminations and wrong results, it is imperative to provide a complete input deck and the
output files. The "input deck" encompasses not only the main input file with keyword settings, but
also all of the external files referenced inside unless they are available under the official `data`
directory, so that the job can be actually run and tested on the developers' side. Suspected wrong
results should have the precise location in the output and the reference expectation pointed out.

```{note}
The input file does not have to use the intended chemical structure and composition in the original
encounter. For the [minimal reproducer](https://en.wikipedia.org/wiki/Minimal_reproducible_example),
any simplified system is fine and the accuracy-controlling parameters can be tuned down, as long as
the input can reliably trigger the problem. Not only would this reduce the demand on computational
resources while reproducing, but also confidential research information would not be disclosed.
```

Lastly, please kindly understand that, despite the CP2K developers having knowledge about the
algorithm infrastructures and program implementations, they may not be suitable for answering all of
the questions arising from practice, especially those pertaining to niche research areas where
apprehending the science and acquiring the skills will require much more extensive academic training
than learning to use a program. The best party to consult for guidance of this type would be the
tutor, advisor, experienced colleagues or collaborators in real life, and when attempting to
reproduce reported findings, the original authors. This is not denying any personal potential to
teach oneself at no cost, but rather hinting the necessity of communicating with the right
professional people which does not have substitutes.

## What can I do for the community?

Potential forms of contribution, apart from engaging in the discussions, include:

- Participating the project development as instructed on [](../development/onboarding.md);
- Enriching the [cp2k-examples](https://github.com/cp2k/cp2k-examples) repository with example
  inputs, outputs, pre- and post-analysis scripts. Interpretation and discussion of the results from
  the program to complete the workflow would be nice to have.

It is also strongly advised to share the input files as well as structures as supplementary
materials in a publication. This will not only help other curious readers see the full potential of
CP2K in terms of scientific and engineering applications, but also bridge the gap between
theoretical configurations and input setup syntax.
