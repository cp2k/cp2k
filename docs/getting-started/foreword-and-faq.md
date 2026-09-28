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
features in the `QUICKSTEP` module of CP2K, as elaborated on [](../methods/dft/k-points). This is an
very active field of development, with many new implementations and performance enhancement becoming
available only recently; so, stay tuned and do not get distracted by outdated unofficial accounts.

## What is the software ecosystem of CP2K like?

There is a large arsenal of third-party auxiliary programs and libraries, interfaced with CP2K in
one form or another, that offer supplementary utilities and can be readily integrated into the
computation and pre- and post-processing workflow. See for instance:

- The aforementioned [](#Iannuzzi2026) that mentions in section 11.1.3 a number of such tools;
- The list of [tools for simplifying your life with CP2K](https://www.cp2k.org/tools) on cp2k.org;
- And the `technologies` section of this manual.

## Where can I reach out to the community for asking questions, making contributions, etc.?

Please find the `SUPPORT.md` document for discussion venues, recommended practice for requesting
help, and possible ways of contribution. These contents are previously on this page but have been
refactored for brevity.
