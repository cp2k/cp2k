# Community support

The CP2K community of end users, developers, and other friends and colleagues uses several venues
for discussions, on which you can ask for help, request for features, share feedbacks and ideas and
tools, and learn about new releases and upcoming events.

1. The official CP2K User Forum hosted on Google Groups:
   - Access from <https://groups.google.com/group/cp2k>.
   - Sign in Google Groups with a Google account, apply to join and wait for approval.
   - This is the main place for most of the daily discussions in the community; developers are ready
     to offer support on the *technical* side.
   - Additionally, a read-only [mirror](https://lists.cp2k.org/listinfo/cp2k-user) and
     [archives](https://lists.cp2k.org/archives/cp2k-user/) for the discussion threads are available
     and does not need a Google account to view.
1. The GitHub repositories of the CP2K organization:
   - In the main `cp2k/cp2k` repository there are [issues](https://github.com/cp2k/cp2k/issues/) and
     [discussions](https://github.com/cp2k/cp2k/discussions).
   - These are for topics most relevant to the program development and code implementation, such as
     *reproducible* bug reports, well-defined feature requests and revisions to the documentation or
     manual. If uncertain, use the User Forum above first and let developers decide whether to bring
     up the matter to GitHub.
1. Other third-party generic venues are also noted to have a section about CP2K usage:
   - [Matter Modeling Stack Exchange](https://mattermodeling.stackexchange.com/questions/tagged/cp2k)
     has a tag for CP2K.
   - For Chinese users, there is also a CP2K category in the First-principles subforum of the
     [Computational Chemistry Commune](http://bbs.keinsci.com/forum-105-1.html?typeid=42).

## I want to ask a question

The general etiquette for requesting tech support online has been summarized nicely by Eric S.
Raymond's [How To Ask Questions The Smart Way](http://www.catb.org/~esr/faqs/smart-questions.html).
(**Disclaimer**: this link does not imply any connection between the original author and the CP2K
developers, nor does it suggest that the original author may be contacted for assistance.)

Before submitting a question, please refer to the existing resources in the very first place and see
if it has already been addressed somewhere:

- The official [CP2K manual](https://manual.cp2k.org/) of method documentations and input references
  have been renovated and expanded with more coverage. Be aware that there are distinctive sets of
  input references, with [the "trunk" one](https://manual.cp2k.org/trunk/CP2K_INPUT.html) for the
  latest *development* version; if instead a past stable release is used, select the corresponding
  one from [manual by version](https://manual.cp2k.org/trunk/versions.html) before reading.
  - The [troubleshooting](https://manual.cp2k.org/trunk/getting-started/troubleshooting.html) page
    has a catalog of frequently encountered warning and error messages and is worth checking too.
- The venues mentioned above may also have questions from similar scenarios submitted by others in
  the past that have been answered, either with detailed suggestions or just some general tips on
  how to diagnose and adjust.

> [!CAUTION]
>
> For the time being, it is not recommended to seek for unofficial CP2K-specific suggestions from
> generic large language models (LLM). Even if they have been trained on a refined and verified
> corpus of CP2K materials one day, they can still hallucinate and generate superficially convincing
> but factually incorrect responses, which are distracting and potentially misleading future users.
> Unless capable of validate and willing to take the responsibility for any content produced by
> artificial intelligence, as with human authors of a formal academic publication, do not bother
> mentioning anything directly from AI in the discussion at all. Consider the more traditional tools
> of a search engine or a translator.

> [!NOTE]
>
> Please kindly understand that, despite the CP2K developers having knowledge about the algorithm
> infrastructures and program implementations, they may not be suitable for answering every kind of
> questions arising from practice, in particular those pertaining to niche application-oriented
> research areas where apprehending the science and acquiring the skills will require much more
> extensive academic training than learning to use a program. Oftentimes the best parties to consult
> for in-depth guidance would be the tutor, advisor, experienced colleagues or collaborators in real
> life, and when attempting to reproduce some findings reported in literature, the original authors.
> Chatting for free with volunteers via public online venues cannot substitute formally contacting
> and communicating with the specialists that are appropriately professional and credible.

If it has been determined that a new question is necessary, please compose and phrase it with a
sufficient level of details, accuracy, and clarity. Approach the process in the same way as making a
presentation to general audience, or even writing the "Methods" section in a formal English academic
publication; this includes giving explanations to uncommon acronyms (say, the abbreviated name of a
specific class of materials, or anything else that is not on our page for acknowledged
[Acronyms](https://manual.cp2k.org/trunk/acronyms.html)) and traceable citations (with publication
title and DOI link, instead of merely showing a screenshot or a paragraph of copy-pasted text). Here
is a checklist for the details that are usually expected to be included in the description:

- The release or git version and the installation method of CP2K, especially the compiler version;
- Custom revisions to the source code, if any;
- For problems related to installation, the hardware specification, the configuration for linked
  libraries, and the distribution source and means of preparation (say, with some package managers,
  environment-controlling modules, or just a build from source);
- For bad performance, the resource allocation of CPU/GPU cores and RAM, the method of execution
  (say, parallel or serial, process-level MPI or thread-level OpenMP or a hybrid parallel run), the
  core binding pattern and whether any job schedulers/management system are in use, and the final
  TIMING report at the end of an output log.
- For error terminations and wrong results, an expected behavior or reference value, as well as a
  *complete* input deck and the output files, encompassing not only the main input file with keyword
  settings, but also all of the external files referenced inside unless they are available under the
  official `data` directory, so that the job can be actually run and tested on the developers' side.

> [!TIP]
>
> - It is encouraged to try out the latest development version from the master branch of the GitHub
>   repository whenever situation permits, as it is likely to have some patches resolving the known
>   problems, and if not, works on which will benefit the next release version. This is not denying
>   the support for legacy versions, but encouraging updates to at least decently recent releases.
> - The input file does not have to use the intended chemical structure, composition, box size, and
>   number of CPU cores and amount of memory allocated in the original encounter. For the
>   [minimal reproducer](https://en.wikipedia.org/wiki/Minimal_reproducible_example), any simplified
>   system is fine and the accuracy-controlling parameters can be tuned down, as long as the input
>   reliably triggers the problem. Refer to the input files under the `tests` directory for typical
>   systems used in development testing. Not only would this reduce the demand on computational
>   resources while reproducing, but also confidential research information would not be disclosed.

## I want to contribute something

Thank you for your interest and efforts. So, potential forms of contribution include:

- Engaging in the discussions in the venues mentioned above;
- Participating the project development, as instructed on [Onboarding](CONTRIBUTING.md);
- Enriching the [cp2k-examples](https://github.com/cp2k/cp2k-examples) repository with example
  inputs, outputs, pre- and post-analysis scripts. Interpretation and discussion of the results from
  the program to complete the workflow would be welcome.
- Also, sharing some representative input files as well as structures as supplementary materials in
  a publication will not only help other curious readers see the full potential of CP2K in terms of
  scientific and engineering applications but also bridge the gap between theories and input syntax.

Please note that the main `cp2k/cp2k` repository is distributed under the [GPL-2.0 license](LICENSE)
in its entirety.
