# Contributing to tblite

First off, thank you for considering contributing to `tblite`.
Please take a moment to review this guidelines to make the contribution process
simple and effective for all involved.

Respecting these guidelines helps communicate that you respect the time of
the developers who manage and develop this open source project.
In return, they should return this respect by addressing your problem,
evaluating changes, and helping you handle your pull requests.


## Reporting a Bug

A bug is a *demonstratable problem* caused by the code in this repository.
Good bug reports are extremely valuable for us - thank you!

Before opening a bug report:

1. Check if the issue has already been reported.
2. Check if it still is an issue or has already been fixed?
   Try to reproduce it with the latest version from the `master` branch.
3. Isolate the problem and create a reduced test case.

A good bug report should not leave others needing to chase you up for more
information. So please try to be as detailed as possible in your report,
answer at least these questions:

1. Which version of `tblite` are you using? The current version is always
   a subject to change, so be more specific.
2. What is your environment (your laptop, the cluster of the university)?
3. What steps will reproduce the issue?
   We have to reproduce the issue, so we need all the input files.
4. What would be the expected outcome?
5. What did you see instead?

All these details will help people to fix any potential bugs.


## Suggesting a New Feature

Feature requests are welcome. But take a moment to find out if your idea fits
the scope and goals of the project. It is up to you to provide a strong
argument to convince the project's developers of the benefits of this feature.
Please provide as much detail and context as possible.


## Implementing a New Feature

Contributions are welcome via Github pull requests.

- Each pull request should implement *one* feature or fix *one* bug.
  If you want to add or fix more than one thing, submit more than one
  pull request.
- Do not commit changes to files that are irrelevant to your feature or
  bugfix (*e.g.* `.gitignore`).
- Be willing to accept criticism and work on improving your code.
- Make sure the code compiles and the tests run successful on more than
  your local machine (*e.g.* on the cluster of your university).
- Please sign-off your commits


### For New Contributors

If you never created a pull request before, welcome :tada:.
You can learn how from [this great tutorial](https://app.egghead.io/courses/how-to-contribute-to-an-open-source-project-on-github)

Don't know where to start?
You can start by looking through these [issues](https://github.com/awvwgk/tblitbe/issues).


## Sign Your Work

The sign-off is a simple line at the end of the explanation for a commit. All 
commits needs to be signed. Your signature certifies that you wrote the patch or
otherwise have the right to contribute the material. The rules are pretty simple,
if you can certify the below (from [developercertificate.org](https://developercertificate.org/)):

```
Developer Certificate of Origin
Version 1.1

Copyright (C) 2004, 2006 The Linux Foundation and its contributors.
1 Letterman Drive
Suite D4700
San Francisco, CA, 94129

Everyone is permitted to copy and distribute verbatim copies of this
license document, but changing it is not allowed.

Developer's Certificate of Origin 1.1

By making a contribution to this project, I certify that:

(a) The contribution was created in whole or in part by me and I
    have the right to submit it under the open source license
    indicated in the file; or

(b) The contribution is based upon previous work that, to the best
    of my knowledge, is covered under an appropriate open source
    license and I have the right under that license to submit that
    work with modifications, whether created in whole or in part
    by me, under the same open source license (unless I am
    permitted to submit under a different license), as indicated
    in the file; or

(c) The contribution was provided directly to me by some other
    person who certified (a), (b) or (c) and I have not modified
    it.

(d) I understand and agree that this project and the contribution
    are public and that a record of the contribution (including all
    personal information I submit with it, including my sign-off) is
    maintained indefinitely and may be redistributed consistent with
    this project or the open source license(s) involved.
```

Then you just add a line to every git commit message:

    Signed-off-by: Joe Smith <joe.smith@example.com>

Use your real name (sorry, no pseudonyms or anonymous contributions.)

If you set your `user.name` and `user.email` git configs, you can sign your
commit automatically with `git commit -s`.
