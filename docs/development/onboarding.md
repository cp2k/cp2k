# Onboarding

CP2K invites the community to contribute to its development! Documentation improvements, bug fixes,
performance enhancements, portability issues, new features, new methods ... you are encouraged to
contribute so that the community as a whole can benefit.

CP2K is a large project, there is no way to study the full code, and start when that is done! Start
working on small patches first, that are easy to code, to test and to integrate. Small patches are
easier to review and thus will be more quickly merged to the Git master branch. As experience with
the code grows, tackling larger projects becomes realistic. Developers who contribute regularly can
join the [developers team](https://github.com/orgs/cp2k/teams/cp2k-developers) on Github.

## Prepare the patch

- Fork the CP2K Repository on GitHub via "Fork" on https://github.com/cp2k/cp2k
  - Even as a member of the development team you will not be able to push directly to the `master`
    branch of the https://github.com/cp2k/cp2k repository. This is by intention and we would like to
    ask you to send changes as Pull Requests instead.
- On your machine, clone the repository via
  `git clone --recursive https://github.com/YOURNAME/cp2k.git`
- Change into the created cp2k directory: `cd cp2k`
- Start working in a new branch: `git checkout -b my-new-feature`
- Make your changes to the code
- Add new and changed files: `git add ...`
- Commit the changes: `git commit`
  - Please follow the [guidelines](https://github.com/cp2k/cp2k/wiki/Contribution-Guidelines)
- Do the first push to your fork, give the remote branch the same name as the local one:
  `git push -u origin my-new-feature`
- Do some more work, then repeat point 6. and 7.
- Push your new changes to the remote repository via `git push` (note: the
  `-u origin my-new-feature` does not have to be repeated)
- Use the GitHub interface at https://github.com/YOURNAME/cp2k to create a pull-request

## Update your copy of the master/rebase your patch

- To update the 'master' of your fork to the same state as the 'master' of the CP2K repository:
  - tell your local git repository once about the remote:
    `git remote add upstream https://github.com/cp2k/cp2k.git`
  - make sure you are on the right branch: `git checkout master`
  - `git rebase` your current branch on top of the cp2k/cp2k master:
    `git pull --rebase upstream master`
- To update a branch with patches onto the updated master branch:
  - `git checkout my-new-feature`
  - `git rebase master`, this may generate rebase/merge-conflicts you should resolve now. If you got
    lost, you can always use `git rebase --abort` (with Git 1.8 and newer) to revert the attempted
    rebase
  - after a rebase of a branch with commits which was already pushed to a remote, you have to
    force-push: `git push --force`

## Submit the patch

Following these guidelines will avoid common mistakes and make it easier to integrate patches. It
usually takes less than one hour:

- `./make_pretty.sh` to auto-format the code (variables).
- Compile the code (see [](../getting-started/build-from-source) or
  [](../getting-started/build-with-spack)).
- Prepare and add testcases suitable for regtesting the code (i.e. that run quickly through all new
  code paths)
- Run these testcases by hand. Be sure to use bounds checking, and valgrind to check for undefined
  variables or memory leaks (see debug page)
- Run the full regression testing suite, be sure the new testcases run correctly, and no other tests
  are incorrectly affected.
- pay attention to:
  - Are new files, including tests, visible ? Use `git add` as needed.
  - Are new test directories registered in `tests/TEST_DIRS`?
  - Are only the intended files and code modified ? Use `git revert` as needed
  - Does the code contain stray write statements or debug info ?
  - Is all new code sufficiently documented and explained ?
  - Are the input keywords clearly described ?
  - Are proper citations added to the bibliography ?
- If any of the above steps required changes to the code, go back to the first step, otherwise go to
  the next step.
- Go to https://github.com/YOURNAME/cp2k, GitHub will usually notify you directly that you have a
  branch from which you could create a Pull Request (PR) and since it is a fork of another
  repository, it will suggest to make the PR against that original repository (for instance the CP2K
  master branch)
  - You can check the status of your PR at
    [https://github.com/cp2k/cp2k/pulls](https://github.com/cp2k/cp2k/pulls)
  - We use `git rebase` to keep the
    [Git history linear](https://github.com/cp2k/cp2k/wiki/CP2K-CI#git-history), meaning that in
    your final Pull Request there can't be any merge commits
- The PR will trigger the
  [CP2K Continuous Integration (CI) system](https://github.com/cp2k/cp2k/wiki/CP2K-CI) to check
  conventions, compatibility, and correctness
  - In the case of success, the PR will be merged by one of the CP2K administrators
  - In the case of error, please check what's wrong in the CI logs, fix it in your branch, and
    commit/push again. The CI will automatically rerun on the new version of the code (no need to
    close the PR, it will be automatically updated!). Repeat until the CI runs fine.
  - The administrators may start additional tests if necessary.
- You will be notified if any of the tests on our [Dashboard](https://dashboard.cp2k.org/index.html)
  breaks. Submit a new PR to fix a bug found revealed by the dashboard caused by your PR. Ask the
  admin team for help if necessary.
