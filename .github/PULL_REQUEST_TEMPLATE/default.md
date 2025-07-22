## Workflow for merging PRs

Please read these instructions carefully and follow the steps below before 
requesting a review by the maintainers. This way we can ensure a smoother 
review process and your changes will be merged sooner.

Additionally, if this is the first PR you open in EXTPAR make sure to read the 
[*"Information for EXTPAR Developers"*](https://c2sm.github.io/extpar/development/) 
section in the documentation.

### Checklist

- [ ] Provide a detailed description of your changes.
- [ ] If you implemented a new feature:
  - [ ] Document it in the correct Markdown file(s) under the 
[`docs/`](https://github.com/C2SM/extpar/tree/master/docs) directory.
  - [ ] [Add a new test](https://c2sm.github.io/extpar/testing/#add-a-new-test) 
or make sure your changes are already tested.
- [ ] Your code follows the [style guidelines](https://c2sm.github.io/extpar/development/#coding-rules-and-best-practices).
- [ ] Your changes only touch the files/lines relevant for you.
- [ ] All four required checks pass (see 
[*"Testing and debugging"*](#testing-and-debugging) for more details).
- [ ] No conflicts with the base branch.

If all the points above are satisfied you can ask for a review by selecting 
`stelliom` as a reviewer.

For any questions please ping `@stelliom` on this PR.

### Testing and debugging

The most important test for PRs is the one labeled 
*"EXTPAR Testsuite on Jenkins"*. This checks that the results of all testcases 
(described by the namelists in 
[`test/testsuite/data`](https://github.com/C2SM/extpar/tree/master/test/testsuite/data)) 
did not change compared to the references.

You can launch the testsuite by writing `launch jenkins` as a comment in the 
PR that you want to test. Once completed, the result of the testsuite will be 
shown on the PR (failure or success).

If you need more details on the testsuite results (e.g., if you are trying to 
debug an error or you are simply unsure why the tests fail) you can launch the 
testsuite with `launch jenkins(debug)`. This will run the tests as usual, but, 
once completed, you will be given a URL (via a comment on the PR) to access the 
logfiles and namelists of all tests that were run.