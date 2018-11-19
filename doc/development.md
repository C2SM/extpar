# Information for EXTPAR developers

## Git and Github
The Extpar code is developed using the git version control system and the Github web interface. 
Outstanding bugs and requested features are tracked using the Issues section of the Github repository.  Additionally, automated testing of newly developed features is integrated into the Github interface using the Jenkins CI tool.  

## Development workflow
The development policy is borrowed from the Fieldextra COSMO code repository 
maintained by Jean-Marie Bettems, and was inspired
by the document http://nvie.com/posts/a-successful-git-branching-model

### Main branches
The master repository holds two main branches with an infinite lifetime
* master 
* develop

The **master** branch only contains code which are released versions. 
All commits on the master branch are tagged (git tag -a vX.Y.Z).
Only the core development team is allowed to modify the master branch.

The **develop** branch is used to collect new features for the next release. All commits
on the develop branch should pass the tests of the technical testsuite. Only the core
development team is allowed to modify the develop branch.

### Supporting branches
Any new code development should be done in a **topic** branch. Topic branches are merged
back into develop branch by opening a pull request. Code must be peer reviewed by the
source code administrator.

A **release** branch supports the preparation of a new production release. It is branched
off from the develop branch and merged into the master branch. It is named
"release_vX.Y.Z", where vX.Y.Z is the name of the new release.

Supporting branches are removed once successfully merged in one of the main branch.

### Testing new developments
Once a developer has finished developing a new feature or bug fix, they should make a 
pull request on the Github repository from their topic branch into the develop branch.  
Then, they should write the following comment into the pull request conversation: "launch jenkins"
This will start the automated testing, and the code will be compiled and tested on Kesch and mistral.

If the tests fail, then the developer should fix the issues and resubmit the testing on Jenkins.  
Once all of the tests are passing, then they should notify the source code administrator that the pull
request is ready for review and merging into the develop branch.  

## Coding rules and best practices

1. All features available in Fortran 2008 as far as supported by Intel,
GCC, and NAG are allowed.

2. Use up to the allowed 132 character per line, but not more. Note
that this includes comments.

3. Indentation rules

| Code feature  | Num. of indentation characters |
| ------------- |-------------| 
| program indentation      | 2 |
| type definition          | 2 |
| do loops                 | 2 |
| if constructs            | 2 |
| continuation             | 5 (with leading &) |
| all directives           | 0 |

4. Always use IMPLICIT NONE and PRIVATE/PUBLIC once only in modules header.

5. Do not add add USE statements after CONTAINS.

6. Fortran keywords should be in capital letters with the exception of len,
in, out, and inout.

7. Do not use tabs, deprecated, or obsolete features.

8. Do not overspecify declarations - especially if standard types are expected.
