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

### Developments with new input data sets
Any new Extpar code that is accompanied by a new input data file or files should be added with a simultaneous pull request
in both this code repository (for the code changes) and in the [extpar-input-data repository](https://gitlab.dkrz.de/extpar-data/extpar-input-data)
(for the addition of the input data files). The **topic** branch for both pull requests should have the same name, in order to enable
the synchronization of the code and input-data repositories.

### Testing new developments
Once a developer has finished developing a new feature or bug fix, they should make a 
pull request on the Github repository from their topic branch into the develop branch.  
Then, they should write the following comment into the pull request conversation: "launch jenkins"
This will start the automated testing, and the code will be compiled and tested on Kesch and mistral.

If the tests fail, then the developer should fix the issues and resubmit the testing on Jenkins.  
Once all of the tests are passing, then they should notify the source code administrator that the pull
request is ready for review and merging into the develop branch.  

## Fortran Code

### Logging
In case you want to add some additional prints in Extpar, please use the logger described below.

CALL the built-in logger-functions in order to print messages or variables in the specific logfile of each Extpar executable.
The logger has three different levels of messages to print:

1. logging%info(your_message): info-prints for better orientation during code execution, variables or other stuff.

2. logging%warning(your_message): warnings, like wrong namelist-inputs, unsupported NetCDF versions or problems with some data points.

3. logging%error(your_message, __FILE__,__LINE__): errors that occur during I/O, allocation, that requires an abort of Extpar. 

As your_message needs to be a sequence of characters, use
WRITE(message_text,*)var_x, 'is now', var_y
and then CALL logging%inf0(message_text) to print the values of var_x and var_y to the logfile.

For quick debugging-prints ONLY use
WRITE(logging%fileunit,*)var_x, 'is now', var_y 

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

## Python Code

### Structure
The organization of the Python programmes is slightly different to the Fortran ones, because Python does not need to be compiled prior to execution.   
The main python scripts *extpar_alb_to_buffer.py*, *extpar_cru_to_buffer.py*, *extpar_ndvi_to_buffer.py* and *emiss_to_buffer.py* can be treated like the Fortran binaries and copied to the run-directory. Make sure the *namelist.py* is also present at the run-directory.  

All self-written Python-modules are stored in [lib](../python/lib) and do not need to be copied to the respective run-directory, rather the environment variable **PYTHONPATH** needs to be set the following:  

unset PYTHONPATH  
export PYTHONPATH=*absolute_path_to_python/lib*

### Logging
In case you want to add some additional prints in Extpar, please use the logger described below.

CALL the built-in logger-functions in order to print messages or variables in the specific logfile of each Extpar executable.
The logger can print variables as well as strings. Use formatted strings (f') in case you want to combine variables and strings.
The logger has 4 different levels of messages to print:

1. loggin.debug(your_message): Mean/Max/Min of variables needed for development, more detailed information about code execution.

2. logging.info(your_message): info-prints for better orientation during code execution, variables or other stuff.

3. logging.warning(your_message): warnings, like wrong namelist-inputs, unsupported NetCDF versions or problems with some data points.

4. logging.error(your_message): errors that occur during I/O, allocation or wrong namelist parameters, that requires an abort of Extpar. The programm does not stop automatically after the call of logging.error, so a sys.exit(1) follows the logging.error()

Default logging level is info, so only messages from logging.info(), logging.warning() and logging.error() are written to the logfile. Adjust the level of the logger right at the beginning of each Python executable to level=logging.DEBUG to also print logging.debug().

### Coding rules and best practices

The Python code needs to fulfill the [Pep8 coding standard](https://www.python.org/dev/peps/pep-0008/).
The testsuite provides a test to check these requirements. To check if your code aligns with the Pep8 coding standard
execute the [pep8_checker](../python/pep8_checker.sh).  
The most important coding rules and best practices are the following:

1. put (short) docstrings at the beginning of each function or class,  
also keep the content list in each module file up to date.

2. Check **ALL** namelist parameters from *namelist.py* for correctness before they are used in the code.

3. Limit the number of characters per line to 79.

4. Always use 4 indentation characters.
