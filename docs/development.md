# Information for EXTPAR Developers

## Git and Github

The EXTPAR code is developed using the Git version control system and the [Github web interface :material-open-in-new:](https://github.com/C2SM/EXTPAR){:target="_blank"}. 
Outstanding bugs and requested features are tracked using the [Issues :material-open-in-new:](https://github.com/C2SM/EXTPAR/issues){:target="_blank"} section of the Github repository.  Additionally, automated testing of newly developed features is integrated into the Github interface using the Jenkins CI tool.  

### Main branches

The **master** branch only contains code which are released versions. 
All commits on the master branch are tagged (`git tag -a vX.Y.Z`).
Only the core development team is allowed to modify the master branch.

### Supporting branches

Any new code development should be done in a **topic** branch. Topic branches are merged
back into master by opening a pull request. Code must be peer reviewed by the
source code administrator.

Supporting branches are removed once successfully merged in one of the master branch.

### Developments with new input data sets

Any new EXTPAR code that is accompanied by a new input data file or files should be added with a simultaneous pull request
in both this code repository (for the code changes) and in the [extpar-input-data repository :material-open-in-new:](https://gitlab.dkrz.de/extpar-data/extpar-input-data){:target="_blank"}
(for the addition of the input data files). The **topic** branch for both pull requests should have the same name, in order to enable
the synchronization of the code and input-data repositories.

### Testing new developments

Once a developer has finished developing a new feature or bug fix, they should make a 
pull request on the Github repository from their topic branch into the rc-branch.  
Then, they should write the following comment into the pull request conversation: 

```
launch jenkins
```

This will start the automated testing, and the code will be compiled and tested on co2 (ETH) and Levante (DKRZ).

If the tests fail, then the developer should fix the issues and resubmit the testing on Jenkins.  
Once all of the tests are passing, then they should notify the source code administrator that the pull
request is ready for review and merging into the rc-branch.  

## Fortran Code

### Logging

In case you want to add some additional prints in EXTPAR, please use the logger described below.

CALL the built-in logger-functions in order to print messages or variables in the specific logfile of each EXTPAR executable.
The logger has three different levels of messages to print:

1. `logging%info(your_message)`: info-prints for better orientation during code execution, variables or other stuff.

2. `logging%warning(your_message)`: warnings, like wrong namelist-inputs, unsupported NetCDF versions or problems with some data points.

3. `logging%error(your_message, __FILE__,__LINE__)`: errors that occur during I/O, allocation, that requires an abort of EXTPAR. 

As your_message needs to be a sequence of characters, use

```fortran
WRITE(message_text,*)var_x, 'is now', var_y
```

and then 

```fortran
CALL logging%inf0(message_text) 
```

to print the values of `var_x` and `var_y` to the logfile.

For quick debugging-prints ONLY use

```fortran
WRITE(logging%fileunit,*)var_x, 'is now', var_y 
```

### Coding Rules and Best Practices

1. All features available in Fortran 2008 as far as supported by Intel,
GCC, and NAG are allowed.

2. Use up to the allowed 132 character per line, but not more. Note
that this includes comments.

3. Indentation rules:

    | Code feature  | Num. of indentation characters |
    | ------------- |-------------| 
    | program indentation      | 2 |
    | type definition          | 2 |
    | do loops                 | 2 |
    | if constructs            | 2 |
    | continuation             | 5 (with leading &) |
    | all directives           | 0 |

4. Always use `IMPLICIT NONE` and `PRIVATE`/`PUBLIC` once only in modules header.

5. Do not add `USE` statements after `CONTAINS`.

6. Fortran keywords should be in capital letters with the exception of `len`,
`in`, `out`, and `inout`.

7. Do not use tabs, deprecated, or obsolete features.

8. Do not overspecify declarations - especially if standard types are expected.

## Python Code

### Logging

In case you want to add some additional prints in EXTPAR, please use the logger described below.

CALL the built-in logger-functions in order to print messages or variables in the specific logfile of each EXTPAR executable.
The logger can print variables as well as strings. Use formatted strings (`f'`) in case you want to combine variables and strings.
The logger has 4 different levels of messages to print:

1. `logging.debug(your_message)`: Mean/Max/Min of variables needed for development, more detailed information about code execution.

2. `logging.info(your_message)`: info-prints for better orientation during code execution, variables or other stuff.

3. `logging.warning(your_message)`: warnings, like wrong namelist-inputs, unsupported NetCDF versions or problems with some data points.

4. `logging.error(your_message)`: errors that occur during I/O, allocation or wrong namelist parameters, that requires an abort of EXTPAR. The programm does not stop automatically after the call of logging.error, so a `raise` follows the logging.error()

Default logging level is info, so only messages from `logging.info()`, `logging.warning()` and `logging.error()` are written to the logfile. Adjust the level of the logger right at the beginning of each Python executable to `level=logging.DEBUG` to also print `logging.debug()`.

### Coding rules and best practices

The Python code needs to fulfill the [Pep8 coding standard :material-open-in-new:](https://www.python.org/dev/peps/pep-0008/){:target="_blank"}.
A GitHub action automatically formats Python code for you.
