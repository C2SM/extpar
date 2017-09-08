#!/usr/bin/env python2

"""
COSMO TECHNICAL TESTSUITE

This script runs a set of tests defined in testlist.xml and checks the results
for correctness using a set of checkers which can be defined for each test.

For further help and command line arguments execute this script without arguments.
"""

# built-in modules
import os, sys, string, struct
import optparse as OP
import xml.etree.ElementTree as XML
import logging as LG
import ConfigParser
import ast

# private modules
sys.path.append(os.path.join(os.path.dirname(__file__), "./tools")) # this is the generic folder for subroutines
from ts_error import StopError, SkipError
from ts_utilities import system_command, change_dir
import ts_logger as LG
from ts_test import Test
from default_values import DefaultValues

# information
__author__     = "Nicolo Lardelli, Xavier Lapillonne, Oliver Fuhrer, Santiago Moreno"
__copyright__  = "Copyright 2012-2017, COSMO Consortium"
__license__    = "MIT"
__version__    = "2.1"
__date__       = "May 15, 2017"
__email__      = "cosmo-wg6@cosmo.org"
__maintainer__ = "xavier.lapillonne@meteoswiss.ch"


def parse_config_file(filename):
    config = ConfigParser.RawConfigParser()

    # create empty conf object
    conf = type('configuration', (), {})()

    # save base directory where testsuite is executed
    conf.basedir = os.getcwd()

    try:
        config.read(filename)
        conf.l_files     = ast.literal_eval(config.get('ts_config','l_files'))
        conf.nl_ts_switch= config.get('ts_config','nl_ts_switch')
        conf.config_nl   = config.get('ts_config','config_nl')
        conf.par_file    = config.get('ts_config','par_file')
        conf.dt_file     = config.get('ts_config','dt_file')
        conf.res_file    = config.get('ts_config','res_file')
        conf.pert_avail  = config.get('ts_config','pert_avail')
        conf.yufile      = config.get('ts_config','yufile')
        conf.dual_params = ast.literal_eval(config.get('ts_config','dual_params'))

    except Exception as e:
        print('Error while reading config file '+filename+':')
        print(e)
        raise # this exits with full traceback
        
    return conf 


def parse_cmdline():
    """parse command line options"""

    # the parser is initialized with its description and its epilog
    parser = OP.OptionParser(description=
                "Description: this script run a series of tests defined in testlist.xml. For "+
                "each test a set of checks are carried out, see checkers/README for more "+
                "information",
              epilog=
                "Example: ./testsuite.py -n 16 --color -f --exe=cosmo --mpicmd=\'aprun -n\' -v 1")
    
    # defines the number of processor, the number of processors is not a test specific option,
    # it overrides the values present in the namelist INPUT_ORG
    parser.set_defaults(nprocs=DefaultValues.nprocs)
    parser.add_option("-n",type="int",dest="nprocs",
               help=("number of processors (nprocx*nprocy+nprocio) to use [default=%d]" % DefaultValues.nprocs))
    
    # defines the number of I/O processors, not a test specific option
    parser.set_defaults(nprocio=DefaultValues.nprocio)
    parser.add_option("--nprocio",type="int",dest="nprocio",
               help="set number of asynchronous IO processor, [default=<from namelist>]")
 
    # defines the behavior of testsuite after fail or crash
    parser.add_option("-f","--force",action="store_true",dest="force",default=False,
               help="do not stop upon error")
    
    # set the level of verbosity of the standard output
    parser.set_defaults(v_level=DefaultValues.v_level)
    parser.add_option("-v",type="int",dest="v_level",help=("verbosity level 0 to 3 [default=%d]" % DefaultValues.v_level))

    # specifies the syntax of the mpi command
    parser.set_defaults(mpicmd=DefaultValues.mpicmd)
    parser.add_option("--mpicmd",dest="mpicmd",type="string",
               help=("MPI run command (e.g. \"mpirun -n\") [default=\"%s\"]" % DefaultValues.mpicmd))

    # defines the executable name, this overides the definition in testlist.xml if any
    parser.set_defaults(exe=DefaultValues.exe)
    parser.add_option("--exe",dest="exe",type="string",
               help="Executable file, [default=<from testlist.xml>]")

    # defines if the output will be colored or not
    parser.add_option("--color",action="store_true",dest="color",default=False,
               help="Select colored output [default=False]")

    # overrides the hstop/nstop options in the namelist and execute the cosmo runs with a given number of steps
    parser.set_defaults(steps=DefaultValues.steps)
    parser.add_option("--steps",dest="steps",type="int",action="store",
               help="Run only specified number of timesteps [default=<from namelist>]")

    # defines if a wrapper for job submission has to be written or not
    parser.add_option("-w","--wrapper",action="store_true",dest="use_wrappers",default=False,
               help="Use wrapper instead of executable for mpicmd [default=False]")

    # defines the filename for the redirected standard output
    parser.set_defaults(stdout=DefaultValues.stdout)
    parser.add_option("-o",dest="stdout",type="string",action="store",
               help="Redirect standard output to selected file [default=<stdout>]")
    
    # defines the behaviour of the redirected standard output, if appended or overwritten
    parser.add_option("-a","--append",action="store_true",default=False,dest="outappend",
               help="Appends standard output if redirection selected [default=False]")
    
    # only one test is executed
    parser.add_option("--only",dest="only",type="string",action="store",
               help="Run only one test define as type,name (e.g. --only=cosmo7,test_1)")

    # update namelist (no run). This is useful to quickly change all namelist at once 
    parser.add_option("--update-namelist",dest="upnamelist",action="store_true",default=False,
               help="Use Testsuite to update namelists (no tests executed)")

    # force bit-reproducible results
    parser.add_option("--force-match",dest="forcematch",action="store_true",default=(DefaultValues.forcematch == 1),
               help="Force bit-reproducible results")

    parser.add_option("--tune-thresholds",dest="tune_thresholds",action="store_true",default=DefaultValues.tune_thresholds,
               help="Change thresholds to always at least return OK")

    parser.add_option("--update-thresholds",dest="update_thresholds",action="store_true",default=DefaultValues.update_thresholds,
               help="Update the thresholds")

    parser.add_option("--tuning-iterations",dest="tuning_iterations",action="store",default=DefaultValues.tuning_iterations,
               help="Defines how many times the tuning gets executed")

    parser.add_option("--reset-thresholds",dest="reset_thresholds",action="store_true",default=DefaultValues.reset_thresholds,
               help="Set all thresholds to 0.0 before tuning")

    # update namelist (no run). This is useful to quickly change all namelist at once 
    parser.add_option("--update-yufiles",dest="upyufiles",action="store_true",default=False,
               help="Define new references (no tests executed)")

    # specifies the namelist file
    parser.set_defaults(testlist=DefaultValues.testlist)
    parser.add_option("-l","--testlist",dest="testlist",type="string",action="store",default=DefaultValues.testlist,
               help=("Select the testlist file [default=%s]" % DefaultValues.testlist))

    # timeout value for individual tests
    parser.set_defaults(timeout=DefaultValues.timeout)
    parser.add_option("-t","--timeout",dest="timeout",type="int",action="store",default=None,
               help=("Timeout in s for each test [default=%s]" % DefaultValues.timeout))

    # working directory
    parser.add_option("--workdir",dest="workdir",type="string",action="store",default="./work",
               help="Working directory [default=./work]")

    # specifies the tolerance file name for the tolerance checker
    parser.set_defaults(tolerance=DefaultValues.tolerance)
    parser.add_option("--tolerance",dest="tolerance",type="string",action="store",default=DefaultValues.tolerance,
               help=("Select the tolerance file name [default=%s]" % DefaultValues.tolerance))

    # parse
    try:
        (options,args)=parser.parse_args()
    except (OP.OptionError,TypeError):
        sys.exit("problem parsing command line arguments (check ./testsuite.py -h for valid arguments)")

    # convert forcematch options from logical to integer
    if options.forcematch:
        options.forcematch = 1
    else:
        options.forcematch = 0

    return options


def parse_xmlfile(filename, logger):

    try: 
        xmltree = XML.parse(filename)
    except Exception as e:
        logger.error('Error while reading xml file '+filename+':')
        logger.error(e)
        sys.exit(1) # this exits without traceback
        #raise # this exits with full traceback

    return xmltree.getroot()


def setup_logger(options):

    # instantiate logger class
    logger = LG.Logger(options.stdout, options.outappend, options.color)

    # set verbosity level
    if options.v_level <= 0:
      logger.setLevel(LG.ERROR)
    elif options.v_level == 1:
      logger.setLevel(LG.WARNING)
    elif options.v_level == 2:
      logger.setLevel(LG.INFO)
    elif options.v_level >= 3:
      logger.setLevel(LG.DEBUG)

    return logger
    

def main():
    """read configuration and then execute tests"""

    # definition of structure carrying global configuration
    # search for config file in current path, otherwise takes
    # default configuration file in testsuite source directory
    if os.path.isfile("./testsuite_config.cfg"):
        config_file = "./testsuite_config.cfg"
    elif os.path.isfile(os.path.join(os.path.dirname(__file__),"./testsuite_config.cfg")):
        config_file = os.path.join(os.path.dirname(__file__),"./testsuite_config.cfg")
    else:
        #logger not initialize at this stage, use print and exit
        print("Error: Missing configuration file testsuite_config.cfg")
        sys.exit(1)

    conf = parse_config_file(config_file)
    
    # parse command line arguments
    options = parse_cmdline()
        
    # redirect standard output (if required)
    logger = setup_logger(options)

    # hello world!
    logger.important('TESTSUITE '+__version__)

    # parse the .xml file which contains the test definitions
    logger.info('Parsing XML ('+options.testlist+')')
    root = parse_xmlfile(options.testlist, logger)

    # generate work directory
    status = system_command('/bin/mkdir -p '+options.workdir+'/', logger, throw_exception=False)
    if status:
      exit(status)

    # loops over all the tests
    stop = False
    for child in root.findall("test"):

        # create test object
        mytest = Test(child, options, conf, logger)
        
        if mytest.run_test():
            # run test
            try:

                # if upyufiles=True, no model run.
                if options.upyufiles:
                    logger.important('Update YU* files mode, no run')
                    mytest.update_yufiles()
                #
                elif options.update_thresholds:
                    logger.important('Updating the thresholds on the current runs')
                    mytest.options.tune_thresholds = True
                    mytest.log_file = 'exe.log'
                    mytest.check()
                # if upnamelist=True, no model run.
                elif options.upnamelist:
                    logger.important('Update namelist mode, no run')
                    mytest.prepare() # prepare test directory and update namelists
                    mytest.update_namelist() #copy back namelist in typedir
                else:
                    if(mytest.options.tune_thresholds):
                        mytest.options.pert = 0
                        for i in range(int(mytest.options.tuning_iterations)):
                            mytest.prepare() # prepare test directory and update namelists
                            logger.important("Iteration number {0}".format(i+1))
                            mytest.prerun() # last preparations (dependencies must have finished)
                            mytest.start()  # start test
                            mytest.wait()   # wait for completion of test
                            mytest.check()  # call checkers for this test
                            mytest.options.reset_thresholds = False
                            # 1: Perturb only in the first timestep
                            # 2: Perturb in every iteration
                            mytest.options.pert = 2
                    else:
                        mytest.options.pert = 0
                        mytest.prepare() # prepare test directory and update namelists
                        mytest.prerun() # last preparations (dependencies must have finished)
                        mytest.start()  # start test
                        mytest.wait()   # wait for completion of test
                        mytest.check()  # call checkers for this test

            except SkipError as smessage:
                mytest.result = 15 # SKIP
                logger.warning(smessage)

            except StopError as emessage:
                if str(emessage).strip():
                    logger.error(emessage)
                    if not options.force:
                        stop = True

            # write result
            mytest.write_result()

            # return into the base directory after each test
            status = change_dir(conf.basedir, logger)

            # exit if required
            if stop:
                break

    # end of testsuite std output
    logger.important('FINISHED')


if __name__ == "__main__":
    main()

