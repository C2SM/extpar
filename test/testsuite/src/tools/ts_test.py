#!/usr/bin/env python

"""
COSMO TECHNICAL TESTSUITE

This module implements the Test baseclass which can be used to setup, run and
evaluate a single test.
"""

# built-in modules
import os, sys, copy, math, re, glob

# private modules
from ts_error import StopError, SkipError
from ts_utilities import dir_path, status_str, pretty_status_str, system_command, change_dir, write_environ
from ts_fortran_nl import get_param, replace_param

# information
__author__     = "Nicolo Lardelli, Xavier Lapillonne, Oliver Fuhrer, Santiago Moreno"
__email__      = "cosmo-wg6@cosmo.org"
__maintainer__ = "xavier.lapillonne@meteoswiss.ch"


class Test:
    """Class representing a test and allows setting up, running and evaluating a test"""

    def __init__(self, node, options, conf, logger):
        # store private information for this test
        self.node = node                # xml node
        self.name = node.attrib['name'] # assign the name to the test
        self.type = node.attrib['type'] # set test type
        self.description = node.findtext('description') # storage of the description string
        self.options = copy.copy(options) # test options
        self.conf = copy.copy(conf) # storage of the auxiliary parameters
        self.logger = logger            # store logger
        self.result = 30                # default to CRASH

        # define prerun actions
        if node.findtext('prerun'):
          self.prerun_actions = node.findtext('prerun').split(',')
        else:
          self.prerun_actions = []

        # setup of directory paths
        self.basedir = dir_path(conf.basedir)
        self.inputdir = dir_path(self.basedir + 'data/' + self.type) # set path for input directory

        self.rundir = dir_path(self.basedir) + dir_path(options.workdir) + dir_path(self.type) + dir_path(self.name)
        if node.findtext('namelistdir'):
            self.namelistdir = dir_path(self.basedir + 'data/' + node.findtext('namelistdir'))
        else:
            #default is  self.type/self.name
            self.namelistdir = dir_path(self.basedir + 'data/' + self.type + '/' + self.name)


        if node.findtext('refoutdir'):
            # set reference output directory, two possibility
            # 1. relative to current path, if name starts with ../
            # 2. relative to the main directory is starts with a string
            refoutdir = dir_path(node.findtext('refoutdir'))
            pattern = '[.][.][/](.*)'
            matchobj = re.match(pattern,refoutdir)
            if matchobj:
                self.refoutdir = self.rundir+refoutdir
            else:
                self.refoutdir = self.basedir+'data/'+refoutdir
        else:
            #default is  self.namelistdir
            self.refoutdir = self.namelistdir

        # set dependecy directory
        depend = self.node.findtext("depend")
        if depend != None:
            depend = depend.strip(' ')
            pattern = '[.][.][/].*'
            matchobj = re.match(pattern,depend)
            # relative to rundir
            if matchobj:
                self.dependdir = depend
            else:
                self.dependdir = self.basedir + '/' + depend
        else:
            self.dependdir = None

        # set executable
        if self.options.exe is not None:
            self.executable = self.options.exe
        elif node.findtext("executable"):
            self.executable = self.node.findtext("executable")
        else:
            raise SkipError('An executable must be defined in the command line or in testlist.xml')

        # overide nprocs if define in xml
        if node.findtext("nprocs"):
            self.nprocs=int(self.node.findtext("nprocs"))
        else:
            self.nprocs=self.options.nprocs

        # set tolerance folder name (used by tolerance checker)
        self.tolerance=self.options.tolerance


    def run_test(self):
        """ check whether this test should be carried out in case "only" option is used"""
        if self.options.only is not None:
            if self.options.only=='%s,%s' %(self.type,self.name):
                return True
            else:
                return False
        else:
            return True

    def __prepare_print(self):
        """print infos about the upcoming test at the start"""
        self.logger.info('')
        self.logger.important('TEST {0}/{1}: {2}'.format(self.type, self.name, self.description))

    def prepare(self):
        """prepare test directory and namelists for this test"""

        self.__setup_directory()

        self.__setup_executable()

        self.__adapt_namelists()

        if 'extpar' not in self.executable:

            self.__set_parallelization()

            self.__set_timesteps()

            self.__set_pert()

        self.__prepare_print()


    def prerun(self):
        """check dependencies and perform any prerun actions"""

        # change to run directory
        status = change_dir(self.rundir, self.logger)

        # check whether dependencies have run
        if self.dependdir != None:
            try:
                f = open(self.dependdir + '/' + self.conf.res_file, "r")
                dresult = f.readline()

                if dresult == 'CRASH':
                    raise SkipError('Required test %s has crashed' %(self.dependdir))
               # if dresult == 'SKIP':
               #     raise SkipError('Required test %s has been skipped' %(self.dependdir))
            except IOError:
                raise SkipError('No result file %s for required test %s' \
                                %(self.conf.res_file,self.dependdir))

        # defines the procedures in case of restart test
        if 'restart' in self.prerun_actions:

            # test with restart are not allowed with a too short number of timesteps
            if self.options.steps is not None and self.options.steps<60:
                raise SkipError('Restart is not compatible with short tests')

            # copy restart file
            status = system_command('/bin/cp '+self.dependdir+'/output/lr* '+self.rundir+'output/', self.logger, throw_exception=False)
            if status:
                raise SkipError('Problem with restart file from '+self.dependdir)


    def start(self):
        """launch test"""

        self.logger.info('Starting test')

        # change to run directory for this test
        status = change_dir(self.rundir, self.logger)

        # generate launch command
        self.log_file = 'exe.log'
        redirect_output = '> %s 2>&1' %(self.log_file)

        if self.options.mpicmd == '':
            run_cmd = ''
        else:
            if '&NTASKS' in self.options.mpicmd:
                #special case when n nodes cannot be given as last argument of
                #mpicmd command, e.g. with mpirun_rsh
                run_cmd = self.options.mpicmd.replace('&NTASKS',str(self.nprocs))
            else:
                run_cmd = self.options.mpicmd + ' ' + str(self.nprocs)

        # writes the wrapper script in case a wrapper run of testsuite is required
        if self.options.use_wrappers:
            f = open('wrapper.sh','w')
            f.write('#!/bin/sh\n')
            f.write('./'+self.executable+redirect_output+'\n')
            f.close()
            status = os.chmod('wrapper.sh',0o755)
            if status:
                raise StopError('Problem changing permissions on wrapper.sh')
            run_cmd = run_cmd + ' ./' + 'wrapper.sh'
        else:
            run_cmd=run_cmd + ' ./' + self.executable + ' ' + redirect_output
        
        if 'extpar' in self.executable:
            run_cmd='./'+self.executable

        # displays the run command
        self.logger.info('Executing: '+run_cmd)

        # executes the run command
        status = system_command(run_cmd, self.logger, issue_error=False, timeout=self.options.timeout)


    def wait(self):
        """wait for completion of test"""

        # currently assumes that tests are run sequentially and all work is done in the start() method
        self.logger.info('Test finished')


    def check(self):
        """perform checks"""

        # change to base directory
        status = change_dir(self.basedir, self.logger)

        # scan for the checker within the xml tree
        checkerlist = []
        checker_nodes = self.node.findall("checker")
        for el in checker_nodes:
            checkerlist.append(el.text)   

        # assignement of the environment variables for checker
        write_environ(self)

        # traversing of the checkerlist
        summary_list = []
        for checker in checkerlist:

            self.logger.debug(checker+' START')

            # run checker and save result
            checker_result,soutput = system_command(os.path.join(os.path.dirname(__file__), "../checkers/")+checker, self.logger, \
                                                      return_output=True,throw_exception=False, \
                                                      issue_error=False)
            # print checker output
            for line in soutput.split('\n'):
                if not line=='':
                    self.logger.chckinfo(line)


            summary_list.append(checker_result)
                
            # display the subsummary for the checkers
            self.logger.result(1, checker_result, checker)
            
            self.logger.debug(checker+' END')

        # compute summary result    
        self.result = max(summary_list)

        # in case of fail or crash signal a stop to the testsuite
        if self.result >= 20:
            raise StopError


    def write_result(self):
        """print result of current test to stdout as well as result file"""

        # print the final result 
        self.logger.result(0, self.result, 'RESULT %s/%s: %s' %(self.type,self.name,self.description))

        # change to run directory for this test
        status = change_dir(self.rundir, self.logger)

        # write in a file (this is used for test dependency)
        f = open(self.conf.res_file, "w")
        f.write(status_str(self.result))
        f.close()
       
       

    def update_namelist(self):
        """copy back modified namelists into data folder"""

        pattern = '(.*)/'+self.type+'/'+self.name+'(.*)'
        text = self.namelistdir

        # checks if is the test is a titular test
        if re.match(pattern,text):
            status = change_dir(self.rundir, self.logger)
            self.logger.important('Updating namelist data/'+self.type+'/'+self.name)
            cmd = 'cp INPUT* '+self.namelistdir
            self.logger.debug('Executing: '+cmd)
            status = system_command(cmd, self.logger)
            self.result = 0 # MATCH
        else:
            raise SkipError('No test repository ' +'data/'+self.type+'/'+self.name)

    def update_yufiles(self):
        """copy back YU* files into the data folder
        Should be used to generate new reference files"""
        

        pattern = '(.*)/'+self.type+'/'+self.name+'(.*)'
        # only update those test where self.refoutdir=self.namelistdir
        text = self.namelistdir

        # checks if is the test is a titular test
        if re.match(pattern,text):
            status = change_dir(self.rundir, self.logger)
            if not os.path.exists(self.conf.yufile):
                raise SkipError('No file ' +self.conf.yufile+' in '+self.rundir)

            self.logger.info('Updating exe.log YU* ' + self.namelistdir)
            cmd = 'cp exe.log YU* '+self.namelistdir
            self.logger.debug('Executing: '+cmd)
            status = system_command(cmd, self.logger)
            self.result = 0 # MATCH
        else:
            raise SkipError('No test repository ' +'data/'+self.type+'/'+self.name)



    def __setup_directory(self):
        """generate test directory including all required links and sub-directories"""

        self.logger.info('Creating directory for '+self.name)

        node = self.node

        # create run directory and move there
        status = system_command('/bin/mkdir -p '+self.rundir, self.logger)
        status = change_dir(self.rundir, self.logger)

        # removal of all the possible pre-existing files
        status = system_command('/bin/rm -r -f *', self.logger)
        
        # explicit copy of the namelists (copy is required since we will apply the change_par)
        status = system_command('/bin/cp -f '+self.namelistdir+'INPUT* .', self.logger)

        # explicit copy of the namelists for the Python modules of Extpar
        status = system_command('/bin/cp -f '+self.namelistdir+'*.py .', self.logger)

        # copy the icon grid if exists
        if not glob.glob(dir_path(self.namelistdir) + 'icon_grid_*')==[]:
            status = system_command('/bin/cp -f '+self.namelistdir+'icon_grid_* ./', self.logger)

        # copy of the climatological fields needed for the icon test if exists
        if not glob.glob(dir_path(self.namelistdir) + 'clim_*')==[]:
            status = system_command('/bin/cp -f '+self.namelistdir+'clim_* ./', self.logger)

        # copy of the auxiliary input parameters if exists
        if not glob.glob(os.path.join(dir_path(self.inputdir)+'in_aux/', '*'))==[]:
            status = system_command('/bin/cp -f -r '+dir_path(self.inputdir)+'in_aux/* ./', self.logger)

        # linking input binary fields
        status = system_command('/bin/ln -s '+dir_path(self.inputdir)+'input .', self.logger)
        # generation of the output folder
        status = system_command('/bin/mkdir -p output', self.logger)


    def __setup_executable(self):

        # choose the executable
        if self.executable == None:
            raise SkipError('No executable defined in argument list or xml file')
       
        self.logger.info('Fetching executable '+self.basedir+self.executable)

        # copy of the executable
        if 'extpar' in self.executable:
            if not os.path.exists(self.basedir+'bin/'+self.executable):
               raise SkipError('Executable '+self.basedir+'bin/'+self.executable+' does not exist')
        else:
            if not os.path.exists(self.basedir+self.executable):
               raise SkipError('Executable '+self.basedir+self.executable+' does not exist')
        if 'extpar' in self.executable:
            status = system_command('/bin/cp '+self.basedir+'bin/* .', self.logger)
        else:
            status = system_command('/bin/cp '+self.basedir+self.executable+' .', self.logger)
        

    def __adapt_namelists(self):

        self.logger.info('Modify namelists (according to XML specification)')
        
        # change_par from the xml file
        changepar_list = self.node.findall("changepar")

        # iterate over all changepars
        for chpar in changepar_list:
            filename = chpar.attrib['file']
            if filename is None:
                self.logger.error('changepar encountered without file attribute')
                continue

            filename = str(filename)
            newparname = str(chpar.attrib['name'])
            
            # look if optional attribute occurrence exists
            try:
                occurrence = int(chpar.attrib['occurrence'])
            except:
                occurrence = 1

            # check if parameter is one of the dual parameters
            parname = newparname
            for (param1, param2) in self.conf.dual_params:
                if param1 == parname:
                    if get_param(filename,parname,occurrence=occurrence) == '':
                        parname = param2
                if param2 == parname:
                    if get_param(filename,parname,occurrence=occurrence) == '':     
                        parname = param1

            value = chpar.text
            modstring = newparname + '=' + str(value)
            replace_param(filename, parname, modstring, occurrence=occurrence)

            # if nprocio has been overwritten, remove configuration (XML should have precedence)
            if parname == 'nprocio':
                self.options.nprocio = None
    
    def __set_pert(self):
        """set perturbation in the parameter file to true or false depending on the given option"""
        if self.conf.pert_avail == 'True':
             pert = self.options.pert
             replace_param(self.conf.par_file,'itype_pert',' itype_pert = %i' %pert)

    def __set_parallelization(self):

        self.logger.info('Set domain decomposition and number of I/O PEs')

        ### extract number of I/O processors
        if self.options.nprocio is not None:
            nprocio = self.options.nprocio
        else:
            if get_param(self.conf.par_file,'num_iope_percomm') != '':
               num_iope_percomm = int(get_param(self.conf.par_file,'num_iope_percomm'))
               num_asynio_comm = int(get_param(self.conf.par_file,'num_asynio_comm'))
               nprocio = num_asynio_comm * num_iope_percomm
            else:
               nprocio = int(get_param(self.conf.par_file,'nprocio'))

        # sets the number of I/O processors
        if get_param(self.conf.par_file,'num_iope_percomm') != '':
           if nprocio == 0:
              replace_param(self.conf.par_file,'num_iope_percomm',' num_iope_percomm= 0')
           else:   
              replace_param(self.conf.par_file,'num_iope_percomm',' num_iope_percomm= 1')
           replace_param(self.conf.par_file,'num_asynio_comm',' num_asynio_comm= %i' %nprocio)
        else:
           replace_param(self.conf.par_file,'nprocio',' nprocio= %i' %nprocio)

        # generates the parallelist
        parlist = []
        parlist = self.set_parallelization(self.nprocs,nprocio)

        # select the parallelization
        ap = int(self.node.findtext("autoparallel"))
        if ap > len(parlist):
            raise SkipError('The selected autoparallel number is too large for the given number of processor (not enough decompositions available)')

        # writes the new MPI decomposition
        nprocx = parlist[ap-1][0]
        nprocy = parlist[ap-1][1]
        replace_param(self.conf.par_file, 'nprocx', ' nprocx= %i' %nprocx)
        replace_param(self.conf.par_file, 'nprocy', ' nprocy= %i' %nprocy)
                                 
        # echo to log
        self.logger.info('Processors distribution set to ' + 
            '(nprocx,nprocy,nprocio)=(%i,%i,%i)' %(nprocx, nprocy, nprocio))


    def __set_timesteps(self):
 
        if self.options.steps is not None:

            self.logger.info('Number of timesteps (nstop=%i)' %(self.options.steps))

            modstring = 'nstop=' + str(self.options.steps)
            parname = 'nstop'
            filename = self.conf.par_file
            if get_param(filename,'nstop') == '':
                parname = 'hstop'
            replace_param(filename, parname, modstring)


    @staticmethod
    def set_parallelization(nprocs,nprocio):
        """return a list of possible tuples of domain decompositions"""

        nxy = int(nprocs)-int(nprocio)
        if nxy < 1:
            raise ValueError('*** ERROR: The number of total processor'\
                                 ' is smaller equal the number of I\O processors')
        parlist = []

        for i in range(nxy,0,-1):
            if nxy%i == 0:
                parlist.append([i,nxy/i])
        
        # sort parlist by aspect ratio of solutions
        parlist = sorted(parlist, key=lambda tuple: aspect_ratio(tuple[0],tuple[1]))

        ## # swap parlist[0] and parlist[-2] so that nxy,1 is not the first parallelization
        ## if len(parlist)>1:
        ##     par_tmp = parlist[0]
        ##     parlist[0] = parlist[-2]
        ##     parlist[-2] = par_tmp

        return parlist    
    

def aspect_ratio(nprocx,nprocy):
    """compute aspect ratio of decomposition"""
    nprocx = abs(float(nprocx))
    nprocy = abs(float(nprocy))
    if (nprocx > nprocy):
        return nprocx/nprocy
    else:
        return nprocy/nprocx


