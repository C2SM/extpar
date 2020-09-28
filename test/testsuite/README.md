  This README is for the Extpar testsuite, which is adapted from the COSMO technical testsuite.


  To run the Extpar testsuite:

  1.  Compile the Extpar code following the instructions in the main directory.

  2.  Download the testsuite data by running the get_data.sh script in the data folder:
      cd data; ./get_data.sh; cd ..

  3.  Copy the Extpar executables to the /bin folder: cp ../../bin/* bin

  4.  If you want to run on a new machine, add all host-dependent paths and variables to the runscripts for ICON and COSMO in the /bin folder.

  5.  To run the testsuite, either make a submit script that calls the testsuite 
      (like the submit.tsa.sh script) or call the testsuite directly from the command line:
      ./src/testsuite.py -v 1 --exe=your_run_script -o testsuite.out
      This will generate a testsuite.out file which contains the results of the testsuite.  
      More information about the testsuite can be found in the src/doc folder.  
      The testsuite script options can be found with the help options:
      ./src/testsuite.py -h

  To add a test to the Extpar testsuite:

  1. Modify the testlist.xml file to add the new test.  Alternatively, you could also add a new 
     testlist XML file (with a new name).  The testlist which is run can be chosen from the testsuite
     command line.  

  2. Make a folder in the data folder for the new test containing the INPUT_* files and the namelist.py for the Python-CDO modules.  

  3. Send me (katherine.osterried@env.ethz.ch) any binary reference files to upload to the ftp site. 
     I will also modify the get_data.sh script accordingly.  


