#!/usr/bin/env python2

"""
COSMO TECHNICAL TESTSUITE

General purpose script to compare two YUCHKDAT output files
"""

# built-in modules
import os, sys, string, bisect

# information
__author__     = "Xavier Lapillonne"
__email__      = "cosmo-wg6@cosmo.org"
__maintainer__ = "xavier.lapillonne@meteoswiss.ch"


def cmp_(file1,file2, \
             v_level=0,minval=1e-15, \
             nts=[10,100,200], \
             tol_list=[1e-15,1e-15,1e-15]):

    # compare two YUCHKDAT file1, file2 with tolerance tol_*
    # Values smaller than minval are not considered.
    # Different threshold are used for different set depending on the time.
    # A set is identified by a line starting with the keyword "Check".
    # If the line contains the keyword "step", than the threshold for the corresponding step is used.
    # If the line does not contain the keyword "step", than the threshold for time 0 is used.
    # A valid line contains 10 (from 0 to 9) columns, with column 3,6,9 containing real numbers
    # v_level:verbose level
    #        -1 -> no print out
    #         0 -> max diff 
    #         1 -> show all lines with differences above tol_
    # if minval is -1 compares absolute differences
    # the comparison is only done for overlapping time steps


    # Valid line contains 10 (from 0 to 9) elements, with column 3,6,9 containing real numbers
    # see isValidLine method for details
    # Note : only numbers at position RealPos are compared
    ValidLineSize=10
    RealPos=[3,6,9]

    setStartKeyword='Check'
    setStepKeyword='step:'  #assumes that 'step:' is followed by a number

   # check file existence
    if not(os.path.exists(file1)):
        print 'File '+file1+' does not exist'
        return -1
    elif not(os.path.exists(file2)):
        print 'File '+file2+' does not exist'
        return -1

    
    # open file
    data1=open(file1).readlines()
    data2=open(file2).readlines()

    # variables initialisation
    error_count  = 0       #number of error detected
    nvalid_line = 0

    tol = tol_list[0]      # tolerence for t
   
    maxdiff=0.       #maximal diff
    maxdiff_line=0
    maxdiff_step=0

    print_header = True
    header = '  Errors above threshold :\n' + \
             '  var        ee    lev       min      imin   jmin         max      imax   jmax           mean          step       error'
    

    if v_level==0:
        if minval==-1:
            print_out='Absolute error:\n'
        else:
            print_out='Relative error:\n'

    if v_level>0:
        if minval==-1:
            print 'Comparing absolute differences ...'
        else:
            print 'Comparing relative differences, min. value is %1.0e ...' %(minval)
    
        
    # first on file1
    nline=min(len(data1),len(data2));

    for i in range(nline):
        l1=data1[i].split()
        l2=data2[i].split()

        #set theashold if new set
        if l1 and (setStartKeyword in l1[0]):
            #use step if defined
            if setStepKeyword in l1[-2]:
                try:
                    step=int(l1[-1])
                except ValueError:
                    print '!! Warning : comp_yuchkdat, format not recognized'
                    step=0 # use default step 0
            else:
                step=0

            #get threashold index
            thInd=bisect.bisect_left(nts,step)
            if thInd >= len(tol_list):
                if v_level>0:
                    print '!!WARNING step=%i > nts[end]=%i' %(step,nts[-1])
                    print '!!You may want to check tolerance threshold'
                thInd=len(tol_list)-1
            
            #get actual threashold
            tol=tol_list[thInd]
            
            

        #check that l1 and l2 are valid
        if (isValidLine(l1,ValidLineSize,RealPos) and
            isValidLine(l2,ValidLineSize,RealPos)):

            nvalid_line+=1    
        
                
            #----------------------------------------------------------------------------------------------------    
            # Comparing lines
        
            varname=l1[0]
            varname2=l2[0]
                
            #check that it is the same variable in both file
            if varname.strip()!=varname2.strip():
                print '!! Error: Variables differ'
                print ' %s at line %i in file %s' %(varname,i,file1)
                print ' %s at line %i in file %s' %(varname2,i,file2)
                error_count+=1
                return -1

            #compare numerical values on this line
            for id in RealPos:  
                n1=float(l1[id])
                n2=float(l2[id])

            #absolute diffference
                if minval==-1:   
                    ldiff=abs(n1-n2)
            #relative diffference
                elif (abs(n1)>minval or abs(n2)>minval):   
                    ldiff=abs((n1-n2)/(max(abs(n1),abs(n2))))
                else:
                    ldiff=0
 
            #check if larger than tol        
            if (ldiff>tol):
                error_count+=1
                if (ldiff > maxdiff):
                    maxdiff=ldiff
                    maxdiff_line=i+1
                    maxdiff_step=step
                    line1=data1[i]
                    line2=data2[i]
                    

                # print line 
                if (v_level==1):
                    if print_header:
                        print header
                        print_header=False

                    print '>' + data1[i].rstrip()+ '     %i      ' %(step)
                    print '<' + data2[i].rstrip()+ '     %i        %2.1e \n' %(step,ldiff)

        else: #not a valid line
            previousLineWasValid=False
 
    #end of loop over lines

    #print if error detected for verbose 0
    if (v_level==0) and (error_count>0):
        print 'Errors above threshold: %i , max diff  %e at line %i, step %i' %(error_count,maxdiff,maxdiff_line,maxdiff_step)
        print header
        print line1
        print line2

    if v_level>0 and error_count==0:
        print 'no difference above threshold'

    #check there there vas at leaste one valid line
    if nvalid_line==0:
       print '!!Waring: there was no valid line, file cannot be compared'
       return -1

    return error_count

#----------------------------------------------------------------------------
# Local functions

# Valid line contains n elements, real numbers at positions real_list
def isValidLine(Line,n,rList):
    if len(Line)==n:
        test=True
        for id in rList:
            test=(test and isReal(Line[id]))
    else:
        test=False

    return test
    
            

# test a string and return true if it can be converted
# to float and if it is not a integer : must contain a "."
#
def isReal(string):
    try:
        a=float(string)
        if "." in string:
            test=True
        else:
            test=False
    except ValueError:
        test=False
    
    return test

        
    

#-----------------------------------
#execute as a script 
if __name__ == "__main__":

    if len(sys.argv)==3:
        cmp_(sys.argv[1],sys.argv[2])
    elif len(sys.argv)==4:
        cmp_(sys.argv[1],sys.argv[2],int(sys.argv[3]))
    elif len(sys.argv)==5:
        cmp_(sys.argv[1],sys.argv[2],int(sys.argv[3]), \
                 float(sys.argv[4]))
    elif len(sys.argv)==6:
        cmp_(sys.argv[1],sys.argv[2],int(sys.argv[3]), \
                 float(sys.argv[4]), \
                 [float(el) for el in sys.argv[5].split(',')])
    elif len(sys.argv)==7:
        cmp_(sys.argv[1],sys.argv[2],int(sys.argv[3]), \
                 float(sys.argv[4]), \
                 [float(el) for el in sys.argv[5].split(',')], \
                 [float(el) for el in sys.argv[6].split(',')])

    else:
        print '''USAGE : compare_yuchkdat.py file1 file2 [v_level=0 minval=1e-15 
                     nts=10,100,200
                     tol=1e-15,1e-15,1e-15]

DEFINITION :      compare two YUCHKDAT file1, file2 with tolerance tol_*
     Values smaller than minval are not considered.
     Different threshold are used for different set
     A set is defined as contiguous valid lines 
     A valid line contains 10 (from 0 to 9) columns, with column 3,6,9 
     containing real numbers
     Assumes set are separated by at leas one invalid line (ex: blank line)
     v_level:verbose level
            -1 -> no print out
             0 -> max diff 
             1 -> show all lines with differences above tol_
     if minval is -1 compares absolute differences
     the comparison is only done for overlapping time steps'''
    

