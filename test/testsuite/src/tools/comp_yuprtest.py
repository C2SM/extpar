#!/usr/bin/env python2

"""
COSMO TECHNICAL TESTSUITE

General purpose script to compare two YUPRTEST output files
"""

# built-in modules
import os, sys, string

# information
__author__     = "Xavier Lapillonne, Nicolo Lardelli"
__email__      = "cosmo-wg6@cosmo.org"
__maintainer__ = "xavier.lapillonne@meteoswiss.ch"


def cmp_(file1,file2, \
             v_level=0,minval=1e-15, \
             nts=[10,100,200], \
             tol_ts=[1e-15,1e-15,1e-15], \
             tol_as=[1e-15,1e-15,1e-15]):

    # compare two YUPRTEST file1, file2 with tolerance tol_*
    # Values smaller than minval are not considered.
    # Different threshold are used for different time step.
    # v_level:verbose level
    #        -1 -> no print out
    #         0 -> max diff over all variables for each time step, short
    #         1 -> max diff over all variables for each time step, print lines
    #         2 -> show all lines with differences above tol_
    # if minval is set to -1 compares absolute differences
    # the comparison is only done for overlapping time steps
    
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
    error_count = 0       #number of error detected
    tol_t = tol_ts[0]     #set tolerence for t
    tol_a = tol_as[0]     #set tolerence for all other variables
   
    maxdiff_t=0.     #maximal diff per time step for t or p
    line1_t=''       #line file1 with max diff
    line2_t=''       #line file2 with max diff
    lnum_t=0


    maxdiff_a=0.      #maximal diff per time step for all variables
    line1_a=''        #line file1 with max diff
    line2_a=''        #line file2 with max diff
    lnum_a=0 
  

    print_header = True
    header = ' '
    
    comment_type='#' #comment at the begining wil be skiped
    
    
    ntstep=0
    lerror_t= False
    lerror_a = False

    if v_level==0:
        print_out='./tools/comp_yuprtest.py ' + file1 + ' ' + file2 + ' ' + str(v_level) + ' ' + str(minval) + \
                  ' ' + ','.join([str(x) for x in nts]) + ' ' + ','.join([str(x) for x in tol_ts]) + \
                  ' ' + ','.join([str(x) for x in tol_as]) + '\n'
        if minval==-1:
            print_out+='Absolute error:\n'
        else:
            print_out+='Relative error:\n'
        print_out+='   nt    max_all         t     Test \n'

    if v_level>0:
        if minval==-1:
            print 'Comparing absolute differences ...'
        else:
            print 'Comparing relative differences, min. value is %1.0e ...' %(minval)
    

    # check that files are not empty
    if len(data1)==0:
        print 'file ' + file1 + ' is empty!'
        return -1
    if len(data1)<=4:
        print 'file ' + file1 + ' contains only header!'
        return -1

    if len(data2)==0:
        print 'file ' + file2 + ' is empty!'
        return -1
    if len(data2)<=4:
        print 'file ' + file2 + ' contains only header!'
        return -1

    # set file counter
    i1=0
    i2=0

    # remove the headers part (all lines starting with comment_type)
    # first on file1
    while True:
        line=data1[i1].split()
        if line[0]!=comment_type:
            break
        i1=i1+1

    # then on file2
    while True:
        line=data2[i2].split()
        if line[0]!=comment_type:
            break
        i2=i2+1

    # Set the file counters to identical time step
    while True:
        #check eof
        if (i1>=len(data1)) or (i2>=len(data2)):
            print 'Files %s and %s do not have overlapping time steps and can not be compared.' %(file1,file2)
            return -1
            
        l1=data1[i1].split()
        l2=data2[i2].split()
        nt1=int(l1[1])
        nt2=int(l2[1])
        
        if nt1<nt2:
            i1=i1+1
        elif nt1>nt2:
            i2=i2+1
        elif nt1==nt2:
            break
        else:
            break
        

    ntstep=nt1
    leof=False

    #----------------------------------------------------------------------------------------------------
    #loop over file lines
    while True:
       
        #check eof
        if (i1>=len(data1)) or (i2>=len(data2)):
            leof=True
        #read file
        else:
            l1=data1[i1].split()
            l2=data2[i2].split()

        #----------------------------------------------------------------------------------------------------
        #prepare printout if new time step or eof
        if (int(l1[1]) != ntstep ) or leof:
            if (v_level==0):
                if (lerror_t) or (lerror_a):
                    print_out+='%4i     %1.2e     %1.2e     FAILED \n' %(ntstep,maxdiff_a,maxdiff_t)
                else:
                    print_out+='%4i     %1.2e     %1.2e     OK     \n' %(ntstep,maxdiff_a,maxdiff_t)

            #print if verbose=1 and error at this step
            if (lerror_t) and (v_level==1):
                if print_header:
                    print header
                    print_header=False
                print 'nt=%i, max rel. er. t,p: %1.1e above threshold %1.1e, at line %i' %(ntstep,maxdiff_t,tol_t,lnum_t)
                print '>'+ line1_t.rstrip()
                print '<'+ line2_t
            if (lerror_a) and (v_level==1):
                if print_header:
                    print header
                    print_header=False    
                print 'nt=%i, max rel. er. all: %1.1e above threshold %1.1e, at line %i' %(ntstep,maxdiff_a,tol_a,lnum_a)
                print '>'+ line1_a.rstrip()
                print '<'+ line2_a

                
            #exit loop if eof
            if leof:
                break
                
            #set step and reset local error counter
            ntstep=int(l1[1])
            lerror_t=False       
            lerror_a =False   
            maxdiff_t=0.
            maxdiff_a=0.
               

            # Set threshold
            #
            # tol[i] is set for t=[nts[i] nts[i+1]]

            for i in range(len(nts)):
                # update threshold if step larger than nts
                if ntstep >= nts[i]:
                    #only update if i<=len(tol_ts)
                    if (i<len(tol_ts)-1): 
                        tol_t=tol_ts[i+1]
                        tol_a=tol_as[i+1]
                        
                
        #----------------------------------------------------------------------------------------------------    
        # Comparing lines
        
        varname=l1[0]
        varname2=l2[0]
                
        #check that it is the same variable in both file
        if varname.strip()!=varname2.strip():
            print '!! Error: Variables differ'
            print ' %s at line %i in file %s' %(varname,i1+1,file1)
            print ' %s at line %i in file %s' %(varname2,i2+1,file2)
            error_count+=1
            lerror_t=True      
            lerror_a =True 
            return -1

        #check that it is the same time step
        if int(l1[1])!=int(l2[1]):
            print '!! Error: Time steps differ'
            print ' nt=%s at line %i in file %s' %(l1[1],i1+1,file1)
            print ' nt=%s at line %i in file %s' %(l2[1],i2+1,file2)
            error_count+=1
            lerror_t=True      
            lerror_a =True
            return -1

        #compare numerical values on this line
        for j in range(len(l1)):                    
            pr_line=True
            if is_num(l1[j]):
                n1=float(l1[j])
                n2=float(l2[j])

                #absolute diffference
                if minval==-1  and is_not_int(n1):   #note: int are not considered (min-max index)   
                    ldiff=abs(n1-n2)
                #relative diffference
                elif abs(n1)>minval  and is_not_int(n1):   #note: int are not considered (min-max index)                          
                    ldiff=abs((n1-n2)/n1)
                else:
                    ldiff=0
                
                #Use tol_t threshold for temperature field
                if (varname in ['T']):
                    # save max
                    if ldiff > maxdiff_t:
                        maxdiff_t=ldiff
                        line1_t=data1[i1]
                        line2_t=data2[i2]
                        lnum_t=i1+1

                    #check if larger than tol
                    if ldiff > tol_t:
                        error_count+=1
                        lerror_t=True
                        # print line 
                        if (v_level==2 and pr_line):
                            if print_header:
                                print header
                                print_header=False
                            if pr_line:
                                print '>' + data1[i1].rstrip()
                                print '<' + data2[i2]
                                pr_line=False 

                #Use tol_a threshold for all other fields
                else:
                    # save max
                    if ldiff > maxdiff_a:
                        maxdiff_a=ldiff
                        line1_a=data1[i1]
                        line2_a=data2[i2]
                        lnum_a=i1+1
                                
                    #check if larger than tol
                    if ldiff > tol_a:
                        error_count+=1
                        lerror_a=True
                        # print line 
                        if (v_level==2 and pr_line):
                            if print_header:
                                print header
                                print_header=False
                            if pr_line:
                                print '>' + data1[i1].rstrip()
                                print '<' + data2[i2]
                                pr_line=False 


        # moves forward the 2 file counters
        i1=i1+1
        i2=i2+1


    #print if error detected for verbose 0
    if (v_level==0) and (error_count>0):
        print print_out

    if v_level>0 and error_count==0:
        print 'no difference above threshold'

    return error_count

#----------------------------------------------------------------------------
# Local functions
def is_num(x):
    test=True
    try:
        a=float(x)
    except ValueError:
        test=False
        
    return test

def is_not_int(x):
    #note: if x=0.0 is return True for this test
    test=False
    try:
        a=int(x)
        if x!=0.0:
            if a/x!=1:
                test=True
        else:
            test=True
            
            
    except ValueError:
        test=True
        
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
    elif len(sys.argv)==8:
        cmp_(sys.argv[1],sys.argv[2],int(sys.argv[3]), \
                 float(sys.argv[4]), \
                 [float(el) for el in sys.argv[5].split(',')], \
                 [float(el) for el in sys.argv[6].split(',')], \
                 [float(el) for el in sys.argv[7].split(',')])
    else:
        print '''USAGE : compare_files.py file1 file2 [v_level=0 minval=1e-15 
                     nts=10,100,200
                     tol_ts=1e-15,1e-15,1e-15
                     tol_as=1e-15,1e-15,1e-15 ]
DEFINITION : Compare relative differences between two YUPRTEST file1, file2 
with tolerance tol_*. Values smaller than minval are not considered.
Different thresholds are used for different time steps. Thresholds tol_ts are used for 
the temperature field while tol_as are used for all other fields.
if minval set to -1 compares absolute instead of relative differences. '''
    

