#!usr/bin/env python2

"""
COSMO TECHNICAL TESTSUITE

Various tools used to read and modify fortran namelists
"""

# built-in modules
import os, sys, string, re
from ts_error import StopError, SkipError

# information
__author__     = "Xavier Lapillonne, Nicolo Lardelli"
__email__      = "cosmo-wg6@cosmo.org"
__maintainer__ = "xavier.lapillonne@meteoswiss.ch"


namelist_pattern=re.compile(
    r""" ( (?P<varname> [a-zA-Z]\w*)[ ]* = [ ]*                    # this reads the variable name part, it has to start
                                                                   # with a letter and can have one single space before '='
      (?P<arg>                                                     # defines the general argument list
       (?P<tal> (([ ]? , [ ]?)? (['].*?[']))+ )  |                 # defines the text argument list
       (?P<nal> (([ ]? , [ ]?)? ([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?))+)  |  # defines the numerical argument list
       (?P<bla> [.] [a-zA-Z]+ [.]) )                               # defines the bool-like argument
      [ ]? ,?) """, re.VERBOSE )


def get_param(filename, param, ignore_comments=True, occurrence=1):
    """retrieve a parameter from a Fortran namelist file"""

    comt = '!' # the comment character
    foundparam = False
    zcount = 0
 
    try:
        data = open(filename).readlines()
    except :
        raise SkipError('get_param: Error while opening '+filename)

    for line in data:

        # remove comment if required
        if ignore_comments:
            i_cmt = line.find(comt)
            if i_cmt > -1: line = line[0:i_cmt]
        else:
            line = line.replace(comt,' ')
            
        # search for  patterns, i.e. assignements
        matchobj = re.finditer(namelist_pattern,line)
        if matchobj is None:
            continue
        for assignement in matchobj:
            if assignement.group('varname') == param:
                zcount += 1
                if zcount == occurrence:
                    return assignement.group('arg')

    return ''     


def replace_param(filename, param, newparamstr, occurrence=1):
    """replace a namelist parameter in a Fortran namelist file"""

    comt = '!' # the comment character
    istatus = -1
    ierr = 0
    occ_counter = 0
    possible = False



    # check that newparamstr is of the form  param2 = val2    
    if not '=' in newparamstr:
        raise SkipError('replace_param: newparamstr should'
                        'be a string of the form "param2 = val2" ')
 
    # open file if exists
    try:
        data = open(filename).readlines()
    except:
        raise SkipError('replace_param: Error while opening '+filename)

    new_data = []
    
    for line in data:

        # remove end of line charcater
        line = line.strip('\n')
        # separates comments from rest of the line
        i_cmt = line.find(comt)
        if i_cmt == -1: i_cmt = len(line) # if no comments on this line
        line_c = line[i_cmt:]
        line = line[:i_cmt]
        
        # search for attention-relevant patterns, i.e. assignements
        matchobj = re.finditer(namelist_pattern,line)
        if matchobj is None:
            continue
        for assignement in matchobj:
            if assignement.group('varname') == param:
                # change the entige group with the newparamstring
                occ_counter += 1
                
                if occ_counter == occurrence:
                    line = line.replace(assignement.group(),newparamstr+',')
                    modification = True
                    possible = True
        
        new_data.append(line+line_c+'\n')
                             
    fout = open(filename,'w')
    fout.write(string.join(new_data,''))
    fout.close

    # if no return has been so far executed the code will have no significance any more
    if not possible:
        raise SkipError('replace_param: unable to successfully find parameter '+param+' into '+newparamstr)

    if not modification:
        raise SkipError('replace_param: unable to successfully modify the '+occurrence+
                        'th occurrence of parameter '+param+' into '+newparamstr)


