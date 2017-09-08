#!/usr/bin/env python

import os, re, sys, unittest


# information
__author__      = "Santiago Moreno"
__email__       = "cosmo-wg6@cosmo.org"
__maintainer__  = "xavier.lapillonne@meteoswiss.ch"

# get name of myself
myname = os.path.basename(__file__)
header = myname + ': '

class Pattern:
    """
    The Pattern class is the parent class of all patterns. 

    Every child of Pattern needs to define what happens when it gets matched 
    and what happens when it doesn't get matched.
    e.g.:
    
    class OccurrencePattern(Pattern):
        # the pattern has to return 0,10,20 as the testsuite works with these numbers
        def pattern_match(self):
            print ("pattern got matched!")
            return self._check_ok()

        def pattern_no_match(self):
            print ("pattern not found")
            return self._check_failed()
    
        pattern_hits = 0
        first_line_hit = 0
    """
    def __init__(self, description, search_pattern):
        """
        description     will be used in logging as the name of the pattern
        search_pattern   expects a regex pattern which will be used to check if the warning exists
        """
        self.description = description
        self.search_pattern = search_pattern

    def match(self, text):
        """
        searches the given text for the defined search_pattern.
        match() returns true if it gets found and false if it doesn't.
        instead of searching for one hit this function gets all the hits
        and saves the count for later use in logging.
        match() also counts the linenumber of the first hit for the purpose of logging
        """
        regex = re.compile(self.search_pattern)
        hits = 0
        hit_once = False
        first_line_hit = 0
        if isinstance(text, str):
            text = text.split("\n")
        for line in text:
            if not hit_once:
                first_line_hit += 1
            if(regex.search(line)):
                hit_once = True
                hits += 1
        self.first_line_hit = first_line_hit
        self.pattern_hits = hits
        return hit_once
    def _check_match(self):
        return 0

    def _check_ok(self):
        return 10

    def _check_failed(self):
        return 20

    def _check_crash(self):
        return 30

    def check(self, text, verbose):
        self.verbose = verbose
        if self.match(text):
            return self.pattern_match()
        else:
            return self.pattern_no_match()


class OccurrencePattern(Pattern):
    """will return FAIL if it doesn't get matched, MATCH if it does"""
    def pattern_match(self):
        if self.verbose > 2:
            print(header + "The Pattern " + self.description + " has been found " 
            + str(self.pattern_hits) + " times. First occurrence on line " + str(self.first_line_hit))
        return self._check_match()

    def pattern_no_match(self):
        if self.verbose > 0:
            print(header + "no occurrence of " + self.description + " has been found")
        return self._check_failed()


class ErrorPattern(Pattern):
    """will return FAIL if it gets matched, MATCH if it doesn't"""
    def pattern_match(self):
        if self.verbose > 0:
            print(header + "The Pattern " + self.description + " has been found "
            + str(self.pattern_hits) + " times. First occurrence on line " + str(self.first_line_hit))
        return self._check_failed()

    def pattern_no_match(self):
        return self._check_match()


class WarningPattern(Pattern):
    """will print a warning if this is found in the log, return ok code"""
    def pattern_match(self):
        if self.verbose > 1:
            print(header + "The Pattern " + self.description + " has been found " 
            + str(self.pattern_hits) + " times. First occurrence on line " + str(self.first_line_hit))
        return self._check_match()

    def pattern_no_match(self):
        return self._check_match()

class OccurrenceCrashPattern(OccurrencePattern):
    """will return crash code if this pattern is not found"""
    def _check_failed(self):
        return 30

class FileChecker:
    """
    The FileChecker class goes through a given list of patterns stored in pattern_list
    to then output the most severe error returned by said patterns
    """
    def __init__(self):
        self.pattern_list = []
    
    def read_file(self, filepath):
        """reads the given filepath and saves the resulting
        string in the variable text for later use"""
        try:
            file = open(filepath, 'r')
            text = file.readlines()
            file.close()
            return text

        except:
            if verbose:
                print(header + 'failed to open ' + filepath)
                return 30  # CRASH

    def add_pattern(self, pattern):
        """appends a Pattern object to the pattern list"""
        self.pattern_list.append(pattern)

    def add_pattern_list(self, pattern_list):
        """extends the pattern_list with the given list"""
        self.pattern_list.extend(pattern_list)

    def check_patterns(self, verbose, text):
        """
        go through the pattern list execute the check method if the result ist more severe than 
        the currently highest one make it the current highest one
        """
        result = 0
        for pattern in self.pattern_list:
            pattern_result = pattern.check(text, verbose)
            if result < pattern_result:
                result = pattern_result
        return result

    def check(self, logfile, verbose):
        text = self.read_file(logfile)
        return self.check_patterns(verbose, text)
     
#-----------------------unit tests------------------------------


class PatternTests(unittest.TestCase):
    warningpattern = WarningPattern("WARNING pattern", "WARNING")
    errorpattern = ErrorPattern("ERROR pattern", "ERROR")
    occurrencepattern = OccurrencePattern("OCCURRENCE pattern", "OCCURRENCE")
    verbosity_level = 0
    def  test_warning_pattern(self):
        self.assertEquals(self.warningpattern.check("WARNING", self.verbosity_level ), 0)
        self.assertEquals(self.warningpattern.check("not an expected String", self.verbosity_level), 0)
        
    def  test_error_pattern(self):
        self.assertEquals(self.errorpattern.check("ERROR", self.verbosity_level) , 20)
        self.assertEquals(self.errorpattern.check("not an expected String", self.verbosity_level), 0)
        
    def  test_occurrence_pattern(self):
        self.assertEquals(self.occurrencepattern.check("OCCURRENCE", self.verbosity_level), 0)
        self.assertEquals(self.occurrencepattern.check("not an expected String", self.verbosity_level), 20)

        
class FileCheckerTest(unittest.TestCase):
    warningpattern = WarningPattern("WARNING pattern", "WARNING")
    errorpattern = ErrorPattern("ERROR pattern", "ERROR")
    occurrencepattern = OccurrencePattern("OCCURRENCE pattern", "OCCURRENCE")
    verbosity_level = 0
    test_patterns = [warningpattern, errorpattern, occurrencepattern]
    filechecker = FileChecker()
    filechecker.add_pattern_list(test_patterns)

    def  test_match(self):
        self.assertEqual(self.filechecker.check_patterns
                (self.verbosity_level, "qwerqwerqwerOCCURRENCEwwegwerg"), 0)
        
    def  test_ok(self):
        self.assertEqual(self.filechecker.check_patterns
                (self.verbosity_level, "qVQWEvOCCURRENCEqwerqwfWARNING"), 0)
        self.assertEqual(self.filechecker.check_patterns
                (self.verbosity_level, "WETOUWEBIBWEVBIWARNINGqVQWEvOCCURRENCEqwerqwf"), 0)
        
    def  test_fail(self):
        self.assertEqual(self.filechecker.check_patterns
                (self.verbosity_level, "qwerqwWARNINGerqwerOCCURRENCEwwegERRORwerg"), 20)
        self.assertEqual(self.filechecker.check_patterns
                (self.verbosity_level, "qwerqERRORwerqwerOCCURRENCEwwegWARNINGwerg"), 20)
        self.assertEqual(self.filechecker.check_patterns
                (self.verbosity_level, "WARNING"), 20)
    def test_multiline(self):
        text = """
            + + + + + + + + + + + + + + + +
            + RUNNING IN DOUBLE PRECISION +
            + + + + + + + + + + + + + + + +

            SETUP OF THE LM
            INITIALIZATIONS
            Info about KIND-parameters:   iintegers / MPI_INT =            4  1275069467
                                          int_ga    / MPI_INT =            4  1275069467
                                          INPUT OF THE NAMELISTS
                           *** NOTE: Old 10 OCCURRENCE digit date format is used

                            ==== Code information used to build this binary ====
                            Binary name ....: cosmo

                            Library name ......: cosmo5.0n1
                            Tag name ..........: git@github.com:Santiago-Moreno/cosmo-pompa.git
                            Revision number ...: 5eabd73  
                            Checkin-Date ......: 2016-02-16 15:15:29 +0100
                            Put into production: (missing)
                            Checkout-Date .....: (missing)
                            Code is modified ..: true
                            Compile-Date ......: Wed Mar  9 11:15:40 CET 2016
                            Compiled by .......: msa
                            Compiled on .......: lema
        """ 
        
        self.assertEqual(self.filechecker.check_patterns
                (self.verbosity_level, text), 0)
        
        self.assertEqual(self.filechecker.pattern_list[2].pattern_hits, 1)
        self.assertEqual(self.filechecker.pattern_list[2].first_line_hit, 11)

if __name__ == '__main__':
    unittest.main()
