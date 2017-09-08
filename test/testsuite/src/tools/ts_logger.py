#!/usr/bin/env python2

"""
COSMO TECHNICAL TESTSUITE

This module implements a logger class which handles messages from the testsuite
to the user. It is based on the standard Python logger class. The following logging
levels are supported: DEBUG, INFO, WARNING, ERROR
"""

# built-in modules
import sys, string
import logging as LG

# private modules
from ts_utilities import status_str, pretty_status_str

# information
__author__     = "Oliver Fuhrer"
__email__      = "cosmo-wg6@cosmo.org"
__maintainer__ = "xavier.lapillonne@meteoswiss.ch"


# status column length
STAT_COLUMN = 12

# logging levels
DEBUG       = LG.DEBUG     # 10
INFO        = LG.INFO      # 20
WARNING     = LG.WARNING   # 30
CHCKINFO    = 35           # 35
IMPORTANT = 35             # 35
ERROR     = LG.ERROR     # 40

# formatting
FORMAT = {
    DEBUG     : " "*(STAT_COLUMN+2-7) + "[ DBG ] %(msg)s",
    INFO      : "%(msg)s",
    WARNING   : "*** WARNING: %(msg)s",
    CHCKINFO  : "%(msg)s",
    IMPORTANT : "%(msg)s",
    ERROR     : "*** ERROR: %(msg)s"
}


class MyFormatter(LG.Formatter):
    """Custom formatter which allows different formatting for different levels"""
    
    def __init__(self, fmt="%(levelno)s: %(msg)s"):
        #super(MyFormatter,self).__init__(fmt)
        LG.Formatter.__init__(self, fmt)

    def format(self, record):

        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._fmt

        # Replace the original format with one customized by logging level
        self._fmt = FORMAT.get(record.levelno)

        # Call the original formatter class to do the grunt work
        #result = super(MyFormatter, self).format(record)
        result = LG.Formatter.format(self, record)

        # Restore the original format configured by the user
        self._fmt = format_orig

        return result


class Logger:
    """simple logger which is used throughout the testsuite"""

    color = True

    def __init__(self, filename, append=False, color=False):
        self.filename = filename
        if append:
            self.mode = 'a'
        else:
            self.mode = 'w'
        self.logger = LG.getLogger('testsuite')
        if filename:
            self.handler = LG.FileHandler(filename,mode=self.mode,delay=False)
        else:
            self.handler = LG.StreamHandler(sys.stdout)
        self.formatter = MyFormatter()
        self.handler.setFormatter(self.formatter)
        self.logger.addHandler(self.handler)
        self.logger.setLevel(INFO)
  
    def __del__(self):
        LG.shutdown()
  
    def setLevel(self, lvl):
        self.logger.setLevel(lvl)
  
    def getLevel(self):
        return self.logger.getEffectiveLevel(lvl)

    def log(self, lvl, msg, *args, **kwargs):
        self.logger.log(lvl, msg, *args, **kwargs)
  
    def debug(self, msg, *args, **kwargs):
        self.log(DEBUG, msg, *args, **kwargs)
  
    def info(self, msg, *args, **kwargs):
        prefix = ' '*(STAT_COLUMN+2) + ' '
        self.log(INFO, prefix + msg, *args, **kwargs)
  
    def warning(self, msg, *args, **kwargs):
        self.log(WARNING, msg, *args, **kwargs)

    def chckinfo(self, msg, *args, **kwargs):
        prefix = ' '*(STAT_COLUMN+4) + ' '
        self.log(CHCKINFO, prefix + msg, *args, **kwargs)
  
    def important(self, msg, *args, **kwargs):
        prefix = '[' + '-'*STAT_COLUMN + '] '
        self.log(IMPORTANT, prefix + msg, *args, **kwargs)
  
    def result(self, indent, status, msg, *args, **kwargs):
        slen = len(status_str(status))
        status = pretty_status_str(status, Logger.color, indent==0)
        pad = STAT_COLUMN - slen - 2*indent
        prefix=' '*(2*indent) + '[' + ' '*(pad/2) + status + ' '*(pad-pad/2) + '] '
        self.log(IMPORTANT, prefix + msg, *args, **kwargs)
  
    def error(self, msg, *args, **kwargs):
        self.log(ERROR, msg, *args, **kwargs)
  
    def flush(self):
        self.handler.flush()

