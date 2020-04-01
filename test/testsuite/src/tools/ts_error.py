#!/usr/bin/env python

"""
COSMO TECHNICAL TESTSUITE

Error classes used for exception handling
"""

# built-in modules
import re

# information
__author__     = "Nicolo Lardelli, Xavier Lapillonne"
__email__      = "cosmo-wg6@cosmo.org"
__maintainer__ = "xavier.lapillonne@meteoswiss.ch"


class SkipError(RuntimeError):
    """Exception when a test has to be skipped (not applicable or unmatched requirement"""
    pass

class StopError(RuntimeError):
    """Exception when a test has encountered an error"""
    pass


