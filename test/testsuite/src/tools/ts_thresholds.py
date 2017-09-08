#!/usr/bin/env python

"""
COSMO TECHNICAL TESTSUITE

Python module to handle testsuite thresholds. Thresholds can be read/written
from/to a Python dictionary or a file of the form...

 minval = 1e-12
  steps =          3          8         20         60
      * =   1.00e-13   1.00e-10   1.00e-06   1.00e-02
      T =   1.00e-11   1.00e-08   1.00e-05   1.00e+00

...where "*" are the default thresholds, which can be overridden for any
variable via custom thresholds. The minimum value (minval) is the
magnitude below which numbers should no longer be compared against
thresholds.

t = Thresholds("THRESHOLDS")         # initialize object from file (from dict works also)
print(t)                             # convert to string and echo to stdout
t == t                               # one can compare if thresholds are equal
t.to_dict()                          # convert thresholds to dictionary
t.to_file('THRESH')                  # write to file
print(t.variables)                   # print names of variables with special thresholds
t.steps                              # list containing steps with defined thresholds
t.add_variable('TKE')                # add a variable with specific thresholds
                                     # (inhertis from the default thresholds)
t.remove_variable('TKE')             # remove a variable with specific thresholds
t.add_step(40)                       # add a timestep with a defined threshold
                                     # (initial value is interpolated from existing)
t.removeStep(40)                     # remove a timestep with a defined threshold
t['T']                               # return thresholds for a specific variable
t['TKE'] = [1e-12, 1e-9, 1e-05, 0.1] # set thresholds for a specific variable
del(t['TKE'])                        # remove thresholds for a specific variable
t.minval = 1.0e-7                    # values below this minimum are not considered
t.mode = 'log'                       # set interpolation mode to logarithmic
t.get_threshold('T',17)              # retrieve a threshold for a certain timestep
                                     # (values in between defined timesteps are interpolated)
t.digits = 2                         # number of digits to retain for thresholds
t.increase_factor = 2.0              # increase factor when setting thresholds
t.update_threshold('T',17,1.8e-6)    # update threshold with a specific value

"""

# built-in modules
import unittest
import tempfile
import re, os
import math

# information
__author__      = "Oliver Fuhrer, Santiago Moreno"
__email__       = "cosmo-wg6@cosmo.org"
__maintainer__  = "oliver.fuhrer@meteoswiss.ch"

# get name of myself
myname = os.path.basename(__file__)
header = myname + ': '

class Thresholds(object):
    """class to return continuous thresholds from discrete set of values of
       timesteps where thresholds can be dependent on variable"""

    def __init__(self, t):
        """initialize threshold object from dict or file"""

        # minmum value below which values are not considered
        self._minval = None 
        # steps where thresholds change
        self._steps = []
        # list of default thresholds
        self._default = []
        # list of variables for which specific thresholds are defined
        self._variables = []
        # list of lists which contains thresholds
        self._thresholds = []
        # interpolation mode (const, linear or log)
        self._mode = 'const'
        # Create non-existing variables when during tuning
        self._create_nonexisting_variables = False
        # number of digits to retain for thresholds
        self._digits = 2
        # factor to increase values for threshold setting
        self._increase_factor = 10.0
        # Excluded variables
        self._excluded_variables = [ "CHKDAT" ]
        
        if t:
            if isinstance(t, dict):
                self.from_dict(t)
            elif isinstance(t, str):
                if re.match('[\n ]', t):
                    self.from_str(t)
                else:
                    self.from_file(t)
            else:
                raise NotImplementedError('initializer form ' + type(t) + ' not implemented yet')

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.to_dict() == other.to_dict()
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __getitem__(self, index):
        return self.__get_threshold_values(index)

    def __setitem__(self, index, value):
        self.add_variable(index)
        self.__set_threshold_values(index, value)

    def __delitem__(self, index):
        self.remove_variable(index)

    def __str__(self):
        s = "   minval = " + str(self._minval) + "\n"
        s += "    steps = " + ' '.join(["{0:10d}".format(x) for x in self._steps]) + "\n"
        s += "        * = " + ' '.join(["{0:10.2e}".format(x) for x in self._default]) + "\n"
        v, t = self.__compress()
        for i in range(len(v)):
            s += "%9s = " % v[i] + \
                 ' '.join(["{0:10.2e}".format(x) for x in t[i]]) + "\n"
        return s.rstrip()

    @property
    def minval(self):
        return self._minval

    @minval.setter
    def minval(self, value):
        self._minval = value

    @property
    def variables(self):
        return self._variables

    @property
    def steps(self):
        return self._steps

    @property
    def mode(self):
        return self._mode

    @mode.setter
    def mode(self, value):
        if value in ['const', 'linear', 'log']:
            self._mode = value
        else:
            raise ValueError('Illegal mode specified')

    @property
    def digits(self):
        return self._digits

    @digits.setter
    def digits(self, value):
        if value < 1 or value > 10:
            raise ValueError('Illegal number of digits specified')
        self._digits = value

    @property
    def increase_factor(self):
        return self._increase_factor

    @increase_factor.setter
    def increase_factor(self, value):
        if value < 1.0 or value > 10.0:
            raise ValueError('Illegal increase_factor value specified')
        self._increase_factor = value

    def __compress(self):
        v = []
        t = []
        # generate list of variables with unique thresholds
        for i in range(len(self._variables)):
            try:
                j = t.index(self._thresholds[i])
                v[j] = v[j] + ',' + self._variables[i]
            except:
                v.append(self._variables[i])
                t.append(self._thresholds[i])
        # sort variables with identical thresholds alphabetically
        for i in range(len(v)):
            v[i] = ','.join(sorted(v[i].split(',')))
        return v, t

    def from_str(self, string):
        """parse thresholds from string"""
        self._variables = []
        for line in string.split('\n'):
            if not line:
                continue
            (key, value) = line.split('=')
            key = key.strip()
            if key == 'minval':
                self._minval = float(value)
            elif key == 'steps':
                self._steps = [int(x) for x in value.split()]
            elif key == '*':
                self._default = [float(x) for x in value.split()]
            else:
                for y in key.split(','):
                    self._variables.append(y)
                    self._thresholds.append([float(x) for x in value.split()])

    def to_dict(self):
        """return a dictionary with thresholds"""
        d = {}
        d["minval"] = self._minval
        d["steps"] = self._steps
        d["*"] = self._default
        v, t = self.__compress()
        for x, y in zip(v, t):
            d[x] = y
        return d

    def from_dict(self, d):
        """parse a dictionary"""
        self._minval = d["minval"]
        self._default = d["*"]
        self._steps = d["steps"]
        self._variables = []
        self._thresholds = []
        vars = [x for x in d.keys() if x not in ["minval", "*", "steps"]]
        for x in vars:
            for y in x.split(','):
                self._variables.append(x)
                self._thresholds.append(d[x])

    def to_file(self, filename):
        """write thresholds to a file"""
        f = open(filename, 'w')
        f.write(str(self))
        f.close()

    def from_file(self, filename):
        """read thresholds to a file"""
        self._variables = []
        f = open(filename, 'r')
        data = f.read()
        f.close()
        self.from_str(data)

    def __get_index_from_var(self, variable):
        i = [i for i, x in enumerate(self._variables) if x == variable]
        if i:
            return i[0]
        else:
            return None

    def __get_threshold_values(self, variable):
        """returns a list of thresholds values for a variable"""
        i = self.__get_index_from_var(variable)
        if i is not None:
            return self._thresholds[i]
        else:
            return self._default

    def _set_thresholds_to_zero(self):
        """Set the thresholds to zero, with the exception of the excluded variables"""
        for i, threshold in enumerate(self._thresholds):
            for j in range(len(threshold)):
                # Only set non-excluded variables to zero
                if self._variables[i] not in self._excluded_variables:
                    self._thresholds[i][j] = 0.0
        for i in range(len(self._default)):
            self._default[i] = 0.0

    def __set_threshold_values(self, variable, thresholds):
        """Set a list of thresholds values for a variable"""
        if len(thresholds) != len(self.steps):
            raise ValueError('Must supply one threshold for each step')
        i = self.__get_index_from_var(variable)
        if i is not None:
            self._thresholds[i] = thresholds
        else:
            self._default = thresholds

    def __get_indices(self, step):
        """find inidices to use for interpolation at step"""

        imin = [i for i, x in enumerate(self._steps) if step >= x]
        imax = [i for i, x in enumerate(self._steps) if step <= x]
        if imax:
            imax = imax[0]
        else:
            imax = imin[-1]
        if imin:
            imin = imin[-1]
        else:
            imin = imax
        return imin, imax

    def __interpolate_threshold(self, y1, y2, x):
        """interpolate treshold between two values"""
        if self._mode == 'const':
            return y2
        elif self._mode == 'linear':
            return y1 * (1.0 - x) + y2 * x
        elif self._mode == 'log':
            return math.exp(math.log(y1) * (1.0 - x) + math.log(y2) * x)
        else:
            raise ValueError('Illegal mode encountered')

    def __compute_threshold(self, value):
        """compute a threshold for a specific value"""
        x = value * self._increase_factor
        shift = -int(math.floor(math.log10(x)))
        factor = math.pow(10.0, shift)
        thresh = math.ceil(x * factor) / factor
        return thresh

    def get_threshold(self, variable, step):
        """return a threshold for a certain variable and step"""
        t = self.__get_threshold_values(variable)
        imin, imax = self.__get_indices(step)
        stepmin = self._steps[imin]
        stepmax = self._steps[imax]
        if stepmin != stepmax:
            return self.__interpolate_threshold(t[imin], t[imax], float(step - stepmin) / (stepmax - stepmin))
        else:
            return t[imax]

    def update_threshold(self, variable, step, value):
        """update threshold values to accomodate current value"""
        t = self.get_threshold(variable, step)
        x = self.__compute_threshold(value)
        if value > t:
            print (header + " thresholds had to be changed at: var= " + str(variable)
                    + " step = " + str(step) + " value = " + str(value))
            if variable not in self._thresholds and self._create_nonexisting_variables:
                self.add_variable(variable)
            t = self.__get_threshold_values(variable)
            imin, imax = self.__get_indices(step)
            stepmin = self._steps[imin]
            stepmax = self._steps[imax]
            if stepmin != stepmax:
                position = float(step - stepmin) / float(stepmax - stepmin)
                tnew = self.__interpolate_threshold(t[imin], value, position)
                t[imax] = max(t[imax], tnew)
            else:
                t[imax] = value
            t[imax] = self.__compute_threshold(t[imax])
            for i in range(imax, len(t)):
                t[i] = max(t[i], t[imax])
            self.__set_threshold_values(variable, t)

    def update_default_thresholds(self, default_variable='*'):
        """Update the default threshold, usually * with the maximum of the variables of all 
        the timesteps"""
        # Loop over all variables:
        default_values = self.__get_threshold_values(default_variable)
        for var in self._variables:
            # Exclude the variable if it is excluded or the default variable
            if var in self._excluded_variables or var == default_variable:
                continue
            var_values = self.__get_threshold_values(var)
            for step_idx in range(len(self._steps)):
                default_values[step_idx] = max(default_values[step_idx], var_values[step_idx])
        self.__set_threshold_values(default_variable, default_values)

    def add_variable(self, variable):
        """insert a new variable where special thresholds are defined"""
        if variable in self._variables:
            return
        self._variables.append(variable)
        self._thresholds.append(list(self._default))

    def remove_variable(self, variable):
        """remove a variable where special thresholds are defined"""
        try:
            i = self._variables.index(variable)
        except:
            return
        del (self._variables[i])
        del (self._thresholds[i])

    def __insert_step(self, step):
        """insert a new step and set thresholds to None"""
        imin, imax = self.__get_indices(step)
        if imin == imax and imin == 0:
            i = 0
        elif imin == imax and imax == len(self._steps) - 1:
            i = imax + 1
        else:
            i = imax
        self._steps.insert(i, step)
        self._default.insert(i, self._default[imin])
        for j in range(len(self._thresholds)):
            self._thresholds[j].insert(i, self._thresholds[j][imin])

    def add_step(self, step):
        """insert a new step and interpolate thresholds"""
        if step in self._steps:
            return
        # backup increase_factor value
        increase_factor = self._increase_factor
        self._increase_factor = 1.0
        # compute thresholds at new step
        x = []
        for variable in ['*'] + self._variables:
            x.append(self.get_threshold(variable, step))
        # add new step
        self.__insert_step(step)
        # update thresholds that have been previously computed
        for variable in ['*'] + self._variables:
            self.update_threshold(variable, step, x.pop(0))
        # revert back increase_factor value
        self._increase_factor = increase_factor

    def removeStep(self, step):
        """remove a step where thresholds are defined"""
        try:
            i = self._steps.index(step)
        except:
            return
        del (self._steps[i])
        del (self._default[i])
        for j in range(len(self._variables)):
            del (self._thresholds[j][i])


class Test(unittest.TestCase):
    def setUp(self):
        self._s = """
 minval = 1e-12
  steps =          3          8         20         60
      * =   1.00e-13   1.00e-10   1.00e-06   1.00e-02
      T =   1.00e-11   1.00e-08   1.00e-05   1.00e+00
"""
        self._d = {
            "minval": 1.0e-12,
            "steps": [3, 8, 20, 60],
            "*": [1.0e-13, 1.0e-10, 1.0e-06, 1.0e-02],
            "T": [1.0e-11, 1.0e-08, 1.0e-05, 1.0e-00]
        }
        f = tempfile.NamedTemporaryFile(mode='w+b', delete=False)
        self._filename = f.name
        f.write(self._s)
        f.close()

    def tearDown(self):
        pass

    def test_basic_setup(self):
        t1 = Thresholds(self._s)  # from string
        t2 = Thresholds(self._d)  # from dictionary
        self.assertEqual(t1, t2)
        t3 = Thresholds(self._filename)  # from file
        self.assertEqual(t1, t3)

    def test_output(self):
        t1 = Thresholds(self._s)
        t2 = Thresholds(self._s)
        # test string output and input
        t2.from_str(str(t2))
        self.assertEqual(t1, t2)
        # test dictionary output and input
        t2.from_dict(t2.to_dict())
        self.assertEqual(t1, t2)
        # test file output and input
        t2.to_file(self._filename)
        t2.from_file(self._filename)
        self.assertEqual(t1, t2)

    def test_thresholds(self):
        t = Thresholds(self._s)
        self.assertAlmostEqual(t.get_threshold('PP', 8), 1.0e-10)
        self.assertAlmostEqual(t.get_threshold('T', 8), 1.0e-8)
        self.assertAlmostEqual(t.get_threshold('PP', 1), 1.0e-13)
        self.assertAlmostEqual(t.get_threshold('T', 1), 1.0e-11)
        t1 = t.get_threshold('PP', 8)
        t3 = t.get_threshold('PP', 20)
        t.mode = 'const'
        t2 = t.get_threshold('PP', 11)
        self.assertAlmostEqual(t2, t3)
        t.mode = 'linear'
        t2 = t.get_threshold('PP', 11)
        self.assertAlmostEqual(t2, float(11 - 8) / (20 - 8) * (t3 - t1) + t1)
        t.mode = 'log'
        t2 = t.get_threshold('PP', 11)
        self.assertAlmostEqual(t2, math.exp(float(11 - 8) / (20 - 8) * (math.log(t3) - math.log(t1)) + math.log(t1)))

    def test_updating(self):
        t = Thresholds(self._s)
        t.digits = 4
        t.increase_factor = 4.0
        t.update_threshold('PP', 8, 0.5e-7)
        self.assertAlmostEqual(t.get_threshold('PP', 8), 2.0e-7)
        t.update_threshold('T', 8, 0.5e-9)
        self.assertAlmostEqual(t.get_threshold('T', 8), 2.0e-9)
        t.mode = 'linear'
        t.update_threshold('PP', 40, 1.00e-02)
        self.assertAlmostEqual(t.get_threshold('PP', 60), 8.0e-2)

    def test_modifying(self):
        t = Thresholds(self._s)
        t.mode = 'linear'
        t.add_step(1)
        self.assertAlmostEqual(t.get_threshold('PP', 1), t.get_threshold('PP', 3))
        self.assertAlmostEqual(t.get_threshold('PP', 1), t.get_threshold('PP', 2))
        t.add_step(80)
        self.assertAlmostEqual(t.get_threshold('PP', 60), t.get_threshold('PP', 80))
        self.assertAlmostEqual(t.get_threshold('PP', 60), t.get_threshold('PP', 70))
        t.add_step(10)
        self.assertAlmostEqual(t.get_threshold('PP', 10), 1.7e-7)
        t.removeStep(1)
        t.removeStep(10)
        t.removeStep(80)
        self.assertEqual(t, Thresholds(self._s))
        t.add_variable('QV')
        t.increase_factor = 1.0
        t.update_threshold('QV', 3, 1.0e-12)
        t.update_threshold('QV', 8, 1.0e-9)
        t.update_threshold('QV', 20, 5.0e-5)
        t.update_threshold('QV', 60, 1.0e-1)
        self.assertAlmostEqual(t.get_threshold('QV', 15), 2.9167083333333337e-05)
        t.remove_variable('QV')
        self.assertEqual(t, Thresholds(self._s))


if __name__ == "__main__":
    unittest.main()
