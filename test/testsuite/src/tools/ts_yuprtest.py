#!/usr/bin/env python
"""
COSMO TECHNICAL TESTSUITE

Python class to wrap around YUPRTEST files and allow the owner
to use it as an iterator in order to read the data line-by-line

c = Compare('YUPRTEST.ref', 'YUPRTEST.ref', 'THRESHOLDS')
c.compare_data()
c.thresholds.mode = "const"
c.thresholds.increase_factor = 1.0
c.update_thresholds()
c.compare_data()
c.print_results()
"""

# built-in modules
import unittest
from tempfile import NamedTemporaryFile
from itertools import izip
import os, math

# other modules
from ts_thresholds import Thresholds

# information
__author__ = "Oliver Fuhrer, Santiago Moreno"
__email__ = "cosmo-wg6@cosmo.org"
__maintainer__ = "oliver.fuhrer@meteoswiss.ch"

# get name of myself
myname = os.path.basename(__file__)
header = myname + ': '


def column(matrix, i):
    """return a specific column from a list of lists (matrix)"""
    return [row[i] for row in matrix]


class Yuprtest(object):
    """class to wrap around a YUPRTEST file"""

    def __init__(self, filename):
        self._filename = filename  # name of associated YUPRTEST file
        self._file = None  # file handle
        self._raw = []  # raw data
        self._data = []  # processed data
        self._headerlines = 0  # number of header lines
        self._lineno = 0  # current line number (for iterator)
        self.__read_data()

    def __iter__(self):
        return self

    def __read_data(self):
        """read data in file and do some basic processing"""
        lineno = 0
        header = True
        self._file = open(self._filename)  # open file
        for line in self._file:  # read data and parse line
            lineno += 1
            data = self.__parse_line(line, lineno)
            if data:
                header = False
                self._raw.append(line.strip())
                self._data.append(data)
            else:
                if not header:
                    raise IOError('Parse error on line ' + str(lineno))

    def __parse_line(self, line, lineno):
        data = line.strip().split()
        if (len(data) == 0):  # check for zero length lines
            self._headerlines += 1
            return None
        if (data[0][0] == "#"):  # remove comment lines
            self._headerlines += 1
            return None
        if (len(data) != 10):
            raise ValueError('Strange record found in ' + self.filename_ +
                             ' on line number ' + str(lineno))
        data.pop(8)  # remove j-position of maximum
        data.pop(7)  # remove i-position of maximum
        data.pop(5)  # remove j-position of minimum
        data.pop(4)  # remove i-position of minimum
        # entry: var step level minval maxval meanval
        data = [
            data[0],
            int(data[1]),
            int(data[2]),
            float(data[3]),
            float(data[4]),
            float(data[5])
        ]
        return data

    def __check_timesteps(self):
        """Check that entries per variable is multiple of # timesteps and that
           all timesteps are present"""
        steps = column(self._data, 1)
        counts = [steps.count(step) for step in set(steps)]
        if not len(set(counts)) == 1:
            raise ValueError('Variables in file contain different timesteps')

    def __str__(self):
        """pretty print thresholds"""
        return (self._data)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self._data == other._data
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    @property
    def variables(self):
        return sorted(list(set(column(self._data, 0))))

    @property
    def steps(self):
        return sorted(list(set(column(self._data, 1))))

    @property
    def levels(self):
        return sorted(list(set(column(self._data, 2))))

    @property
    def data(self):
        return self._data

    def __next__(self):
        return self.next()

    def next(self):
        if self._lineno < len(self._data):
            self._lineno += 1
            return self._data[self._lineno - 1]
        else:
            raise StopIteration()

    def getline(self):
        for i in self._data:
            yield i

    def skip(self, step):
        pass


class YuprLine(object):

    def __init__(self, status, diff, var, step, level, thresh, pos):
        self._status = status
        self._diff = diff
        self._var = var
        self._step = step
        self._level = level
        self._thresh = thresh
        self._pos = pos

    @property
    def pos(self):
        return self._pos

    @property
    def thresh(self):
        return self._thresh

    @property
    def level(self):
        return self._level

    @property
    def step(self):
        return self._step

    @property
    def var(self):
        return self._var

    @property
    def diff(self):
        return self._diff

    @property
    def status(self):
        return self._status


class Compare(object):
    """class to compare two YUPRTEST files using a thresholds object"""

    def __init__(self, filename1, filename2, thresh):
        self._filename1 = filename1
        self._filename2 = filename2
        self._yu1 = Yuprtest(filename1)
        self._yu2 = Yuprtest(filename2)
        self._threshold = Thresholds(thresh)

    @property
    def thresholds(self):
        return self._threshold

    def __compute_difference(self, ref, value):
        """calculate difference between to threshold values depending on the given minval"""
        diff = 0.0
        if self._threshold.minval < 0.0:
            diff = abs(value - ref)  # absolute difference
        elif abs(ref) > self._threshold.minval:
            diff = abs((value - ref) / ref)  # relative difference
        return diff

    def __compare_values(self, var, timestep, ref, value):
        """gets the difference of two values based on which a status code is returned"""
        diff = self.__compute_difference(ref, value)
        if self._mode == "compare":
            thresh = self._threshold.get_threshold(var, timestep)
            if diff == 0.0:
                return 0, diff, thresh  # MATCH
            if diff <= thresh:
                return 1, diff, thresh  # OK
            else:
                return 2, diff, thresh  # FAIL
        elif self._mode == "update":
            if (diff > 0.0):
                self._threshold.update_threshold(var, timestep, diff)
            return 0, None, None

    def __compare_entry(self, ref, data):
        """compares min, max and mean of two YUPRTEST files and
        then returns everything needed to evaluate if the test succeeded"""
        (var, step, level, minval1, maxval1, meanval1) = ref
        (var2, step2, level2, minval2, maxval2, meanval2) = data
        if (var, step, level) != (var2, step2, level2):
            raise ValueError(
                'Non-matching data entries cannot be compared on line' +
                self._lineno)
        (status1, diff1,
         thresh) = self.__compare_values(var, step, minval1, minval2)
        (status2, diff2,
         thresh) = self.__compare_values(var, step, maxval1, maxval2)
        (status3, diff3,
         thresh) = self.__compare_values(var, step, meanval1, meanval2)
        status = max([status1, status2, status3])
        diff = max([diff1, diff2, diff3])
        pos = ["minimum", "maximum",
               "mean"][[diff1, diff2, diff3].index(min([diff1, diff2, diff3]))]
        return YuprLine(status, diff, var, step, level, thresh, pos)

    def __update_status(self, yupr_line):
        rel = yupr_line.diff
        if yupr_line.thresh > 0.0:
            rel = yupr_line.diff / yupr_line.thresh
        step = str(yupr_line.step)
        var = yupr_line.var

        # update status for this timestep
        if not step in self._status.keys():
            self._status[step] = yupr_line.status
        self._status[step] = max(self._status[step], yupr_line.status)
        # update maxdiff for this timestep and variable
        if not step in self._maxdiff.keys():
            self._maxdiff[step] = {}
        if not var in self._threshold.variables:
            var = [yupr_line.var, "*"]
        else:
            var = [yupr_line.var]
        for x in var:
            if not x in self._maxdiff[step].keys():
                self._maxdiff[step][x] = [0, -float('Inf'), -float('Inf')]
            if yupr_line.status >= self._maxdiff[step][x][
                    0] and rel > self._maxdiff[step][x][1]:
                self._maxdiff[step][x] = [
                    yupr_line.status, rel, yupr_line.diff, yupr_line.level,
                    yupr_line.thresh, yupr_line.pos
                ]

    def print_results(self):
        """print results"""
        print('%5s  %s  %7s  %s' % ('nt', '  '.join([
            '%9s' % x for x in self._threshold.variables + ['other']
        ]), 'status', 'reason'))
        for step in [str(x) for x in self._yu1.steps]:
            vals = [
                self._maxdiff[step][x][2]
                for x in self._threshold.variables + ['*']
            ]
            stat = ["MATCH", "OK", "FAIL"][self._status[step]]
            reason = 'none'
            if stat == "FAIL":
                # generate reason
                for var in self._maxdiff[step].keys():
                    if var == "*":
                        continue
                    m = self._maxdiff[step][var]
                    if 'mmax' in locals():
                        if m[0] >= mmax[0] and m[1] > mmax[1]:
                            mmax = m + [var]
                    else:
                        mmax = m + [var]
                if 'mmax' in locals() and mmax[0] > 1:
                    reason = str(mmax[5]) + ' of ' + str(mmax[-1]) + ' on level ' + str(mmax[3]) \
                            + ' (%9.2e' % mmax[2] + ' > ' + str(mmax[4]) + ')'
            print('%5d  %s  %7s  %s' %
                  (int(step), '  '.join(['%9.2e' % x
                                         for x in vals]), stat, reason))

    def compare_data(self):
        """compare two yu files line by line and return the highest error"""
        self._mode = "compare"
        self._lineno = 0
        self._maxdiff = {}
        self._status = {}
        for x, y in izip(self._yu1.getline(), self._yu2.getline()):
            yupr_line = self.__compare_entry(x, y)
            self.__update_status(yupr_line)
            self._lineno += 1
        stat = max(self._status.values())
        # fix thresholds variables which were not encountered
        for step in [str(x) for x in self._yu1.steps]:
            if step not in self._maxdiff.keys():
                print(
                    "WARNING: Not enough reference data, comparison only until max Time steps reference."
                )
                print(" Time steps reference:  {steps}".format(
                    steps=str(self._yu2.steps)))
                print(" Time steps comparison: {steps}".format(
                    steps=str(self._yu1.steps)))
                break
            for var in self._threshold.variables + ["*"]:
                if not var in self._maxdiff[step].keys():
                    self._maxdiff[step][var] = [
                        float('NaN'),
                        float('NaN'),
                        float('NaN'),
                        float('NaN')
                    ]
        return stat

    def reset_thresholds(self):
        """Reset the thresholds in a loaded file"""
        self._threshold._set_thresholds_to_zero()

    def update_thresholds(self):
        """Updates the thresholds of the corresponding threshold file"""
        # Note: The effect update is done in the __compare_values
        self._mode = "update"
        self._lineno = 0
        for x, y in izip(self._yu1.getline(), self._yu2.getline()):
            self.__compare_entry(x, y)
            self._lineno += 1

        # Set the default threshold to the maximum of all the variables
        self._threshold.update_default_thresholds()

    def write_threshold_to_file(self, file_location):
        self._threshold.to_file(file_location)


class Test(unittest.TestCase):

    def setUp(self):
        self._s = """
#    Experiment:  COSMO-Model        Number:    907   SUN 22.01.2012  00:00:00 UTC +          0 H  (SUN 22.01.2012  00:00:00 UTC)
#    ie_tot =   80   je_tot =   70   ke =   60
#
#    var    nt  lev                         min imin jmin                         max imax jmax                        mean
       U     0    1   3.5660915374755859375E+00    9    1   1.5401052474975585938E+01   77   22   8.8077939717442301770E+00
       U     0    2   2.9315328598022460938E+00   32   12   1.3403596878051757812E+01   77   22   8.8707148721225124177E+00
       U     0    3  -1.3106002807617187500E+00   35    8   1.3883832931518554688E+01   79   70   6.5442065665390272144E+00
       U     0    4  -1.6478111598711162600E+00   39    6   1.4910346031188964844E+01   80   70   6.3009189136502516959E+00
       U     0    5   1.2602137718068048257E+00   45    5   1.8625530242919921875E+01   52   68   9.7023587573987377652E+00
       V     0    1  -2.5782619961955823840E+01   75   42  -7.3655395507812500000E+00   52   28  -1.7147092250796632129E+01
       V     0    2  -2.5669403076171875000E+01   78   38  -8.3037261962890625000E+00   35   11  -1.7240270624892069407E+01
       V     0    3  -2.5015502929687500000E+01   50   70  -5.9637451171875000000E+00   41    9  -1.6173610956714036035E+01
       V     0    4  -2.5917953491210937500E+01   45   61  -4.2434954891346858830E+00   42    5  -1.7750392476685011189E+01
       V     0    5  -2.4031661987304687500E+01   78   54  -5.0756072998046875000E+00   46    2  -1.5936751603319409654E+01
     TKE     0    1   1.0000000000000000208E-02   37    1   1.6331740078276316130E-01   71   34   5.1854393734733487953E-02
       U    10    1   3.5660915374755859375E+01    9    1   1.5401052474975585938E+01   77   22   8.8077939717442301770E+01
       U    10    2   2.9315328598022460938E+01   32   12   1.3403596878051757812E+02   77   22   8.8707148721225124177E+01
       U    10    3  -1.3106002807617187500E+01   35    8   1.3883832931518554688E+02   79   70   6.5442065665390272144E+01
       U    10    4  -1.6478111598711162600E+01   39    6   1.4910346031188964844E+02   80   70   6.3009189136502516959E+01
       U    10    5   1.2602137758068048257E+01   45    5   1.8625530242919921875E+02   52   68   9.7023587573987377652E+01
       V    10    1  -2.5782619961955823840E+02   75   42  -7.3655395507812500000E+01   52   28  -1.7147092250796632129E+02
       V    10    2  -2.5669403076171875000E+02   78   38  -8.3037261962890625000E+01   35   11  -1.7240270624892069407E+02
       V    10    3  -2.5015502929687500000E+02   50   70  -5.9637451171875000000E+01   41    9  -1.6173610956714036035E+02
       V    10    4  -2.5917953491210937500E+02   45   61  -4.2434954891346858830E+01   42    5  -1.7750392476685011189E+02
       V    10    5  -2.4031661987304687500E+02   78   54  -5.0756072998046875000E+01   46    2  -1.5936751603319409654E+02
     TKE    10    1   1.0000000000000000208E-03   37    1   1.6331740078276316130E-02   71   34   5.1854393734733487953E-03
"""
        f = NamedTemporaryFile(mode='w+b', delete=False)
        self._filename1 = f.name
        f.write(self._s)
        f.close()
        self._s = """
#    Experiment:  COSMO-Model        Number:    907   SUN 22.01.2012  00:00:00 UTC +          0 H  (SUN 22.01.2012  00:00:00 UTC)
#    ie_tot =   80   je_tot =   70   ke =   60
#
#    var    nt  lev                         min imin jmin                         max imax jmax                        mean
       U     0    1   3.5660915374755859335E+00    9    1   1.5401052474975585918E+01   77   22   8.8077939717442301170E+00
       U     0    2   2.9315328598022460948E+00   32   12   1.3403596878051757812E+01   77   22   8.8707148721225124177E+00
       U     0    3  -1.3106002807617187550E+00   35    8   1.3883832931518554618E+01   79   70   6.5442065665390272144E+00
       U     0    4  -1.6478111598711162620E+00   39    6   1.4910346031188964814E+01   80   70   6.3009189136502516159E+00
       U     0    5   1.2602137758068048247E+00   45    5   1.8625530242919921815E+01   52   68   9.7023587573987377152E+00
       V     0    1  -2.5782619961955823850E+01   75   42  -7.3655395507812500010E+00   52   28  -1.7147092250796632129E+01
       V     0    2  -2.5669403076171875030E+01   78   38  -8.3037261962890625010E+00   35   11  -1.7240270624892069417E+01
       V     0    3  -2.5015502929687500030E+01   50   70  -5.9637451171875000010E+00   41    9  -1.6173610956714036015E+01
       V     0    4  -2.5917953491210937560E+01   45   61  -4.2434954891346858810E+00   42    5  -1.7750392476685011119E+01
       V     0    5  -2.4031661987304687570E+01   78   54  -5.0756072998046875010E+00   46    2  -1.5936751603319409614E+01
     TKE     0    1   1.0000000000000020208E-02   37    1   1.6331740078276316130E-01   71   34   5.1854393734733487253E-02
       U    10    1   3.5660915374755839375E+01    9    1   1.5401052474975515938E+01   77   22   8.8077939717442301370E+01
       U    10    2   2.9315328598022440938E+01   32   12   1.3403596878051717812E+02   77   22   8.8707148721225124377E+01
       U    10    3  -1.3106002807617157500E+01   35    8   1.3883832931518514688E+02   79   70   6.5442065665390272344E+01
       U    10    4  -1.6478111598711162600E+01   39    6   1.4910346031188914844E+02   80   70   6.3009189136502516359E+01
       U    10    5   1.2602137758068078257E+01   45    5   1.8625530242919911875E+02   52   68   9.7023587573987377652E+01
       V    10    1  -2.5782619961955883840E+02   75   42  -7.3655395507812510000E+01   52   28  -1.7147092250796332129E+02
       V    10    2  -2.5669403076171855000E+02   78   38  -8.3037261962890615000E+01   35   11  -1.7240270624892369407E+02
       V    10    3  -2.5015502929687540000E+02   50   70  -5.9637451171875010000E+01   41    9  -1.6173610956714336035E+02
       V    10    4  -2.5917953491210937500E+02   45   61  -4.2434954891346818830E+01   42    5  -1.7750392476685311189E+02
       V    10    5  -2.4031661987304627500E+02   78   54  -5.0756072998046815000E+01   46    2  -1.5936751603319309654E+02
     TKE    10    1   1.0000000000010000208E-03   37    1   1.6331740078271316130E-02   71   34   5.1854393734333487953E-03
"""
        f = NamedTemporaryFile(mode='w+b', delete=False)
        self._filename2 = f.name
        f.write(self._s)
        f.close()
        self._t = """
 minval = 1e-12
  steps =          3          8         20         60
      * =   1.00e-15   1.00e-15   1.00e-15   1.00e-15
    TKE =   1.00e-15   1.00e-15   1.00e-15   1.00e-15
"""

    def tear_down(self):
        pass

    def test_basic_setup(self):
        y1 = Yuprtest(self._filename1)
        y2 = Yuprtest(self._filename1)
        self.assertEqual(y1, y2)
        c = Compare(self._filename1, self._filename1, self._t)

    def test_properties(self):
        y = Yuprtest(self._filename1)
        self.assertEqual(y.variables, [
            "TKE",
            "U",
            "V",
        ])
        self.assertEqual(y.steps, [0, 10])
        self.assertEqual(y.levels, [1, 2, 3, 4, 5])

    def test_iterator(self):
        y = Yuprtest(self._filename1)
        count = 0
        for data in y:
            count += 1
            if count == 4:
                ref = [
                    'U', 0, 4, -1.6478111598711163, 14.910346031188965,
                    6.300918913650252
                ]
                self.assertEqual(data, ref)
        self.assertEqual(count, 22)

    def test_generator(self):
        y = Yuprtest(self._filename1)
        count = 0
        for data in y.getline():
            count += 1
            if count == 4:
                ref = [
                    'U', 0, 4, -1.6478111598711163, 14.910346031188965,
                    6.300918913650252
                ]
                self.assertEqual(data, ref)
        self.assertEqual(count, 22)

    def test_compare(self):
        c = Compare(self._filename1, self._filename2, self._t)
        c.compare_data()
        c.print_results()

    def test_update(self):
        for mode in ["const", "linear", "log"]:
            c = Compare(self._filename1, self._filename2, self._t)
            c.thresholds.mode = mode
            c.thresholds.increase_factor = 2.0
            self.assertEqual(c.thresholds[0][0], 1.0000000000000001e-15)
            c.update_thresholds()
            self.assertEqual(c.compare_data(), 1)
            self.assertEqual(c.thresholds[0][0], 6.9999999999999998e-09)
            c.print_results()


if __name__ == "__main__":
    unittest.main()
