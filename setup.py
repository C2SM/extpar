"""Setup script for the extpar python bindings.
"""

import versioneer
from setuptools import setup
import glob

setup(
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    scripts=glob.glob("bin/*.py"),
    data_files=[(
        # non-python executables
        "bin",
        glob.glob("bin/*.exe") + glob.glob("bin/*.sh") + ["bin/sfmakedepend"],
    )],
)
