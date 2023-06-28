# ALWAYS INCLUDE THIS CELL TO DOCUMENT PACKAGES USED WHEN RUNNING THIS NOTEBOOK
# Record the python and OS version information for this notebook
#
# import main python libraries
import time as pytime
import numpy as np
import pandas as pd
#
# python visualization packages
import hvplot.pandas
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.animation import FFMpegWriter
from matplotlib import __version__
#
# set path to pyko files
import sys
sys.path.append("../")
import pyko
from pyko import OutputClass as OutputClass
from pyko import RunClass as RunClass
#
# io tools
import pickle
import yaml
import pint
#
# gather version information
from platform import version
print('Platform: ',version())
del version
from platform import python_version
print('python version: ',python_version())
del python_version
print('matplotlib version: ', __version__)
del __version__
from hvplot import __version__
print('hvplot version: ', __version__)
del __version__
print('numpy version: ', np.__version__)
print('pandas version: ', pd.__version__)
print('pickle version: ',pickle.format_version)
print('yaml version: ', yaml.__version__)
print('pint version: ', pint.__version__)
print('pyko version: ', pyko.__version__)
print('')

# speed up loading the unit registry by making a local cache
#ureg = pint.UnitRegistry(cache_folder=":auto:")  
ureg = pint.UnitRegistry()
Q_ = ureg.Quantity
#
# if compiling the fortran version for comparison
import os
print(os.popen('which gfortran').read())
print(os.popen('gfortran --version').read())
