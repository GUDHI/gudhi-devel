'''
Helper to preload windows dlls to prevent dll not found errors.
Once a DLL is preloaded, its namespace is made available to any subsequent DLL.
'''
import os
from ctypes import WinDLL
import glob
if os.name == 'nt':
    # convention for storing / loading the DLL from gudhi/.libs/, if present
    try:
        basedir = os.path.dirname(__file__)
    except:
        pass
    else:
        libs_dir = os.path.abspath(os.path.join(basedir, '.libs'))
        if os.path.isdir(libs_dir):
            for filename in glob.glob(os.path.join(libs_dir, '*dll')):
                WinDLL(os.path.abspath(filename))
