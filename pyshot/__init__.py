"""
pyShot
======

Provides

1. Available subpackages
------------------------
load    Functions for loading seismic data into memory.
pick    Functions for picking seismic events in a seismic signal.
cryo    Functions for processing glacial seismic data.
model   Functions for modeling seismic data.
"""

from . import load
#from . import pick
#from . import cryo
#from . import model

# If you want to make certain functions or classes directly accessible
# from pyshot import *
__all__ = ['load', 'pick', 'cryo', 'model']

# Optionally, you can also import specific functions to make them directly accessible
# from .load import some_load_function
# from .pick import some_pick_function
# __all__ += ['some_load_function', 'some_pick_function']