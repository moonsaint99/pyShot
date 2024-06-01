<h1>pyShot</h1>

pyShot is a python package including tools for processing active source seismic shot records, including loading, picking, reflectivity analysis, and Wiechert-Herglotz inversion. The flagship tool in the pyShot toolbox is the GUI-based picker, which allows the user to use their mouse cursor to hover over seismic shot records and make amplitude picks that snap to peaks. The picker tool is built on matplotlib and ObsPy, making it easily extensible using existing functionalities of these two established python packages. Picks are exported to CSV format, which allows them to be assimilated straightforwardly into any data analysis pipeline.

Dependencies:
* [NumPy](https://numpy.org/)
* [SciPy](https://scipy.org/)
* [matplotlib](https://matplotlib.org/)
* [ObsPy](https://docs.obspy.org/)

pyShot is a python 3 module divided into several packages: 
* `load` for loading and preparing active source seismic data for picking
* `pick` for the picking tool and associated functionality
* `cryo` for functions particularly applicable to glacial active source seismology
* `model` for seismic reflectivity modeling tools
