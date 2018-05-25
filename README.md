# SMCM-C 
The Coastal version of the Stochastic Multcloud model

This is the python code of the coastal version of the stochastic multi-cloud-model (SMCM-C).
A in depth describtion of the model can be found in the scientific journal article [Coastal Tropical Convection in a Stochastic Modeling Framework](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017MS001048). 
 
 The code should be running with both python versions 2 and 3. 
 
 ## Pre-Requisits
 You need to have a working [python](http://python.org) 2 or 3 version installed on your computer to use the model. Additionaly you need to have the following external python modules installed.
 
 * [numpy](http://www.numpy.org)
 * [netCDF4](https://unidata.github.io/netcdf4-python)
 * [pandas](http://pandas.pydata.org)
 * [pytz](https://pythonhosted.org/pytz)
 * [matplotlib](https://matplotlib.org) and [basemap](https://matplotlib.org/basemap) for visualisation
 * [timezonefinder](https://pypi.org/project/timezonefinder/#description)
 * [geocoder](https://pypi.org/project/geocoder/#description)
 * [SALib](https://github.com/SALib/SALib)
 
 *SALib* and *geocoder* are only necessary if you want to run a Variance Based Sensitivity Test (SALib) or apply the model to realy world situations.
 
 You can either install the above mentioned additional python packages with the os repository manager (apt, pacman, emerge, port) or via [pip](https://pip.pypa.io/en/stable/installing).

## Running the Model
### The Important Modules
The heart of the model contains two Classes:

 | SMCM (in the model module)      |      Coarsgraining (in the coarsgraining module)|
 | ------                          | ------                                          |
 | Defines the transition-rates    |       Does the coarsgraining and nearest neighbor interaction|
 
 The *Coarsgraining* class inherits from *SMCM* so all you need to do in order to run the model is to create an instance of the 
*Coarsgraning* module. Take also a look at the attached ipython-notebook for details.
### The Model Configuration
The model is configured by the text file *constants.config*. The configuration is saved in the following format:
```
key = value
```
Characters starting with '#' are ignored. You can add as many new keys and values as you like. The *Config* class form the configdir module will automatically detect them and pass the to the SMCM class where they will be made available. The following procedure can be seen as an example.

Suppose you want to add a new variable to the model. The easiest way to achieve that is to add a new entry into the *constants.config* file:
```bash
new_variable = value # A meaningful describtion of what this variable does
```

*New_variable* is then available to the SMCM and Coarsgraining class as a new class instance. If you created the Coarsgraining class with:
```python
>>> CG = Coarsgraining('constants.config', 0.1, 0.4)
```

You can access *new_variable* via:
```python
>>> CG.new_variable
```
### Running an Example 
A short example of how to run the cloud model is available in the [test_run](https://github.com/antarcticrainforest/SMCM-C/blob/master/test_run.ipynb) ipython-notebook

### Running Experiments
The *run_cloudmodel* module contains a some methods to run various experiments, like a sensetivity study. To make use of all features of the modle it is highly recommended to have [mpi4py](http://www.mpi4py.scipy.org/docs) module installed.
