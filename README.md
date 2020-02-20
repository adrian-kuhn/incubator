Incubator
==============
Python Command line application (CLI) containing a genetic algorithm to find best parameters for arbitrary assignments.
Used to evaluate and parametrize algorithms to detect single trees in LiDAR data.
Realized as a master thesis for the canton of Lucerne to obtain the UNIGIS Master of Science (MSc) degree at the Paris-Lodron University of Salzburg.


![Final workflow](doc/incubator.png)

The `Incubator` is build with classes in Python 3.6 and can be used for the evolutionary
development of any algorithm. Starting point is the incubator, which can be started (`breed`) 
with a concrete instance of an `Individual` (e.g. `WatershedAdaptive`), an arbitrary number of `assignments` and various process parameters.
The assignments in this framework were used for single 
tree extraction from LiDAR data in six test areas. Each algorithm is derived from the base class `Individual` 
and thus contains a static `GENE_POOL` and a number of `Gene` packed 
in a `Chromosome`. The incubator initializes the first generation of individuals 
and allows all individuals to solve the task either in serial or parallel processing. 
The target (the reference) and the result of a task are point data sets (`PointCloud`), which 
can calculate the accuracy by means of a similarity function. From this metric the fitness of an algorithm can be calculated. As in other research on single tree extraction, 
the F1-Score is used to minimize false positives and false negatives equally. 
The higher the F1-score, the more precisely an algorithm has solved the problem. 
The incubator ensures that the individuals develop within a generation through an evolutionary 
step (`evolve`) consisting of selection, chromosome crossover, mating and random mutation.

Purpose
-------
The framework was build to test different algorithms for single tree detection from LiDAR data. 
Every algorithm started with a random set of parameters and had to solve the problem.
With the incubator, the best performing parameters per algorithm could be found.
Finally the best performing algorithm was used to detect single trees in the canton of Lucerne.
(See: https://github.com/adrian-kuhn/tree-detection)

Feel free to implement other algorithms with other assignments and new fitness functions
to use the framework for any genetic process.

Precondition and data
---------------------
To reproduce the genetic algorithm for tree detection, the algorithm needs access to the tiled LiDAR point cloud and the tiled DOM and DTM.
(See the settings file `settings/settings.json`) The data should be tiled with **equal raster schema**.

If you want to evolve algorithms based on R code (e.g. Li2012), ensure to install R on your machine.
R code will be run with the library `rpy2` (See: https://rpy2.readthedocs.io/en/latest/)

Installing
----------
The follwing steps describes the installing process in the IT environment of the canton of Lucerne.
1. Setup a new conda environment and install all required packages according the `setup.py` file.
2. Start the Python 3 conda environment with `"C:\Program Files\ArcGIS\Pro\bin\Python\scripts\proenv.bat"`
3. Insert `pip install git+https://app.geo.lu.ch/gitea/GEO/incubator.git` to install the newest release
4. Use credentials for Gitea

Running
-------
The following steps describes the procedure to run the incubator on Windows:
1. Configure all parameters in the settings file `settings/settings.json`
2. Configure the path to the conda environment with all installed libraries in `run.bat`
2. Start the evolutionary process by running the script `run.bat`


Uninstalling
------------
If you want to uninstall incubator (e.g. in case of an update), use the following steps:
1. Check or set in ArcGIS Pro the environment to the incubator environment
2. Start the Python 3 conda environment with `"C:\Program Files\ArcGIS\Pro\bin\Python\scripts\proenv.bat"`
3. Insert `pip uninstall incubator` to uninstall the newest release
4. Answer with `y` to complete the uninstall process

Releases
--------
| Date         | Version | Notes                                                    |
|--------------|---------|----------------------------------------------------------|
| 31.03.2020   | 1.0.0   | Initial release                                          |