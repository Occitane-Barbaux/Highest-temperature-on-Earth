# Code for paper "What is the highest near-surface temperature expected on Earth by 2100?"

This is the companion code for paper "What is the highest near-surface temperature expected on Earth by 2100?" (to be published).


##  Overview

1. Assess observed values and trends in series of annual absolute highest near-surface air temperature on Earth (TXxx) in four datasets.

2. Projections of future TXxx for various emission scenarios + probabilistic assessment of very high values (61.8°C).


## 1. System requirements

### Operating systems
The code has been tested on:
- Ubuntu 22.04.5 LTS

It is expected to run on other modern Linux and Windows systems.

### Software requirements
- Python **3.9–3.13** (Python 3.12.0, tested)
- Climate Data Operators (CDO version 2.0.4 tested)
- Jupyter Notebook

### Python package dependencies
The required Python packages are listed in 'Utilities/requirements.txt' or 'Utilities/environment.yml'

### Hardware requirements
No non-standard hardware is required.  
All analyses can be run on a standard desktop or laptop computer (tested).

### Additional notes
Some Python packages (e.g. 'cartopy') require system-level libraries such as PROJ, GEOS, and ESMF. These dependencies are commonly available via standard package managers and conda distributions.
 
This code relies on [CDO](https://code.mpimet.mpg.de/projects/CDO) and [CmdStan](https://cmdstanpy.readthedocs.io/en/v1.3.0-post2/installation.html#cmdstan-installation) for calculations 

## 2. Installation guide

### Instructions

0. Check for existing installation of [CDO](https://code.mpimet.mpg.de/projects/CDO). (untested, as it's installed by default on every computer here)


1. Clone or download this repository:
   '''bash
   git clone <rhttps://github.com/Occitane-Barbaux/Highest-temperature-on-Earth>
   cd <Highest-temperature-on-Earth>
   '''

2. (Recommended) Create and activate a virtual environment:
   '''bash
    python -m venv venv
    source venv/bin/activate
    '''

3.a Install the required Python packages from your package manager:
   '''bash
   pip install -r Utilities/requirements.txt
   '''

3.b Install [NSSEA](https://github.com/Occitane-Barbaux/NSSEA) and [SDFC](https://github.com/Occitane-Barbaux/SDFC-python) manually from github.

3.c Install [CMDStan](https://cmdstanpy.readthedocs.io/en/v1.3.0-post2/installation.html#cmdstan-installation) by following the specific instructions.

3.d Modify path to your CMDSTan location in *Utilities/Full_run_Model.py*, line 10.

3.e Check package installation in environement by running *Test_Packages_installation.ipynb* once. 

4. Launch Jupyter Notebook:
    '''bash
    jupyter notebook
    '''
    
### Typical installation time

Approximately up to an hour minutes on a normal laptop computer with an existing Python installation.


## 3. Demo

*Rerun_Calculations.ipynb* is provided as a demo of the full algorithm.

It reruns all calculations from raw data to final figure:

1. *Utilities/Data_Preparation_ERA5_025_Restricted_all.py* processes observed and simulated (GCM) data.
2. *Utilities/Full_run_Model.py* applies the algorithm from RR2020 (Simulated prior and 2-steps constraint with observations)
3. *Utilities/Figures.py* produces figures and indexes, including the predictive sampling

Expect a total run time of up to 10 minutes (tested on a laptop)

## 4. Instruction for use

### 4.1 Input data

This study uses:
- Simulated series of annual Tasmean and annual maxima of daily Tasmax from the CMIP6 ensemble. (Extracted from LMD's spirit server with a bash code that took a week, available on demand).
- Observed annual and daily values for Tasmean and Tasmax from ERA5. (Same).
- Observed annual maxima of TXxx from HadEX3. Found [here](https://www.climdex.org/learn/datasets/).
- Observed annual maxima of TXxx from BEST. Found [here, too](https://www.climdex.org/learn/datasets/).
- Observed annual maxima of TXxx from GHCNd. Found [here](https://www.ncei.noaa.gov/data/global-historical-climatology-network-daily/access/).



Due to Github's 100mo file limit, data including Hadex and ERA5 datasets have been excluded and can be found on [10.5281/zenodo.23558536](https://doi.org/10.5281/zenodo.23558535). All other data files are available in folder './Data_origin'.

### 4.2 Organisation

#### Folder ./Utilities
All the defined functions, requirement and STAN code files used in the scripts are saved in the './Utilities' directory.

#### Folder ./Data_origin
All the raw data used in the scripts are saved in the './Data_origin' directory.

Due to Github's 100mo file limit, data including Hadex and ERA5 datasets have been excluded and can be found on [10.5281/zenodo.23558536](https://doi.org/10.5281/zenodo.23558535).

#### Folder ./Figures_Datasets
All the figures produced by *Observations_TXxx.ipynb* are saved in the './Figures_Datasets' directory.

#### Folder ./Rerun_World
All the figures, processed data and intermediates files produced by *Projections_TXxx.ipynb* and *Rerun_Calculations.ipynb* are saved in the './Rerun_World' directory.

### 4.3 Usage of scripts

#### *Projections_TXxx.ipynb*
Calculation of the Instantaneous Return levels and predictive Equivalent Reliability level.

Input: Precalculated projections produced by *Rerun_Calculations.ipynb* + processed data.
##### Beware: To rerun these specific calculations, you need to download the intermediate file *clim_multiscenarios_multisynthesis_allmodels_03261504.nc*  from [Zenodo](https://doi.org/10.5281/zenodo.23558535).  Copy it in folder './Rerun_World/Era5_025/Outputs/World_World/03261504/'.  See Github's file limitations.

Output: Figures 2 and Tables 1, 2 and Supplementary Figures 7, 8 , 9, 10 and 11 and Tables 4, 5

Run time: up to 5 minutes.

#### *Observations_TXxx.ipynb*
Comparison of datasets of the annual absolute highest near-surface air temperature on Earth (TXxx)

Input: Raw observed data
##### Beware: To run these specific calculations, you need to download the data files *ERA5_tasmax_World_allgrid.nc*, *HadEX3_TXx_1901-2018_ADW_61-90_1.25x1.875deg.nc*, *tasmax_World_All_MaxWorld_2016.nc* and *tasmax_World_All_MaxWorld_2017.nc*  from [Zenodo](https://doi.org/10.5281/zenodo.23558535).  Copy it in folder './Data_origin/Other_datasets/' . See Github's file limitations.

Output: Figures 1 and Supplementary 3, 4 and 5.

Run time: up to 10 minutes.

#### *Rerun_Calculations.ipynb*
Redo all projection calculations from raw data.

Input: Raw observed and simulated (GCM) data.

Output: This produces the same figures as *Projections_TXxx.ipynb* (Figures 2 and Tables 1, 2 and Supplementary Figures 7, 8 , 9, 10 and 11 and Tables 4, 5) and adds Supplementary Figure 6 and Table 3 (from data preparation).

Run time: up to 10 minutes.


## 5. Acknowledgments
Portions of the constraint code are based and codeveloped with the Python packages of Yoann Robin and Aurélien Ribes: [NSSEA](https://github.com/yrobink/NSSEA).
Part of this code has been integrated in the following package: [ANKIALE](https://github.com/yrobink/ANKIALE). See [paper](https://gmd.copernicus.org/articles/19/2349/2026/).


## 6. Authors
[Occitane BARBAUX](https://occitane-barbaux.github.io/) (Météo-France, CNRM)

[Aurelien RIBES](https://www.umr-cnrm.fr/spip.php?article23) (Météo-France,CNRM,CNRS,Toulouse University)

[Philippe NAVEAU](https://www.lsce.ipsl.fr/en/pisp/philippe-naveau/) (LSCE,CNRS)

[Nathalie BERTRAND](https://www.irsn.fr/recherche/bureau-dexpertise-hydrogeologie-sur-risques-dinondation-meteorologiques-geotechniques-behrig) (ASNR)

[Margot BADOR](https://scholar.google.com/citations?user=QZ6QvSgAAAAJ&hl=fr)(Cerfacs/CECI,CNRS,Toulouse University)

[Erich FISCHER](https://usys.ethz.ch/en/people/profile.erich-fischer.html)(ETH Zürich Institut für Atmosphäre und Klima)

[Blair TREWIN](https://scholar.google.com/citations?user=R7BQUw8AAAAJ&hl=en)(Australian Bureau of Meteorology)



## 7. Others
"we strongly recommend that you ask at least one colleague who is unfamiliar with your software to install the tool(s), follow the instructions, and provide feedback.": Colleague produced all the figures and is still willing to talk to me after installing all that. Validated.

## 8. License


This software is governed by the CeCILL-C license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL-C
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the 
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also 
therefore means  that it is reserved for developers  and  experienced 
professionals having in-depth computer knowledge. Users are therefore 
encouraged to load and test the software's suitability as regards their 
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the   
same conditions as regards security.                                    

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-C license and that you accept its terms.    



