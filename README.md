# Highest-temperature-on-Earth



##  Content

This is the companion code of paper "What is the highest near-surface temperature expected on Earth by 2100?" (to be published) 

## Main ideas

- *Observations_TXxx.ipynb*: Comparison of datasets of the annual absolute highest near-surface air temperature on Earth (TXxx)
- *Projections_TXxx.ipynb*: Calculation of the Instantaneous Return levels and predictive Equivalent Reliability level


## Necessary packages

This code relies on modified versions of packages [NSSEA](https://github.com/Occitane-Barbaux/NSSEA) and [SDFC](https://github.com/Occitane-Barbaux/SDFC-python) created by [Yoann ROBIN](https://github.com/yrobink/).

### Installation
To create the necessary environement, please:
- Check *cdo* existence.
- Install all python packages in *Utilities/requirements.txt* or *Utilities/environment.yml*.
- Install [NSSEA](https://github.com/Occitane-Barbaux/NSSEA) and [SDFC](https://github.com/Occitane-Barbaux/SDFC-python) manually from github.
- Check package installation in environement by running *Test_Packages_installation.ipynb* once.
- Install [CMDStan](https://cmdstanpy.readthedocs.io/en/v1.3.0-post2/installation.html#cmdstan-installation).
- Modify path to CMDSTan location in *Utilities/Full_run_Model.py*, line 10.

### Missing files
A few data files are missing because they wre over the 100mo github limit. 
You can fin on [10.5281/zenodo.23558536](https://doi.org/10.5281/zenodo.23558535). 


## Authors
[Occitane BARBAUX](https://occitane-barbaux.github.io/) (Météo-France, CNRM)

[Aurelien RIBES](https://www.umr-cnrm.fr/spip.php?article23) (Météo-France,CNRM,CNRS,Toulouse University)

[Philippe NAVEAU](https://www.lsce.ipsl.fr/en/pisp/philippe-naveau/) (LSCE,CNRS)

[Nathalie BERTRAND](https://www.irsn.fr/recherche/bureau-dexpertise-hydrogeologie-sur-risques-dinondation-meteorologiques-geotechniques-behrig) (ASNR)

[Margot BADOR](https://scholar.google.com/citations?user=QZ6QvSgAAAAJ&hl=fr)(Cerfacs/CECI,CNRS,Toulouse University)

[Erich FISCHER](https://usys.ethz.ch/en/people/profile.erich-fischer.html)(ETH Zürich Institut für Atmosphäre und Klima)

[Blair TREWIN](https://scholar.google.com/citations?user=R7BQUw8AAAAJ&hl=en)(Australian Bureau of Meteorology)


