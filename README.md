
## Light Attenuation in the Mesopelagic Zone

Code and workflows used to investigate patterns of light attenuation in the mesopelagic ocean using autonomous glider observations and complementary environmental data.

This repository contains MATLAB and R scripts used for data processing, quality control, statistical analyses, and figure generation supporting the manuscript:

> *What we do in the shadows: using light attenuation to observe mesopelagic ecosystems* by Alice Della Penna, Camrin D. Braun, Joanne O'Callaghan, Christophe Guinet, Moninya Roughan, Christina Schallenberg, and Martin C. Arostegui to be submitted soon.

The objective is to quantify variability in underwater light attenuation and evaluate possible drivers.


- **an01-04** include the scripts for Case study 1 - Glider data from Samoa
- **an05-08** include the scripts for Case study 2 - Triaxus data from the East Australian Current
- **an09-an12** include the scripts for Case study 3 - Southern elephant seals and vessels of opportunities in the Kerguelen region, with the exception of the R codes used to define light attenuation from the Wildlife Computers tags

- **R** codes are used to produce light attenuation profiles 

- The other MATLAB codes (read_xxx.m) are functions called by the scripts in an01-12.

Please reach out to alice.penna@auckland.ac.nz if there are any questions. Please keep in mind that these codes are still in development and are distributed 'as is' with no guarantee or support. Please credit us if you use any of these codes.


## Requirements

### MATLAB

Tested with:
- MATLAB R2023b
- Mapping Toolbox
- Statistics and Machine Learning Toolbox

### R

Tested with:
- R 4.4.0

Required packages:
- tidyverse
- lubridate
- ggplot2
- mgcv
- patchwork
