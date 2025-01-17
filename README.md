# CrypToth
CrypToth can explore hotspots on the protein surface favorable to ligand binding using MSMD simulations with six different probes and then identify hotspots corresponding to cryptic sites by assessing the protein's conformational variability using the topological data analysis. 

Koseki J, et al., CrypToth: Cryptic pocket detection through mixed-solvent molecular dynamics simulations based topological data analysis. bioRxiv. 2024 doi: https://doi.org/10.1101/2024.07.10.602991

![image](https://github.com/user-attachments/assets/ea300d6d-c5cf-43e3-a920-1a10667fcd9b)

To run CrypToth, exprprea (https://github.com/keisuke-yanagisawa/exprorer_msmd), cosmdanalyzer (XXX) and DASI (https://github.com/jkoseki/DAIS) must be installed. 



## Execution protocol of CrypToth
### 1.	Environment Setup
First, you need to create a CrypToth directory as below.

**$ mkdir CrypToth**

The following three systems are required to run CrypToth. These systems are available in each GitHub repository. Install each system in a directory with the name of the system.


**Installation of _exprorer_msmd_**

exprorer_msmd is a system for performing mixed-solvent molecular dynamics (MSMD) simulation using GROMACS automatically.
exprorer_msmd is available in https://github.com/keisuke-yanagisawa/exprorer_msmd


**Installation of _cosmdanalyzer_**

cosmdanalyzer is a system for hotspot detection form output of exprorer_msmd
cosmdanalyzer is available in "GitHub URL".

**Installation of _DAIS including FF score calculator_**
Create a conda environment using conda-TDA.yml. After activating the TDA conda environment, install the bio3d, TDA, readr, data.table, tidyr, stringr, kernlab, tidyverse, dplyr, openxlsx, earth, Rtsne, mclust, gplots, and pheatmap packages in R.

Since the DAIS method can be implemented by executing the R script, the following is not required, but is a procedure to be followed in order to use the uploaded bash script.

1. Create a bin directory directly under the Linux home directory.
2. Save DAIS.bash, AUTO-DAIS.bash, mdcrd2pdb_for-TDA.bash, and PDB_for-TDA.bash in the bin directory and grant execution permission (chmod +x).
3. Create a DAIS-Source directory in the bin directory and save the R script in it.
The following assumes the case for parallelized calculations (R-script_parallel computing).


An example of the directory structure is as follows.

![image](https://github.com/user-attachments/assets/65217c06-7f12-40e3-99c8-e1816f32cf50)


### 2.	Execution of CrypToth
#### 2.1    Detection of hotspots based on MSMD simulation
In CrypToth, it is necessary to perform MSMD simulation using 6 different probes (benzene, isopropanol, phenol, imidazole, acetonitrile, and ethylene glycol) to detect hotpots which are candidates of cryptic sites (Perform MSMD simulations for each of the six types of probes).


##### 2.1.0    Making working directory
You need to create two working directories in the CrypToth directory. For convenience, the PDB ID is used for the names of the working directories.

**$ cd CrypToth**

**$ mkdir 2am9 2am9_WAT**

In the “2am9” directory, you need to create directories in which the results of MSMD simulation are saved for each probe. In the “2am9_WAT” directory, the results of MD simulation in water phase are saved. Here, you create a directory with probe ID name for each probe. Probe IDs are defined as below.

- A00: Benzene
- A01: Isopropanol
- A37: Imidazole
- B71: Acetonitrile
- A20: Phenol
- E20: Ethylene glycol


**$ cd 2am9**

**$ mkdir A00 A01 A37 B71 A20 E20**

*Also create WAT directory which is the directory for MSMD without probe (in water phase).


##### 2.1.1    Performing MSMD simulation using exprorer_msmd
**Input file preparation**<br>
Store the following three files in each probe directory. 

- PDB file: e.g. 2am9.pdb <br>
The input PDB file should be preprocessed as necessary.

- Probe molecule file: e.g. A20.mol2 and A20.pdb <br>
These files are created by performing structure optimization and partial charge assignment for the probe using Gaussian 16 software package. For details, please refer to the GitHub repository of _explorer_msmd_ (https://github.com/keisuke-yanagisawa/exprorer_msmd).

- The YAML file defining the protein, probe molecules, and simulation protocol. <br>
For details about this file, please refer to the GitHub repository of _explorer_msmd_ (https://github.com/keisuke-yanagisawa/exprorer_msmd).


The directory structure is as follows.

![image](https://github.com/user-attachments/assets/072bacf0-0984-4d81-9e3d-0d8a817ae032)

Example of YAML files and probe molecule files for CrypToth are available in This page.


**Running _exprorer_msmd_** <br>
You can execute exprorer_msmd as below.
