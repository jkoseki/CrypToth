# CrypToth
CrypToth can explore hotspots on the protein surface favorable to ligand binding using MSMD simulations with six different probes and then identify hotspots corresponding to cryptic sites by assessing the protein's conformational variability using the topological data analysis. 

Koseki J, et al., CrypToth: Cryptic pocket detection through mixed-solvent molecular dynamics simulations based topological data analysis. bioRxiv. 2024 doi: https://doi.org/10.1101/2024.07.10.602991

![image](https://github.com/user-attachments/assets/ea300d6d-c5cf-43e3-a920-1a10667fcd9b)

To run CrypToth, exprprea (https://github.com/keisuke-yanagisawa/exprorer_msmd), cosmdanalyzer (XXX) and DASI (https://github.com/jkoseki/DAIS) must be installed. 



## Execution protocol of CrypToth
### 1.	Environment Setup
First, you need to create a CrypToth directory as below.

**`$ mkdir CrypToth`**

The following three systems are required to run CrypToth. These systems are available in each GitHub repository. Install each system in a directory with the name of the system.


**Installation of _exprorer_msmd_** <br>
exprorer_msmd is a system for performing mixed-solvent molecular dynamics (MSMD) simulation using GROMACS automatically.
exprorer_msmd is available in https://github.com/keisuke-yanagisawa/exprorer_msmd


**Installation of _cosmdanalyzer_** <br>
cosmdanalyzer is a system for hotspot detection form output of exprorer_msmd
cosmdanalyzer is available in this repository.

**Installation of _DAIS including FF score calculator_** <br>
Create a conda environment using conda-TDA.yml. After activating the TDA conda environment, install the bio3d, TDA, readr, data.table, tidyr, stringr, kernlab, tidyverse, dplyr, openxlsx, earth, Rtsne, mclust, gplots, and pheatmap packages in R.

Since the DAIS method can be implemented by executing the R script, the following is not required, but is a procedure to be followed in order to use the uploaded bash script.

1. Create a bin directory directly under the Linux home directory.
2. Save DAIS.bash, AUTO-DAIS.bash, mdcrd2pdb_for-TDA.bash, and PDB_for-TDA.bash in the bin directory and grant execution permission (chmod +x).
3. Create a DAIS-Source directory in the bin directory and save the R script in it.
4. Set path in the directry, ~/bin.

An example of the directory structure is as follows.

![image](https://github.com/user-attachments/assets/65217c06-7f12-40e3-99c8-e1816f32cf50)


### 2.	Execution of CrypToth
#### 2.1    Detection of hotspots based on MSMD simulation
In CrypToth, it is necessary to perform MSMD simulation using 6 different probes (benzene, isopropanol, phenol, imidazole, acetonitrile, and ethylene glycol) to detect hotpots which are candidates of cryptic sites (Perform MSMD simulations for each of the six types of probes).


##### 2.1.0    Making working directory
You need to create two working directories in the CrypToth directory. For convenience, the PDB ID is used for the names of the working directories.

**$ cd CrypToth**

**$ mkdir 2am9 2am9_WAT**

In the “2am9” directory, you need to create directories in which the results of MSMD simulation are saved for each probe. 
In the “2am9_WAT” directory, the results of MD simulation in water phase are saved. Here, you create a directory with probe ID name for each probe. Probe IDs are defined as below.

- A00: Benzene
- A01: Isopropanol
- A20: Phenol
- A37: Imidazole
- B71: Acetonitrile
- E20: Ethylene glycol


**$ cd 2am9**

**$ mkdir A00 A01 A20 A37 B71 E20**

*Also create WAT directory which is the directory for MSMD without probe (in water phase).


##### 2.1.1    Performing MSMD simulation using exprorer_msmd
**Input file preparation**<br>
Store the following three files in each probe directory. 

- PDB file: e.g. 2am9.pdb <br>
The input PDB file should be preprocessed as necessary.

- Probe molecule file: e.g. A20.mol2 and A20.pdb <br>
These files are created by performing structure optimization and partial charge assignment for the probe using Gaussian 16 software package.
For details, please refer to the GitHub repository of _explorer_msmd_ (https://github.com/keisuke-yanagisawa/exprorer_msmd).

- The YAML file defining the protein, probe molecules, and simulation protocol. <br>
For details about this file, please refer to the GitHub repository of _explorer_msmd_ (https://github.com/keisuke-yanagisawa/exprorer_msmd).


The directory structure is as follows.

![image](https://github.com/user-attachments/assets/072bacf0-0984-4d81-9e3d-0d8a817ae032)

Example of YAML files and probe molecule files for CrypToth are available in This page.


**Running _exprorer_msmd_** <br>
You can execute exprorer_msmd as below.

**$ cd ../exprorer_msmd**

**$ ./exprorer_msmd ../2am9/A20/msmd_protocol_A20.yaml**

In this step, 20 runs of 40 ns MSMD simulation were executed. After the execution, results are saved into “output” directory in the A20 directory.

Then voxel file in OpenDX format which is necessary to calculation of probe occupancy for hotspot detection is generated based on trajectories of 20 runs of MSMD simulation. The voxel file can be generated using the following command.

**$ ./protein_hotspot ../ /2am9/A20 /msmd_protocol.yaml**

maxPMAP_2am9_A20_nVH.dx file is generated in the “output” directory.


**Running exprorer_msmd without probe molecules** <br>
In CrypToth, MD simulation in water phase is also necessary to DAIS analysis. For the MD simulation, protein structure obtained in trajectory at 20 ns of the MSMD simulation ifs used as input PDB file (initial structure). 
For DAIS analysis, 5 runs of the MD simulation should be performed for each probe. The input PDB files for the 5 MD simulations can be prepare as below.


- Open the each “2am9_A20_woWAT_500ps.pdb” files in the “system0”-“system4” directory of each “output” directory.
- Then “MODEL 1” structure is picked out (extract only atoms of the protein molecule) from each 2am9_A20_woWAT_500ps.pdb and save it as the input PDB files.
- For convenience, the name of the pdb files are set to “2am9_A20_0.pdb, 2am9_A20_1.pdb, 2am9_A20_2.pdb, 2am9_A20_3.pdb and 2am9_A20_4.pdb”, respectively.
- make A20_0, A20_1, A20_2, A20_3 and A20_4 directories in WAT directory and store each file in correspondence directories.


**$ cd 2am9_WAT**

**$ mkdir A20_0 A20_1 A20_2 A20_3 A20_4**

YAML file is also necessary. Example of YAML files for the MD simulation (e.g. msmd_protocol_WAT_A20.yaml) is available in this page. It needs to prepare dummy probe files since exprorer_msmd is used for the MD simulation. For convenience, A20.mol2 and A20.pdb are used as dummy probe files. Those files are also stored in WAT directory.

![image](https://github.com/user-attachments/assets/4da39b78-4366-4fcf-aff9-749c105147d0)

execute exprorer_msmd without probe molecules as below.

**$ ./exprorer_msmd ../2am9_WAT/A20_0/msmd_protocol_WAT.yaml**

“2am9_A20_woWAT_10ps.pdb” file is generated in the “system0” directory of the “output” directory of A20_0 directories. This file is necessary to DAIS analysis.

Repeat the same process for the remaining five probes (A00, A01, A37, B71 and E20).


##### 2.1.2    Detection of hotspots based on the results of _exprorer_msmd_
_cosmdanalyzer_ can generate hotspot files showing hotspot position and amino acids contacting hotspots based on the maxPMAP_2am9_probe ID_nVH.dx obtained from each MSMD simulations with the probe (probe ID is A00, A01, A20, A37, B71 and E20, respectively). A setting file (setting.toml) is necessary to execute cosmdanalyzer. Example of setting.toml is available in this page. Then the setting.toml file is stored into cosmdanalyzer directory.


**Running cosmdanalyzer** <br>
You can execute exprorer_msmd as below.

**$ cd ../cosmdanalyzer**

**$ python cosmdanalyzer.py -s setting.toml ../output/ ../2am9/ -v**

After execution of cosmdanalyzer, spot_probe.toml file and A00, A01, A20, A37, B71 and E20 directories are generated in the “output” directory. 
The correspondence between spots and probes is described in the spot_probe.toml file. In the A00 directory, A00_out.pdb file and spots directory are generated. 
A00_out.pdb is a file that shows the location of the hotspot in the protein structure in final frame of the MSMD simulation with the probe. 
In the spots directory, spot PDB files (spots1.pdb, spots2.pdb, …) show atoms of amino acids contacting to each hotspot. 
The structures of proteins in the final frame of MSMD simulation do not differ greatly among probes, so the A00 spots data is used as input data of DAIS analysis.


#### 2.2    Ranking of hotspots for detection of Crypitc site
##### 2.2.1    Preparation for DAIS
Before running DAIS, the molecular dynamics simulations to be compared must be completed.

Create a directory to be used for calculation, and save the structure that will become Contl. and the structure that will become Case(s) in the directory.

After converting the trajectory to a single PDB, delete the first two structures where the initial and final structures are stored, then execute the following command to convert the PDB for DAIS.

**$ PDB_for-TDA.bash**

The following command can then be executed to create a directory structure in which DAIS calculations can be performed.

**$ DAIS.bash**

By executing these, the following directory structure is created.
For example, an example DAIS directory structure for the target protein 1FVR using the A00 solvent molecule is shown below.

![image](https://github.com/user-attachments/assets/11d5d4d2-9667-4cad-ad49-e4a5adc3204a)

This completes the preparation.


##### 2.2.2    Perform DAIS calculations
The following the R script can be used to calculate the FF Scores.

**$ Rscript Estimate-Average-Scores.R**

For each DAIS result, Spot_Scores.csv is written.
There are five Tot.Counts for each Probe molecule, and the average of the counts at each execution is taken, and then the mean and standard deviation of the five count means are calculated.

Then, the enforcement times that are more than a standard deviation away from the average of the mean values are excluded as outliers, and the final FF Score can be obtained by taking the average of the FF Score at each spot.
