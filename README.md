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

Installation of exprorer_msmd
exprorer_msmd is a system for performing mixed-solvent molecular dynamics (MSMD) simulation using GROMACS automatically.
exprorer_msmd is available in https://github.com/keisuke-yanagisawa/exprorer_msmd

Installation of cosmdanalyzer
cosmdanalyzer is a system for hotspot detection form output of exprorer_msmd
cosmdanalyzer is available in GitHub URL.

