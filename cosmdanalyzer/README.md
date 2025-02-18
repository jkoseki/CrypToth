# Manual for Hotspot Analysis Software Based on Co-solvent Molecular Dynamics Calculation Results
> Chie Motono and Kenichiro Imai, AIST

# File List
* cosmdanalyzer/ : Source Code
    * cosmdanalyzer.py : Execution script
    * pyproject.toml : Package management file for poetry
    * setting.toml : Setting file (default values)
* src/ : Main source code
* data/ : Static data required at runtime
* sample_input/ : Sample input
    * multi : Small sample for Multiple probes

# Installation
## Using poetry
Requires an environment with Python 3.10 or higher and poetry installed.
Execute the following command in the directory containing pyproject.toml:

~~~~~~~~~~~~~~~~
poetry install --no-dev
~~~~~~~~~~~~~~~~
To run the program, use:

~~~~~~~~~~~~~~~~
poetry run python ./cosmdanalyzer.py [options]  
~~~~~~~~~~~~~~~~

## Installing dependencies manually
Requires an environment with Python 3.10 or higher. Download cosmdanalyzer and install the necessary Python packages with the following command:

~~~~~~~~~~~~~~~~
pip3 install rdkit griddataformats plotly tomli
~~~~~~~~~~~~~~~~

# Usage
## Example commands using sample input
### Running in a local environment
Execute the following command in the directory containing cosmdanalyzer.py.

For multiple probes:

~~~~~~~~~~~~~~~~
python cosmdanalyzer.py -s setting.toml ../out/ ../sample_input/multi/ -v
~~~~~~~~~~~~~~~~

For a single probe:

~~~~~~~~~~~~~~~~
python cosmdanalyzer.py -s setting.toml ../out/ ../sample_input/single/ -v
~~~~~~~~~~~~~~~~

The results will be output to the directory ../out/.


## Options

~~~~~~~~~~~~~~~~
usage: cosmdanalyzer.py [-h] [-s SETTING]  [--output_detail] out_dir input_dir

~~~~~~~~~~~~~~~~

* positional arguments:
    * out_dir: Output directory
    * input_dir: Input directory

* options:
    * -h, --help: show this help message and exit
    * -s SETTING, --setting SETTING: Path to the configuration file
    * --output_detail: Not required for CrypToth execution. Outputs detailed information on score trends.
    * -v, --verbose: Displays detailed processing information in standard
    * --fpocket_info FPOCKET_INFO: Not required for CrypToth execution. Path to fpocket output (xxxx_info.txt). Only available if fpocket is executable.
    * --fpocket_pdb FPOCKET_PDB: Not required for CrypToth execution. Path to fpocket output PDB file. Only available if fpocket is executable.
    
## Input Directory
A directory containing xxx_nVH.dx files is recognized as a single system. The system directories should be structured as follows:

For multiple Probe inputs:
~~~~~~~~~~~~~~~~
maxPMAP_xxx_probe ID_nVH.dx 
system00/PMAP_xxx_nVH.dx
system01/PMAP_xxx_nVH.dx
...
~~~~~~~~~~~~~~~~

For a single Probe input:
~~~~~~~~~~~~~~~~
PMAP_xxx_nVH.dx
...
~~~~~~~~~~~~~~~~

* xxx_nVH.dx : Co-solvent occupancy file
* xxx_position_check2.pdb : Trajectory PDB
* fpocket_out : fpocket output directory. Only available if fpocket is executable.

If multiple xxx_nVH.dx files exist, composite calculations of multiple probe spots will be performed. However, xxx_nVH.dx files located in subdirectories of the system directory will be ignored.

## Setting File
The setting file is written in TOML format.
Default configuration file:

~~~~~~~~~~~~~~~~
[clustering]
# algorithm = "single_linkage" or "dbscan" or "mean_shift"
algorithm = "dbscan" 
occupancy = 0.0004
extend    = 0.0
spot_marge_rate = 0.2

[clustering.single_linkage]
threshold = 3.0

[clustering.dbscan]
epsilon = 3.0
min_pts = 7

[clustering.mean_shift]
bandwidth = 3.0

[score]
temperature = 300.0
solvent_radius = 1.4
resolution = 256
fpocket_threshold = 0.0

[score.weight]
gfe            = 1.0
size           = 1.0
protrusion     = 1.0
convexity      = 1.0
compactness    = 1.0
hydrophobicity = 1.0
charge_density = 1.0
flexibility    = 1.0
fpocket        = 1.0
~~~~~~~~~~~~~~~~

* clustering : Hotspot同定のためのクラスタリングの設定値
    * algorithm : The clustering algorithm for hotspots. Choose from "single-linkage," "DBSCAN," or "mean-shift." 
    * occupancy : Threshold for voxel selection. Voxels where (probe occupancy probability / number of heavy atoms in the probe) >= occupancy are selected.
    * extend : Expands the hotspot voxels by the specified Å.
    * spot_marge_rate : [0.0, 1.0] Spots that overlap above the specified ratio in multiple-probe input are treated as the same hotspot.
    
* score : Settings for Hotspot Score Calculation
    * temperature : Absolute temperature (K) used when creating the input trajectory.
    * solvent_radius : Solvent radius.
    * resolution : Approximates tWeights for each score component.
    * fpocket_threshold : [0.0, 1.0] If an fpocket output pocket overlaps with a hotspot above the specified ratio, the hotspot score is assigned accordingly. Only applicable when fpocket is executable.

## Output
When using a system directory (or its parent directory if output), the following files will be generated in the output directory:

* all_info.txt : Scores for all probes (frame count-weighted average)
* spot_probe.toml : Correspondence table between hotspots and probes
* basename (probeID)/ : Directory for each probe
    * basename_info.txt : Score file
    * basename.pml : PyMOL input file
    * basename_out.pdb : PDB file with hotspot spatial configuration added as atom name APOL in the final frame
    * spots : PDB file recording patch atoms (atoms of protein surface) corresponding to each hotspot
    * score_detail.csv : Trends of each score (mean, variance, min, max, first quartile, median, third quartile), only output when the --output_detail option is specified.

### Score File for Each Patch (Hotspot and Corresponding Protein Surface Feature Analysis)
Example output:

~~~~~~~~~~~~~~~~
Patch 1
    score :     394.091
    gfe :   1.337
    size :  388.927
    protrusion :    1.000
    convexity :     0.022
    compactness :   3.307
    hydrophobicity :    -0.930
    charge_density :    -0.015
    flexibility :   0.442
    fpocket :   n/a

Patch 2
...
~~~~~~~~~~~~~~~~

Scores for each patch are output.

### PyMOL Input File
Loading the output basename.pml in PyMOL will visualize the 3D structure of the hotspots.


### Hotspot and Probe Correspondence File
Multiple probes corresponding to a hotspot are recorded in the format:
Hotspot name = Probe1 ID, Probe2 ID, ...

~~~~~~~~~~~~~~~~
Patch 0 = A00
Patch 1 = A01, A00
...
~~~~~~~~~~~~~~~~

### score_detail.csv
Column information:

* spot : Hotspot number
* score_type : Type of score represented in the row
* mean : Mean
* variance : Variance
* 0.0 : Minimum value
* 0.25 : First quartile
* 0.5 : Median
* 0.75 : Third quartile
* 1.0 : Maximum value

Output example

~~~~~~~~~~~~~~~~
spot,score_type,mean,variance,0.0,0.25,0.50,0.75,1.0
1,size,388.9274,19.73281,379.7122,385.784,388.1892,391.6664,398.7536
1,protrusion,1,0,1,1,1,1,1
1,convexity,0.02166458,5.456679e-06,0.01765554,0.02062114,0.02174732,0.02299604,0.0277241
1,compactness,3.307329,0.07657145,2.908801,3.11562,3.170944,3.561596,3.787375
1,charge_density,-0.01541033,4.736121e-06,-0.01944439,-0.01679969,-0.01558867,-0.01391103,-0.01187199
2,size,524.7099,17.29361,519.2607,521.8222,524.3902,526.5038,533.8469
...
~~~~~~~~~~~~~~~~
