###Getting started

PopNet is available at: https://github.com/ParkinsonLab/PopNet

Requirements:
	Operating System:
		Ubuntu 14.04 or equivalent Linux-based OS

	Dependencies:
		Python 3.6.5 or later
		Python packages:
			- Numpy
			- Matplotlib
			- Scikit-learn 
			- Pandas
		MCL (available from http://micans.org/mcl/)
	
	Visualization:
		Cytoscape 2.7.0 or later (available from www.cytoscape.org/)
		enhancedGraphics plugin for Cytoscape (available from http://apps.cytoscape.org/apps/enhancedgraphics)

Quick Start:
	Clone contents to local folder
	Edit config file
	cd /path/to/popnet/folder
	python PopNet.py /path/to/config/file

###Introduction

PopNet is a set of Python scripts that generates a XGMML file to be visualized in Cytoscape as well as a json file that can be viewed on popnetd3. 
The scripts are run by executing the runner script, PopNet.py, with the path to a configuration file as argument. All configurable
options of the program are contained in the configuration file. The output XGMML file need to be manually loaded
into cytoscape for visualization. Information from several steps of the program are written to other files in the
output folder, which may be useful for extracting particular bits of information, optimizing parameters, or debugging

###Input Requirements

PopNet now expects a tab separated SNP table, formatted as:

\#CHROM	POS	SAMPLE1	SAMPLE2	...
XX_CHRI	1	A	T	...
XX_CHRI	2	G	C	...

PopNet expects that all individual genomes to be aligned to a common reference. In other words, that the coordinates
in each genome refers to the same location. In addition, all chromosome names should be in the format:
	XX_ChrI
	XX_ChrXIV

###Configuration File

PopNet requires a configuration file that includes the output directory, input path, clustering parameters, segment length,
and whether to use the autogroup feature.

Requirements:
output_directory=/path/to/outputdirectory
input_path=/path/to/input/file
reference=name_of_reference #this is relevant only if the first sample in your SNP table is called REF. Otherwise put 'None'
section_length=(integer) #between 1000 and 1000000, inclusive
S1_iVal=(float) #between 1 and 20. Should be 8 unless you're sure.
S1_piVal=(float) #between 1 and 20. Should be 19 unless you're sure.
S2_iVal=(float) #between 1 and 20. Default is 4.
S2_piVal(float) #between 1 and 20. Default is 1.5.
autogrouop=(True/False) #the autogroup feature will ignore your S2_iVal/piVal and determine it based on Silhouette coefficient. 

###Output

Output files are generated at the location specified in the config file, and are always named the same way.

The key file for network visualization is located at:
	/path/to/results/cytoscape/cytoscapeGenome.xgmml

Other potentially helpful output include:
	Heatmaps.pdf               Metrics to help you decide on the Inflation (I) and Pre-inflation (pI) parameters for 
	                           secondary clustering

	log.txt                    A log file that includes all the run parameters

	Genome_nexus.nex           A neighbor-net of the population

	groups.txt.mci             The raw matrix used in secondary clustering. You can quickly try out the effects of
	                           different I and pI values.

	persistentResult.txt       The clustering results of each individual chromosome segment during primary clustering

	results.txt                A tab-delimited file containing all the SNPs used

	/cytoscape/tabNetwork.tsv  A tab-delimited file containing all chromosome paintings in the network. Useful for
	                           focusing on specific locations on the genome.	

The remaining files are either intermediate files or for debug purposes.

###Visualization

	Start Cytoscape
	Install the enhancedGraphics plugin if not installed
	Select 'Import Network from File'
	Find and Select 'cytoscapeGenome.xgmml' described in the Output section
	Select Layout -> Profuse Force Directed Layout or another of your choice
	In the Control Panel -> Style -> Properties -> Paint -> Custom Paint 1 select Image/Chart 1
	Image/Chart 1 will now be in the list of properties
	Under Image/chart 1, set Column = Gradient, Mapping Type = Passthrough Mapping
	The chromosome paintings will now be visible
	Adjust other visual properties as needed


###Additional Diagnostics

The NodeSummary.py script, included in the PopNet directory, is able to generate stacked bar graphs (similar those seen in supplemental figure 2) to aid in the determination of parameters such as section length and gap penalty. Currently, the script has not been optimized for user experience, and requires some direct editing by the user.

The first step is to run PopNet once under each of the condition being compared (i.e. with section length = 2000, 4000, 6000.. etc), and **saving the resulting xgmml as IXXPIXXSXXXX.xgmml**, where the first two 'XX' are the I and PI values used multiplied by 10, and the 'XXXX' following S is the value of the variable parameter. (i.e. if section length is being varied, and the run is done with I = 4, PI = 1.5, Section length = 8000, the file should be saved as I40PI15S8000.xgmml). All the numerical values need to be integers. If the file is not in this format it would not be recognized by NodeSummary.py.  

Place the generated xgmml files into a folder. This is the directory to be specified in the NodeSummary.py. Open NodeSummary.py, and go to the main function at the bottom. The parameters to be specified include the directory containing the input files, title of the graph, axis titles, the output file's name, and the bins. The bins control the size of the features represented by each stack on the stacked bar graph. If all the features are fall into the same stack (e.g. because they are too large or small), the bins can be adjusted to offer better resolution. Please note that the script only looks at recombinant features (i.e. regions where a sample has inherited genes from an ancestry other than its own). 

Run the script to generate the graph. The output file will be placed in the directory of the input files. 

###Examples

Two example datasets are provided, one for yeast in .SNPs format and one for toxoplasma in tabular format. Note that either format can be used for any species, the examples simply give two different species in two formats.

The two example configuration scripts can be used to analyze the example datasets. 

To run the yeast example:

First, edit the Example_Config_Yeast file to change the base_directory and output_directory to the absolute path of /examples/Nucmer and /examples/Nucmer/output on your computer. The prefilled default serves to illustrate the format, and will not work.

Then:
	cd PopNet
	python scripts/FullRunner.py Example_Config_Yeast.txt
	Open Cytoscape
	Select 'Import Network from File'
	Find and Select /examples/Nucmer/output/cytoscape/cytoscapeGenome.xgmml
	Follow subsequent steps described in Visualization to visualize

A similar procedure can be used to run the Toxoplasma dataset using Example_Config_Toxo.txt and /examples/Tabular/Toxo20.txt

Pregenerated results are already placed in the coresponding /output folders for reference.
