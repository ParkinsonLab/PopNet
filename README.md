### Welcome

I've recently moved the nn branch to master. There has been a significant upgrades to the software since publication, so please read through this readme file. 

I encourage you to use the browser-based version at [www.compsysbio.org/popnetd3](www.compsysbio.org/popnetd3). It runs this program in the backend on our server, thereby saving you both the trouble and the resources. If you have any questions regarding either this program or the browser-based version, please contact me at [javi.zhang@mail.utoronto.ca](mailto:javi.zhang@mail.utoronto.ca). 

Another advantage of the browser version is that the visualizer can handle more samples. I find that cytoscape is slow at rendering chromosome paintings at more than 20 samples, while the online visualizer can handle up to 200. 

Finally, 700 samples with 150k positions each appear to be the upper limit for a typical desktop with 32G of RAM. If your dataset is bigger, I would recommend using a high-performance computing resource. 


### Getting started

Try the new browser-based PopNetD3 at [www.compsysbio.org/popnetd3](www.compsysbio.org/popnetd3)

Submit jobs to be run on our server, and view your results right in the browser!

Requirements:
```
Operating System:
	Ubuntu or equivalent Linux-based OS

Dependencies:
	Python 3.6.5 or later
	Python packages:
		- Numpy
		- Matplotlib
		- Scikit-learn 
		- Pandas
	MCL (available from [http://micans.org/mcl/](http://micans.org/mcl))

Visualization:
	Cytoscape 2.7.0 or later (available from www.cytoscape.org/)
	enhancedGraphics plugin for Cytoscape (available from http://apps.cytoscape.org/apps/enhancedgraphics)
```

Quick Start:
```
Clone contents to local folder
Edit config file
cd /path/to/popnet/folder
python PopNet.py /path/to/config/file
```
### Introduction

PopNet is a set of Python scripts that generates a XGMML file to be visualized in Cytoscape as well as a json file that can be viewed on popnetd3. 
The scripts are run by executing the runner script, PopNet.py, with the path to a configuration file as argument. All configurable
options of the program are contained in the configuration file. The output XGMML file need to be manually loaded
into cytoscape for visualization. Information from several steps of the program are written to other files in the
output folder, which may be useful for extracting particular bits of information, optimizing parameters, or debugging

### Input Requirements

PopNet now expects a tab separated SNP table, formatted as:

\#CHROM	POS	SAMPLE1	SAMPLE2	...  
XX_CHRI	1	A	T	...  
XX_CHRI	2	G	C	...  

PopNet expects that all individual genomes to be aligned to a common reference. In other words, that the coordinates
in each genome refers to the same location. In addition, all chromosome names should be in the format:
	XX_ChrI
	XX_ChrXIV

### Configuration File

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

### Output

Output files are generated at the location specified in the config file, and are always named the same way.

The key file for network visualization is located at:
	/path/to/results/name.xgmml

where name is the name of your input file. 

Other potentially helpful output include:
	Heatmaps.pdf               Metrics to help you decide on the Inflation (I) and Pre-inflation (pI) parameters for 
	                           secondary clustering. Only created if autogroup = True in config.  

	log.txt                    A log file that includes all the run parameters  

	overall_similarity.tsv     The raw matrix used in secondary clustering. Derived from the co-clustring frequency of samples.  

	pclusters.txt              The clustering results of each individual chromosome segment during primary clustering  

	results.txt                A tab-delimited file containing all the SNPs used  

	chromosome_painting.tsv    A tab-delimited file containing all chromosome paintings in the network.  

The remaining files are either intermediate files or for debug purposes.

### Visualization
```
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
```

### Examples

The examples include a data set, toxo20.txt, and a config file, example_config.txt. 

To run the yeast example:

First, edit the example_config.txt file to change the input_path and output_directory to the absolute path of /examples/toxo20.txt and a desired output location on your computer. The prefilled default serves to illustrate the format, and will not work.

Then:
```
	cd PopNet
	python3 PopNet.py examples/example_config.txt
	Open Cytoscape
	Select 'Import Network from File'
	Find and Select the xgmml file generated at your output location
	Follow subsequent steps described in Visualization to visualize
```
