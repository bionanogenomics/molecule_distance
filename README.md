<img src="/images/Bionano-Logo.png" alt= “Bionano” width="5%" height="5%" title="GA"/>

# Molecule Distance Script #
The molecule distance workflow supplements Bionano's Guided Assembly and Rare Variant Analysis pipelines as a downstream analysis tool to identify mosaicism in repeat expansion genes or unstable insertions between 2 DLE1 molecule labels.

This script calculates and visualizes the distance between two labels for all molecules (the exact labels chosen are a input parameter). There are 3 main steps:
1.  Parse input files from `alignmolvref/merge/*` to get complete molecule distance dataset.
2.  Trim complete dataset to remove outliers beyond 3 s.d.
3.  Plot complete and trimmed data:

## Set-up and Usage

This script accepts both Rare Variant Analysis (RVA) and Guided Assembly (GA) files as input,specifically intermediary output files from `alignmolvref/merge/*`. However, pre-processing is recommended for RVA to dramatically reduce runtime.

Skip to **Launching Molecule Distance Scripts** for quick-run.

Full command examples in ``` test/test_run_molecule_distance_full_script_commands.sh ```


### Set-up/dependencies ###

* This script can only be run on machine with Linux-based operating system because it relies on several Linux tools & Python.
* There is an optional script "mol_dist_reqs" that will install the required R packages for the command line interface. This can be installed onto the local computer directly or through running this script.
* Additional R packages:

```
mclust
prabclus
ggplot2
```

### Local Guided Assembly (GA) via Command Line ###
It is recommended to launch assemblies from molecule .bnx files via command line, since assembly intermediate files -`alignmolvref` step - are needed as input for molecule distance script. Assemblies downloaded from Access will likely not contain `alignmolvref` intermediate files. The output file is similar to an Enfocus Fragile X assembly, with the difference that it estimates repeat expansion distances along the locus specified by the .csv file instead of _FMR1_. More information on Solve (3.7.2) Local Guided Assembly can be found in [local_guided_assembly repository](https://bitbucket.org/bionanoclinicalaffairs/local_guided_assembly/src/master/).


### Rare Variant Analysis (RVA) Input Pre-processing ###
The molecule distance script also accepts other assemblies' (RVA and de novo) `alignmolvref` files as input. However, unlike GA and de novo, RVA output files (.xmap & .cmap) are not partitioned by chromosome. Therefore, a pre-processing step is recommended to:
	1. Select .xmap rows with `RefContigID` matching chromosome of interest
	2. Filter r.cmap & q.cmap files for rows matching filtered .xmap file RefContigID & QryContigID columns, respectively.

Launch
```
python filter_map_files.py chr \
/path/to/*_r.cmap \
/path/to/*.xmap \
/path/to/*_q.cmap \
/path/to/output/name_for_plot_files
```

For example, the following command specifies filtering of molecules from chrX of the human genome.
```
# filter q.cmap & r.cmap by xmap IDs relevant to chr 23
python /home/users6/sshukor/molecule_distance_scripts/filter_map_files.py 23 \
/home/output/contigs/alignmolvref/merge/exp_alignmolvref_r.cmap \
/home/output/contigs/alignmolvref/merge/exp_alignmolvref.xmap \
/home/output/contigs/alignmolvref/merge/exp_alignmolvref_q.cmap \
/home/users6/sshukor/nih_niaid_mol_distance/2006039/2006039_chr23 \
```

This script outputs filtered r.cmap, .xmap, and q.cmap files to be used as input for molecule distance script below.


### Launching Molecule Distance Scripts ###

After activating `R` conda environment. From the command line, the general usage is:
```
Rscript molecule_distance_full.R \
/path/to/r_cmap \
/path/to/xmap \
/path/to/q_cmap \
start_label end_label path/to/output/name_for_plot_files
```

Example command:
```
# launch molecule distance script on guided assembly alignmolvref files from chr 4, calculating molecule distances containing reference labels 7722 and 7725
Rscript /home/users6/sshukor/molecule_distance_scripts/molecule_distance_full.R \
/home/.../RFC1_16_guided_assembly/output/contigs/alignmolvref/merge/exp_refineFinal1_contig4_r.cmap \
/home/.../RFC1_16_guided_assembly/output/contigs/alignmolvref/merge/exp_refineFinal1_contig4.xmap \
/home/.../RFC1_16_guided_assembly/output/contigs/alignmolvref/merge/exp_refineFinal1_contig4_q.cmap \
7722 7725 /home/.../RFC1_16_guided_assembly/mol_distance \
```

Alternatively, if the .xmap and .cmap files have been processed by `molecule_distance_full.R`, the following script takes in `*_complete_data.csv`, and does steps 2. to trim dataset for outliers and 3. data plotting (e.g. delta from referece, histograms, auto-clustering). This script offers the flexibility to customize plots without executing resource intensive step 1. to parse .cmap & .xmap input files.
```
Rscript /archived/molecule_distance_delta_ref.R \
/path/to/out_handle_complete_data.csv \
/path/to/output/name_for_plot_files
```

Example command:
```
Rscript /home/users6/sshukor/molecule_distance_scripts/archived/molecule_distance_delta_ref.R \
/home/.../2006039/2006039_nih_niaid_complete_data.csv 13021 13022 \
/home/.../2006039/2006039_output
```


### Script Output ###
The script outputs multiple files:

* Output csv of molecule distances and summary statistics.
* Output csv of molecule distances trimmed for outliers beyond 3 s.d.
* Generate various plots (Examples below).
	
#### Complete data histogram ####

<img src="/tests/test_molecule_distance_full/guided_output/ga_test_histogram.png" alt= “GA” width="1%" height="1%" title="histogram"/>

#### Complete data violin plots ####

<img src="/tests/test_molecule_distance_full/guided_output/ga_test_violin_plot.png" alt= “GA” width="5%" height="5%" title="violin plots"/>

#### Complete data barplot ####
Distance between specified start and end labels for all molecules

<img src="/tests/test_molecule_distance_full/guided_output/ga_test_complete_distance_barplot.png" alt= “GA” width="5%" height="5%" title="complete data barplot"/>

#### Trimmed data barplot ####
Trimmed data excludes datapoints beyond 3 standard deviation from molecule distance mean.

<img src="/tests/test_molecule_distance_full/guided_output/ga_test_trimmed_distance_barplot.png" alt= “GA” width="5%" height="5%" title="trimmed data barplot"/>

#### Trimmed data delta barplot  #### 

<img src="/tests/test_molecule_distance_full/guided_output/ga_test_trimmed_delta_barplot.png" alt= “GA” width="5%" height="5%" title="trimmed data delta plot"/>

#### Auto Clustering ####
Uses the Mclust package to determine clusters for the trimmed molecule distances using a maximum-likelihood method. The example below shows clustered molecule distances between labels 26242 and 26247 in chrX.

<img src="/tests/test_molecule_distance_full/guided_output/ga_test_auto_clustering.png" alt= “GA” width="5%" height="5%" title="GMM plot"/>

#### Mapping Molecules to Allele ####
Guided assembly was launched on .bnx molecules, and its outputs were used as input for both Enfocus and molecule distance. Therefore, molecules can be mapped between independent molecule distance and Enfocus assembly runs. `mol_distance_to_enfocus_map.ipynb` does exactly this, with the addition of assigning alleles to each molecule. However, assigning alleles to Enfocus maps is a manual process, and depends on repeat expansion lengths independently estimated using de novo and Enfocus. Once alleles were assigned to Enfocus maps, then the notebook script is able to complete the final step to assign alleles to molecule distance output. Refer to  `mol_distance_to_enfocus_map.ipynb` Jupyter notebook for more information.

### Interpretation of Results ###

Bar plot showing difference in molecule distances vs reference label distances. Bars with y>0 are longer than reference, and the second right half of the x-axis suggests unstable expansion between the specified labels.

<img src="/tests/test_molecule_distance_full/guided_output/ga_test_auto_clustering.png" alt= “GA” width="5%" height="5%" title="GA"/>

## Contact ###

Molecule distance script is authored by Joyce Lee, Syukri Shukor, and Jillian Burke. For any questions, please reach out to Syukri Shukor (sshukor@bionano.com), Andy Pang (apang@bionano.com), or support@bionano.com for questions and issues.