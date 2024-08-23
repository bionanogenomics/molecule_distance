<img src="/images/Bionano-Logo.png" alt= “Bionano” width="5%" height="5%" title="GA"/>

# Molecule Distance Script #
The molecule distance workflow supplements Bionano's Guided Assembly and Rare Variant Analysis pipelines as a downstream analysis tool to identify mosaicism in repeat expansion genes or unstable insertions between 2 molecule labels.

This script calculates and visualizes the distance between two labels for all molecules (the exact labels chosen are a input parameter). There are 3 main steps:

1. Parse input files from `alignmolvref/merge/*` to get complete molecule distance dataset.
2. Trim complete dataset to remove outliers beyond 3 s.d.
3. Plot complete and trimmed data:

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
The molecule distance script also accepts other assemblies' (RVA and de novo) `alignmolvref` files as input. However, unlike GA and de novo, RVA output files (.xmap & .cmap) are not partitioned by chromosome. Directly running molecule distance script on RVA output will result in long runtimes. Therefore, a pre-processing step using `filter_map_files.py` is recommended to: 

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
# filter q.cmap & r.cmap by xmap IDs relevant to chr X (23)
python ~/.../molecule_distance_scripts/filter_map_files.py 23 \
~/.../output/contigs/alignmolvref/merge/exp_alignmolvref_r.cmap \
~/.../output/contigs/alignmolvref/merge/exp_alignmolvref.xmap \
~/.../output/contigs/alignmolvref/merge/exp_alignmolvref_q.cmap \
~/.../output/sample_1/sample_1_chr23 \
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
Rscript ~/.../molecule_distance_scripts/molecule_distance_full.R \
~/.../RFC1_16_guided_assembly/output/contigs/alignmolvref/merge/exp_refineFinal1_contig4_r.cmap \
~/.../RFC1_16_guided_assembly/output/contigs/alignmolvref/merge/exp_refineFinal1_contig4.xmap \
~/.../RFC1_16_guided_assembly/output/contigs/alignmolvref/merge/exp_refineFinal1_contig4_q.cmap \
7722 7725 ~/.../RFC1_16_guided_assembly/mol_distance \
```

Alternatively, if the .xmap and .cmap files have been processed by `molecule_distance_full.R`, the following script takes in `*_complete_data.csv`, and does only steps: 

2. to trim dataset for outliers and 
3. data plotting (e.g. delta from reference, histograms, auto-clustering). 

This script offers the flexibility to customize plots without executing resource intensive step 1. to parse .cmap & .xmap input files.

```
Rscript /archived/molecule_distance_delta_ref.R \
/path/to/out_handle_complete_data.csv \
/path/to/output/name_for_plot_files
```

Example command:
```
Rscript ~/.../molecule_distance_scripts/archived/molecule_distance_delta_ref.R \
~/.../sample_1/sample_1_complete_data.csv 13021 13022 \
~/.../sample_1/sample_1_output
```


### Script Output ###
The script outputs multiple files:

* Output csv of molecule distances and summary statistics.
* Output csv of molecule distances trimmed for outliers beyond 3 s.d.
* Generate various plots (Examples below).

#### Complete data violin plots ####

<img src="/tests/test_molecule_distance_full/guided_output/ga_test_violin_plot.png" alt= “GA” width="50%" height="50%" title="violin plots"/>

#### Complete data barplot ####
Distance between specified start and end labels for all molecules

<img src="/tests/test_molecule_distance_full/guided_output/ga_test_complete_distance_barplot.png" alt= “GA” width="50%" height="50%" title="complete data barplot"/>

#### Trimmed data barplot ####
Trimmed data excludes datapoints beyond 3 standard deviation from molecule distance mean.

<img src="/tests/test_molecule_distance_full/guided_output/ga_test_trimmed_distance_barplot.png" alt= “GA” width="50%" height="50%" title="trimmed data barplot"/>

#### Trimmed data delta barplot ####
The y-axis plots delta, which is the distances between molecules' labels subtracted by the distance between reference labels.

<img src="/tests/test_molecule_distance_full/guided_output/ga_test_trimmed_delta_barplot.png" alt= “GA” width="50%" height="50%" title="trimmed data delta plot"/>

#### Auto Clustering ####
Uses the Mclust package to determine clusters for molecule distances using a maximum-likelihood method. The example below shows clustered molecule distances between labels 26242 and 26247 in chrX.

<img src="/tests/test_molecule_distance_full/guided_output/ga_test_auto_clustering.png" alt= “GA” width="50%" height="50%" title="GMM plot"/>

### Interpretation of Results ###

Bar plot showing molecule distances between specified labels. There are two clusters observed in blue and red. The blue cluster refer to the normal allele, while red cluster on the right half of the x-axis suggests unstable expansion.

<img src="/tests/test_molecule_distance_full/guided_output/ga_test_auto_clustering.png" alt= “GA” width="50%" height="50%" title="GA"/>

## Contact ###

Molecule distance script is authored by Joyce Lee, Syukri Shukor, Jillian Burke, and Andy Pang. For any questions, please reach out to Syukri Shukor (sshukor@bionano.com) or Andy Pang (apang@bionano.com) for questions and issues.
