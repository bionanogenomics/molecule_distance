#!/bin/bash

# Commands to run molecule_distance_full script multicontigs
# the molecule distance full script is configured for both RVA & Guided assembly alignmolvref/merge output files (xmap & cmap files)

# Activate conda
conda activate R \


# Guided Assembly output
Rscript /home/users6/sshukor/molecule_distance_scripts/molecule_distance_full.R \
/home/IrysView1/users9/jburke/repeat_expansion_project/guided_assembly_enfocus_output/RUMC_CNBP_08_*_pipelines/output/contigs/alignmolvref/merge/exp_refineFinal1_contig3_r.cmap \
/home/IrysView1/users9/jburke/repeat_expansion_project/guided_assembly_enfocus_output/RUMC_CNBP_08_*_pipelines/output/contigs/alignmolvref/merge/exp_refineFinal1_contig3.xmap \
/home/IrysView1/users9/jburke/repeat_expansion_project/guided_assembly_enfocus_output/RUMC_CNBP_08_*_pipelines/output/contigs/alignmolvref/merge/exp_refineFinal1_contig3_q.cmap 26242 26247 \
/home/users6/sshukor/molecule_distance_scripts/tests/test_molecule_distance_full/guided_output/ga_test \

Rscript /home/users6/sshukor/molecule_distance_scripts/archived/molecule_distance_delta_ref.R \
/home/users6/sshukor/molecule_distance_scripts/tests/test_molecule_distance_full/guided_output/ga_test_complete_data.csv 26242 26247 \
/home/users6/sshukor/molecule_distance_scripts/tests/test_molecule_distance_full/guided_delta_ref_output/ga_test_delta_ref



# RVA output
# Requires filtering of xmap & q.cmap files as output files are not partitioned by chr
# select xmap rows that match X chromosome
grep -v "^#" exp_alignmolvref.xmap | awk '{if($3 == 23) print}' >> /home/users6/sshukor/nih_niaid_mol_distance/1941659/chr23.xmap \

# filter q.cmap file for rows matching existing xmapIDs
python /home/users6/sshukor/nih_niaid_mol_distance/scripts/filter_cmap.py \
/home/users6/sshukor/nih_niaid_mol_distance/1941659/chr23.xmap \
/home/users3/apang/human/NIH-NIAID/1941659/alignmolvref/merge/exp_alignmolvref_q.cmap \
/home/users6/sshukor/nih_niaid_mol_distance/1941659 \

Rscript /home/users6/sshukor/molecule_distance_scripts/molecule_distance_full.R \
/home/users3/apang/human/NIH-NIAID/1941659/alignmolvref/merge/exp_alignmolvref_r.cmap \
/home/users6/sshukor/nih_niaid_mol_distance/1941659/chr23.xmap \
/home/users6/sshukor/nih_niaid_mol_distance/1941659/chr23_q.cmap 13021 13022 \
/home/users6/sshukor/molecule_distance_scripts/tests/test_molecule_distance_full/rva_output/rva_test \

# Deactivate conda
conda deactivate
