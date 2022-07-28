#!/bin/bash
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=M.I.vanGinneken@umcutrecht.nl



output_dir="patient2"
haplohiv_folder="/hpc/dla_lti/dvanginneken/HaploHIV_Daphne/haplohiv2/patient2_out/"
table="/hpc/dla_lti/dvanginneken/MergedTimepoints_PredictHaplo/Mastertable_timepoints.tsv"

mkdir -p ${output_dir}/consensus
if [ "$(ls -A ${output_dir}/consensus)" ]; then
    echo "    Already grepped consensus name"
else
    for sample in `ls ${haplohiv_folder}/ec`
    do
        Rscript ConsensusFiles.R $table $haplohiv_folder $output_dir $sample
    done
fi


snakemake -j 8 --config output_dir=$output_dir haplohiv_folder=$haplohiv_folder table=$table --

