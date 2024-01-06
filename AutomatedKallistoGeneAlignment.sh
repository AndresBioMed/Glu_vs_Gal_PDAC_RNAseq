#!/bin/bash

echo "----Welcome to Automated Kallisto Gene alignment [AKG]----"
echo "-- Made by Andrés Gordo, 2023-2024 --"
echo "This script will align all your samples in a new folder using Kallisto, and checking their quality with FastQC and MultiQC."


# Check if mamba is installed
if ! command -v mamba &> /dev/null; then
    echo "mamba is not installed. Installing mamba..."
    conda install mamba
    exit 1
fi

# Check if kallisto=0.48 is installed
if ! mamba list | grep -q "kallisto=0.48"; then
    echo "kallisto=0.48 is not installed. Installing kallisto=0.48..."
    mamba install -c bioconda "kallisto=0.48"
fi

# Check if fastqc is installed
if ! mamba list | grep -q "fastqc"; then
    echo "fastqc is not installed. Installing fastqc..."
    mamba install -c bioconda "fastqc"
fi

# Check if multiqc is installed
if ! mamba list | grep -q "multiqc"; then
    echo "multiqc is not installed. Installing multiqc..."
    mamba install -c bioconda "multiqc"
fi

echo "All required packages are installed."

# create the folders for the output
mkdir -p ~/new_AKG
mkdir -p ~/new_AKG/fastqc
mkdir -p ~/new_AKG/kallisto
mkdir -p ~/new_AKG/index

# Prompt the user for the input folder
read -p "Enter the absolute path (use realpath) to the input folder (containing *.gz files): " input_folder

# Prompt the user for the output folder
read -p "Enter the absolute path to the reference genome file (containing a *.fa file): " genome_file

# Prompt the user for the threads to be used
read -p "Enter the number of threads available in your machine: " threads

echo "Input folder:  $input_folder"
echo "Genome file: $genome_file"

# Analyze .gz files with FastQC
echo "--> AKG will now analyze the quality of your samples with FastQC"
cd $input_folder || exit
fastqc *.gz -t $threads

# Move files to the fastqc folder
cd $input_folder || exit
mv *fastqc* ~/new_AKG/fastqc
echo "--> AKG has finished analyzing the quality of your samples with FastQC"
echo "--> AKG will now create an index based on your reference genome"
# create the index for the reference genome
cd ~/new_AKG/index ||  exit
kallisto index -i Homo_sapiens.GRCh38.cdna.all.index $genome_file
echo "--> AKG has finished the index"

#Get the index realpath
index_file=$(realpath "$(find ~/new_AKG/index -type f -name "*.index" -print -quit)")
echo "Using index file: $index_file"

cd ~/new_AKG/kallisto || exit

# Iterate through all ".gz" files in the input folder
for input_file in "$input_folder"/*.gz; do
    # Extract the base name of the input file (without the path and extension)
    base_name=$(basename -s .fastq.gz "$input_file")

    # Create an output folder for the sample
    sample_output_folder="$base_name"
    mkdir -p "$sample_output_folder"
    echo "-> The sample $base_name is being aligned now by Kallisto"

    # Run kallisto quant for each input file
    kallisto quant -i "$index_file" -o "$sample_output_folder" -t "$threads" --single -l 250 -s 30 "$input_file" > "$sample_output_folder/$base_name.log" 2>&1
    echo "-> Kallisto has finished aligning $base_name, a log file has also been produced"
done

echo "--> AKG has now Finished processing all your samples"
echo "Summarising results via MultiQ"

cd ~/new_AKG || exit
multiqc -d .
echo "AKG has finished, the final report has been produced alongside with the pseudoalignments"
