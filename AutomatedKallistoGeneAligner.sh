#!/bin/bash

echo "----Welcome to Automated Kallisto Gene alignment [AKG]----"
echo "-- Made by Andrés Gordo, 2023-2024 --"
echo "This script will align all your samples in a new folder using Kallisto, and checking their quality with FastQC and MultiQC. Optionally, it can also use trim_galore for quality trimming before FastQC."

# Check if mamba is installed
if ! command -v mamba &> /dev/null; then
    echo "mamba is not installed. Installing mamba..."
    conda install mamba
    exit 1
fi

# Check if kallisto is installed
if mamba list | grep -q "kallisto"; then
    installed_version=$(mamba list | awk '/kallisto/ {print $2}')
    if [ "$installed_version" != "0.48.0" ]; then
        echo "kallisto is installed, but not version 0.48. Installing kallisto=0.48..."
        mamba install -c bioconda "kallisto=0.48"
    else
        echo "kallisto=0.48 is already installed."
    fi
else
    echo "kallisto is not installed. Installing kallisto=0.48..."
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

# Check if trim_galore is installed
if ! mamba list | grep -q "trim-galore"; then
    echo "trim_galore is not installed. Installing trim_galore..."
    mamba install -c bioconda "trim_galore"
fi

echo "All required packages are installed."




# Prompt the user for the input folder
read -p "Name you new RNAseq experiment: " experiment_name
read -p "Enter the absolute path (use realpath) to the input folder (containing *.gz files): " input_folder
# Prompt user for FASTQ file extension
echo "Enter the extension of the FASTQ files (e.g., .fq.gz or .fastq.gz):"
read -r file_extension

read -p "Do you want to use trim_galore for quality trimming before FastQC? (y/n): " use_trim

if [ "$use_trim" = "y" ]; then
    read -p "Do you want to delete the original files after trimming? (y/n): " delete_original
fi
# Prompt the user to choose between single-end and paired-end data
read -p "Enter 's' for single-end data or 'p' for paired-end data: " data_type

# Create the folders for the output
mkdir -p ~/"$experiment_name"
mkdir -p ~/"$experiment_name"/fastqc
mkdir -p ~/"$experiment_name"/kallisto

if [ "$data_type" = "p" ]; then
    read -p "Enter the extension for the first paired-end sample (e.g., _f1): " ext1
    read -p "Enter the extension for the second paired-end sample (e.g., _r2): " ext2
fi

# Prompt the user for the threads to be used
read -p "Enter the number of threads available in your machine: " threads

# Check if the user already has an index file
read -p "Do you already have an index file? (y/n): " has_index

if [ "$has_index" = "y" ]; then
    # Prompt the user for the index file path
    read -p "Enter the absolute path to the index file: " index_file
else
    # Prompt the user for the output folder
    read -p "Enter the absolute path to the reference genome file (containing a *.fa file): " genome_file
    
    mkdir -p ~/"$experiment_name"/index
    echo "--> AKG will now create an index based on your reference genome"

    # Create the index for the reference genome
    cd ~/"$experiment_name"/index || exit
    kallisto index -i Homo_sapiens.GRCh38.cdna.all.index "$genome_file"
    echo "--> AKG has finished creating the index"

    # Get the index realpath
    index_file=$(realpath "$(find ~/"$experiment_name"/index -type f -name "*.index" -print -quit)")
fi

echo "Using index file: $index_file"
echo "Input folder:  $input_folder"



# Trim_galore and move files to trim_galore_output folder
if [ "$use_trim" = "y" ]; then
    mkdir -p ~/"$experiment_name"/trim_galore_output
    echo "--> AKG will now use trim_galore for quality trimming before FastQC"
    cd "$input_folder" || exit

    if [ "$data_type" = "s" ]; then
        # Single-end data
        for input_file in "$input_folder"/*"$file_extension"; do
            base_name=$(basename "$input_file" "$file_extension")
            echo "$base_name"
            trim_galore -o ~/"$experiment_name"/trim_galore_output -j 4 "$input_file"
            mv ~/"$experiment_name"/trim_galore_output/"$base_name"_trimmed.fq.gz ~/"$experiment_name"/trim_galore_output/"$base_name""$file_extension"
            if [ "$delete_original" = "y" ]; then
                rm "$input_file"
            fi
        done
    elif [ "$data_type" = "p" ]; then
        # Paired-end data
        for input_file1 in "$input_folder"/*"$ext1""$file_extension"; do
            base_name=$(basename -s "$ext1$file_extension" "$input_file1")
            echo "$base_name"
            input_file2="$input_folder/$base_name$ext2$file_extension"
            
            if [ -e "$input_file2" ]; then
                trim_galore --paired -o ~/"$experiment_name"/trim_galore_output -j 4 "$input_file1" "$input_file2"
                mv ~/"$experiment_name"/trim_galore_output/"$base_name""$ext1"_val_1.fq.gz ~/"$experiment_name"/trim_galore_output/"$base_name""$ext1""$file_extension"
                mv ~/"$experiment_name"/trim_galore_output/"$base_name""$ext2"_val_2.fq.gz ~/"$experiment_name"/trim_galore_output/"$base_name""$ext2""$file_extension"
                if [ "$delete_original" = "y" ]; then
                    rm "$input_file1" "$input_file2"
                fi

            else
                echo "Warning: Paired-end file for $base_name not found. Skipping..."
            fi
        done
    else
        echo "Invalid data type. Please choose 's' for single-end or 'p' for paired-end data."
        exit 1
    fi

    input_folder=~/"$experiment_name"/trim_galore_output
    echo "--> AKG has finished trimming with trim_galore"
fi


# Analyze .gz files with FastQC
echo "--> AKG will now analyze the quality of your samples with FastQC"
cd $input_folder || exit
fastqc *.gz -t $threads

# Move files to the fastqc folder
cd $input_folder || exit
mv *fastqc* ~/"$experiment_name"/fastqc
echo "--> AKG has finished analyzing the quality of your samples with FastQC"

cd ~/"$experiment_name"/kallisto || exit

# Iterate through all ".gz" files in the input folder
if [ "$data_type" = "s" ]; then
    # Single-end data
    for input_file in "$input_folder"/*.gz; do
        base_name=$(basename -s "$file_extension" "$input_file")
        sample_output_folder="$base_name"
        mkdir -p "$sample_output_folder"
        echo "-> The sample $base_name is being aligned now by Kallisto"
        kallisto quant -i "$index_file" -o "$sample_output_folder" -t "$threads" --single -l 250 -s 30 "$input_file" > "$sample_output_folder/$base_name.log" 2>&1
        echo "-> Kallisto has finished aligning $base_name, a log file has also been produced"
    done
elif [ "$data_type" = "p" ]; then
    # Paired-end data

    for input_file1 in "$input_folder"/*"$ext1""$file_extension"; do
        base_name=$(basename -s "$ext1""$file_extension" "$input_file1")
        input_file2="$input_folder/$base_name$ext2""$file_extension"
        
        if [ -e "$input_file2" ]; then
            sample_output_folder="$base_name"
            mkdir -p "$sample_output_folder"
            echo "-> The sample $base_name is being aligned now by Kallisto"
            kallisto quant -i "$index_file" -o "$sample_output_folder" -t "$threads" "$input_file1" "$input_file2" > "$sample_output_folder/$base_name.log" 2>&1
            echo "-> Kallisto has finished aligning $base_name, a log file has also been produced"
        else
            echo "Warning: Paired-end file for $base_name not found. Skipping..."
        fi
    done
else
    echo "Invalid data type. Please choose 's' for single-end or 'p' for paired-end data."
    exit 1
fi

echo "--> AKG has now Finished processing all your samples"
echo "Summarising results via MultiQC"

cd ~/"$experiment_name" || exit
multiqc -d .
echo "AKG has finished, the final report has been produced alongside with the pseudoalignments"
