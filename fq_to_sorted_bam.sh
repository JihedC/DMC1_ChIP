#!/bin/bash

# CVL: lets make sure this script can take the first part of the input filename (e.g. "Chip1_S_DA10") as argument

# CVL: First, we create a bash function that prints how this script should be used (if it is called with --help):
usage() {
    echo -e "\n Usage: $0 -p <prefix>

    -p,--prefix     Filenames to be analyzed will be ${prefix}_1.fq.gz and ${prefix}_2.fq.gz
    "
    exit 1;
}

# CVL: this loops picks up the arguments one by one and puts them in variables.
while [ $# -ne 0 ]; do
    case $1 in
        -p | --prefix)
            shift
            PREFIX="$1"
            ;;
        -h | --help)
            usage
            ;;
        --)
            ;;
        *)
            usage
            ;;
    esac
    shift
done

# CVL: Check that all mandatory arguments have been passed to the script, otherwise, print the usage:
# CVL: Checks if PREFIX is set:
if [[ -z "$TMPDIR"/"${PREFIX}" ]]; then
  # -z string flag is TRUE is length of true is non-zero
    echo -e "\n No prefix was supplied when calling $0"
    usage
fi

# CVL: Now we continue with the actual script
# CVL: NOTE, I did a find-and-replace operation on all your ${i} and replaced them with ${PREFIX}. Check if that produces the right filenames

### Modules required:
  # module bowtie/2.2.4
  # module samtools/0.1.19

### Files required:
  # tomato_3.00_inx bowtie2 indexed reference genome
  # Fastq files (6 couples of files as inputs):
    #Chip1_S_DA10_1.fq.gz Chip1_S_DA10_2.fq.gz
    #Chip2_S_DA11_1.fq.gz Chip2_S_DA11_2.fq.gz
    #Chip3_S_DA12_1.fq.gz Chip3_S_DA12_2.fq.gz
    #Chip4_S_DA19_1.fq.gz Chip4_S_DA19_2.fq.gz
    #Chip5_S_DA20_1.fq.gz Chip5_S_DA20_2.fq.gz
    #Chip6_S_DA21_1.fq.gz Chip6_S_DA21_2.fq.gz
  # Scripts in the pipeline:
      # /home/jihed/scripts/foo.py
      # /home/jihed/scripts/multi_unique_extract_pairend.r
  # copy all files to "$TMPDIR" in job script

### Output files :
  # "$TMPDIR"/"${PREFIX}"_lowmiss_unique_both_sort.bam
  # "$TMPDIR"/"${PREFIX}"_lowmiss_unique_both_sort.bam

echo "Mapping bowtie2 sample "$TMPDIR"/"${PREFIX}" to tomato genome..."
# CVL: Not sure what bowtie2 expects in terms of input for -1 and -2, but realize that you pass RELATIVE paths now (i.e. only the filename).
# CVL: It may be better to pass a full path, as then the functionality is independend on where this bash script is executed.
# CVL: My suggestion would actually be to make sure you also copy these input files to the TMPDIR, and then pass the arguments as e.g. -1 "${TMPDIR}"/"${PREFIX}"_1.fq.gz
# CVL: The same goes for all your intermediately created files, like the "$TMPDIR"/"${PREFIX}".bam, I would write these on scratch by replacing that part with "${TMPDIR}"/"${PREFIX}".bam.
# CVL: That goes for all the intermediate files throughout your code: most if these are written to the current directory, which means wherever you 'cd'-ed before you executed the script. It is more robust to make sure that all these files are written under "${TMPDIR}".
bowtie2 --very-sensitive --no-discordant --no-mixed -x "$TMPDIR"/tomato_inx_3.00 -1 "$TMPDIR"/"${PREFIX}"_1.fq.gz -2 "$TMPDIR"/"${PREFIX}"_2.fq.gz | samtools view -h -b -F 12 > "$TMPDIR"/"${PREFIX}".bam
# -1 : forward read
# -2 : reverse read
#-p 16 -k 10
echo "Removing unmapped reads from bam..."
samtools view -h -F 4 "$TMPDIR"/"${PREFIX}".bam > "$TMPDIR"/"${PREFIX}"_mapped.bam
# -h : include header in SAM output
# -F : only include reads with none of the FLAGS in INT present \
# here INT = 4, removes all the reads unmapped
echo "Overwriting previous bam file..."
samtools view -h -o "$TMPDIR"/"${PREFIX}".sam "$TMPDIR"/"${PREFIX}"_mapped.bam
# overwrite the previous bam file with filtered sam### allow a maximum of 2 mismatches in alignment ([^0-9] matches characters not in the range of 0 to 9)
echo "Filtering for maximum 2 mismatches in alignment..."
samtools view -Sh "$TMPDIR"/"${PREFIX}".sam | grep -e "^@" -e "XM:i:[012][^0-9]" > "$TMPDIR"/"${PREFIX}"_lowmiss.sam
# [012][^0-9] : 0 or 1 or 2 mismatches, excluding all number that could come after with the [^0-9]
# 'XM:' is an element of the sam file indicating the mismatches number
sed -n 1,15p "$TMPDIR"/"${PREFIX}"_lowmiss.sam > "$TMPDIR"/"${PREFIX}"_lowmiss_multi_header.sam
#Don't know what is -n 1,9p, may be generates only headers that is used later
#1,9p seems to indicate that the header would contain 9 lines, which is not enough for tomato, needed to chnage it to 15
echo "Including only the reads that are properly mapped in pair..."
samtools view -S -f 0x02 "$TMPDIR"/"${PREFIX}"_lowmiss.sam | grep -v "XS:i:" |python "$TMPDIR"/foo.py > "$TMPDIR"/"${PREFIX}"_lowmiss_unique.txt
# '-f 0x02' only include the reads that are properly mapped in pair
# 'grep -v "XS:i:' search that are not matching with "XS:i:" (unique mapping)
# do not forget to add the command 'python' before the path to the script
# both script need to be chmod +x before running the pipeline
cat "$TMPDIR"/"${PREFIX}"_lowmiss_multi_header.sam "$TMPDIR"/"${PREFIX}"_lowmiss_unique.txt > "$TMPDIR"/"${PREFIX}"_lowmiss_unique_ori.sam
#concatenate both files

# bowtie2 MAPQ scores >=42 correspond to uniquely mapping reads
samtools view -h -q 42 "${PREFIX}"_lowmiss_unique_ori.sam > "${PREFIX}"_lowmiss_unique.sam
samtools view -bS -o "${PREFIX}"_lowmiss_unique.bam "${PREFIX}"_lowmiss_unique.sam
samtools sort  "${PREFIX}"_lowmiss_unique.bam -o "${PREFIX}"_lowmiss_unique_sort.bam
samtools index "${PREFIX}"_lowmiss_unique_sort.bam

echo "Combining unique and multiple mapping reads..."
### identify and filter multiply mapping reads, and combine with uniquely mapping reads ("both")
samtools view -S -f 0x02 "$TMPDIR"/"${PREFIX}"_lowmiss.sam | grep "XS:i:"|python "$TMPDIR"/foo.py > "$TMPDIR"/"${PREFIX}"_lowmiss_multi.txt
# '-f 0x02' only include the reads that are properly mapped in pair
# 'grep -v "XS:i:' search that are matching with "XS:i:" (multiple mapping)
cat "$TMPDIR"/"${PREFIX}"_lowmiss_multi_header.sam "$TMPDIR"/"${PREFIX}"_lowmiss_multi.txt > "$TMPDIR"/"${PREFIX}"_lowmiss_multi.sam
# concatenate only header with the paired-reads aligning multiple times
samtools view -h -q 10 "$TMPDIR"/"${PREFIX}"_lowmiss_multi.sam > "$TMPDIR"/"${PREFIX}"_lowmiss_multi_fq10.sam
# -q 10 include only reads with mapping quality >= 10, why is the quality only 10 ?
samtools view -S -f 0x02 "$TMPDIR"/"${PREFIX}"_lowmiss_multi_fq10.sam | grep "XS:i:" | python "$TMPDIR"/foo.py > "$TMPDIR"/"${PREFIX}"_lowmiss_multi_fq10.txt
# same as before but this time only for read quality of 10
wc -l "$TMPDIR"/"${PREFIX}"_lowmiss_multi_fq10.txt >> "$TMPDIR"/"${PREFIX}"_mdim.stats 2>&1
echo "Running R script..."
Rscript "$TMPDIR"/multi_unique_extract_pairend.r "$TMPDIR"/"${PREFIX}"mdim.stats "$TMPDIR"/"${PREFIX}"_lowmiss_multi_fq10.txt "$TMPDIR"/"${PREFIX}"_MU.RData "$TMPDIR"/"${PREFIX}"_lowmiss_multi_fq10_unique.txt
#the pipeline here takes the multiple aligned reads and pipe it into Rscript to use
#the total mapping reads are locaed in the *_lowmiss.sam
#If XS is present the read is mapping multiple time

#R script rank the mutliple aligning reads and either choose the top one or random reads
cat "$TMPDIR"/"${PREFIX}"_lowmiss_multi_header.sam  "$TMPDIR"/"${PREFIX}"_lowmiss_multi_fq10_unique.txt > "$TMPDIR"/"${PREFIX}"_lowmiss_multi_fq10_unique.sam

rm "$TMPDIR"/"${PREFIX}"_lowmiss_multi.txt "$TMPDIR"/"${PREFIX}"_lowmiss_multi_fq10.txt "$TMPDIR"/"${PREFIX}"_lowmiss_multi_fq10_unique.txt
# remove unnecessary files
samtools view -bS -o "$TMPDIR"/"${PREFIX}"_lowmiss_multi_fq10_unique.bam "$TMPDIR"/"${PREFIX}"_lowmiss_multi_fq10_unique.sam
samtools sort "$TMPDIR"/"${PREFIX}"_lowmiss_multi_fq10_unique.bam -o "$TMPDIR"/"${PREFIX}"_lowmiss_multi_fq10_unique_sort.bam
samtools index "$TMPDIR"/"${PREFIX}"_lowmiss_multi_fq10_unique_sort.bam
samtools merge "$TMPDIR"/"${PREFIX}"_lowmiss_unique_both.bam "$TMPDIR"/"${PREFIX}"_lowmiss_unique_sort.bam "$TMPDIR"/"${PREFIX}"_lowmiss_multi_fq10_unique_sort.bam
# merge 2 sorted bam files
samtools sort "$TMPDIR"/"${PREFIX}"_lowmiss_unique_both.bam -o "$TMPDIR"/"${PREFIX}"_lowmiss_unique_both_sort.bam
samtools index "$TMPDIR"/"${PREFIX}"_lowmiss_unique_both_sort.bam
echo "Removing unnecessary files ..."
### remove no longer required .txt and .sam files
rm "$TMPDIR"/"${PREFIX}"_mapped_lowmiss_unique.txt
rm "$TMPDIR"/"${PREFIX}".sam "$TMPDIR"/"${PREFIX}"_mapped.bam "$TMPDIR"/"${PREFIX}"_lowmiss.sam "$TMPDIR"/"${PREFIX}"_lowmiss_multi_header.sam "$TMPDIR"/"${PREFIX}"_lowmiss_unique_ori.sam "$TMPDIR"/"${PREFIX}"_lowmiss_unique.sam "$TMPDIR"/"${PREFIX}"_lowmiss_multi.sam "$TMPDIR"/"${PREFIX}"_lowmiss_multi_fq10.sam "$TMPDIR"/"${PREFIX}"_lowmiss_multi_fq10_unique.sam
# As the files on the scratch I don't think this step is necessary, the job script should have a line to recover the files I want right?
# CVL: yes, you need to make a segment in the job script that copies back your output files to a directory in your home.
# CVL: also, it is good practice to cleanup the tempdir once you are done: rm -rf "${TMPDIR}"/* (in the job script). Lisa does this automatically after every job however, so it's not strictly needed.
