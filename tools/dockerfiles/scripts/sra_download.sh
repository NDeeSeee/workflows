#!/bin/bash

SRA_IDS=()
PARAMS=()

for i in "$@"; do
    if [[ "$i" = "--split-files" ]] || [[ "$i" = "--split-3" ]]; then
        echo "Adding fastq-dump downloading parameter $i"
        PARAMS+=($i)
    else
        echo "Adding SRR identifier to download $i"
        SRA_IDS+=($i)
    fi
done;

echo "### Single FASTQ files statistics" > single_fastq_stats.md
echo "### Merged FASTQ files statistics" > merged_fastq_stats.md
echo "### SRA metadata" > srr_metadata.md
echo "### Errors and Warnings" > debug.md


DATA_TYPE="TenX"                  # setting it to TenX by default

for SRA in ${SRA_IDS[@]}; do

    esummary -db sra -id $SRA -mode json|python3 -m json.tool > ${SRA}_meta.json
    esummary -db sra -id $SRA -mode xml > ${SRA}_meta.xml
    echo "#### ${SRA}" >> srr_metadata.md
    echo "**Title:** `esummary -db sra -id $SRA -mode xml | xtract -pattern Summary -element Title`" >> srr_metadata.md
    echo "**Bioproject:** `esummary -db sra -id $SRA -mode xml | xtract -pattern DocumentSummary -element Bioproject`" >> srr_metadata.md
    echo "**Biosample:** `esummary -db sra -id $SRA -mode xml | xtract -pattern DocumentSummary -element Biosample`" >> srr_metadata.md
    echo "**Study:** `esummary -db sra -id $SRA -mode xml | xtract -pattern DocumentSummary -element Study@name`" >> srr_metadata.md
    echo "**Spots:** `esummary -db sra -id $SRA -mode xml | xtract -pattern DocumentSummary -element Statistics@total_spots`" >> srr_metadata.md
    echo "**Bases:** `esummary -db sra -id $SRA -mode xml | xtract -pattern DocumentSummary -element Statistics@total_bases`" >> srr_metadata.md
    echo "**Size:** `esummary -db sra -id $SRA -mode xml | xtract -pattern DocumentSummary -element Statistics@total_size`" >> srr_metadata.md

    echo "Attempting to prefetch $SRA as TenX type"
    prefetch --type TenX $SRA > /dev/null 2>&1
    EXIT_CODE=$?
    if [[ $EXIT_CODE -ne 0 ]]
    then
        echo "Failed to prefetch $SRA as TenX type with exit code $EXIT_CODE"
        echo "Attempting to prefetch $SRA as sra type"
        prefetch --type sra $SRA > /dev/null 2>&1
        EXIT_CODE=$?
        if [[ $EXIT_CODE -ne 0 ]]
        then
            echo "- Failed to prefetch $SRA as sra type with exit code $EXIT_CODE. Exiting" >> debug.md
            break
        fi
        DATA_TYPE="sra"           # changing data type to "sra" as from now we assume that all SRR should have sra type
        echo "Extracting downloaded $SRA with ${PARAMS[@]} parameters"
        fastq-dump --gzip --log-level info ${PARAMS[@]} $SRA/$SRA.sra > /dev/null 2>&1
        j=1
        for FASTQ in $SRA*.gz; do
            echo "#### `basename $FASTQ`" >> single_fastq_stats.md
            echo "**`zcat $FASTQ | wc -l`** lines, **`stat -c%s $FASTQ`** bytes, top **5** reads" >> single_fastq_stats.md
            echo "\`\`\`" >> single_fastq_stats.md
            echo "`zcat $FASTQ | head -n 20`" >> single_fastq_stats.md
            echo "\`\`\`" >> single_fastq_stats.md
            echo "Adding $FASTQ to read_$j.fastq.gz"
            cat $FASTQ >> read_$j.fastq.gz
            rm -f $FASTQ
            (( j++ ))
        done;
        rm -rf $SRA
    else
        if [[ DATA_TYPE = "sra" ]]
        then
            echo "- Previous SRR was downloaded as sra data type. Current $SRA has TenX type. Exiting" >> debug.md
            break
        fi
        cellranger bamtofastq $SRA/*.bam extracted_fastq
        cat extracted_fastq/*/*_L00*_R1_00*.fastq.gz > read_1.fastq.gz
        cat extracted_fastq/*/*_L00*_R2_00*.fastq.gz > read_2.fastq.gz
        rm -rf $SRA extracted_fastq
        if [[ ${#SRA_IDS[@]} -ne 1 ]]
        then
            echo "- Merging multiple SRR identifiers extacted to BAM is not correct. Skipping remaining SRR" >> debug.md
            break
        fi
    fi
done;

for MERGED in read*.gz; do
    echo "#### `basename $MERGED`" >> merged_fastq_stats.md
    echo "**`zcat $MERGED | wc -l`** lines, **`stat -c%s $MERGED`** bytes" >> merged_fastq_stats.md
done;

cat srr_metadata.md single_fastq_stats.md merged_fastq_stats.md debug.md > report.md
rm -f srr_metadata.md single_fastq_stats.md merged_fastq_stats.md debug.md