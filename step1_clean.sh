### Required Arguments ###
fastq1=$1
fastq2=$2

### Optional Arguments ###
outdir=$3
if [ ! -n "$outdir" ]; then
    outdir=`pwd`;
fi
fastp=$4
adapter=$5
if [ ! -n "$adapter" ]; then
    adapter="Illumina";
fi

if [ $adapter = "Illumina" ]; then
    forward="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
    reverse="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";

elif [ $adapter = "MGI" ]; then
    forward="AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA";
    reverse="AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG";

else
    forward=$6;
    reverse=$7;
    if [ -n "$forward" ] && [ -n "$reverse" ]; then
        adapter="UserDefined: $adapter";
    else
        echo "ERROR: Invalid adapter option: ${adapter}.";
        echo "For user defined adapter, arguments 6 and 7 should be given."
        exit;
    fi
fi

### Preperation and Software ###
prefix=$(echo `basename $fastq1` | awk -F . '{print $1}')

############
### Main ###
############

echo "[Adapter Used]: $adapter"
echo "[Adapter Used]: Forward=$forward"
echo "[Adapter Used]: Reverse=$reverse"

if [ ! -s "${outdir}/${prefix}.clean_R1.fq.gz" ] && \
   [ ! -s "${outdir}/${prefix}.clean_R2.fq.gz" ] && \
   [ ! -s "${outdir}/${prefix}.fastp.finish" ]; then
    echo "[FASTP] Start: `date`"
    time $fastp \
        --in1=$fastq1 --out1=${outdir}/${prefix}.clean_R1.fq.gz \
        --in2=$fastq2 --out2=${outdir}/${prefix}.clean_R2.fq.gz \
        --unpaired1=${outdir}/${prefix}.unpaired_R1.fq.gz \
        --unpaired2=${outdir}/${prefix}.unpaired_R2.fq.gz \
        --failed_out=${outdir}/${prefix}.filter_failed.fq.gz \
        --html=${outdir}/${prefix}.report.html \
        --json=${outdir}/${prefix}.report.json \
        --qualified_quality_phred=12 \
        --unqualified_percent_limit=50 \
        --n_base_limit=10 \
        --adapter_sequence=$forward \
        --adapter_sequence_r2=$reverse \
        --disable_trim_poly_g \
        --thread=16
    echo "[FASTP] End: `date`"
    echo "Practice_makes_perfect" > ${outdir}/${prefix}.fastp.finish

else
    info1="${outdir}/${prefix}.clean_R1.fq.gz exists"
    info2="${outdir}/${prefix}.clean_R2.fq.gz exists"
    info3="${outdir}/${prefix}.fastp.finish exists"
    echo "[Skip] fastp: ${info1}, ${info2} and ${info3}."
    exit
fi
