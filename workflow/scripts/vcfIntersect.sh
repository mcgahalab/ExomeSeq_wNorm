# UNSAGE: vcfCompare.sh output_dir Tumor_name Normal_name file1.vcf file2.vcf file3.vcf...

## Using bcftools, best explanation here: https://sourceforge.net/p/vcftools/mailman/message/32239335/
# This will create three VCFs, each containing only sites which appear in
# two or more files. The old vcf-isec works the way you expected, it would
# produce a single file, possibly containing a mixture of records from all
# three input files. This usually is not a desired behaviour as the VCFs
# may contain different sample columns, therefore this was discontinued.

module load vcftools/0.1.15
module load tabix/0.2.6
module load igenome-human/hg38
module load gatk/3.8
module load java/8
module load samtools/1.9

OUTPUT_DIR=$1
T=$2
N=$3
VCF_LIST=${@:4}

BASE_DIR=$(pwd)
cd $OUTPUT_DIR

for V in $(echo $VCF_LIST);do

  BASE=${V##*/}

  ## Sort VCFs
  vcf-sort -c ${BASE_DIR}/${V} > $BASE$'_sorted.vcf'

  ## Left align indels
  java -jar $gatk_dir/GenomeAnalysisTK.jar \
     --unsafe \
     -T LeftAlignAndTrimVariants \
     -R $REF \
     --variant $BASE$'_sorted.vcf' \
     -o $BASE$'_sorted_left_aligned.vcf' \
     --dontTrimAlleles

  bgzip -c $BASE$'_sorted_left_aligned.vcf' > $BASE$'.gz'
  tabix -p vcf $BASE$'.gz'

  V_LIST=${V_LIST}${BASE}'.gz '

done

bcftools isec -p $T --nfiles +2 $V_LIST
