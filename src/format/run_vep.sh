FILE=$1
OUTPATH=$2
NAMEFILE=$(basename $FILE)
zcat $FILE | awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' | vep -i STDIN --assembly GRCh38 --no_stats --cache --symbol --protein --canonical --offline --dir /workspace/datasets/vep --format vcf --custom /workspace/datasets/gnomad/gnomad.genomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz,gnomADg,vcf,exact,0,A --output_file $OUTPATH$NAMEFILE.vep.bed