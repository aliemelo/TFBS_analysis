#!/bin/bash

while getopts i:o:b:g:c:f:d: flag
do
case "${flag}" in
                i) genelist=${OPTARG}
                        ;;
                o) outfolder=${OPTARG}
                         ;;
                b) annotation=${OPTARG}
                         ;;
                g) genome=${OPTARG}
                         ;;
                c) chrsize=${OPTARG}
                         ;;
                f) feat=${OPTARG}
                         ;;
                *) echo "Invalid option: -$flag" ;;
        esac
done

# get current directory
cwd=$(pwd)
# create output directory
echo "### Creating output directory ###"
echo ""
mkdir $cwd/$outfolder

# the mikado genes have a problem with the name, so we need to fix it first
sed '/^mikado/s/.\{1\}$/.&/' $genelist | sed '/^mikado/s/.\{6\}/&.scga7/' > $cwd/$outfolder/genelist.corrected

# get the annotation for the target genes
echo "### Getting annotation for the target genes ###"
echo ""
grep -Fw -f $cwd/$outfolder/genelist.corrected $annotation | awk -v var="$feat" '$3 == var {print}' > $cwd/$outfolder/genes_entries.gff
python3 /projects/aliemelo/resources/scripts/create_start_bed.py -f $cwd/$outfolder/genes_entries.gff

# check which genes don't have a match in the annotation file
grep -Fw -f $cwd/$outfolder/genelist.corrected $annotation | awk -v var="$feat" '$3 == var {print}' | cut -f1 | grep -vFf - $cwd/$outfolder/genelist.corrected > $cwd/$outfolder/unmatched_genes

# if no matches found, stop script
if ! [ -s $cwd/$outfolder/genes_entries.gff ]; then
    echo "No matches found in annotation file"
    exit 0
fi

# run TSSFinder to find TSS positions
echo "### Running TSSFinder ###"
echo ""
tssfinder \
--model /projects/mauro/TSSFinder/articleData/athaliana/models/athaliana.0 \
--start $cwd/$outfolder/genes_entries.gff.start.bed \
--genome $genome \
--output $cwd/$outfolder/ \
--max_seq_size 2000 > $cwd/$outfolder/log.tssfinder.2000.txt

# define the promoter coordinates
echo "### Defining the promoter region ###"
echo ""
bedtools flank -s -i $cwd/$outfolder/out.tss.bed -g $chrsize -l 800 -r 0 > $cwd/$outfolder/promoters_800nt.bed

# get the promoter sequences
bedtools getfasta -name -s -fi $genome -bed $cwd/$outfolder/promoters_800nt.bed -fo $cwd/$outfolder/promoters_800nt.fa.tmp

# remove sequences shorter than 100 nt
awk -v RS=">" -v FS="\n" '{for(i=2;i<NF;i++) {l+=length($i)}; if(l>100) printf ">%s", $0}' $cwd/$outfolder/promoters_800nt.fa.tmp > $cwd/$outfolder/promoters_800nt.fa.tmp.tmp

# shorten ID of sequences
sed 's/;Name=.*//' $cwd/$outfolder/promoters_800nt.fa.tmp.tmp | sed 's/;Parent=.*//' | sed 's/ID=.*mikado.//' | sed 's/ID=.*model.//' | sed 's/ID=.*TU.//' > $cwd/$outfolder/promoters_800nt.fa.tmp.tmp.tmp
mv $cwd/$outfolder/promoters_800nt.fa.tmp.tmp.tmp $cwd/$outfolder/promoters_800nt.fa
rm $cwd/$outfolder/*tmp*

# obtain background model
echo "### Creating background model ###"
echo ""
fasta-get-markov $cwd/$outfolder/promoters_800nt.fa $cwd/$outfolder/background_model.meme_format

# scan sequences
echo ""
echo "### Scanning promoter sequences with FIMO ###"
echo ""
fimo \
--bfile $cwd/$outfolder/background_model.meme_format \
--max-stored-scores 500000 \
--thresh 1e-6 \
--o $cwd/$outfolder/fimo.out \
/projects/aliemelo/resources/databases/JASPAR2018_CORE_plants_non-redundant_pfms_meme.txt \
$cwd/$outfolder/promoters_800nt.fa > $cwd/$outfolder/log.fimo.txt

# filter results for q-value <= 0.05
echo "### Filtering results ###"
echo ""
awk '($9 <= 0.05) {print}' $cwd/$outfolder/fimo.out/fimo.tsv > $cwd/$outfolder/fimo.out/fimo.qvalue0.05.tsv

python3 /projects/aliemelo/resources/scripts/fimo_tsv_to_gff.py --tsv $cwd/$outfolder/fimo.out/fimo.qvalue0.05.tsv
sort -k1,1 -k4,4n $cwd/$outfolder/fimo.out/fimo.qvalue0.05.gff > $cwd/$outfolder/fimo.out/ordered.fimo_results.gff
bedtools merge -i $cwd/$outfolder/fimo.out/ordered.fimo_results.gff -s -c 9 -o distinct -d -6 > $cwd/$outfolder/fimo.out/fimo.qvalue0.05.merged_overlaps.gff
python3 /projects/aliemelo/resources/scripts/merged_bed_to_csv.py --bed $cwd/$outfolder/fimo.out/fimo.qvalue0.05.merged_overlaps.gff
python3 /projects/aliemelo/resources/scripts/fimo_prep.py --file $cwd/$outfolder/fimo.out/fimo.qvalue0.05.merged_overlaps.csv --out $cwd/$outfolder/fimo.out/mapping.results --alt T

# remove intermediary files
rm $cwd/$outfolder/fimo.out/ordered.fimo_results.gff


# check classes of TFs mapped to the target genes
cut -f1 $cwd/$outfolder/fimo.out/fimo.qvalue0.05.tsv | sort | uniq | sed '/^#/d' | sed '/^[[:space:]]*$/d' > $cwd/$outfolder/fimo.out/jaspar_ids

grep -F -f $cwd/$outfolder/fimo.out/jaspar_ids /projects/aliemelo/resources/databases/jaspar_core_plants_classes.csv > $cwd/$outfolder/fimo.out/classes

grep -oFf $cwd/$outfolder/fimo.out/jaspar_ids /projects/aliemelo/resources/databases/jaspar_core_plants_classes.csv | grep -vFf - $cwd/$outfolder/fimo.out/jaspar_ids > $cwd/$outfolder/fimo.out/unmatched_tfs

echo "### The script has finished ###"
echo ""
