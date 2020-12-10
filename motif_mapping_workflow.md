# 1. Dependencies and required files
In order to execute this workflow, you'll need to have installed:
1. [bedtools](https://github.com/arq5x/bedtools2)
2. [TSSFinder](https://github.com/tssfinder/tssfinder.github.io)
3. [FIMO from the MEME Suite](http://meme-suite.org/doc/fimo.html)
4. [fasta-get-markov](http://meme-suite.org/doc/fasta-get-markov.html) or [CreateBackgroundModel](http://bioinformatics.intec.ugent.be/MotifSuite/standalones.php) to create the background model
   * Alternatively, you can provide a previously computed background model (it must be in the MEME background format)

You'll also need a file containing your chromossome names and their respective sizes:
```bash
cat <YOURFASTAFILE> | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > chrom_size.tsv
```

The genome sequence (fasta) and annotation (gff) are also needed. We analysed the following organisms:
1. *Arabidopsis thaliana* (TAIR10): [genome](https://tinyurl.com/y63nexnj "Ensembl's FTP page") | [gff](https://tinyurl.com/y45ztpdq "Ensembl's FTP page")
2. *Saccharum spontaneum* AP85-411: [genome and gff](http://www.life.illinois.edu/ming/downloads/Spontaneum_genome/)
3. *Saccharum* hybrid cultivar SP80-3280: [genome](https://www.ncbi.nlm.nih.gov/nuccore/QPEU01000000) | gff

**Note:** upload the models and scripts to the repository and add links to them.

**Other files:**
1. Models for TSSFinder:
   * *A. thaliana*: `/projects/mauro/TSSFinder/articleData/athaliana/models/athaliana.0`
   * *O. sativa*: `/projects/mauro/TSSFinder/articleData/osativa/models/osativa.0`
2. Bedfile containing the start positions of the target genes
```bash
# filter gff file to get only the gene entries
awk '$3 == "gene" {print}' <GFF_FILE> > genes_entries.gff

# create a bedfile with the start coordinates of genes
# this creates a file called <filename>.start.bed
python3 /projects/aliemelo/resources/scripts/create_start_bed.py -f genes_entries.gff
```
3. JASPAR CORE Plants PFMs: [Plants PFMs (non-redundant) single batch file in MEME format](http://jaspar.genereg.net/download/CORE/JASPAR2020_CORE_plants_non-redundant_pfms_meme.txt)

# 2. Preparing your promoter sequences
Here we defined the promoter as the region of 800 nt upstream of the TSS.

## 2.1 Prediction of the TSS positions (TSSFinder)
```bash
nohup time tssfinder \
--model <MODELORGANISM> \
--start genes_entries.gff.start.bed \
--genome <GENOMEFASTAFILE> \
--output ./ \
--max_seq_size 2000 > log.tssfinder.2000.txt &
```
**Note:** the output is a file called `out.tss.bed`.

## 2.2 Delimitation of the promoter sequences
Define the promoter coordinates
```bash
bedtools flank -s -i out.tss.bed -g chrom_size.tsv -l 800 -r 0 > promoters_800nt.bed
```
Get the promoter sequences
```bash
bedtools getfasta -name -s -fi <GENOME> -bed promoters_800nt.bed -fo promoters_800nt.fa.tmp
```
Remove sequences shorter than 100 nt
```bash
awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 100 {print ">"$0}' promoters_800nt.fa.tmp > promoters_800nt.fa
```

> :heavy_exclamation_mark: If your FASTA IDs are too long, you might need to shorten them to only show the gene's name (as opposed to the whole description of that gene).
```bash
sed 's/;Name=.*//' promoters_800nt.fa | sed 's/ID=//' > promoters_800nt.fa.new
rm promoters_800nt.fa
mv promoters_800nt.fa.new promoters_800nt.fa
```

**Optional steps:**

Verify the lengths of your promoter sequences:
```bash
cat promoters_800nt.fa | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > promoters_length
```

Add the length of each sequence to its ID:
```bash
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' promoters_800nt.fa | awk -F '\t' '{printf("%s_%d\n%s\n",$1,length($2),$2);}' > promoters_800nt.fa.new
rm promoters_800nt.fa
mv promoters_800nt.fa.new promoters_800nt.fa
```

# 3. Scan promoters for motif occurences
## 3.1 Obtain the background model for your sequences
Here we will be using the CreateBackgroundModel software followed by the conversion of its output to the MEME background format
```bash
nohup time CreateBackgroundModel \
-f promoters_800nt.fa \
-b background_model > log.create_background_model.txt &

# convert to MEME format
python3 SCRIPT -f background_model
```
## 3.2 Scan sequences
We used the JASPAR CORE Plants
```bash
nohup time fimo \
--bfile background_model.meme_format \
--max-stored-scores 500000 \
--thresh 1e-6 \
/projects/aliemelo/resources/databases/JASPAR2018_CORE_plants_non-redundant_pfms_meme.txt \
promoters_800nt.fa > log.fimo.txt &
```

# 4. Assess results
Filter results for q-values <= 0.01
```bash
cd fimo.out/
awk '($9 <= 0.01) {print}' fimo.tsv > fimo.qvalue0.01.tsv
```
Many of the mapped PWMs overlap each other, so we can merge those sites into one locus
```bash
# convert the tsv file into a gff file
python3 /projects/aliemelo/resources/scripts/fimo_tsv_to_gff.py --tsv fimo.qvalue0.01.tsv

# sort the gff file by chromosome and start position
sort -k1,1 -k4,4n fimo.qvalue0.01.tsv_converted_to_gff > ordered.fimo_results.gff

# merge overlapping sites
bedtools merge -i ordered.fimo_results.gff -s -c 9 -o distinct > fimo_results_merged_overlaps.gff

# transform merged gff into a csv file (to be opened as a spreadsheet)
python3 /projects/aliemelo/resources/scripts/merged_bed_to_csv.py --bed fimo_results_merged_overlaps.gff

# make a table with genes names and their corresponding mapped PWMs
python3 /projects/aliemelo/resources/scripts/fimo_prep.py --file fimo_results_merged_overlaps.gff.csv --alt T
```



