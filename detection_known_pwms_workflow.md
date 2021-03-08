# Dependencies and required files
In order to execute this workflow, you'll need to have installed:
1. [bedtools](https://github.com/arq5x/bedtools2)
2. [TSSFinder](https://github.com/tssfinder/tssfinder.github.io)
3. [FIMO from the MEME Suite](http://meme-suite.org/doc/fimo.html)
4. [fasta-get-markov](http://meme-suite.org/doc/fasta-get-markov.html)
   * Alternatively, you can provide a previously computed background model (it must be in the MEME background format).

**Other files:**
1. Models for TSSFinder:
   * [A. thaliana](resources/TSSFinder_models/athaliana)
   * [O. sativa](resources/TSSFinder_models/osativa)
   * Models for different organisms are available in the [TSSFinder download page](https://tssfinder.github.io/download.html)
2. JASPAR CORE Plants PFMs: [Plants PFMs (non-redundant) single batch file in MEME format](http://jaspar.genereg.net/download/CORE/JASPAR2020_CORE_plants_non-redundant_pfms_meme.txt)
3. *Saccharum* hybrid cultivar SP80-3280: [genome](https://www.ncbi.nlm.nih.gov/nuccore/QPEU01000000) | gff
4. File containing your chromosomes names and their respective sizes, separated by tab. You can obtain this file by running this command line:
```bash
cat <YOURFASTAFILE> | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > chrom_size.tsv
```

# Case study: co-expressed genes
We're gonna apply the workflow to [this set of co-expressed genes](test_sample_files/test_sample_list.txt) in *Saccharum* hybrid cultivar SP80-3280. Here we will break down each step of the workflow with its respective command line call. 

## 1. Preparing your promoter sequences
We defined the promoter as the region of 800 bp upstream of the TSS. First, we predict the TSS position using TSSFinder. Note that we need the start positions of each gene in BED format. The file used for this example can be found [here](test_sample_files/genes_entries.gff.start.bed).

### 1.1 Prediction of the TSS positions (TSSFinder)
```bash
nohup time tssfinder \
--model TSSFinder_models/athaliana \
--start genes_entries.gff.start.bed \
--genome <GENOMEFASTAFILE> \
--output ./ \
--max_seq_size 2000 > log.tssfinder.2000.txt &
```
**Note:** the output is a file called `out.tss.bed`.

### 1.2 Delimitation of the promoter sequences
1. Define the promoter coordinates:
```bash
bedtools flank -s -i out.tss.bed -g chrom_size.tsv -l 800 -r 0 > promoters_800nt.bed
```
2. Get the promoter sequences:
```bash
bedtools getfasta -name -s -fi <GENOME> -bed promoters_800nt.bed -fo promoters_800nt.fa.tmp
```
3. Remove sequences shorter than 100 nt:
```bash
awk -v RS=">" -v FS="\n" '{for(i=2;i<NF;i++) {l+=length($i)}; if(l>100) printf ">%s", $0}' promoters_800nt.fa.tmp > promoters_800nt.fa
```
Alternatively, you can use seqtk:
```bash
seqtk seq -L 100 promoters_800nt.fa.tmp > promoters_800nt.fa
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

## 2. Scan promoters for motif occurences
### 2.1 Obtain the background model for your sequences
Here we will be using the fasta-get-markov script from the MEME Suite:
```bash
fasta-get-markov promoters_800nt.fa background_model.meme_format
```

### 2.2 Scan sequences
We used the JASPAR CORE Plants as the database. We altered the options --max-stored-scores to 500000 and --thresh to 1e-6
```bash
nohup time fimo \
--bfile background_model.meme_format \
--max-stored-scores 500000 \
--thresh 1e-6 \
JASPAR2018_CORE_plants_non-redundant_pfms_meme.txt \
promoters_800nt.fa > log.fimo.txt &
```

## 3. Assess results
Filter results for q-values <= 0.01:
```bash
cd fimo.out/
awk '($9 <= 0.01) {print}' fimo.tsv > fimo.qvalue0.01.tsv
```
Your final result is the file `fimo.qvalue0.01.tsv`. The next steps are optional, and you are encouraged to explore the results as you see fit.

**Optional steps**
* Many of the mapped PWMs overlap each other, so we can merge those sites into one locus:
```bash
# convert the tsv file into a gff file
python3 fimo_tsv_to_gff.py --tsv fimo.qvalue0.01.tsv

# sort the gff file by chromosome and start position
sort -k1,1 -k4,4n fimo.qvalue0.01.tsv_converted_to_gff > ordered.fimo_results.gff

# merge overlapping sites
bedtools merge -i ordered.fimo_results.gff -s -c 9 -o distinct > fimo_results_merged_overlaps.gff
rm ordered.fimo_results.gff
```
File with merged loci: `fimo_results_merged_overlaps.gff`.

* Get a list of the genes with its corresponding mapped PWMs:
```bash
# transform merged gff into a csv file (to be opened as a spreadsheet)
python3 merged_bed_to_csv.py --bed fimo_results_merged_overlaps.gff

# make a table with genes names and their corresponding mapped PWMs
python3 fimo_prep.py --file fimo_results_merged_overlaps.gff.csv --out mapping.results --alt T
```
File with list of mapped PWMs: `mapping.results`.

* The coordinates in the output GFF files are in relation to the input sequences, that is, the promoter sequences used for the scanning step. To convert those coordinates to the where they would be in the chromosome, you can use this script:
```bash
python3 convert_genomic_position_v2.py \
-f fimo.qvalue0.01.gff \
--bed /promoters_800nt.bed
```
Output file: `mapping.results`.

Check the class for each of the mapped PWM:
```bash
cut -f1 fimo.qvalue0.01.tsv | sort | uniq | sed '/^#/d' | sed '/^[[:space:]]*$/d' > jaspar_ids

grep -F -f jaspar_ids jaspar_core_plants_classes.csv > classes
```
Output file: `classes`.
