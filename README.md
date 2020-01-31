# 'Host-dependent induction of disease tolerance to infection by tetracycline antibiotics' - Colaço *et al.*, 2019
## RNA-Seq Analysis

### Lung Dataset
- Samples from *Mus musculus*, strain C57BL/6J;
- 4 Groups:
    - 5 Non-Infected + Injected with PBS
    - 5 Non-Infected + Injected with 1.75 ug/g Doxy
    - 5 Infected with *E.coli* + Injected with PBS
    - 5 Infected with *E.coli* + Injected with 1.75 ug/g Doxy

#### Quality Control
This quality control is done using the program **Fastqc** (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/); considering this program produces a report per file, we can also use a program to merge all the reports into one - **MultiQC** (https://multiqc.info/).

```
mkdir fastqc_reports
mkdir mQC_reports

fastqc raw/*.txt.gz -o fastqc_reports -t 5
```
#### Alignment
For the alignment steps, we're going to use **STAR Alignment** (https://github.com/alexdobin/STAR)

###### 1) Create the Reference Genome Indexes

First of all, create a new folder to store the genome indexes and also create a shortcut for it:
```
mkdir gen_index
export gen_index=./gen_index
```

Then, you can assess it by just using:

```
cd $gen_index
```

To align our samples, we need the reference genome indexes, that are just the corresponding reference genome (.fasta) and its corresponding annotation (.gtf). To get both of these files, we go to the Ensembl website (https://www.ensembl.org/info/data/ftp/index.html) and get the corresponding files:

```
wget ftp://ftp.ensembl.org/pub/release-97/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz

wget ftp://ftp.ensembl.org/pub/release-97/gtf/mus_musculus/Mus_musculus.GRCm38.97.gtf.gz
gunzip Mus_musculus.GRCm38.97.gtf.gz
```
We are going to use the most recent version of the annotation (**97 release**). Now that we have the files, we proceed to use STAR with the option of  `genomeGenerate`

```
STAR --runThreadN 10 --runMode genomeGenerate --genomeDir $gen_index --genomeFastaFiles $gen_index/Mus_musculus.GRCm38.dna.primary_assembly.fa --sjdbGTFfile $gen_index/Mus_musculus.GRCm38.97.gtf 

```
Options explained:

- `--runThreadN 10` Nr. of Cores used
- `--runMode genomeGenerate` Argument for the program to know what is going to run
- `--genomeDir $gen_index/` Path to the genome indexes
- `--genomeFastaFiles $gen_index/Mus_musculus.GRCm38.dna.primary_assembly.fa` Path to the fasta file of the reference genome
- `--sjdbGTFfile $gen_index/Mus_musculus.GRCm38.95.gtf` Path to the gtf file of the reference genome

Afterwards, you can proceed to the alignment step. For this, we will perform a cycle but, before, we need to create the several required folders to store the necessary files:
```
mkdir aligned
mkdir qualimap_aligned

```
Now, we go to the directory where the merged files are and proceed with the cycle:
```
cd raw/
export gen_index=../../gen_index #this step helps in case of using a different window
```
```
for f in *.txt.gz;
do STAR --genomeDir $gen_index --readFilesIn $f --readFilesCommand zcat --sjdbGTFfile $gen_index/Mus_musculus.GRCm38.97.gtf --quantMode GeneCounts --runThreadN 24 --outFileNamePrefix ../aligned/"$f"_ --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx;
qualimap rnaseq -bam ../aligned/"$f"_*.bam -gtf $gen_index/Mus_musculus.GRCm38.97.gtf -outdir ../qualimap_aligned/$f --java-mem-size=40G;
rm -rf ../aligned/$f*.bam; done
```

This cycle will perform the alignment and perform the corresponding report for the file of the alignment.
