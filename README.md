# 'INSERT PAPER TITLE HERE' - Colaço et al., 2019
## RNA-Seq Analysis

### Lung Dataset
- Samples from *Mus musculus*, strain C57BL/6J;
- 4 Groups:
    5 Non-Infected + Injected with PBS
    5 Non-Infected + Injected with 1.75 ug/g Doxy
    5 Infected with *E.coli* + Injected with PBS
    5 Infected with *E.coli* + Injected with 1.75 ug/g Doxy

#### Quality Control
This quality control is done using the program **Fastqc** (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/); considering this program produces a report per file, we can also use a program to merge all the reports into one - **MultiQC** (https://multiqc.info/).

``mkdir fastqc_reports
mkdir mQC_reports

fastqc raw/*.txt.gz -o fastqc_reports -t 5``
