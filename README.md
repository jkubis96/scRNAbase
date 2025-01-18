### scRNAbase: Pre-Processing Toolkit for scRNA-Seq

<br />



<p align="right">
<img  src="https://github.com/jkubis96/Logos/blob/main/logos/jbs_current.png?raw=true" alt="drawing" width="180" />
</p>


### Author: Jakub Kubiś

<div align="left">
 Institute of Bioorganic Chemistry<br />
 Polish Academy of Sciences<br />
 Laboratory of Single Cell Analyses<br />

<br />

<p align="left">
<img  src="fig/lsca.png" alt="drawing" width="250" />
</p>
</div>


<br />


## Description

A comprehensive toolkit for the basic analysis of single-cell RNA sequencing (scRNA-seq) data without cell barcodes, relying exclusively on unique molecular identifiers (UMIs). Developed at the Laboratory of Single-Cell Analyses, Institute of Bioorganic Chemistry, PAS, in Poznań, Poland.


Core functionalities include:
- Genome Downloading
- Reads Quality Control - [fastp](https://github.com/OpenGene/fastp)
- UMI Selection - [UMI-tools](https://github.com/CGATOxford/UMI-tools)
- Genome Annotation - [STAR](https://github.com/alexdobin/STAR)
- Annotation Adjustment (names repairing, UTR extending) - [GTF-tool](https://github.com/jkubis96/GTF-tool)
- Reads Mapping - [STAR](https://github.com/alexdobin/STAR)
- BAM Indexing - [samtools](https://github.com/samtools/samtools)
- Features selection - [featureCounts](https://github.com/ShiLab-Bioinformatics/subread)

scRNAbase enables researchers to perform essential pre-processing analyses on scRNAseq datasets in a simple and efficient way.



#### Installation - docker pull

```
docker pull jkubis96/scrna_base:1.0.5
```

<br />



#### Genome directories hierarchy

genome/ ................................................... Main directory<br />
├── Homo_sapiens/ ......................................... Subdirectory for Homo sapiens genome<br />
│   └── index/ .............................................Subfolder for indexes and annotations<br />
│       ├── 100/ ...........................................Index for read length 100<br />
│       └── 150/ ...........................................Index for read length 150<br />
├── Mus_musculus/ ..........................................Subdirectory for Mus musculus genome<br />
│   └── index/ .............................................Subfolder for indexes and annotations<br />
│       ├── 100/ ...........................................Index for read length 100<br />
│       └── 150/ ...........................................Index for read length 150<br />

<br />


#### Genome Downloading

genome_downloading - download genome and annotation files<br />

Options for genome_downloading:<br />
  --genome_link (URL) - URL to download the genome file (required)<br />
  --annotation_link (URL) - URL to download the annotation file (required)<br />
  --species (string) - name of the species (e.g., Homo_sapiens) (required)<br />

<br />


```
cd genome
docker run -it -v $(pwd):/app/genome jkubis96/scrna_base:1.0.5 bash -c "cd /app/genome && bash"

Example:

IndexPip genome_downloading --genome_link URL --annotation_link URL --species Homo_sapiens
```

<br />


#### Genome Annotation

 genome_indexing - create reference files and index the genome


Options for genome_indexing:<br />
  --species (string) - name of the species (e.g., Homo_sapiens) (required)<br />
  --reads_length (int) - length of reads for STAR index (required)<br />
  --CPU (int) - number of threads to use (optional, default: number of CPU cores - 2)<br />
  --memory (int) - amount of memory for STAR (optional, default: all available RAM)<br />

Additional parameters:<br />
  --optimize (bool) - run GTF/GFF3 file adjustment (optional, default: TRUE)<br />
  --extend (bool) - extend parameter for genome prep (optional, default: FALSE)<br />
  --five_prime_utr (int) - length of 5' UTR (optional, default: 400)<br />
  --three_prime_utr (int) - length of 3' UTR (optional, default: 1000)<br />
  --coding_elements (list) - list of coding annotation elements (optional, default: EXON,CDS,TRANSCRIPT,MRNA)<br />
  --space (int) - minimal differential factor for separating features [genes] (optional, default: 100000)<br />

<br />


```
cd genome
docker run -it -v $(pwd):/app/genome jkubis96/scrna_base:1.0.5 bash -c "cd /app/genome && bash"

Example:

IndexPip genome_indexing --species Homo_sapiens --reads_length 100 --optimize TRUE --extend TRUE
```



<br />


#### Fastq - experimental data directories

data/ ......................................................Main directory for experimental FASTQ data<br />
├── name.1_R1.fastq.gz .....................................Read 1 of sample 1 in FASTQ format (compressed)<br />
├── name.1_R2.fastq.gz .....................................Read 2 of sample 1 in FASTQ format (compressed)<br />
├── name.2_R1.fastq.gz .....................................Read 1 of sample 2 in FASTQ format (compressed)<br />
├── name.2_R2.fastq.gz .....................................Read 2 of sample 2 in FASTQ format (compressed)<br />
├── results/ ...............................................Directory for analysis results<br />
│   └── name.1.fastp_report.html ...........................Quality control report for sample 1 (generated by fastp)<br />
│   └── name.2.fastp_report.html ...........................Quality control report for sample 2 (generated by fastp)<br />
│   └── matrices/ ..........................................Subdirectory for count matrices<br />
│       ├── gene_id_exon_name.1_genes_count_matrix .........Count matrix  exon sample 1<br />
│       ├── gene_id_UTR_name.1_genes_count_matrix ..........Count matrix  UTR sample 1<br />
│       ├── gene_id_exon_name.2_genes_count_matrix .........Count matrix  exon sample 2<br />
│       ├── gene_id_UTR_name.2_genes_count_matrix ..........Count matrix  UTR sample 2<br />


<br />

#### Reads Mapping

matrices_creating - function processes FASTQ files to create gene and transcript count matrices. It includes steps such as read trimming, UMI extraction, mapping, deduplication, and matrix creation.<br />

Options:<br />
  --umi_length (int) - length of the Unique Molecular Identifier (UMI). Required<br />
  --reads_length (int) - read length for the input sequencing reads. Required<br />
  --species (string) - species name (e.g., Homo_sapiens, Mus_musculus). Required<br />
  --qc_reads (TRUE|FALSE) - whether to perform quality control of reads (TRUE by default)<br />
  --multi (TRUE|FALSE) - count reads mapped to multiple locations (FALSE by default)<br />
  --annotation_names (string) - names of annotations to include, provided as a list. Defaults to 'gene_id,gene_name'<br />
  --annotation_side (string) - specifies the side of the annotation, e.g., 'exon,intron,UTR,five_prime_UTR,three_prime_UTR,transcript'<br />
  --CPU (int) - number of CPU threads to use. Defaults to (total CPU cores - 2)<br />



Note:<br />
  - All required arguments (--umi_length, --reads_length, --species) must be provided.<br />
  - Ensure that the genome directory and annotation files are correctly set up using the IndexingPip function.<br />

<br />

```
cd data
docker run -it -v $(pwd):/data -v genome:/app/genome jkubis96/scrna_base:1.0.5

Example:

CountFeatures matrices_creating --umi_length 12 --reads_length 100 --species Homo_sapiens --qc_reads TRUE --annotation_names gene_id,gene_name --annotation_side exon,intron,UTR,five_prime_UTR,three_prime_UTR,transcript
```
<br />

<br />


#### Have fun JBS©