These are scripts I used in calling, filtering, and analysing structural variants in SB clones.

<br>

**Basecall and demuliplex**

Raw data (pod5) were re-basecalled using dorado v0.7.2 with sup and modbase (5mC_5hmC and 6mA) models.

```
module purge
module load Dorado/0.7.2

dorado basecaller --device cuda:all --recursive --emit-moves \
--kit-name SQK-NBD114-96 \
--sample-sheet sample_sheet.csv \
sup,5mC_5hmC,6mA \
pod5_dir/ > /file.bam 
```

<br>

Most runs use the 96 barcode ligation kit (SQK-NBD114-96). if using other kit, change the --kit-name, e.g. SQK-RBK114.96 if using 96 barcode rapid kit.

Next, demultiplex the bam file using dorado demux

```
module purge
module load Dorado/0.7.2

dorado demux --sample-sheet sample_sheet.csv \
--output-dir output_DIR/ \
--no-classify file.bam        #from dorado basecaller above \
--threads 8
```

<br>

**Combine bam, convert to fastq, filter, and align to reference**

If the same library was run across multiple runs/flowcells, combine the bam files using samtools cat

```
module purge
module load SAMtools/1.19-GCC-12.3.0

samtools cat \
file1.bam \
file2.bam \
file3.bam \
-o output_DIR/.combined.bam \
--threads 8
```

We can run in loop for multiple samples, if combine after demultiplexing.

```
module purge
module load SAMtools/1.19-GCC-12.3.0

cd output_DIR            #DIR containing the bam files
file_names=(*.bam)
sample_names=("${file_names[@]%.*}")          #create the list of sample names

for sample_name in "${sample_names[@]}"
do
samtools cat \
$DIR1/${sample_name}.bam \
$DIR2/${sample_name}.bam \
$DIR3/{sample_name}.bam \
-o output_DIR/${sample_name}.combined.bam \
--threads 8
done
```

<br>

Next, convert the bam file into fastq, include some filtering (length, quality) if we want to, and map the fastq to reference. you can run in loop as above for multiple samples

```
module purge
module load SAMtools/1.19-GCC-12.3.0
module load minimap2/2.28-GCC-12.3.0 
module load chopper/0.5.0-GCC-11.3.0

samtools fastq -T '*' \
file.combined.bam | chopper -q 10 | gzip > file.combined.fastq.gz

minimap2 -ax lr:hq -y --secondary=no -t 16 \
ref.fasta \
file.combined.fastq.gz | samtools sort -@8 -o \ file.combined.aligned.sorted.bam
samtools index file.combined.aligned.sorted.bam
```

<br>

**Call variants**

Call structural variants using sniffles. Default mapq score is 20, which means that ony reads with mapq > 20 will be used in variant calling. As we are mapping to diploid genomes, there are regions with mapq = 0 or 1 as they mapped to both chromosomes equally well, so we have tried including all reads in this variant calling (by including --mapq 0)

```
module purge
module load Sniffles/2.3.3-gimkl-2022a-Python-3.10.5

sniffles --input file.combined.aligned.sorted.bam \
--snf file.snf \
--vcf file.vcf \
--tandem-repeats TE_annotation.bed \
--output-rnames           #to output the read names in the vcf file \
--threads 8 \
--reference ref.fasta \
--allow-overwrite \
#--mapq 0
```

<br>

**Filter variants**

We have tried different filtering criteria. 
1. Filter out all variants shared with SB1031, as these are highly likely artefacts, as our genome is built based on SB1031 genome. For this, we use truvari bench which compare a sample vcf with SB1031 vcf, and output shared and unshared variants.

    Truvari required vcf to be in bgzip and indexed (using tabix). Truvari is not installed on NeSI, so installed our own in a conda environment.

```
    module load Miniconda3/4.9.2
    source activate /nesi/project/uoo03533/env/truvari
    module load tabix/0.2.6-GCCcore-9.2.0

    #bgzip and index
    bgzip -c file.vcf > file.vcf.gz
    tabix -p vcf file.vcf.gz

    #run truvari bench
    truvari bench -b SB1031.vcf.gz \        #reference vcf to compare to
    -c $file.vcf.gz \                       #sample vcf to compare
    -o DIR_truvari/sample_name/             #output directory

    #copy and re-name the fp.vcf.gz, which contains variants present in comparison (sample vcf), but not reference vcf, and decompress
    cp DIR_truvari/sample_name/fp.vcf.gz \             
    DIR_truvari/sample_name/filtered.vcf.gz
    bgzip -d DIR_truvari/sample_name/filtered.vcf.gz
```

These are the output details of truvari bench:
![truvari_bench](https://github.com/user-attachments/assets/bf5fd6f8-7fad-4842-a0c6-163f8385936a)



<br>

2. Filter based on QC matrices in sniffles output using bcftools
    Example below filter variants that are "PASS", "PRECISE", and with allele frequency >=0.4
    
    ```
    module purge
    module load  BCFtools/1.19-GCC-11.3.0
    
    bcftools view -i 'FILTER="PASS"' file.vcf | \
    grep -v 'IMPRECISE' | \
    bcftools view -i 'INFO/AF>=0.4' > file_filtered.vcf
    ```

<br>

**Annotate variants**

Next, we can annotate variants with gene information using VEP. all files needed to be in compressed form

```
module purge
module load VEP/107.0-GCC-11.3.0-Perl-5.34.1

vep --offline \
-i file.filtered.vcf.gz \
--gtf gene_annotation.gtf.gz \
--fasta ref.fasta.gz \
-o file.filtered.annotated.vcf \
--vcf \                              #output as vcf file
--force_overwrite \
--stats_text \                       #write the stats in text file, useful to plot stats
--stats_html
```

<br>

**Compare variants**

We can comparestructural variants among different samples, among same samples from different run or with different coverage, etc, using tools such as intervene, or SURVIVOR.

1. Intervene
   
   Intervene can create Venn, Upset, or heatmap diagram, and also generate set of variants shared or unshared among the samples. Intervene can do up to 6-way comparison for Venn diagram. 
   
   ```
    module load Miniconda3/4.9.2
    source activate /nesi/project/brins03581/env/intervene
    export PYTHONUSERBASE=/nesi/nobackup/brins03581/Cen/python_userbase/intervene
    export PATH=$PATH:/nesi/nobackup/brins03581/Cen/python_userbase/intervene/bin
    conda install -c conda-forge pandas=1.5.3
    
    intervene venn -i file1.vcf.gz file2.vcf.gz file3.vcf.gz \
    -o output_DIR/ \
    --save-overlaps \
    --project project_name \
    --bedtools-options header \
    --title title_of_svg_file
    
    ```

2. SURVIVOR
   
   ```
    module purge
    module load SURVIVOR/1.0.7

    #get all the vcf we want to compare in one directory
    cd DIR\
    ls *vcf > samplefile

    # can check the content of the sample file:
    file1.vcf
    file2.vcf
    file3.vcf, etc

    #First, run SURVIVOR merge as below
    #SURVIVOR merge
    #File with VCF names and paths
    #max distance between breakpoints (0-1 percent of length, 1- number of bp) 
    #Minimum number of supporting caller
    #Take the type into account (1==yes, else no)
    #Take the strands of SVs into account (1==yes, else no)
    #Estimate distance based on the size of SV (1==yes, else no).
    #Minimum size of SVs to be taken into account.
    #Output VCF filename

    SURVIVOR merge samplefile 1000 1 1 1 0 30 output.merged.vcf

    #run SURVIVOR genComp
    #Merged Vcf file
    #Normalize output (1==yes, else no)
    #Output: pariwise overlap matrix

    SURVIVOR genComp /output.merged.vcf \
    0 matrix.txt

    #run this perl script to get the matrix overlap file, which we can use in R to create Venn or upset plot
    perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' output.merged.vcf | \
    sed -e 's/\(.\)/\1 /g' > matrix.overlap.txt
    
    
    #in R
    
    setwd ("working_Dir")
    
    library(readxl)
    library(ggVennDiagram)
    library(UpSetR)
    
    matrix=read.table("matrix.oerlap.txt", header=T)
    
    #plot the Venn and upset diagram
    #Convert dataframe into list
    set_list <- list(
      X1 = which(matrix$X1 == 1),
      X2 = which(matrix$X2 == 1),
      X3  = which(matrix$X3 == 1),
      X4  = which(matrix$X4 == 1),
      X5 = which(Smatrix$X5 == 1)
    )

    ggVennDiagram(set_list) +
      scale_fill_gradient(low = "white", high = "blue") +
      ggtitle("Title") +
      theme(plot.title = element_text(size = 32, face = "bold",
      hjust=0.5, margin=margin(0,0,30,0)))

    upset(matrix, 
      sets = c("X1", "X2", "X3", "X4", "X5"),
      sets.bar.color = "blue",
      order.by = "freq",
      mainbar.y.label = "Intersection Size",
      sets.x.label = "Set Size",
      text.scale = 1.5)
    ```

 <br>
 
 
 
 
    














