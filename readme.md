# EGRVA: An effective gene-based rare variant association analysis pipeline for case–control studies of disease

As a complementary pipeline for rare variant analysis on GWAS chips, EGRVA is straightforward and cost-efficient.

## 1.Pre-requisite software

Table 1. List of pre-requisite software and the available information.

| No.  | Software  | Version | Availability                                                 |
| :--- | --------- | ------- | ------------------------------------------------------------ |
| 1    | Michigan  | 1.2.4   | https://imputationserver.sph.umich.edu/index.html            |
| 2    | PLINK     | 1.9     | http://www.cog-genomics.org/plink2/                          |
| 3    | LiftOver  |         | https://genome.ucsc.edu/cgi-bin/hgLiftOver                   |
| 4    | SHAPEIT   | 2.r900  | https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html |
| 5    | HRC check | 4.2.11  | https://www.well.ox.ac.uk/~wrayner/tools/                    |
| 6    | BCFTOOLS  | 1.9     | http://samtools.github.io/bcftools/bcftools.html             |
| 7    | ANNOVAR   |         | http://annovar.openbioinformatics.org/en/latest/             |

Make sure the pre-requisite software listed in Table 1  has been installed. 

## 2.EGRVA pipeline

### Step 2.1 Four data pre-processing steps

1. Pre-QC executed by PLINK.

   ```
   $Plink --bfile file1 --keep-allele-order --hwe 0.00001 --geno 0.05 --maf 0.005 --max-maf 0.01 --mind 0.05 --make-bed --out file2
   ```

2. Use LiftOver to convert the genome coordinates into hg19.

   ```python
   python liftOverPlink.py -m file.map -p file.ped -o plink_ped -c hg18ToHg19.over.chain.gz -e liftOver
   ```

   NOTE:

   "liftOverPlink.py" and "hg18ToHg19.over.chain.gz"  can be downloaded  from https://github.com/ruijiali/EGRVA. If your dataset build is hg19, skip this step.

3. Use SHAPEIT2 as the phasing tool.

4. Check the PLINK bim file against the HRC reference SNP list with the script .

   ```
   $PLINK --bfile file --freq
   ```

   ```perl
   $perl HRC-1000G-check-bim.pl -b file.bim -f plink.frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.vcf -h
   ```

   NOTE:

   "HRC-1000G-check-bim.pl"  can be downloaded  from https://github.com/ruijiali/EGRVA.

### Step 2.2 Genotype imputation(Minimac4)

Sign up the Michigan Imputation service (https://imputationserver.sph.umich.edu/index.html) and upload vcf files separated by chromosomes. Generally, the filling operation time will not be too long.

NOTE:

We recommend HRC r1.1 2016 for Reference Panel, 0.3 for rsq Filter, and Population according to your samples.

### Step 2.3 Post-processing steps and QC.

1. Filter SNPs with R2 < 0.8.

   ```
   for i in {1..22}; do bcftools view -i 'R2>.8' -Oz chr${i}.dose.vcf.gz > file1; done
   ```

2. Changing the ID to '%CHROM:%POS:%REF:%ALT'.

   ```
   for i in {1..22}; do bcftools norm -Ou -m+any file1| bcftools annotate --output-type b --output file2  -I '%CHROM:%POS:%REF:%ALT'; done
   ```

3. Cut

   ```
   do bcftools view file2| cut -f3  | awk '{a[$0]++; if(a[$0]==2) print; if (a[$0]>=2) print}'; 
   ```

4. Convert bcf format to plink format.

   ```
   $PLINK --bcf file2--keep-allele-order --allow-extra-chr 0 --const-fid --split-x b37 no-fail --vcf-idspace-to _ --make-bed --out file4 
   ```

5. Check for duplicates and deletions.

   Run "duplicate.R" in R software, R reports "No duplicate SNPs".

   NOTE:

   "duplicate.R"  can be downloaded  from https://github.com/ruijiali/EGRVA.

6. Merge chromosomes.

   ```
   $PLINK --bfile file4 --merge-list mergelist.txt --biallelic-only --make-bed --out file5
   ```

7. Post-QC executed by PLINK.

   ```
   $PLINK --bfile file5 --keep-allele-order --hwe 0.00001 --geno 0.05 --maf 0.005 --max-maf 0.01 --make-bed --mind 0.05 --noweb --out file6
   ```

8. Remove multiallelic SNPs.

   ```
   $Plink --bfile file6 --snps-only --make-bed --out file7
   ```

### Step 2.4 Gene-based functional annotation.

1. Convert plink to VCF4 file.

   ```
   $PLINK --bfile file7 --recode vcf-iid --out file8
   ```

2. ANNOVAR annotation

   ```
   perl convert2annovar.pl -format vcf4old file8 > file9;
   perl file9 humandb/ -buildver hg19 -out file10 -remove -protocol related database -operation g,f
   ```

### Step 2.5 Statistical analysis

1. Count allele frequency seperately by case and control.

   ```
   $Plink --bfile file10 --keep case_samole_file -extract exonic_file  --recode vcf-iid --out  file11_case
   $Plink --bfile file10 --keep control_samole_file -extract exonic_file  --recode vcf-iid --out  file11_control
   ```

2. Combine to get the  file"12.csv"

3. Run "mirage.R" in R software.

   NOTE:

   "mirage.R" can be downloaded  from https://github.com/ruijiali/EGRVA. MIRAGE is available at https://xinhe-lab.github.io/mirage.

### Step 2.6 Bioinformatic analyses

Table 2. List of bioinformatics analysis tools and the available information.

| No.  | Tools                                            | Availability                       |
| :--- | ------------------------------------------------ | ---------------------------------- |
| 1    | NCBI database                                    | https://www.ncbi.nlm.nih.gov/gene/ |
| 2    | GeneCards database                               | https://www.genecards.org/         |
| 3    | Uniprot                                          | https://www.uniprot.org/           |
| 4    | Human Protein Atlas database                     | https://www.proteinatlas.org/      |
| 5    | online Mendelian Inheritance in the Man database | https://www.omim.org/              |
| 6    | human disease database                           | https://www.malacards.org/         |
| 7    | PathCards                                        | https://pathcards.genecards.org/   |