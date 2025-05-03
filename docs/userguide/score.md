# Calculate Statistics

`sai` provides the `score` command to calculate the number of uniquely shared sites (U statistic) and the quantile of the allele frequencies in such sites (Q statistic).

To see available options, we can use the following command:

```
sai score -h
```

This will show information for each argument:

| Argument | Description |
| - | - |
| -h, --help    | show this help message and exit |
| --vcf         | Path to the VCF file containing variant data |
| --chr-name    | Chromosome name to analyze from the VCF file. |
| --ref         | Path to the file with reference population identifiers. |
| --tgt         | Path to the file with target population identifiers. |
| --src         | Path to the file with source population identifiers. |
| --win-len     | Length of each genomic window in base pairs. Default: 50,000. |
| --win-step    | Step size in base pairs between consecutive windows. Default: 10,000. |
| --anc-alleles | Path to the BED file with ancestral allele information. If ancestral allele information is not provided, filtering will be performed for each variant based on whether the allele frequency of any allele (assuming biallelic) meets the specified condition during the calculation of the statistics. Default: None. |
| --ploidy      | Ploidy values for reference, target, and source populations (in that order). Default: 2 2 2. |
| --w           | Frequency threshold for variants in the reference population; only variants with frequencies below this threshold are included in the analysis. Default: 0.01. |
| --y           | List of allele frequency conditions for the source populations. Each value must be in the form =X, >X, <X, >=X, or <=X (e.g., =0.7, >0.8, <0.1, >=0.5, <=0.2). The number of values must match the number of source populations in the file specified by `--src`; the order of the allele frequency conditions should also correspond to the order of source populations in that file. Default: =1. |
| --output      | Output file path for saving results. |
| --stat        | Type of statistic to compute: UXX or QXX, where XX is a percentage-like index indicating a threshold in the target population. For example, `U50` means the allele frequency is greater than 0.5, and `Q95` means the allele frequency is greater than or equal to the 95th percentile among sites meeting the specified conditions. |

## Input files

`sai` requires a biallelic VCF file containing all populations, including the reference population (assumed to have no introgression), the target population (assumed to have received introgressed material), and the source population (assumed to have provided the introgressed material). It also assumes that all variants observed in these populations are already included in the input VCF file. The program does not fill in missing genotypes (e.g., treating unobserved genotypes as homozygous reference, as done in [MaLAdapt](https://github.com/xzhang-popgen/maladapt/blob/fba58dac6278d1b4cc573c78ab7c56911b4072b0/empirical/1empirical_compute-features.py#L318-L326)). An example VCF file can be found [here](https://github.com/xin-huang/sai/blob/main/examples/data/1KG.nea_den.chr9.example.vcf.gz).

Three tab-delimited files are required to specify the population information for individuals in the reference, target, and source populations. An example file is shown below and can also be found [here](https://github.com/xin-huang/sai/blob/main/examples/data/1KG.nea_den.samples.txt).

```
NEA	AltaiNeandertal
DEN	Denisova
```

This file specifies the individuals in the source populations. In this example, there are two source populations, NEA and DEN (first column), and the second column lists the corresponding sample names matching those in the VCF header.

## Output files

An example output of the U statistic calculation is shown below.

| Chrom | Start | End | Ref | Tgt | Src | N(Variants) | U50(w<0.01,y=(=1.0,=1.0)) | Candidate |
| - | - | - | - | - | - | - | - | - |
| 9 | 16400001 | 16440000 | AFR+EAS | EUR | NEA,DEN | 1082 | 0 | NA |
| 9 | 16440001 | 16480000 | AFR+EAS | EUR | NEA,DEN | 983  | 0 | NA | 
| 9 | 16480001 | 16520000 | AFR+EAS | EUR | NEA,DEN | 985  | 0 | NA |
| 9 | 16520001 | 16560000 | AFR+EAS | EUR | NEA,DEN | 966  | 0 | NA |
| 9 | 16560001 | 16600000 | AFR+EAS | EUR | NEA,DEN | 1029 | 0 | NA |
| 9 | 16600001 | 16640000 | AFR+EAS | EUR | NEA,DEN | 1063 | 0 | NA |
| 9 | 16640001 | 16680000 | AFR+EAS | EUR | NEA,DEN | 1153 | 0 | NA |
| 9 | 16680001 | 16720000 | AFR+EAS | EUR | NEA,DEN | 999  | 1 | 9:16715657 |
| 9 | 16720001 | 16760000 | AFR+EAS | EUR | NEA,DEN | 1007 | 2 | 9:16758710,9:16759357 |
| 9 | 16760001 | 16800000 | AFR+EAS | EUR | NEA,DEN | 1060 | 1 | 9:16770229 |
| 9 | 16800001 | 16840000 | AFR+EAS | EUR | NEA,DEN | 960  | 1 | 9:16800789 |
| 9 | 16840001 | 16880000 | AFR+EAS | EUR | NEA,DEN | 923  | 0 | NA |
| 9 | 16880001 | 16920000 | AFR+EAS | EUR | NEA,DEN | 512  | 0 | NA |

Each row in the output file corresponds to a genomic window, summarizing the result for that region. The columns are as follows:

| Column | Description |
| - | - |
| Chrom | Chromosome number where the window is located. |
| Start | Start position (1-based, inclusive) of the window. |
| End   | End position (inclusive) of the window. |
| Ref   | Reference population. |
| Tgt   | Target population. |
| Src   | Source population. |
| N(Variants) | Number of variants found in this genomic window. Sites with missing genotypes in any population are excluded. If ancestral allele information is provided, variants whose ancestral allele differ from both the reference and alternative allele in the VCF file are also excluded. |
| UXX or QXX | Number of uniquely shared sites or the quantile statistic in this window. |
| Candidate | Variant(s) (`chrom:pos`) that meet the specified condition. For the U statistic, this refers to variant(s) whose allele frequency in the reference population is below a specified threshold, allele frequency in the target population exceeds a specified threshold, and allele frequency in the source population meets a defined condition (e.g., fixed). For the Q statistic, this refers to variant(s) that satisfy the following three conditions: (1) the allele frequency in the reference population is below a specified threshold; (2) the allele frequency in the source population meets a defined condition (e.g., fixed); and (3) among all variants meeting (1) and (2), the allele frequency in the target population is greater than or equal to a specified percentile (e.g., the 95th percentile) within this subset of variants in the window. |

## Examples

The following example estimates the U statistic using biallelic single nucleotide polymorphisms (SNPs) in the region `chr9:16400001-16900000` from [the 1000 Genomes Project](https://ftp.ncbi.nih.gov/1000genomes/ftp/release/20130502/). The Neanderthal genome was obtained from [here](http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Altai/), and the Denisovan genome from [here](http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Denisova/). The reference population includes all African populations (ESN, GWD, LWK, MSL, YRI), excluding the admixed African populations (ACB and ASW), and all East Asian populations (CDX, CHB, CHS, JPT, KHV). The target population consists of all European populations (CEU, FIN, GBR, IBS, TSI). The source populations comprise one Neanderthal individual (AltaiNeandertal) and one Denisovan individual (Denisova).

To detect adaptively introgressed variants from **both Neanderthals and Denisovans**, we can use the following command:

```
sai score --vcf examples/data/1KG.nea_den.chr9.example.vcf.gz \
          --ref examples/data/1KG.ref.samples.txt \
          --tgt examples/data/1KG.tgt.samples.txt \
          --src examples/data/1KG.nea_den.samples.txt \
          --anc-alleles examples/data/hg19.chr9.example.anc.alleles.bed \
          --chr-name 9 --stat U50 --w 0.01 --y =1 =1 \
          --win-len 40000 --win-step 40000 \
          --output examples/results/1KG.nea_den.chr9.example.both.U50.scores.tsv
```

We provide a BED file to specify the ancestral allele for each variant, based on alignments from [Ensembl](https://ftp.ensembl.org/pub/release-75/fasta/ancestral_alleles/).
In this example, the `--w`, and `--y` arguments define the filtering conditions: the allele frequency of a variant must be less than 0.01 in the reference population and fixed (`=1`) in both Neanderthal and Denisovan genomes. The `50` in `--stat U50` indicates that the allele frequency must be greater than 0.5 in the target population. The output is shown above and can be found [here](https://github.com/xin-huang/sai/blob/main/examples/results/1KG.nea_den.chr9.example.both.U50.scores.tsv).

To estimate the Q95 statistic for the same data, simply change `--stat U50` to `--stat Q95`:

```
sai score --vcf examples/data/1KG.nea_den.chr9.example.vcf.gz \
          --ref examples/data/1KG.ref.samples.txt \
          --tgt examples/data/1KG.tgt.samples.txt \
          --src examples/data/1KG.nea_den.samples.txt \
          --anc-alleles examples/data/hg19.chr9.example.anc.alleles.bed \
          --chr-name 9 --stat Q95 --w 0.01 --y =1 =1 \
          --win-len 40000 --win-step 40000 \
          --output examples/results/1KG.nea_den.chr9.example.both.Q95.scores.tsv
``` 

The output can be found [here](https://github.com/xin-huang/sai/blob/main/examples/results/1KG.nea_den.chr9.example.both.Q95.scores.tsv).

To detect adaptively introgressed variants from **Neanderthal introgression only**, we can use:

```
sai score --vcf examples/data/1KG.nea_den.chr9.example.vcf.gz \
          --ref examples/data/1KG.ref.samples.txt \
          --tgt examples/data/1KG.tgt.samples.txt \
          --src examples/data/1KG.nea_den.samples.txt \
          --anc-alleles examples/data/hg19.chr9.example.anc.alleles.bed \
          --chr-name 9 --stat U50 --w 0.01 --y =1 =0 \
          --win-len 40000 --win-step 40000 \
          --output examples/results/1KG.nea_den.chr9.example.nea.specific.U50.scores.tsv
```

Here, we set `--y =1 =0` to require that the variant be fixed in the Neanderthal individual and absent in the Denisovan. This ensures that only Neanderthal-specific variants are considered as potential candidates for introgression. The output can be found [here](https://github.com/xin-huang/sai/blob/main/examples/results/1KG.nea_den.chr9.example.nea.specific.U50.scores.tsv).

To detect adaptively introgressed variants from **Denisovan introgression only**, we can use:

```
sai score --vcf examples/data/1KG.nea_den.chr9.example.vcf.gz \
          --ref examples/data/1KG.ref.samples.txt \
          --tgt examples/data/1KG.tgt.samples.txt \
          --src examples/data/1KG.nea_den.samples.txt \
          --anc-alleles examples/data/hg19.chr9.example.anc.alleles.bed \
          --chr-name 9 --stat U50 --w 0.01 --y =0 =1 \
          --win-len 40000 --win-step 40000 \
          --output examples/results/1KG.nea_den.chr9.example.den.specific.U50.scores.tsv
```

In this case, `--y =0 =1` requires that the variant be fixed in the Denisovan individual and absent in the Neanderthal, thereby identifying variants potentially introgressed from Denisovans only. The output can be found [here](https://github.com/xin-huang/sai/blob/main/examples/results/1KG.nea_den.chr9.example.den.specific.U50.scores.tsv).

Any number of source populations is supported, including one, two, or more. However, note that the order of values provided to `--y` must match the order of source populations as specified in the file given by `--src`. For example, in `examples/data/1KG.nea_den.samples.txt`, the first source population is Neanderthal and the second is Denisovan. Therefore, `--y =1 =0` specifies that the variant must be fixed in Neanderthals and absent in Denisovans, while `--y =0 =1` specifies the reverse.

Providing `--y` values in the wrong order (e.g., reversing the values without changing the order of the populations in `--src`) will result in incorrect interpretation of source-specific introgression signals.
