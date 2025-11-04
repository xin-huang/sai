# Calculate Statistics

`sai` provides the `score` command to calculate different statistics.

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
| --win-len     | Length of each genomic window in base pairs. Default: 50,000. |
| --win-step    | Step size in base pairs between consecutive windows. Default: 10,000. |
| --anc-alleles | Path to the BED file with ancestral allele information. If ancestral allele information is not provided, filtering will be performed for each variant based on whether the allele frequency of any allele (assuming biallelic) meets the specified condition during the calculation of the statistics. Default: None. |
| --output      | Output file path for saving results. |
| --config      | Path to the YAML configuration file specifying the statistics to compute, ploidy settings, and population group file paths. |

## Input files

`sai` requires a biallelic VCF file containing all populations, including the reference population (assumed to have no introgression), the target population (assumed to have received introgressed material), and the source population (assumed to have provided the introgressed material). It also assumes that all variants observed in these populations are already included in the input VCF file. The program does not fill in missing genotypes (e.g., treating unobserved genotypes as homozygous reference, as done in [MaLAdapt](https://github.com/xzhang-popgen/maladapt/blob/fba58dac6278d1b4cc573c78ab7c56911b4072b0/empirical/1empirical_compute-features.py#L318-L326)). An example VCF file can be found [here](https://github.com/xin-huang/sai/blob/main/examples/data/1KG.nea_den.chr9.example.vcf.gz).

Three tab-delimited files are required to specify the population information for individuals in the reference, target, and source populations. An example file is shown below and can also be found [here](https://github.com/xin-huang/sai/blob/main/examples/data/1KG.nea_den.samples.txt).

```
NEA	AltaiNeandertal
DEN	Denisova
```

This file specifies the individuals in the source populations. In this example, there are two source populations, NEA and DEN (first column), and the second column lists the corresponding sample names matching those in the VCF header.

A configuration file in [YAML format](https://yaml.org/) is required to specify the statistics to compute, ploidy settings, and population group file paths. An example can be found [here](https://github.com/xin-huang/sai/blob/main/examples/data/1KG.nea_den.chr9.example.both.config.yaml). In the example below, the `statistics` section lists which statistics will be computed, while entries such as $U$ and $Q$ specify reference (`ref`), target (`tgt`), and source (`src`) populations with filtering conditions. The `ploidies` section defines the ploidy level for each population, typically 2 for diploids. The populations section provides file paths to sample lists that assign individuals to each population.

```
statistics:
  Danc: True
  DD: True
  Dplus: True
  df: True
  fd: True
  U:
    ref: 
      AFR+EAS: 0.01
    tgt: 
      EUR: 0.5
    src:
      NEA: "=1"
      DEN: "=1"
  Q:
    ref: 
      AFR+EAS: 0.01
    tgt: 
      EUR: 0.95
    src:
      NEA: "=1"
      DEN: "=1"

ploidies:
  ref:
    AFR+EAS: 2
  tgt:
    EUR: 2
  src:
    NEA: 2
    DEN: 2

populations:
  ref: "examples/data/1KG.ref.samples.txt"
  tgt: "examples/data/1KG.tgt.samples.txt"
  src: "examples/data/1KG.nea_den.samples.txt"
```

## Output files

An example output is shown below.

| Chrom | Start | End | Ref | Tgt | Src | Outgroup | N(Variants) | Danc.NEA | Danc.DEN | DD.NEA | DD.DEN | Dplus.NEA | Dplus.DEN | df.NEA | df.DEN | fd.NEA | fd.DEN | U | Q |
| -     | -     | -   | -   | -   | -   | -        | -           | -        | -        | -      | -      | -         | -         | -      | -      | -      | -      | - | - |
| 9 | 16400001 | 16440000 | AFR+EAS | EUR | NEA,DEN | NA | 1082 | -0.12355029565477203 | -0.09559566661094993 | -7.478853072043918 | -4.739448231247437 | -0.14512264357342178 | -0.09212831462235394 | -0.021122428026740657 | -0.005965214141977462 | -0.03412663364767112 | -0.011254553402083362 | 0 | nan |
| 9 | 16440001 | 16480000 | AFR+EAS | EUR | NEA,DEN | NA | 983 | 0.23853287087522393 | 0.09720568452718485 | 15.824720328820732 | 2.6519198302249976 | 0.24254139814192216 | 0.04198624490861517 | 0.06679895300702358 | -0.027612395797831944 | 0.09662599519731511 | -0.06072062935032934 | 0 | nan |
| 9 | 16480001 | 16520000 | AFR+EAS | EUR | NEA,DEN | NA | 985 | 0.11191758509638665 | 0.1371820707366662 | 4.845944176212569 | 6.488856543280008 | 0.08062408365129729 | 0.10781131967371557 | -0.0059518440634471394 | 0.006339448215182395 | -0.009777268672625093 | 0.011579864693126562 | 0 | nan |
| 9 | 16520001 | 16560000 | AFR+EAS | EUR | NEA,DEN | NA | 966 | 0.008890710141582193 | 0.0017874385890340458 | -3.2895208116381127 | -3.762141437091735 | -0.07409036881532824 | -0.08473527378923923 | -0.08168801665077327 | -0.0882036818431493 | -0.10989219575082872 | -0.11581875355584634 | 0 | nan |
| 9 | 16560001 | 16600000 | AFR+EAS | EUR | NEA,DEN | NA | 1029 | 0.0024719610094924406 | 0.013019347645432937 | -3.0949185837356765 | -2.3295899207926993 | -0.07022371194458008 | -0.05285840228768851 | -0.024174513563260166 | -0.021736551375934953 | -0.05323265434587958 | -0.051743268070762606 | 0 | 0.0008946322067594431 |
| 9 | 16600001 | 16640000 | AFR+EAS | EUR | NEA,DEN | NA | 1063 | 0.022253891584592442 | -0.05846778040121223 | 3.258157404777677 | -3.199939253368683 | 0.07998846049314343 | -0.06970294049702096 | 0.025940599557366877 | -0.0056209360304606635 | 0.10418231535859726 | -0.01922836158621726 | 0 | 0.0 |
| 9 | 16640001 | 16680000 | AFR+EAS | EUR | NEA,DEN | NA | 1153 | 0.10931478521694758 | -0.02727035560668536 | 11.405249061188428 | -2.9954578087033354 | 0.17733519261177239 | -0.04657505369825686 | 0.07092801200718006 | -0.01271958977736748 | 0.14831905506347323 | -0.0603694010038118 | 0 | nan |
| 9 | 16680001 | 16720000 | AFR+EAS | EUR | NEA,DEN | NA | 999 | 0.15113133237922183 | 0.14710421254212552 | 18.309291078923295 | 17.05495597841523 | 0.25761194634556367 | 0.24391495232933705 | 0.1898707100981595 | 0.17809065689988668 | 0.22245899486649723 | 0.3768632297159813 | 1 | 0.6282306163021869 |
| 9 | 16720001 | 16760000 | AFR+EAS | EUR | NEA,DEN | NA | 1007 | 0.4129294073240463 | 0.20572780178794584 | 53.56721772539367 | 29.49165522736596 | 0.49454474017007327 | 0.2692256133787326 | 0.45835432326244324 | 0.2539470613190382 | 0.5958374858338679 | 0.354370074530527 | 2 | 0.6828528827037773 |
| 9 | 16760001 | 16800000 | AFR+EAS | EUR | NEA,DEN | NA | 1060 | 0.40683688578049354 | 0.1438067334121798 | 40.19722340559816 | 16.581497128341063 | 0.5263549295804029 | 0.2153953867869153 | 0.5848435316604244  | 0.21982990433919478 | 0.5444309862037052 | 0.23506396309346742 | 1 | 0.6988071570576541 |
| 9 | 16800001 | 16840000 | AFR+EAS | EUR | NEA,DEN | NA | 960 | 0.40565740562034003 | 0.22930692442169218 | 24.16718735207801 | 13.914000126226764 | 0.5158532156040262 | 0.29611831851629133 | 0.3520651233378468 | 0.1984114833463481 | 0.6411180111472388 | 0.5318202119747509 | 1 | 0.7176938369781312 |
| 9 | 16840001 | 16880000 | AFR+EAS | EUR | NEA,DEN | NA | 923 | 0.21619895332973812 | 0.19275643241641574 | 12.284513159140396 | 10.489264807977534 | 0.22271722145226366 | 0.19016951529508055 | 0.09059952770351988 | 0.058657408877032656 | 0.08230019337435454 | 0.06443318954095166 | 0 | 0.15382703777335985 |
| 9 | 16880001 | 16920000 | AFR+EAS | EUR | NEA,DEN | NA | 512 | 0.3881193782134191 | 0.3589872405156586 | 20.739761825870175 | 19.506952333617342 | 0.48085645624316853 | 0.4522734663012016 | 0.5199698990819357  | 0.5044977777777778 | 0.47736503902045885 | 0.55657933398083 | 0 | nan |

Each row in the output file corresponds to a genomic window, summarizing the result for that region. The columns are as follows:

| Column | Description |
| - | - |
| Chrom | Chromosome number where the window is located. |
| Start | Start position (1-based, inclusive) of the window. |
| End   | End position (inclusive) of the window. |
| Ref   | Reference population. |
| Tgt   | Target population. |
| Src   | Source population. |
| Outgroup | Outgroup population. |
| N(Variants) | Number of variants found in this genomic window. Sites with missing genotypes in any population are excluded. If ancestral allele information is provided, variants whose ancestral allele differ from both the reference and alternative allele in the VCF file are also excluded. |
| Danc.src1 | Value of the $D_{anc}$ statistic using src1 as the source population. |
| Danc.src2 | Value of the $D_{anc}$ statistic using src2 as the source population. |
| DD.src1 | Value of the $D_D$ statistic using src1 as the source population. |
| DD.src2 | Value of the $D_D$ statistic using src2 as the source population. |
| Dplus.src1 | Value of the $D^+$ statistic using src1 as the source population. |
| Dplus.src2 | Value of the $D^+$ statistic using src2 as the source population. |
| df.src1 | Value of the $d_f$ statistic using src1 as the source population. |
| df.src2 | Value of the $d_f$ statistic using src2 as the source population. |
| fd.src1 | Value of the $f_d$ statistic using src1 as the source population. |
| fd.src2 | Value of the $f_d$ statistic using src2 as the source population. |
| U | Value of the $U$ statistic using both src1 and src2 as the source populations. |
| Q | Value of the $Q$ statistic using both src1 and src2 as the source populations. |

## Examples

The following example estimates all the statistics defined by the example configuration above with biallelic single nucleotide polymorphisms (SNPs) in the region `chr9:16400001-16900000` from [the 1000 Genomes Project](https://ftp.ncbi.nih.gov/1000genomes/ftp/release/20130502/). The Neanderthal genome was obtained from [here](http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Altai/), and the Denisovan genome from [here](http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Denisova/). The reference population includes all African populations (ESN, GWD, LWK, MSL, YRI), excluding the admixed African populations (ACB and ASW), and all East Asian populations (CDX, CHB, CHS, JPT, KHV). The target population consists of all European populations (CEU, FIN, GBR, IBS, TSI). The source populations comprise one Neanderthal individual (AltaiNeandertal) and one Denisovan individual (Denisova).

To detect adaptively introgressed variants from **both Neanderthals and Denisovans**, we can use the following command:

```
sai score --vcf examples/data/1KG.nea_den.chr9.example.vcf.gz \
          --anc-alleles examples/data/hg19.chr9.example.anc.alleles.bed \
          --chr-name 9 \
          --win-len 40000 --win-step 40000 \
          --config examples/data/1KG.nea_den.chr9.example.both.config.yaml \
          --output examples/results/both/1KG.nea_den.chr9.example.both.stats.scores.tsv
```

We provide a BED file to specify the ancestral allele for each variant, based on alignments from [Ensembl](https://ftp.ensembl.org/pub/release-75/fasta/ancestral_alleles/).
In this example, for the $U$ statistic, the `ref` entry in the configuration file specifies that the allele frequency must be less than 0.01 in the AFR+EAS population, the `src` entries require the allele to be fixed (`=1`) in both Neanderthal and Denisovan genomes, and the `tgt` entry indicates that the allele frequency must be greater than 0.5 in the EUR population. Similarly, for the $Q$ statistic, the `tgt` entry specifies the quantile of the derived allele frequencies in the EUR population at 0.95. The output is shown above and can be found [here](https://github.com/xin-huang/sai/blob/main/examples/results/both/1KG.nea_den.chr9.example.both.stats.scores.tsv). For $D^+$, $D_{anc}$, $d_f$ and $f_d$, if no outgroup is provided, then the derived allele frequency is assumed to be 1 in the outgroup.

The `src` entries can be set as Nea: "=1", Den: "=0" or the reverse to identify lineage-specific variants:

```
statistics:
  U:
    ref: 
      AFR+EAS: 0.01
    tgt: 
      EUR: 0.5
    src:
      NEA: "=1"
      DEN: "=0"
  Q:
    ref: 
      AFR+EAS: 0.01
    tgt: 
      EUR: 0.95
    src:
      NEA: "=1"
      DEN: "=0"
```

An example configuration file for Neanderthal-specific introgressed variants can be found [here](https://github.com/xin-huang/sai/blob/main/examples/data/1KG.nea_den.chr9.example.nea.specific.config.yaml), with its corresponding output available [here](https://github.com/xin-huang/sai/blob/main/examples/results/nea_specific/1KG.nea_den.chr9.example.nea.specific.stats.scores.tsv). Similarly, an example configuration file for Denisovan-specific introgressed variants is provided [here](https://github.com/xin-huang/sai/blob/main/examples/data/1KG.nea_den.chr9.example.den.specific.config.yaml), with the output available [here](https://github.com/xin-huang/sai/blob/main/examples/results/den_specific/1KG.nea_den.chr9.example.den.specific.stats.scores.tsv). Note that only the $U$ and $Q$ statistics support multi-source introgression filtering.

In addition, log files corresponding to the $U$ or $Q$ statistics are generated to record the positions of variants used in their calculation. An example is shown below:

| Chrom | Start | End | U_SNP |
| - | - | - | - | 
| 9 | 16400001 | 16440000 | NA |
| 9 | 16440001 | 16480000 | NA |
| 9 | 16480001 | 16520000 | NA | 
| 9 | 16520001 | 16560000 | NA | 
| 9 | 16560001 | 16600000 | NA |
| 9 | 16600001 | 16640000 | NA |
| 9 | 16640001 | 16680000 | NA |
| 9 | 16680001 | 16720000 | 9:16715657 |
| 9 | 16720001 | 16760000 | 9:16758710,9:16759357 |
| 9 | 16760001 | 16800000 | 9:16770229 |
| 9 | 16800001 | 16840000 | 9:16800789 |
| 9 | 16840001 | 16880000 | NA |
| 9 | 16880001 | 16920000 | NA |
