# SAI

`sai` is a Python package for **S**tatistics for **A**daptive **I**ntrogression. It detects candidate regions of adaptive introgression from population genomic datasets. Currently, it supports:

- The average difference of sequence divergence ($D_D$ statistic) proposed by [Huang et al. (2025)](https://doi.org/10.1093/molbev/msaf295).
- The $D^+$ and $D_{anc}$ statistics proposed by [Fang et al. (2024)](https://doi.org/10.1371/journal.pgen.1010155).
- The distance fraction ($d_f$ statistic) proposed by [Pfeifer and Kapan (2019)](https://doi.org/10.1186/s12859-019-2747-z).
- The number of uniquely shared sites ($U$ statistic) and the quantile of the derived allele frequencies in such sites ($Q$ statistic) proposed by [Racimo et al. (2017)](https://doi.org/10.1093/molbev/msw216).
- The dynamic estimator of the proportion of introgression ($f_d$ statistic) proposed by [Martin et al. (2015)](https://doi.org/10.1093/molbev/msu269).

`sai` does not require phased data, and supports an arbitrary number of source/donor populations (i.e., populations assumed to provide introgressed material) and arbitrary ploidy.

## Requirements

`sai` works on Linux operating systems and tested with the following:

- matplotlib=3.9.1
- natsort=8.4.0
- numpy=1.26.4
- pandas=2.2.1
- pydantic=2.11.7
- pysam=0.23.0
- python=3.9.19
- pytest=8.1.1
- pytest-cov=6.0.0
- scikit-allel=1.3.7
- scipy=1.12.0

## Installation

Users can first install [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html), and then install `sai` using the following commands:

```
git clone https://github.com/xin-huang/sai
cd sai
mamba env create -f build-env.yaml
mamba activate sai
pip install .
```

Users can also install `sai` from [PYPI](https://pypi.org/):

```
pip install sai-pg
```

## Help

To get help information, users can use:

```         
sai -h
```

This will display information for two commands:

| Command | Description |
| - | - |
| score | Run the score command based on specified parameters |
| outlier | Detect and output outlier rows based on quantile thresholds |

If you need further help, such as such as reporting a bug or suggesting a feature, please open an [issue](https://github.com/xin-huang/sai/issues).
