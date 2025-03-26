# SAI

`sai` is a Python package for **S**tatistics for **A**daptive **I**ntrogression, which detects candidate regions of adaptive introgression from population genomic datasets. Currently, it supports the calculation of the number of uniquely shared sites (U statistic) and the quantile statistic (Q statistic) proposed by [Racimo et al. (2017)](https://doi.org/10.1093/molbev/msw216).

## Requirements

`sai` works on Linux operating systems and tested with the following:

    - natsort=8.4.0
    - numpy=1.26.4
    - pandas=2.2.1
    - pysam=0.23.0
    - python=3.9.19
    - pytest=8.1.1
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

## Help

To get help information, users can use:

```         
sai -h
```

This will display information for three commands:

| Command | Description |
| - | - |
| score | Run the score command based on specified parameters |
| outlier | Detect and output outlier rows based on quantile thresholds |
| plot | Generate a scatter plot of U vs Q statistics |

If you need further help, such as such as reporting a bug or suggesting a feature, please open an [issue](https://github.com/xin-huang/sai/issues).
