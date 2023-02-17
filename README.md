# High-Throughput Cellular Biology

## Getting Started
### Requirements
- `python 3.11`
- `requirements.txt`

**If using a Mac**

- `command line developer tools`

### Installation

**Using pyenv**

```
cd [path_to_your_local_folder]
git clone https://github.com/yvesmartindestaillades/highthroughputcellularbiology
cd highthroughputcellularbiology
python3.11 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
mkdir figs
```

**Using conda**

```
conda create -n [environment_name] python=3.11
conda activate [environment_name]
git clone https://github.com/yvesmartindestaillades/highthroughputcellularbiology
cd highthroughputcellularbiology
conda install --file requirements.txt
mkdir figs
```

### Download the data
```
cd [path_to_your_local_folder]/highthroughputcellularbiology
mkdir data
wget "https://drive.google.com/u/2/uc?id=1BB0kFP8v6oxVoDzSlG_E-c8rxv4jBWSj&amp;export=download&amp;confirm=t&amp;uuid=791ead2c-e9ef-4dbe-b4b5-9e9216c3c47c&amp;at=ACjLJWmLtrTCfcic3eaawMjLUlUw:1674012296979&confirm=t&uuid=9c0d55e2-911e-432d-becd-2057489390b7&at=ACjLJWnkUr3NbyEmgdAQ0Pg5LKNJ:1674012359603" -P data
mv "data/uc?id=1BB0kFP8v6oxVoDzSlG_E-c8rxv4jBWSj&amp;export=download&amp;confirm=t&amp;uuid=791ead2c-e9ef-4dbe-b4b5-9e9216c3c47c&amp;at=ACjLJWmLtrTCfcic3eaawMjLUlUw:1674012296979&confirm=t&uuid=9c0d55e2-911e-432d-becd-2057489390b7&at=ACjLJWnkUr3" data/df.feather
```

You should get the following:
```
/highthroughputcellularbiology
├── data
│   └── df.feather
├── figs
├── notebooks
├── src
├── (venv)
├── requirements.txt
├── LICENSE
├── Makefile
└── README.md
```

## Usage

### Notebook organization

Each notebook is a self-contained script that can be run independently. The notebooks are organized in the following order:

**Computational Pipeline**
1. Demultiplexing: `i_demultiplexing.ipynb`, keys values and plots about the demultiplexing. 
2. DMS-MaPseq: `ii_DMSMaPseq.ipynb`, mutations per read, base coverage, population average.
3. Error estimation: `iv_error_estimation.ipynb`, messy thoughts about the error estimation.

**Results**
1. Barcodes: `iii_barcodes.ipynb`, study of barcode mutations and barcode replicates.
2. Reproducibility: `iv_reproducibility.ipynb`, barcode replicates and biological replicates.
3. Signal is proportional to fraction folded: `v_proportional.ipynb`, mutation rates vs DMS concentration, temperature, and free energy.

**K-fold**
1. K-fold: `kfold.ipynb`, mutation rates heatmap, confidence interval on the k-fold estimation.

### How to run the notebooks
1. Open the notebook in your favorite editor (e.g. VSCode, PyCharm, Jupyter Notebook, etc.)
2. Run the first cell to import the libraries and the functions
3. Some plots maybe require some inputs, e.g a sample, a reference, etc. Just change the value as you need.

Ex: 

```
###### Change sample ######
sample = '01_02_S23_reads'   
###########################    

plot = plots.barcode_comparison_scatter_plot(study, sample)
```
### How to get the list of samples, families, references, sections etc
The list of samples, families, and references are stored in the handbook `notebooks/misc/handbook.ipynb`. Just print the corresponding cell to get the list.

### How to extract the data from the dataframe
DREEM's object `study` provides a filtering function called `get_df`. This function returns a dataframe with the filtered data. Check `notebooks/misc/handbook.ipynb`. for more details.

### How to extract the data from a figure
1. Run the corresponding cell.
2. The data is stored in the variable `plot['data']`.
3. You can save the data in a csv file using the following code:

```
import pandas as pd
pd.DataFrame(plot['data']).to_csv('data.csv', index=False)
```

### How to download the figures
Figures generated in the google colab are automatically outputted in the folder `figs`, and in the subfolder corresponding to the notebook.
