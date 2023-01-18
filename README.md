# High-Throughput Cellular Biology

### Requirements
- `python 3.11`
- `requirements.txt`

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
├── requirements.txt
├── LICENSE
├── Makefile
└── README.md
```

The series of notebooks in the folder `notebooks` shows the plots for Lauren Halger's paper: «high-throughput cellular biology».

The plotting functions and the data procesing methods are stored in the `src` folder.

The data is stored in the folder `data`.

Figures generated in the google colab are automatically outputted in the folder `figs`.
