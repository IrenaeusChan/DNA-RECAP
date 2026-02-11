# DNA-RECAP (DNA Repair Event Classification, Analysis, and Processing)

## Overview
This tool was written to complement CIBAR-seq, a Cas9-Induced Break and Repair assay used to explore how RPE cells handle induced double-strand breaks.
Details for this method and manuscript can be found here: <insert CIBAR-seq manuscript>

## System Requirements
#### OS Requirements
This code is supported for all operating systems that can run python v3.12.6 or above. The code found in this repository has been tested on macOS:
- macOS: Sequoia v15.5
- Linux: Ubuntu v20.04.6

#### Python Packages:
```
args==0.1.0
click==8.0.3
clint==0.5.1
numpy==1.24.3
pandas==2.0.0
python-dateutil==2.8.2
numba==0.59.0
```

## Installation
#### Installing Python
```
wget https://www.python.org/ftp/python/3.12.6/Python-3.12.6.tgz \
    && tar -xzf Python-3.12.6.tgz \
    && cd Python-3.12.6 \
    && ./configure --prefix=/opt/python-3.12.6 \
    && make \
    && make install \
PATH="/opt/python-3.12.6:$PATH"
```

#### Install using PIP
```
sudo pip3 install git+https://github.com/IrenaeusChan/DNA-RECAP.git@main
```

#### Install from Github
```
git clone https://github.com/IrenaeusChan/DNA-RECAP
cd DNA-RECAP
python3 setup.py install
```
Please use sudo if required

#### Install using conda
```
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -e .
```

## Usage

## License

