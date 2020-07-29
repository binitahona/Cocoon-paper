#!/usr/bin/env bash

#credit to Colas, Giacomo, and possibly others

set -e

# Set the conda environment name
export ENVIRONMENT_NAME= hawc_analysis

# Set the number of threads
export N_THREADS=4

# Create our conda environment installing the Fermi ST as well as threeML and the
# externals needed
conda create -y --name $ENVIRONMENT_NAME -c threeml -c conda-forge/label/cf201901 -c fermi fermitools threeml boost=1.63 cmake zeromq cppzmq healpix_cxx=3.31 pytest==3.9.3  matplotlib numba pyyaml==3.13 yaml==0.1.7 fermipy 

# Activate the conda environment
source activate $ENVIRONMENT_NAME

#This version of numpy may not be available via conda.
pip install numpy==1.15.3

pip install naima

# Install root_numpy making sure it is built against the installed version of ROOT
pip uninstall root_numpy
export MACOSX_DEPLOYMENT_TARGET=10.10
pip install --no-binary :all: root_numpy

#another package needed for HAL
pip install reproject

# Install HAL
pip uninstall hawc_hal -y ; pip install git+https://github.com/threeML/hawc_hal.git
 
# Write setup file
# NOTE: variables will be expanded
cat > $HOME/init_conda_hawc.sh <<- EOM
