# How to run biological neural network simulations

## Getting Started

### Prerequisites

Anaconda Python 4.3.13 (Python 2.7)
Neuron 7.4 (https://neuron.yale.edu/neuron/)
Brain Modeling Toolkit (https://alleninstitute.github.io/bmtk/)

### Installation

The first step is to install the Brain Modeling Toolkit Package together with necessary instructions:
https://alleninstitute.github.io/bmtk/installation.html

### Running simulations and visualizing the results

To run simulations:

1) Copy the files from Network_description.txt to the "network" folder.
2) Run simulations in the bash shell: 
'''python run_bionet.py config.json'''
2.1) To run simulations in the parallel mode, use the following instructions:
'''https://alleninstitute.github.io/bmtk/bionet'''
3) Once simulations are finished, run the Jupyter notebook 'analyze_simulation.ipynb' to vizualise the results

## Contributors

* **Anatoly Buchin** - [abuchin](https://github.com/abuchin)

## License

Copyright 2020. Allen Institute. All rights reserved

## Level of support

We are planning on occasional updating this tool with no fixed schedule. Community involvement is encouraged through both issues and pull requests.
