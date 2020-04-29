# How to run biological neural network simulations

## Getting Started

### Prerequisites

Anaconda Python 4.3.13 (Python 2.7)
Neuron 7.4 (https://neuron.yale.edu/neuron/)
Brain Modeling Toolkit: workshop2018 version (https://github.com/AllenInstitute/bmtk/tree/release/workshop2018/docs/tutorial)

### Requirements

Enviroment requirements could be found in requirements folder: requirements.txt and epilepsy_human_dg_network.yml file

To install the enviroment run:
'''conda env create -f epilepsy_human_dg_network.yml'''

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

Copyright © 2020. Allen Institute. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Redistributions and use for commercial purposes are not permitted without the Allen Institute’s written permission. For purposes of this license, commercial purposes are the incorporation of the Allen Institute's software into anything for which you will charge fees or other compensation or use of the software to perform a commercial service for a third party. Contact terms@alleninstitute.org for commercial licensing opportunities.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## Level of support

We are planning on occasional updating this tool with no fixed schedule. Community involvement is encouraged through both issues and pull requests.
