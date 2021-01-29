# Fractional diffusion theory of balanced heterogeneous neural networks

This repository houses the companion code for the article [*Fractional diffusion theory of balanced heterogeneous neural networks*](https://doi.org/10.1103/PhysRevResearch.3.013083).

## Citation
Asem Wardak and Pulin Gong, *Phys. Rev. Research* **3**, 013083 (2021)

## Requirements
The simulations require only Python and NEST, while some analyses additionally require MATLAB and Mathematica; see the file descriptions below.

## Usage
The file `wardakgong2021_sim.py` sets parameters for and simulates the network model functions in `wardakgong2021_networkmodel.py`.
The resultant simulation data is analysed in `wardakgong2021_analysis.py` and `nb/wardakgong2021_networkquantities.nb`.

The files `wardakgong2021_[defs|calcs|plots].nb` contain analytical calculations involving the stationary solution to the fractional Fokker-Planck equation, while `wardakgong2021_montecarlo.nb` contains code used to obtain, via the Monte-Carlo method, the mean-field analytical firing rate distribution.
