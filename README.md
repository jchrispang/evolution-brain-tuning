# evolution-brain-tuning
Code and data for the manuscript "Evolutionary shaping of human brain dynamics"

Using computational modeling to study changes in the tuning of brain dynamics across different species

## File descriptions

1. +models: package folder containing different brain network model functions
2. +utils: package folder containing various utility analysis and visualization functions
3. artworks: folder containing silhouette artworks of different species used in the paper
4. data: folder containing source data to reproduce the main figures of the paper
5. demo_simulation.m: code to demonstrate how to simulate different models
6. generate_paper_figures.m: code to generate the main figures of the paper
7. generate_paper_suppfigures.m: code to generate the supplementary figures of the paper

## Installation

Download the repository and add to your Matlab path

## Dependencies

Some important aspects you need to do before running generate_paper_X.m

1. Download [cbrewer](https://au.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab) and add to Matlab path.
2. Download the latest version (2019) of the [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet/) and add to Matlab path.
3. Download the gifti [library](https://github.com/gllmflndn/gifti) and add to Matlab path.
4. Download Stuart Oldham's [repository](https://github.com/StuartJO/plotSurfaceROIBoundary) for drawing ROI boundaries on a surface and add to Matlab path. 

## Compatibility

The codes have been tested on versions of MATLAB from R2017a to R2020b.

## Note on running generate_paper_X.m

Note that the configuration of your computer (e.g., screen resolution) affects how the figures created by running generate_paper_figures.m and generate_paper_suppfigures.m will look. Hence, they will not 100% visually match the figures in the paper, but the scientific contents are replicated.

## Citation (TO BE UPDATED)

If you use our code in your research, please cite us as follows:

J.C. Pang, J.K. Rilling, J.A. Roberts, M.P. van den Heuvel, L. Cocchi, Evolutionary shaping of human brain dynamics, (2022) (DOI: XXX)

## Further details

Please contact james.pang1@monash.edu if you need any further details.
