# GTSE1-14A Astral Microtubule Analysis

<strong>Destabilization of long astral microtubules via Cdk1-dependent removal of GTSE1 from their plus-ends facilitates prometaphase spindle orientation</strong>  
Divya Singh, Nadine Schmidt, Franziska Müller, Tanja Bange, and Alexander W. Bird  
Under revision at Current Biology, 2020  
[bioRxiv](https://www.biorxiv.org/content/10.1101/2020.05.23.111989v1.full)

This repository contains code for analyzing astral microtubules as presented in the paper.

- Histograms and dot plots of astral length
- Histograms and dot plots of astral counts
- Histograms of microtubules touching the cell cortex

----

- <code>read_comets_file.m</code>: Read a CSV file containing EB1 comet coordinates  
- <code>read_poles_file.m</code>: Read a CSV file containing spindle pole coordinates  
- <code>create_astral_lengths.m</code>: Computes Euler distance spindle pole and EB1 comet and classifies between inner spindle or astral microtubules  
- <code>compare_cell_lines.m</code>: Plots histograms of astral length or number of astrals among different conditions  
- <code>astrals_over_threshold.m</code>: Counts the number of astrals above a particular length threshold
- <code>line_scan_cortical.m</code>: Counts the number of microtubules touching the cell cortex in prometaphase
- <code>line_scan_intensities_plot.m</code>: Plots intensities of EB3 and GTSE1 at microtubule-plus ends in different cell cycle phases

----

