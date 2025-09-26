# zero-phase-lag
To cite this code: https://doi.org/10.5281/zenodo.17209298
Code used for manuscript: https://doi.org/10.1101/2025.01.04.631256

# Folders:
Phase lag matrices
- phase_lag_matrices.m - overall script to take pre-processed EEG data and produce matrices of phase lag between region-pairs in each frequency band
- lag_undirected.m - calculates the mean phase difference between pairs of ROIs
- PLF_lag_padding.m - computes the phase difference between region pairs

Functional Connectivity Methods Code
- functional_networks_AEC_orthAEC.m - takes pre-processed EEG data and produces functional connectivity matrices
- functional_networks_PLV_wPLI.m - takes pre-processed EEG data and produces functional connectivity matrices
- functional_networks_coh_img_coh.m - takes pre-processed EEG data and produces functional connectivity matrices

Quantifying zero-phase-lag functional connectivity
- histogram_lag_every_subject.m - This script produces histograms of phase distribution between -pi and +pi radians. The histogram combines FC across all subjects.
- near_zero_phase_versus_distance.m - Produces table of proportion of near-zero-phase-lag FC vs Euclidean distance per participant in output

Performance as FC derived biomarkers
- adjacency_matrix_ICC.m - produces a figure of test-retest reliability (ICC) for each method in each freq band
- mean_matrix_plot_EEG_DWI.m - plots mean EEG and DWI matrices across participants
- DWI_EEG_correlation.m - correlates DWI matrix with EEG matrices with each method in each frequency band + contains code to obtain WM tract distance and proxy of WM transmission time for post hoc analyses
- select_ROI_cognition.m - takes the whole FC adjacency matrix and only keeps ROIs at are based on spatial working memory literature. ROI numbers and literature references included
- shortest_pathlength_cognition.m - produces normalised, weighted shortest pathlength from adjacency matrices for each FC method, in 4-8 Hz
- mean_strength_age_prediction.m - produces table of mean strength for each FC method, in each frequency band

Post hoc analysis
- strength_time_penalisation_variables.m - produces streamline count, proxy for WM signal transmission time, penalisability of phase delay (penalisation scale), in preparation for post hoc mediation analysis
- lag_FC_strength_relationship.m - looks at the relationship between phase lag and strength of functional connections

Package for diffusion imaging analysis: https://www.mr-startrack.com/ 
