# zero-phase-lag
Code used for manuscript: https://doi.org/10.1101/2025.01.04.631256

# Folders:
Quantifying zero-phase-lag functional connectivity
- histogram_lag_every_subject.m - This script produces histograms of phase distribution between -pi and +pi %radians. The histogram combines FC across all subjects.
- near_zero_phase_versus_distance.m - Produces table of proportion of near-zero-phase-lag FC vs Euclidean distance per participant in output

Performance as FC derived biomarkers
- adjacency_matrix_ICC.m - produces a figure of test-retest reliability (ICC) for each method in each freq band
- mean_matrix_plot_EEG_DWI.m - plots mean EEG and DWI matrices across participants
- DWI_EEG_correlation.m - correlates DWI matrix with EEG matrices with each method in each frequency band + contains code to obtrain WM tract distance and proxy of WM transmission time for post hoc analyses
- select_ROI_cognition.m - takes the whole FC adjacency matrix and only keeps ROIs at are based on spatial working memory literature. ROI numbers and literature references included
- shortest_pathlength_cognition.m - 

Post hoc analysis
- DWI_EEG_correlation.m - correlates DWI matrix with EEG matrices with each method in each frequency band + contains code to obtrain WM tract distance and proxy of WM transmission time for post hoc analyses
