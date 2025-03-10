%% Plot supplement figure 8

addpath("../shared_functions")

% Plot supplement figure 7a 
virmen_cell = importdata("../saved_files/virmen_cell_5.mat");
tbt_cell = importdata("../saved_files/tbt_cell_5.mat");
summary_cell = importdata("../saved_files/t_types_summary_cell_5.mat");
nbins = 50;
linearise_x = true;

ex_offsets = [1.4953,1.4953,1.4953,1.4804,1.4804];

[RMSE_all,corr_vals,p_vals] = plot_biases_and_correlations_06122022(virmen_cell,tbt_cell,summary_cell,nbins,linearise_x,ex_offsets);
