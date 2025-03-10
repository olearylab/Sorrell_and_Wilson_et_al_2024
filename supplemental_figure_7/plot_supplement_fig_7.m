%% Plot supplement figure 7
addpath("../shared_functions")

% Plot supplement figure 7a 
tbt_cell = importdata("../saved_files/tbt_cell_5.mat");
z_binned_res_cell = importdata("z_binned_res_cell_ball_norm.mat");
error_cell = importdata("../figure_5/error_cell_norm.mat");
x_binned_cell = importdata("../figure_5/x_binned_cell.mat");
centres = importdata("centres_50.mat");
centres = centres*0.74;

error_thresh = 1;
bin_space = 5;
num_shuff = 1000;
shuff_limit = 5;
[h_boot,ex_md,ro,ro_p,n_kept,stats_xc] = plot_neuron_residual_correction_correlations_spaced_bins(z_binned_res_cell,error_cell,tbt_cell,x_binned_cell,error_thresh,bin_space,num_shuff,shuff_limit);

% Plot supplement figure 7b
x_binned_cell_from_pix_res = importdata("x_binned_cell_from_pix_res.mat");
plot_ind = 8;
[h_boots_res] = run_many_residuals_correct_corrs(x_binned_cell_from_pix_res,error_cell,tbt_cell,centres,plot_ind,error_thresh);
