%% Plot figure 5

% Set up workspace for plotting
addpath("../shared_functions")

virmen_cell = importdata("../saved_files/virmen_cell_5.mat");
tbt_cell = importdata("../saved_files/tbt_cell_5.mat");
z_binned_cell = importdata("z_binned_cell_ball_norm_cut.mat");
error_cell = importdata("error_cell_norm.mat");
x_binned_cell = importdata("x_binned_cell.mat");
z_correction_corrs_cell = importdata("z_correction_corrs_cell_space5_cut.mat");
z_correction_corrs_p_cell = importdata("z_correction_corrs_p_cell_space5_cut.mat");

nbins = 50;
ex_md = [5,3];
error_thresh = 1;

% Plot figure 5a
bin_space = 5;
[n_ind,n_samps] = plot_neuron_correction_corr_spaced_bins_example(z_binned_cell,x_binned_cell,error_cell,tbt_cell,z_correction_corrs_cell,ex_md,error_thresh,bin_space);

% Get example correlation and p-value for selected neuron
ex_corr = z_correction_corrs_cell{ex_md(1),ex_md(2)}(n_ind);
ex_p = z_correction_corrs_p_cell{ex_md(1),ex_md(2)}(n_ind);

% Plot figure 5b
num_shuff = 1000;
shuff_limit = 5;
[h_boot,ex_md,ro,ro_p,n_kept] = plot_neuron_correction_correlations_spaced_bins(z_binned_cell,error_cell,tbt_cell,x_binned_cell,error_thresh,bin_space,num_shuff,shuff_limit);

all_res_cell = importdata("all_res_cell_corrective_pixels.mat");

% Plot figure 5c
[h_boots] = plot_many_corrective_decoding_check_subs(all_res_cell);

Wout_online_cell = importdata("Wout_online_cell.mat");
bmi_weights_cell = importdata("bmi_weights_cell_corrective_pixels_full.mat");

% Plot figure 5 d, e and f
[dot_prods,dot_prods_yaw,h_boots_w] = plot_va_yaw_weights_compare_subs(Wout_online_cell,bmi_weights_cell);