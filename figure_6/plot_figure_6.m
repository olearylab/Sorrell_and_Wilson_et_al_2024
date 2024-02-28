%% Plot figure 6

% Set up workspace
addpath("../shared_functions")

virmen_cell = importdata("../saved_files/virmen_cell_5.mat");
tbt_cell = importdata("../saved_files/tbt_cell_5.mat");
all_popboot_mean = importdata("all_popboot_means_ball_norm.mat");
all_overall_n_means = importdata("all_overall_n_means_ball_norm.mat");
z_cell_CNN = importdata("../saved_files/z_cell_CNN_5_v2_cut.mat");

accuracies_cell_bvn = importdata("accuracies_cell_neural_50_lrsep_100subs.mat");
accuracies_cell_shuff = importdata("accuracies_cell_neural_50_lr_sep_shuff_true.mat");

sig_results_all = importdata("sig_results_all_cut.mat");
sig_summary_cell = importdata("sig_summary_cell_cut.mat");
all_full_results_both = importdata("all_full_results_both_cut.mat");
full_all_params_means = importdata("full_all_params_means_cut.mat");
sig_check_all = importdata("sig_check_all_cut.mat");

initialise_params;
linearise_x = true;
bootsamps = 100;
nbins = 50;
xnum = 6;

centres = importdata("centres_50.mat");
centres = centres*0.74;

zdata = z_cell_CNN{2,2};
virmen_data = virmen_cell{2,2};
tbt_details = tbt_cell{2,2};
sig_results = []; % Not restricting to neurons with signficant peaks
plot_res = true;
kept_types = [1,3];

% Plotting

% Plot figure 6a
[h_boots_svm] = plot_many_svm_bvn_neural_classifier_results(accuracies_cell_bvn,accuracies_cell_shuff,centres);

% Plot figure 6b
normalise_z = false;
sub_sample = false;
[all_sorted_means,all_sorted_means_joined,I_store,cur_neurons_set] = sort_neurons_plot_07122022(zdata,virmen_data,tbt_details,sig_results,linearise_x,nbins,bootsamps,xnum,sub_sample,normalise_z,plot_res,kept_types);
% Cross validated sorting with half, then plotting for other half
[all_sorted_means,all_sorted_means_joined,I_store,cur_neurons_set] = sort_neurons_plot_crossval(zdata,virmen_data,tbt_details,sig_results,linearise_x,nbins,bootsamps,xnum,sub_sample,normalise_z,plot_res,kept_types);

% Plot figure 6c + supplement 6.1 a 
[h_boots_n] = plot_pop_mean_results_06122022(all_popboot_mean,all_overall_n_means,centres);

% Pick example neurons for plotting
num_examples = 2;
normalise_z = false;
% Get bootstrapped tuning curves for example session
[bootstrap_means, centres] = calc_bootstrap_from_raw_subsample_norm(z_cell_CNN{1,2}, virmen_cell{1,2}, tbt_cell{1,2}, linearise_x, nbins,bootsamps, xnum,sub_sample,normalise_z);
ex_neu = [58,66;164,18;36,110;34,134;17,50;10,24];
ex_dir = [1,1;1,2;1,1;1,1;2,2;1,1];
ex_md = [1,2];

% Plot figure 6 d and e
[h_boots] = plot_param_change_figure_18052023(sig_results_all,sig_summary_cell,all_full_results_both,full_all_params_means,sig_check_all,bootstrap_means,centres,ex_neu,ex_dir,ex_md);

% Plot figure 6f 
z_binned_res_cell = importdata("z_binned_res_cell_ball_norm.mat");
error_cell = importdata("../figure_5/error_cell_norm.mat");
x_binned_cell = importdata("../figure_5/x_binned_cell.mat");
error_thresh = 1;
bin_space = 5;
num_shuff = 1000;
shuff_limit = 5;
[h_boot,ex_md,ro,ro_p,n_kept,stats_xc] = plot_neuron_residual_correction_correlations_spaced_bins(z_binned_res_cell,error_cell,tbt_cell,x_binned_cell,error_thresh,bin_space,num_shuff,shuff_limit);

% Plot figure 6g
x_binned_cell_from_pix_res = importdata("x_binned_cell_from_pix_res.mat");
plot_ind = 8;
[h_boots_res] = run_many_residuals_correct_corrs(x_binned_cell_from_pix_res,error_cell,tbt_cell,centres,plot_ind,error_thresh);

