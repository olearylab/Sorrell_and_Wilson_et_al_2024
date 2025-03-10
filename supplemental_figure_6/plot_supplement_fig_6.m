%% Plot supplemental figure 6

% Set up workspace
addpath("../shared_functions")

virmen_cell = importdata("../saved_files/virmen_cell_5.mat");
tbt_cell = importdata("../saved_files/tbt_cell_5.mat");
all_popboot_mean = importdata("all_popboot_means_ball_norm.mat");
all_overall_n_means = importdata("all_overall_n_means_ball_norm.mat");
z_cell_CNN = importdata("../saved_files/z_cell_CNN_5_v2.mat");

centres = importdata("centres_50.mat");
centres = centres*0.74;

% Plot figure 6c + supplement 6a 
[h_boots_n] = plot_pop_mean_results_06122022(all_popboot_mean,all_overall_n_means,centres);

% plot supplement figure 6b
initialise_params;
av_num = 10;
run_many_assess_mean_time_change_18042023(virmen_cell,z_cell_CNN,tbt_cell,model_params,av_num);

% plot supplement figure 6c
plot_across_sessions_mean_change_18042023(all_overall_n_means);

va_corrs = importdata("va_dec_corrs_29122023_cut.mat");
yaw_corrs = importdata("yaw_dec_corrs_subs_29122023_cut.mat");

% plot supplement figure 6d
[h_boots] = plot_weights_changes_corrs(va_corrs,yaw_corrs);

full_pairs_means_cell = importdata("full_pairs_means_cell_ball_norm_nosub.mat");
% plot supplement figure 6e,f,g
[bootstats_center,bootstats_sem,h_boots_noise] = plot_many_noise_corrs(full_pairs_means_cell);