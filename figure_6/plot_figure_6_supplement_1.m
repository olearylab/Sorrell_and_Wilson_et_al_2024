%% Plot supplemental figure 6.1

% Set up workspace
addpath("../shared_functions")

virmen_cell = importdata("../saved_files/virmen_cell_5.mat");
tbt_cell = importdata("../saved_files/tbt_cell_5.mat");
all_popboot_mean = importdata("all_popboot_means_ball_norm.mat");
all_overall_n_means = importdata("all_overall_n_means_ball_norm.mat");
z_cell_CNN = importdata("../saved_files/z_cell_CNN_5_v2.mat");

centres = importdata("centres_50.mat");
centres = centres*0.74;

% Plot figure 6c + supplement 6.1 a 
[h_boots_n] = plot_pop_mean_results_06122022(all_popboot_mean,all_overall_n_means,centres);

% plot figure 6 supplement 1b
initialise_params;
av_num = 10;
run_many_assess_mean_time_change_18042023(virmen_cell,z_cell_CNN,tbt_cell,model_params,av_num);

% plot figure 6 supplement 1c
plot_across_sessions_mean_change_18042023(all_overall_n_means);
