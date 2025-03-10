%% FIGURE 2
addpath("../shared_functions")

% Prepare workspace for plotting
virmen_cell = importdata("../saved_files/virmen_cell_5.mat");
tbt_cell = importdata("../saved_files/tbt_cell_5.mat");

trial_types_summary_cell = importdata("trial_types_summary_cell_w_train_5.mat");
days_cell{1} = [0,1,2,3,4];
days_cell{2} = [0,1,2,3,4];
days_cell{3} = [0,1,2,3,4];
days_cell{4} = [0,1,3,4];
days_cell{5} = [0,1,2,3,4];

virmen_cell{4,2} = [];

% Plot figure 2c and supplement figure 2b
[summary_results] = plot_BMI_performance_all_mice_w_train_06122022(trial_types_summary_cell,days_cell);

n_trials = 10;
% Plot figure 2d
[all_perf,sig_combined] = plot_early_performance(tbt_cell,n_trials);

summary_cell = importdata("../saved_files/t_types_summary_cell_5.mat");

% Plot figure 2e-f
plot_perf_across_days(summary_cell);