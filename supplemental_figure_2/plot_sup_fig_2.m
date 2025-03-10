% Plot supplement figure 2
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

n_trials = 10;
max_trials = 100;
% Plot supplement figure 2a
[smoothed_perf] = plot_running_performance(tbt_cell,n_trials,max_trials);

% Plot figure 2c and supplement figure 2b
[summary_results] = plot_BMI_performance_all_mice_w_train_06122022(trial_types_summary_cell,days_cell);

% Use 10 trials to determine performance
smooth_n = 10;
% Plot supplement figure 2c
[all_rho,all_p] = trial_performance_progression_plot(tbt_cell,smooth_n);

% adjust smoothing as desired
smooth_n = 5;
% Plot supplement figure 2d
[all_rho,all_p] = trial_length_progression_plot(virmen_cell,tbt_cell,smooth_n);