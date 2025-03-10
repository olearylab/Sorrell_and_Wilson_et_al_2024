%% Supplemental figure 5
addpath("../shared_functions")

% Load data and set parameters
virmen_cell = importdata("../saved_files/virmen_cell_5.mat");
tbt_cell = importdata("../saved_files/tbt_cell_5.mat");
nbins = 50;
linearise_x = true;
num_trials = 2;
offsets = [1.4953,1.4961;1.4953,1.4961;1.4953,1.4961;1.4804,1.4844;1.4804,1.4844];

% Plot figure 5 supplement a and b
plot_single_trials_ball(virmen_cell{1,1},tbt_cell{1,1},nbins,linearise_x,offsets(1,:),num_trials)

% Calculate correlations
[behav_corr_cell,behav_std_corr_cell,behav_full_corr_cell] = run_many_behaviour_correlations(virmen_cell,tbt_cell,offsets);
% Plot figure 5 supplement c
[h_boots] = plot_behav_corrs(behav_corr_cell,behav_std_corr_cell);

% Pitch and Yaw offsets for each mouse
offsets = [1.4953,1.4961;1.4953,1.4961;1.4953,1.4961;1.4804,1.4844;1.4804,1.4844];

ex_days = [2,2,2,3,3];
ex_trials = [1,2,2,2,2];

% Plot
check_spinning_all_w_velocities_v2(virmen_cell,tbt_cell,offsets,ex_days,ex_trials);