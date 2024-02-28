%% FIGURE 3
addpath("../shared_functions")
% prepare workspace for plotting
virmen_cell = importdata("../saved_files/virmen_cell_5.mat");
tbt_cell = importdata("../saved_files/tbt_cell_5.mat");

virmen_data = virmen_cell{1,1};
tbt_details = tbt_cell{1,1};
num_trials = 3;

nbins = 50;
linearise_x = true;
boot_samps = 100;
sub_sample = false;
ex_offset = 1.4804;

virmen_data_2 = virmen_cell{4,1};
tbt_details_2 = tbt_cell{4,1};

virmen_cell{4,2} = [];

% Plot figure 3 a
plot_example_position_trajectories(virmen_data,tbt_details,num_trials)

num_trials = 3; 
% Plot figure 3 b
plot_pva_example_trajectories(virmen_data_2,tbt_details_2,30,linearise_x,ex_offset,num_trials)
% Plot figure 3 c
plot_pva_trajectory_comparisons_0612022(virmen_data_2,tbt_details_2,nbins,linearise_x,boot_samps,sub_sample,ex_offset);

% plot figure 3 d
[h_boots] = trial_length_comparisons_plot_02122022(virmen_cell,tbt_cell);

%%
ex_md = [1,1];
offsets = [1.4953,1.4961;1.4953,1.4961;1.4953,1.4961;1.4804,1.4844;1.4804,1.4844];
linearise_x = true;
virmen_data = virmen_cell{ex_md};
tbt_details = tbt_cell{ex_md};

% Plot figure 3 e
plot_example_trajectories_06122022(virmen_cell, tbt_cell, ex_md, offsets(1,:), linearise_x, nbins);

centres = importdata("centres_50.mat"); % Change if different number of bins used
accuracies_cell_bvn = importdata("accuracies_cell_bvn_lrsep_50_100subs.mat");
accuracies_cell_shuff = importdata("accuracies_cell_bvn_lr_sep_50_shuff_true.mat");

% Plot figure 3 f
[h_boots_class] = plot_many_svm_bvn_classifier_results_25072023(accuracies_cell_bvn,accuracies_cell_shuff,centres);