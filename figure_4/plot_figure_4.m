%% Plot Figure 4

% Set up workspace
addpath("../shared_functions")

virmen_cell = importdata("../saved_files/virmen_cell_5.mat");
tbt_cell = importdata("../saved_files/tbt_cell_5.mat");
virmen_train_cell = importdata("../saved_files/virmen_train_cell_5.mat");
tbt_train_cell = importdata("../saved_files/tbt_train_cell_5.mat");
centres = importdata("centres_50.mat"); % Change if different number of bins used

ex_md = [1,1];
offsets = [1.4953,1.4961;1.4953,1.4961;1.4953,1.4961;1.4804,1.4844;1.4804,1.4844];
linearise_x = true;
nbins = 50;

virmen_data = virmen_cell{ex_md};
tbt_details = tbt_cell{ex_md};

% Calculate mean ball view angle for example session
[mean_binned,std_binned] = calculate_mean_ball_va(virmen_train_cell{1},tbt_train_cell{1},nbins,offsets(1,:));
% Plot figure 4b 
plot_specific_behaviour_examples_norm(virmen_cell{1,2},tbt_cell{1,2},nbins,offsets(1,:),mean_binned,std_binned,centres);

error_thresh = 1;
% Plot figure 4c
[error_cell,x_cell,ro_all,sign_opposite,h_boots] = run_many_check_error_correction_norm(virmen_cell,tbt_cell,virmen_train_cell,tbt_train_cell,nbins,error_thresh);

% Plot figure 4d
[error_cell,x_cell,ro_all,sign_opposite,h_boots_dec] = run_many_check_BMI_correct_errors_norm(virmen_cell,tbt_cell,virmen_train_cell,tbt_train_cell,nbins,error_thresh,true);


