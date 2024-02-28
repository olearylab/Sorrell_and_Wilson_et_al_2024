%% Script for training and testing models

addpath("../shared_functions")

virmen_cell = importdata("../saved_files/virmen_cell_5.mat");
tbt_cell = importdata("../saved_files/tbt_cell_5.mat");
z_cell = importdata("../saved_files/z_cell_CNN_5_v2.mat");

centres = importdata("../saved_files/centres_50.mat"); % Change if different number of bins used

num_subs = 30;
zscore_x = true;
nbins = 50;
run_lambda_calc;
% load lambda_cell into workspace

extra_norm = false;
[all_r2_cell,trial_avs_cell,trial_av_r2_cell,binned_r2_cell,trial_av_cd_r2_cell,trial_av_rmse_cell] = run_many_linear_model_crossvals_varyreg(z_cell,virmen_cell,tbt_cell,num_subs,extra_norm,lambda_cell,zscore_x,nbins);

% load lambda_cell_norm into workspace

extra_norm = true;
[all_r2_cell_norm,trial_avs_cell_norm,trial_av_r2_cell_norm,binned_r2_cell_norm,trial_av_cd_r2_cell_norm,trial_av_rmse_cell_norm] = run_many_linear_model_crossvals_varyreg(z_cell,virmen_cell,tbt_cell,num_subs,extra_norm,lambda_cell_norm,zscore_x,nbins);

save('ballPredictors.mat')