%% Script for training and testing shuffled models

addpath("../shared_functions")

virmen_cell = importdata("../saved_files/virmen_cell_5.mat");
tbt_cell = importdata("../saved_files/tbt_cell_5.mat");
z_cell = importdata("../saved_files/z_cell_CNN_5_v2.mat");

centres = importdata("../saved_files/centres_50.mat"); % Change if different number of bins used

n_shuff = 100;
zscore_x = true;
nbins = 50;
shuff_limit = 300;

% load lambda_cell into workspace

extra_norm = false;
[all_r2_cell_shuff,trial_avs_cell_shuff,trial_av_r2_cell_shuff,binned_r2_cell_shuff,trial_av_cd_r2_cell_shuff,trial_av_rmse_cell_shuff] = run_many_linear_model_crossvals_shuffle(z_cell,virmen_cell,tbt_cell,n_shuff,extra_norm,lambda_cell,zscore_x,nbins,shuff_limit);

% load lambda_cell_norm into workspace

extra_norm = true;
[all_r2_cell_shuff_norm,trial_avs_cell_shuff_norm,trial_av_r2_cell_shuff_norm,binned_r2_cell_shuff_norm,trial_av_cd_r2_cell_shuff_norm,trial_av_rmse_cell_shuff_norm] = run_many_linear_model_crossvals_shuffle(z_cell,virmen_cell,tbt_cell,n_shuff,extra_norm,lambda_cell,zscore_x,nbins,shuff_limit)