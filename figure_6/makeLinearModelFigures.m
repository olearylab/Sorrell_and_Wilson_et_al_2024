load('G:\Dropbox (HMS)\2021 PPC BMI Experiments Shared\BMI_Methods_Paper_Draft\Figure_making\Figure_making_code\linear_models_dan\positionPredictors_10bumps.mat');
%% Script for plotting results
addpath(genpath('../shared_functions'));
centres = importdata("../saved_files/centres_50.mat"); % Change if different number of bins used
ex_md = [1,1];
n_ind = 8;
plot_linear_model_crossvals_examples(trial_av_rmse_cell,trial_av_rmse_cell_norm,trial_avs_cell,trial_avs_cell_norm,ex_md,n_ind,centres);
linkaxes;
clear;

load('G:\Dropbox (HMS)\2021 PPC BMI Experiments Shared\BMI_Methods_Paper_Draft\Figure_making\Figure_making_code\linear_models_dan\positionPredictors_10bumps.mat','trial_av_rmse_cell','trial_av_rmse_cell_norm');
[h_boots_position] = plot_linear_model_scatters(trial_av_rmse_cell,trial_av_rmse_cell_norm);
sgtitle('Position predictors');
set(gcf,'Position',[680,796,1007,302]);
clearvars -except h_boots_position;
load('G:\Dropbox (HMS)\2021 PPC BMI Experiments Shared\BMI_Methods_Paper_Draft\Figure_making\Figure_making_code\linear_models_dan\ballPredictors.mat','trial_av_rmse_cell','trial_av_rmse_cell_norm');
[h_boots_ball] = plot_linear_model_scatters(trial_av_rmse_cell,trial_av_rmse_cell_norm);
sgtitle('Ball predictors');
set(gcf,'Position',[680,796,1007,302]);