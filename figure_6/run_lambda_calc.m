%% Script for creating lambda cell

addpath("../shared_functions")

virmen_cell = importdata("../saved_files/virmen_cell_5.mat");
tbt_cell = importdata("../saved_files/tbt_cell_5.mat");
z_cell = importdata("../saved_files/z_cell_CNN_5_v2.mat");

lambda_vec = [10^-3,10^-2,10^-1,0,1,10,100,1000];
zscore_x = true;
nbins = 50;

accuracy_type = "tav_rmse"; 
% Current alternatives
% accuracy_type = "tav_r2_cd"; % trial averaged coefficient of determination
% accuracy_type = "bin_r2_corrsquare"; % binned corrsquare
% accuracy_type = "tav_r2_corrsquare"; % trial averaged corrsquare
% accuracy_type = "r2_corrsquare"; % all sample corrsquare
%%
% Maintaining mean shifts
extra_norm = false;
[lambda_cell] = run_many_linear_model_crossvals_lambda(z_cell,virmen_cell,tbt_cell,extra_norm,lambda_vec,zscore_x,nbins,accuracy_type);
%%
% removing mean shifts
extra_norm = true;
[lambda_cell_norm] = run_many_linear_model_crossvals_lambda(z_cell,virmen_cell,tbt_cell,extra_norm,lambda_vec,zscore_x,nbins,accuracy_type);