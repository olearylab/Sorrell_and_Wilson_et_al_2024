%% Plot modelling figure - supplement figure 3
% add necessary functions to path
addpath("../shared_functions")

% try
%     addpath('G:\HMS Dropbox\Daniel Wilson\2021 PPC BMI Experiments Shared\BMI_Methods_Paper_Draft\Figure_making\Figure_making_code\shared_functions');
% catch
% end;
initialise_sim_model;
tbt_details = importdata("DW81_20200928_trial_by_trial_details.mat");
virmen_data = importdata("trimmed_virmen_DW81_20200928_train_revised.mat");
alpha = 0; % Purely sensory
[rel_results, va_rel] = basic_va_control_sim_23062023(virmen_data, tbt_details, sim_model,alpha);

alpha = 1; % Purely motor
[rel_results, va_rel] = basic_va_control_sim_23062023(virmen_data, tbt_details, sim_model,alpha);
