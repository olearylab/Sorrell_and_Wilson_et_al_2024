%% Plot cross variable decoding results

% These need to be updated if this script is moved
addpath("../shared_functions")

%% Plotting results of decoding corrective angular velocity using view angle decoder
% Load in results
all_res_mat_cross = importdata("all_res_mat_va_decoding_corrective_vav.mat");
all_res_cell = importdata("../figure_5/all_res_cell_corrective_pixels.mat");

res_ind = 3; % for correlations
% plot
[h_boots] = plot_cross_decoding_results_compare(all_res_mat_cross,all_res_cell,res_ind);