%% Plot cross variable decoding results


%% Plotting results of decoding corrective angular velocity using view angle decoder
% Load in results
all_res_mat = importdata("all_res_mat_va_decoding_corrective_vav.mat");

% plot
plot_cross_decoding_results(all_res_mat);