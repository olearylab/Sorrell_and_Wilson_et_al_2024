%% Plot figure 6 supplement 2

full_pairs_means_cell = importdata("full_pairs_means_cell_ball_norm_nosub.mat");

[bootstats_center,bootstats_sem,h_boots_noise] = plot_many_noise_corrs(full_pairs_means_cell);