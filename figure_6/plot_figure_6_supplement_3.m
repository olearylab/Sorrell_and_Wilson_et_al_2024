%% Plot figure 6 supplement 3

va_corrs = importdata("va_dec_corrs_29122023_cut.mat");
yaw_corrs = importdata("yaw_dec_corrs_subs_29122023_cut.mat");

%
[h_boots] = plot_weights_changes_corrs(va_corrs,yaw_corrs);