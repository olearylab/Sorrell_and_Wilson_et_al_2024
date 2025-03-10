function [] = plot_cross_decoding_results(all_res_mat)
% 21/09/2024

% Funciton for plotting decoding of corrective angular velocity using view
% angle decoder.

% simply plot scatter with mean line of Correlation.

% all_res_mat is num_mice x num_days x 3 (rmse,R2,rho) x 2 (check decoding
% self, cross decoding)

num_mice = size(all_res_mat,1);
num_days = size(all_res_mat,2);

% Scatter plot with mean line
num_sess = num_mice*num_days;
plot_off = linspace(-0.4,0.4,num_sess);

plot_res = squeeze(all_res_mat(:,:,3,2));

figure
scatter(plot_off+1.*ones(1,num_sess),squeeze(plot_res(:))','filled','k')
hold on
plot([1+plot_off(1),1+plot_off(end)],[mean(squeeze(plot_res(:)),'omitnan'),mean(squeeze(plot_res(:)),'omitnan')],'k','LineWidth',2)
title(["Decoding of corrective angular velocities";"using heading direction decoder"])
ylabel("\rho")
xticks(1)
% xticklabels()
% xlabel("Normalised heading deviation magnitude")
ylim([-1,1])
axis('square')

%% Hierarchical bootstrapping
% What kind of statistical test do we run for this.