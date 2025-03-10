function [h_boots] = plot_cross_decoding_results_compare(all_res_mat_cross,all_res_cell,res_ind)
% 30/09/2024

% Funciton for plotting decoding of corrective angular velocity using view
% angle decoder.

% simply plot scatter with mean line of Correlation. Compare to results for
% decoding corrective av with av weights.

% all_res_mat is num_mice x num_days x 3 (rmse,R2,rho) x 2 (check decoding
% self, cross decoding)

% res_ind to determine what results to plot. 2 for R2, 3 for correlation

num_mice = size(all_res_mat_cross,1);
num_days = size(all_res_mat_cross,2);

yaw_ind = 4;

% Store mean RMSE, and R^2
all_means = nan.*ones(num_mice,num_days,3,8);
for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(all_res_cell{m,d})
            
            cur_res = all_res_cell{m,d};
            
            all_means(m,d,:,:) = mean(cur_res(:,:,res_ind,:),1,'omitnan');
           
        end
    end
end

% Scatter plot with mean line
num_sess = num_mice*num_days;
plot_off = linspace(-0.4,0.4,num_sess);

plot_res_cross = squeeze(all_res_mat_cross(:,:,res_ind,2));
% 3 contains results for all corrective samples
plot_res = squeeze(all_means(:,:,3,yaw_ind));

ylabels = ["rmse";"R^2";"\rho"];

figure
scatter(plot_off+1.*ones(1,num_sess),squeeze(plot_res_cross(:))','filled','k')
hold on
plot([1+plot_off(1),1+plot_off(end)],[mean(squeeze(plot_res_cross(:)),'omitnan'),mean(squeeze(plot_res_cross(:)),'omitnan')],'k','LineWidth',2)

scatter(plot_off+2.*ones(1,num_sess),squeeze(plot_res(:))','filled','k')
plot([2+plot_off(1),2+plot_off(end)],[mean(squeeze(plot_res(:)),'omitnan'),mean(squeeze(plot_res(:)),'omitnan')],'k','LineWidth',2)

yline(0,'--','LineWidth',2);

title(["Decoding of"; "corrective angular velocities"])
ylabel(ylabels(res_ind))
xticks([1,2])
xticklabels(["Heading direction";"Angular velocity"])
xlabel("Decoder")
ylim([-1,1])
xlim([0.5,2.5])
axis('square')

%% Hierarchical bootstrapping
% What kind of statistical test do we run for this.

[all_p_boot,all_centres,all_sems] = run_H_boot_ets(plot_res_cross, plot_res,true);

h_boots.all_p_boot = all_p_boot;
h_boots.all_centres = all_centres;
h_boots.all_sems = all_sems;