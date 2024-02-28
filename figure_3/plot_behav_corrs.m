function [h_boots] = plot_behav_corrs(behav_corr_cell,behav_std_corr_cell)
% 25/07/2023

num_mice = size(behav_corr_cell,1);
num_days = size(behav_corr_cell,2);

corrs_mat = nan.*ones(4,4,2,num_mice,num_days);
stds_mat = nan.*ones(4,4,2,num_mice,num_days);

for m = 1:num_mice
    for d = 1:num_days
        if~isempty(behav_corr_cell{m,d})
            for i = 1:2
                corrs_mat(:,:,i,m,d) = behav_corr_cell{m,d,i};
                stds_mat(:,:,i,m,d) = behav_std_corr_cell{m,d,i};
            end
        end
    end
end


%% Plot just across types comparison

figure
for i = 1:2
    within_means = mean([squeeze(corrs_mat(1,1,i,:)),squeeze(corrs_mat(2,2,i,:)),squeeze(corrs_mat(3,3,i,:)),squeeze(corrs_mat(4,4,i,:))],2,'omitnan');
    across_means = mean([squeeze(corrs_mat(1,3,i,:)),squeeze(corrs_mat(2,4,i,:))],2,'omitnan');
    diff_means = mean(across_means-within_means,'omitnan');
    diff_stds = std(across_means-within_means,0,'omitnan');
    errorbar(i,diff_means,diff_stds,'o',"Color",'k',"MarkerFaceColor",'k','LineWidth',2,'MarkerSize',20)
    hold on
end

title(["Correlations in behaviour";"decrease across trial type"])
ylabel(["\Delta Pearson Correlation"])
xticks([1,2])
xlim([0.5,2.5])
xticklabels(["Forward Velocity";"Angular Velocity"])
xtickangle(30)
axis('square')
box off

%% Hierarchical bootstrap
final_means = nan.*ones(num_mice,num_days,2);
for i = 1:2
    within_means = mean([squeeze(corrs_mat(1,1,i,:,:)),squeeze(corrs_mat(2,2,i,:,:)),squeeze(corrs_mat(3,3,i,:,:)),squeeze(corrs_mat(4,4,i,:,:))],2,'omitnan');
    across_means = mean([squeeze(corrs_mat(1,3,i,:,:)),squeeze(corrs_mat(2,4,i,:,:))],2,'omitnan');
    
    cur_means = across_means-within_means;
    for m = 1:num_mice
        final_means(m,1:sum(~isnan(cur_means(m,:))),i) = cur_means(m,~isnan(cur_means(m,:)));
    end
end

% run hierarchical bootstrapping
[all_p_boot,all_centres,all_sems] = run_H_boot_ets(squeeze(final_means(:,:,1)), squeeze(final_means(:,:,2)),false);

h_boots.all_p_boot = all_p_boot;
h_boots.all_centres = all_centres;
h_boots.all_sems = all_sems;
