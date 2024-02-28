function [h_boots] = plot_many_svm_bvn_neural_classifier_results(accuracies_cell_bvn,accuracies_cell_shuff,centres)
% 03/08/2023

% Function for plotting results of svm classification of trials as ball vs 
% bmi using neural data

num_mice = size(accuracies_cell_bvn,1);
num_days = size(accuracies_cell_bvn,2);

means_all = nan.*ones(num_mice,num_days,2,4,length(centres));
means_all_shuff = nan.*ones(num_mice,num_days,2,4,length(centres));

centres = centres*0.74;

nbins = length(centres);

% prepare for summary figure
for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(accuracies_cell_bvn{m,d})
            cur_acc = accuracies_cell_bvn{m,d};
            cur_acc_shuff = accuracies_cell_shuff{m,d};
            
            for i = 1:2
                means_all(m,d,i,1,:) = cur_acc(i,:);
                means_all_shuff(m,d,i,1,:) = mean(cur_acc_shuff(:,i,:));
            end

        end
    end
end

num_sess = 19;
%% Plot
figure
cur_mean = permute(squeeze(mean(squeeze(means_all(:,:,:,1,:)),3,'omitnan')),[3,1,2]);
cur_mean_shuff = permute(squeeze(mean(squeeze(means_all_shuff(:,:,:,1,:)),3,'omitnan')),[3,1,2]);

% H bootstrap means and stds
% Move NaNs to end of rows
orig_cur_means = cur_mean;
for m = 1:num_mice
    d_ind = 0;
    for d = 1:num_days
        if ~isnan(orig_cur_means(1,m,d)) 
            d_ind = d_ind+1;
            cur_mean(:,m,d_ind) = orig_cur_means(:,m,d);

        end
    end
    if d_ind<num_days
        cur_mean(:,m,d_ind+1:num_days) = nan;
    end
end

orig_cur_mean_shuff = cur_mean_shuff;
for m = 1:num_mice
    d_ind = 0;
    for d = 1:num_days
        if ~isnan(orig_cur_mean_shuff(1,m,d)) 
            d_ind = d_ind+1;
            cur_mean_shuff(:,m,d_ind) = orig_cur_mean_shuff(:,m,d);

        end
    end
    if d_ind<num_days
        cur_mean_shuff(:,m,d_ind+1:num_days) = nan;
    end
end

% Run Hierarchical bootstrapping
all_centres = NaN(nbins,2);
all_sems = NaN(nbins,2);
all_p_boot = NaN(nbins,1);
for b = 1:nbins
    [all_p_boot(b),all_centres(b,:),all_sems(b,:)] = run_H_boot_ets(squeeze(cur_mean(b,:,:)), squeeze(cur_mean_shuff(b,:,:)),false);
end

h_boots.all_p_boot = all_p_boot;
h_boots.all_centres = all_centres;
h_boots.all_sems = all_sems;

% boostrap stds
boot_stds = all_sems;
lims_all = zeros(2,nbins,2);
lims_all(1,:,:) = all_centres - boot_stds;
lims_all(2,:,:) = all_centres + boot_stds;

h = fill([centres,fliplr(centres)],[squeeze(lims_all(1,:,1)),fliplr(squeeze(lims_all(2,:,1)))],'k','EdgeColor','none');
set(h,'facealpha',.3)
hold on

plot(centres,all_centres(:,1),'LineWidth',2,'Color','k')

h = fill([centres,fliplr(centres)],[squeeze(lims_all(1,:,2)),fliplr(squeeze(lims_all(2,:,2)))],'k','EdgeColor','none');
set(h,'facealpha',.3)
hold on

plot(centres,all_centres(:,2),'--','LineWidth',2,'Color','k')

xline(200*0.74,'--','Linewidth',2); % Cue end
xline(300*0.74,'--','Linewidth',2); % Turn start
title(["Ball vs BMI Classifier";"Using Neural Data"])
ylim([0,1])
xlabel("Maze position (cm)")

ylabel(["Classification accuracy"])
axis('square')
box off
