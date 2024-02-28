function [] = plot_across_sessions_mean_change_18042023(all_overall_n_means)
% 18/04/2023

% Plot mean changes by day to see if there is a trend over time.

num_mice = size(all_overall_n_means,1);
num_days = size(all_overall_n_means,2);

%% get population means
full_means = nan.*ones(2,num_mice,num_days);
for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(all_overall_n_means{m,d})
            cur_means = all_overall_n_means{m,d};
            cur_means = mean(cur_means,2,'omitnan');
            full_means(:,m,d) = cur_means;         
            
        end
    end
end

%% Plot overall means over time
figure

for m = 1:num_mice
    cur_res = squeeze(full_means(2,m,:))-squeeze(full_means(1,m,:));
    plot(find(~isnan(cur_res)),cur_res(~isnan(cur_res)),'-o','LineWidth',1.5,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    hold on
end
plot(mean(squeeze(full_means(2,:,:))-squeeze(full_means(1,:,:)),1,'omitnan'),'-o','LineWidth',3,'Color','k','MarkerFaceColor','k')
xlabel("Day")
ylabel(["Mean activity difference"; "BMI - Ball (a.u.)"])
title(["Mean activity difference";"across sessions"])
box off
axis('square')
xticks([1,2,3,4])
xlim([0.5,4.5])