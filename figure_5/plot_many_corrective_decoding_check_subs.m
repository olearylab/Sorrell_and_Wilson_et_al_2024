function [h_boots] = plot_many_corrective_decoding_check_subs(all_res_cell)
% 26/06/2023

% Plot results of decoding corrective ball angualr velocities.

num_mice = size(all_res_cell,1);
num_days = size(all_res_cell,2);

% Store mean RMSE, and R^2
all_means = nan.*ones(3,8,num_mice,num_days);
all_means_R2 = nan.*ones(3,8,num_mice,num_days);
for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(all_res_cell{m,d})
            
            cur_res = all_res_cell{m,d};
            % RMSE - legacy
            all_means(:,:,m,d) = mean(cur_res(:,:,1,:),1,'omitnan');
            
            all_means_R2(:,:,m,d) = mean(cur_res(:,:,2,:),1,'omitnan');
           
        end
    end
end

%% Plot

% Add offset of points for plotting
num_sess = num_mice*num_days;
plot_off = linspace(-0.4,0.4,num_sess);

% Corrective angular velocity, R^2, low and high heading deviations 
figure
for i = 1:2
    scatter(plot_off+i.*ones(1,num_sess),squeeze(all_means_R2(i,4,:))','filled','k')
    hold on
    plot([i+plot_off(1),i+plot_off(end)],[mean(squeeze(all_means_R2(i,4,:)),'omitnan'),mean(squeeze(all_means_R2(i,4,:)),'omitnan')],'k','LineWidth',2)
end
title(["Ball angular velocity decoding";"accuracy during corrections"])
ylabel("R^2")
xticks([1,2,3,4])
xticklabels(["|\Delta\theta|\leq1";"|\Delta\theta|>1"])
xlabel("Normalised heading deviation magnitude")
ylim([0,1])
axis('square')


%% Hierarchical bootstrap
% Put nans at end of rows
means_ready = nan.*ones(2,num_mice,num_days);
for m = 1:num_mice
    cur_means = squeeze(all_means_R2(1:2,4,m,:));
    means_ready(:,m,1:sum(~isnan(cur_means(1,:)))) = cur_means(:,~isnan(cur_means(1,:)));
end

[all_p_boot,all_centres,all_sems] = run_H_boot_ets(squeeze(means_ready(1,:,:)), squeeze(means_ready(2,:,:)),false);

h_boots.all_p_boot = all_p_boot;
h_boots.all_centres = all_centres;
h_boots.all_sems = all_sems;
