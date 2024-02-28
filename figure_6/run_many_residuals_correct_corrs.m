function [h_boots] = run_many_residuals_correct_corrs(x_binned_cell,error_cell,tbt_cell,centres,plot_ind,error_thresh)
% 13/09/2023
num_mice = size(x_binned_cell,1);
num_days = size(x_binned_cell,2);
% virmen data indices in x_binned [6,7,15,17,on_ind,res_ind,vav_ind,vav_res_ind];
% plot_ind = 8 to index angular velocity from projected residuals
types_vec = [1,4,7,10];
nbins = 50;

centres = centres*0.74;

num_shuffles = 100;
ro_all = nan.*ones(num_mice,num_days);
ro_bins = nan.*ones(nbins,num_mice,num_days);
ro_bins_shuff = nan.*ones(num_shuffles,nbins,num_mice,num_days);

for m = 1:num_mice
    
    for d = 1:num_days
        if ~isempty(x_binned_cell{m,d})

            tbt_details = tbt_cell{m,d};
            
            % Only use BMI trials
            cur_trials = find(ismember(tbt_details(3,:),types_vec([3,4])));     
            
            % get correlations with corrective movements
            error_mat = error_cell{m,d};
            x_binned = x_binned_cell{m,d};
            
            cur_errors = error_mat(cur_trials,:);
            cur_yaw = squeeze(x_binned(cur_trials,:,plot_ind)); % angular velocity from residual decoded view angle.
            
            ro = corr(cur_errors(abs(cur_errors)>error_thresh),cur_yaw(abs(cur_errors)>error_thresh)); 
            ro_all(m,d) = ro;
            % Correlations in each bin. Only when heading deviations are
            % above threshold
            for b = 1:nbins
                ro_bins(b,m,d) = corr(cur_errors(abs(cur_errors(:,b))>error_thresh,b),cur_yaw(abs(cur_errors(:,b))>error_thresh,b));
            end
            
            % calculate shuffled corelations
            for n = 1:num_shuffles
                cur_yaw_shuff = cur_yaw(:);
                cur_yaw_shuff = cur_yaw_shuff(randperm(length(cur_yaw_shuff)));
                cur_yaw_shuff = reshape(cur_yaw_shuff,size(cur_yaw,1),size(cur_yaw,2));
                for b = 1:nbins
                    ro_bins_shuff(n,b,m,d) = corr(cur_errors(abs(cur_errors(:,b))>error_thresh,b),cur_yaw_shuff(abs(cur_errors(:,b))>error_thresh,b));
                end
            end
             
        end
    end
end


%% H-bootstrapping of correlations

% reorder ro_bins so nans are at the end
orig_ro_bins = ro_bins;
for m = 1:num_mice
    d_ind = 0;
    for d = 1:num_days
        if ~isnan(orig_ro_bins(1,m,d)) 
            d_ind = d_ind+1;
            ro_bins(:,m,d_ind) = orig_ro_bins(:,m,d);

        end
    end
    if d_ind<num_days
        ro_bins(:,m,d_ind+1:num_days) = nan;
    end
end

orig_ro_bins_shuff = ro_bins_shuff;
for m = 1:num_mice
    d_ind = 0;
    for d = 1:num_days
        if ~isnan(orig_ro_bins_shuff(1,1,m,d)) 
            d_ind = d_ind+1;
            ro_bins_shuff(:,:,m,d_ind) = orig_ro_bins_shuff(:,:,m,d);

        end
    end
    if d_ind<num_days
        ro_bins_shuff(:,:,m,d_ind+1:num_days) = nan;
    end
end

% Average across shuffles
ro_bins_shuff = squeeze(mean(ro_bins_shuff));

% Run H-bootstrapping
all_centres = NaN(nbins,2);
all_sems = NaN(nbins,2);
all_p_boot = NaN(nbins,1);
for b = 1:nbins
    [all_p_boot(b),all_centres(b,:),all_sems(b,:)] = run_H_boot_ets(squeeze(ro_bins(b,:,:)), squeeze(ro_bins_shuff(b,:,:)),true);
end

h_boots.all_p_boot = all_p_boot;
h_boots.all_centres = all_centres;
h_boots.all_sems = all_sems;
    
% boostrap stds
boot_stds = all_sems;
lims_all = zeros(2,nbins,2);
lims_all(1,:,:) = all_centres - boot_stds;
lims_all(2,:,:) = all_centres + boot_stds;

figure
h = fill([centres,fliplr(centres)],[squeeze(lims_all(1,:,1)),fliplr(squeeze(lims_all(2,:,1)))],'k','EdgeColor','none');
set(h,'facealpha',.3)
hold on

plot(centres,all_centres(:,1),'LineWidth',2,'Color','k')

h = fill([centres,fliplr(centres)],[squeeze(lims_all(1,:,2)),fliplr(squeeze(lims_all(2,:,2)))],'k','EdgeColor','none');
set(h,'facealpha',.3)
hold on

plot(centres,all_centres(:,2),'--','LineWidth',2,'Color','k')

ylim([-1,1])
title(["Correlation between heading deviation"; "and angular velocity from decoder output"])
xlabel("Linearised position (cm)")
ylabel(["Pearson correlation"])
axis('square')
box off
