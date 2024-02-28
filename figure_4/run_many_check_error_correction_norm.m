function [error_cell,x_cell,ro_all,sign_opposite,h_boots] = run_many_check_error_correction_norm(virmen_cell,tbt_cell,virmen_train_cell,tbt_train_cell,nbins,error_thresh)
% 09/09/2023

num_mice = size(virmen_cell,1);
num_days = size(virmen_cell,2);

types_vec = [1,4,7,10];
num_shuffles = 100;

offsets = [1.4953,1.4961;1.4953,1.4961;1.4953,1.4961;1.4804,1.4844;1.4804,1.4844];

error_cell = cell(num_mice,num_days);
x_cell = cell(num_mice,num_days);

for m = 1:num_mice
    % get mean view angle from training ball trials for specific mouse
    [mean_binned,std_binned] = calculate_mean_ball_va(virmen_train_cell{m},tbt_train_cell{m},nbins,offsets(m,:));
    for d = 1:num_days
        if ~isempty(virmen_cell{m,d})
 
            % calculate binned heading deviation, and other variables for
            % each session
            [error_cell{m,d},x_cell{m,d}] = check_error_correction_norm(virmen_cell{m,d},tbt_cell{m,d},mean_binned,std_binned,nbins,offsets(m,:),false);
            
        end
    end
end

%% Plotting
% Combine left and right
plot_start = 2; % 0 for ball trials, 2 for BMI trials

% Initialise matrices for storing results
ro_all = nan.*ones(num_mice,num_days);
sign_opposite = nan.*ones(num_mice,num_days);
ro_bins = nan.*ones(nbins,num_mice,num_days);
ro_bins_shuff = nan.*ones(num_shuffles,nbins,num_mice,num_days);

for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(virmen_cell{m,d})
            
            error_mat = error_cell{m,d};
            x_binned = x_cell{m,d};
            tbt_details = tbt_cell{m,d};

            cur_trials = find(ismember(tbt_details(3,:),types_vec([plot_start+1,plot_start+2])));
            cur_errors = error_mat(cur_trials,:);
            cur_yaw = squeeze(x_binned(cur_trials,:,3));
            % Calculate correlations
            ro = corr(cur_errors(abs(cur_errors)>error_thresh),cur_yaw(abs(cur_errors)>error_thresh));
            ro_all(m,d) = ro;
            % Calculate binned correlations
            for b = 1:nbins
                ro_bins(b,m,d) = corr(cur_errors(abs(cur_errors(:,b))>error_thresh,b),cur_yaw(abs(cur_errors(:,b))>error_thresh,b));
            end
            
            % Calculate shuffled correlations
            for n = 1:num_shuffles
                cur_yaw_shuff = cur_yaw(:);
                cur_yaw_shuff = cur_yaw_shuff(randperm(length(cur_yaw_shuff)));
                cur_yaw_shuff = reshape(cur_yaw_shuff,size(cur_yaw,1),size(cur_yaw,2));
                for b = 1:nbins
                    ro_bins_shuff(n,b,m,d) = corr(cur_errors(abs(cur_errors(:,b))>error_thresh,b),cur_yaw_shuff(abs(cur_errors(:,b))>error_thresh,b));
                end
            end
            
            % Check for opposite sign fraction - legacy
            error_sign = cur_errors>0;
            yaw_sign = cur_yaw>0;
            sign_opposite(m,d) = sum(error_sign(abs(cur_errors)>error_thresh)~=yaw_sign(abs(cur_errors)>error_thresh))/length(error_sign(abs(cur_errors)>error_thresh));

        end
            
    end
end

%% Get centres
[x_binned,centres] = bin_kin_data(virmen_cell{1,1},6,true,nbins);
centres = centres*0.74;

%% Correlations and H-bootstrapping

% Set in run_H_boot_ets
% boot_samps = 1000;
% num_trials = 4;

% Put NaNs at end of each row for h-bootstrapping
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

% Put NaNs at end of each row for h-bootstrapping
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

% average across shuffles
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

% Plot
figure

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

ylim([-1,1])
title(["Correlation between heading deviation"; "and ball angular velocity"])
xlabel("Linearised Position (cm)")
ylabel(["Pearson Correlation"])
axis('square')
box off
