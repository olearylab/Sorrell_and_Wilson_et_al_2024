function [smoothed_perf] = plot_running_performance(tbt_cell,n_trials,max_trials)
% 10/11/2023

% function for plotting running performance from first to specified max BMI
% trial

% Plot as groups of n_trials
% Also includes a moving average plot.

tbt_cell{4,2} = [];

types_vec = [1,4,7,10];
num_mice = size(tbt_cell,1);
num_days = size(tbt_cell,2);

smoothed_perf = cell(num_mice,1);
smoothed_perf_max = NaN(num_mice,max_trials);
block_perf = cell(num_mice,1);
day1_smoothed_perf = cell(num_mice,1);
av_num = n_trials;

all_bmi_trials = cell(num_mice,1);
all_len = NaN(num_mice,1);
all_len_b = NaN(num_mice,1);
day1_all_len = NaN(num_mice,1);
for m = 1:num_mice
    cur_bmi_trials = [];
    for d = 1:num_days
        if ~isempty(tbt_cell{m,d})
            tbt_details = tbt_cell{m,d};
    
            bmi_trials = tbt_details(3,ismember(tbt_details(1,:),[3,4]));
            cur_bmi_trials = [cur_bmi_trials,bmi_trials];
            if d == 1
                day1_smoothed_perf{m} = movmean(ismember(bmi_trials,types_vec([3,4])),av_num,"Endpoints","discard");
                day1_all_len(m) = length(day1_smoothed_perf{m});
            end
        end
    end
    all_bmi_trials{m} = cur_bmi_trials;
    smoothed_perf{m} = movmean(ismember(cur_bmi_trials,types_vec([3,4])),av_num,"Endpoints","discard");  
    smoothed_perf_max(m,:) = movmean(ismember(cur_bmi_trials(1:max_trials),types_vec([3,4])),av_num);   
    all_len(m) = length(smoothed_perf{m});
    all_len_b(m) = length(all_bmi_trials{m});
end
%%
% Plot moving average
%% Mean and SEM for first 100 trials - with truncations
figure
% Perform bootstrapping for mean and sem
% Could replace with hierarchical, but only one level
boot_samps = 100;
bootstrap_means = zeros(boot_samps,max_trials);
for b = 1:boot_samps
    cur_samples = datasample(1:num_mice,num_mice);
    cur_resamp = smoothed_perf_max(cur_samples,:);
    bootstrap_means(b,:) = mean(cur_resamp,'omitnan');
end

plot_means = mean(bootstrap_means);
plot_sem = std(bootstrap_means);
lims_all = [plot_means-plot_sem;plot_means+plot_sem];

h = fill([1:max_trials,fliplr(1:max_trials)],[squeeze(lims_all(1,:)),fliplr(squeeze(lims_all(2,:)))],'k','EdgeColor','none');
set(h,'facealpha',.3)
hold on
plot(plot_means,'Color','k','LineWidth',2)
% xticklabels([])
yline(0.5,'--','LineWidth',2);
ylim([0,1.1])
yticks([0,0.25,0.5,0.75,1])
box off
% axis('square')
ylabel("Fraction correct")
title(["Moving average performance"; "on first 100 BMI trials"])

%% Blocks of n_trials
all_blocks = NaN(num_mice,ceil(max(all_len_b)./n_trials));
for m = 1:num_mice
    cur_trials = all_bmi_trials{m};
    cur_num = ceil(all_len_b(m)/n_trials);
    for b = 1:cur_num
        if b<cur_num
            all_blocks(m,b) = sum(ismember(cur_trials((b-1)*n_trials+1:b*n_trials),types_vec([3,4])))/n_trials;
        else            
            all_blocks(m,b) = sum(ismember(cur_trials((b-1)*n_trials+1:end),types_vec([3,4])))/(all_len_b(m)-(b-1)*n_trials);
        end
    end
end

%% Blocks of n_trials
cur_blocks = all_blocks(:,1:10);

boot_samps = 100;
bootstrap_means = zeros(boot_samps,10);
for b = 1:boot_samps
    cur_samples = datasample(1:num_mice,num_mice);
    cur_resamp = cur_blocks(cur_samples,:);
    bootstrap_means(b,:) = mean(cur_resamp,'omitnan');
end

plot_means = mean(bootstrap_means);
plot_sem = std(bootstrap_means);
lims_all = [plot_means-plot_sem;plot_means+plot_sem];

figure
h = fill([1:10,fliplr(1:10)],[squeeze(lims_all(1,:)),fliplr(squeeze(lims_all(2,:)))],'k','EdgeColor','none');
set(h,'facealpha',.3)
hold on
plot(plot_means,'Color','k','LineWidth',2)
% xticklabels([])
yline(0.5,'--','LineWidth',2);
ylim([0,1.1])
xlim([0,11])
% xticks([1,5,10])
yticks([0,0.25,0.5,0.75,1])
box off
% axis('square')
ylabel("Fraction correct")
title(["Average performance on first"; "10 groups of 10 BMI trials"])
