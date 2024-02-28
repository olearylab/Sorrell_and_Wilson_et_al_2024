function [h_boots] = trial_length_comparisons_plot(virmen_cell,tbt_cell)
% 13/01/2022

% function for comparing and plotting trials lengths for different trial
% types on BMI days.

num_mice = size(virmen_cell,1);
num_days = size(virmen_cell,2);

colour_vec = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];

%% load and get trial length results
trial_means = cell(num_mice,1);
trial_stds = cell(num_mice,1);
all_trials_m = cell(num_mice,1);
types_vec = [1,4,7,10];
all_trial_lengths = cell(length(types_vec),1);

m_d_trial_lengths = cell(num_mice,num_days,2);
mean_trial_lengths = nan.*ones(num_mice,num_days,2);

m_trial_lengths = cell(num_mice,2);
num_m_ts = nan.*ones(num_mice,2);
for m = 1:num_mice
    cur_trial_means = zeros(length(types_vec),num_days);
    cur_trial_stds = zeros(length(types_vec),num_days);
    cur_trial_lengths = cell(length(types_vec),1);
    for d = 1:num_days
        if ~isempty(virmen_cell{m,d})
            xfull = virmen_cell{m,d}; 
            tbt_details = tbt_cell{m,d};

            ITI = xfull(8,:);
            test_valid = clean_valid_data(ITI);
            trial_num = xfull(12,:);
            trial_lengths = zeros(max(trial_num),1);
            for j = 1:max(trial_num)

                trial_lengths(j) = sum(test_valid(trial_num==j));

            end
            for j = 1:length(types_vec)
                cur_trial_means(j,d) = mean(trial_lengths(tbt_details(3,:)==types_vec(j)));
                cur_trial_stds(j,d) = std(trial_lengths(tbt_details(3,:)==types_vec(j)));
                cur_trial_lengths{j} = [cur_trial_lengths{j};trial_lengths(tbt_details(3,:)==types_vec(j))];
                all_trial_lengths{j} = [all_trial_lengths{j};trial_lengths(tbt_details(3,:)==types_vec(j))];
            end
            
            for j = 1:2
                m_d_trial_lengths{m,d,j} = trial_lengths(ismember(tbt_details(3,:),types_vec([(j-1)*2+1,(j-1)*2+2])));
                mean_trial_lengths(m,d,j) = mean(trial_lengths(ismember(tbt_details(3,:),types_vec([(j-1)*2+1,(j-1)*2+2]))));
                m_trial_lengths{m,j} = [m_trial_lengths{m,j};trial_lengths(ismember(tbt_details(3,:),types_vec([(j-1)*2+1,(j-1)*2+2])))];
            end
        else
            cur_trial_means(:,d) = nan*ones(length(types_vec),1);
            cur_trial_stds(:,d) = nan*ones(length(types_vec),1);
        end
        
    end
    for j = 1:2
        num_m_ts(m,j) = length(m_trial_lengths{m,j});
    end
    
    trial_means{m} = cur_trial_means;
    trial_stds{m} = cur_trial_stds;
    all_trials_m{m} = cur_trial_lengths;
end


% sample period
dt = 1/30;

%% Scatter plot version
circle_size = 100;
mean_trial_lengths = mean_trial_lengths.*dt;
figure
scatter_means = permute(mean_trial_lengths,[3,1,2]);
scatter(scatter_means(1,:),scatter_means(2,:),circle_size,'k','filled')
hold on
plot([10,22],[10,22],'--','color','k','LineWidth',2)
xlim([10,22])
xlim([10,22])
xticks([10,15,20])
yticks([10,15,20])
axis('square')
xlabel("Ball trial lengths (s)")
ylabel("BMI trial lengths (s)")
title("Comparison of Trial Lengths")
%% Hierarchical bootstrap

[p_boot,all_centres,all_sems] = run_H_boot_ets(squeeze(mean_trial_lengths(:,:,1)), squeeze(mean_trial_lengths(:,:,2)),true);

h_boots.all_p_boot = p_boot;
h_boots.all_centres = all_centres;
h_boots.all_sems = all_sems;