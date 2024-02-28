function [all_perf,sig_combined] = plot_early_performance(tbt_cell,n_trials)
% 25/10/2023

% function for plotting performance using only the first n_trials BMI
% trials, and compare to performance on all other BMI trials

tbt_cell{4,2} = []; % excluded session.

types_vec = [1,4,7,10];
num_mice = size(tbt_cell,1);
num_days = size(tbt_cell,2);
all_perf = NaN(num_mice,2);
smoothed_perf = cell(num_mice,1);
av_num = n_trials;

for m = 1:num_mice
    
    tbt_details = tbt_cell{m,1};
    
    bmi_trials = tbt_details(:,ismember(tbt_details(1,:),[3,4]));
    % Get performance on first n trials
    all_perf(m,1) = sum(ismember(bmi_trials(3,1:n_trials),types_vec([3,4])))/n_trials;
    all_perf(m,2) = sum(ismember(bmi_trials(3,n_trials+1:end),types_vec([3,4])))/(size(bmi_trials,2)-n_trials);
    
    smoothed_perf{m} = movmean(ismember(bmi_trials(3,1:n_trials),types_vec([3,4])),av_num);
    
end

% Stats
% Check if early trial performance is above chance when combined across
% mice.

total_correct = sum(all_perf(:,1).*n_trials);
sig_combined = myBinomTest(total_correct,n_trials*num_mice,0.5,'one');

%% All remaining trials performance
rest_perf = NaN(num_mice,1);
for m = 1:num_mice
    tot_correct = 0;
    tot_b = 0;
    for d = 1:num_days
        if ~isempty(tbt_cell{m,d})
            tbt_details = tbt_cell{m,d};
            bmi_trials = tbt_details(:,ismember(tbt_details(1,:),[3,4]));
            if d == 1
                tot_correct = tot_correct + sum(ismember(bmi_trials(3,n_trials+1:end),types_vec([3,4])));
                tot_b = tot_b + size(bmi_trials,2)-n_trials;
            else
                tot_correct = tot_correct + sum(ismember(bmi_trials(3,:),types_vec([3,4])));
                tot_b = tot_b + size(bmi_trials,2);
            end
        end
    end
    rest_perf(m) = tot_correct/tot_b;
end
            
%% Plot simple connected lines - first 10 to all other trials
all_perf(:,2) = rest_perf;

figure
bar(mean(all_perf),'FaceColor',[.7 .7 .7],'LineWidth',2)
xticklabels(["First 10 trials";"All other trials"])
yline(0.5,'--','LineWidth',2);
ylim([0,1.1])
yticks([0,0.25,0.5,0.75,1])
hold on

for m = 1:num_mice
    plot([1,2],all_perf(m,:),'-o','Color',[0,0,0],'MarkerFaceColor',[0,0,0],'MarkerSize',10,'LineWidth',2)
end
box off
xlim([0.5,2.5])
axis('square')
ylabel("Fraction correct")
title(["Performance on first";"10 BMI trials"])