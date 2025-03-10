function [] = plot_va_corrections_time_v2_subplot(virmen_cell, tbt_cell, ex_md,thresh,t_ind,h)
% 18/05/2024

% plot specific examples of trials where mice made view angle corrections
% during the trial.
pos_scale = 0.74;
x_vec = [7];
%% 
% 2 x 2
types_vec = [1,4,7,10];

virmen_data = virmen_cell{ex_md(1),ex_md(2)};
tbt_details = tbt_cell{ex_md(1),ex_md(2)};

% remove ITI
ITI = virmen_data(8,:);
cleaned_valid = clean_valid_data(ITI);
virmen_data = virmen_data(:,cleaned_valid);

va = virmen_data(7,:);
trial_num = virmen_data(12,:);

%% Plot
color_vec = [[0.4660 0.6740 0.1880]];
num_trials = 1;
plot_trials = nan.*ones(2,num_trials);
for i = 1:2
    cur_trials = find(tbt_details(3,:)==types_vec(i+2));
    cur_min_max_va = zeros(length(cur_trials),2);
    for j = 1:length(cur_trials)
        cur_min_max_va(j,1) = min(va(trial_num==cur_trials(j)));
        cur_min_max_va(j,2) = max(va(trial_num==cur_trials(j)));
    end
    if i==1
        t_nums = find(cur_min_max_va(:,1)<(-1*thresh),t_ind);
    else
        t_nums = find(cur_min_max_va(:,2)>thresh,t_ind);
    end
    plot_trials(i,:) = cur_trials(t_nums(t_ind));
end

trial_num = virmen_data(12,:);
% Plot bmi trials
% figure
for i = 1:2  
    subplot(h(i))   
    for j = 1:num_trials
        cur_trial = virmen_data(7,trial_num==plot_trials(i,j));
        t = (1:length(cur_trial))./30;
        plot(t,cur_trial,'Color',color_vec(1,:),'LineWidth',2)
        hold on
    end

    yline(0,'--','LineWidth',2);
    box off
end
ax1 = subplot(1,2,1);
ylabel(["Heading deviation"; "(rad)"]);
title("Left trial")
xlabel("Time (s)")
% axis('square')

ax2 = subplot(1,2,2);
title("Right trial")
xlabel("Time (s)")
% axis('square')

linkaxes([ax1,ax2],'y')