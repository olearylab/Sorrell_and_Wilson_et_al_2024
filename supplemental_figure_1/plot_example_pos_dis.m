function [] = plot_example_pos_dis(results_struct,xfull,plot_trial)
% 10/08/2023

% plot example position decoding highlighting discontinuities
d_col = [0.4940 0.1840 0.5560];

% get valid data indicator
ITI = xfull(8,:);
cleaned_valid = clean_valid_data(ITI);
xfull = xfull(:,cleaned_valid);
xtest_new = results_struct.xtest(cleaned_valid,:);
xprediction_new = results_struct.xprediction(cleaned_valid,:);
trial_num = xfull(12,:);


t = (1:sum(trial_num==plot_trial))/30;

plot(t,0.0074*xtest_new(trial_num==plot_trial,1),'Color',[0.5,0.5,0.5],'LineWidth',2)
hold on 
plot(t,0.0074*xprediction_new(trial_num==plot_trial,1),'Color',d_col,'LineWidth',2)
title("Decoded Y position")
xlabel("Time (s)")
ylabel(["Y position";"(rad/s)"])