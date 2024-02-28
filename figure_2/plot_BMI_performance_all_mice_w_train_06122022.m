function [summary_results] = plot_BMI_performance_all_mice_w_train_06122022(trial_types_summary_cell,days_cell)
% 06/12/2022
% create several performance related plots.

% include training session results for each mouse. Start Days vec at 0.
% First value is for training data.

%%
colour_vec = [[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0]];
ball_colour = [0.4940 0.1840 0.5560];
bmi_colour = [0.4660 0.6740 0.1880];

set(0,'DefaultAxesFontSize',20)
lines_vec = ["-o";"--^";"-o";"--^"];

% get results vector and plot
figure

for m = 1:length(trial_types_summary_cell)
    cur_res = trial_types_summary_cell{m};
    cur_m_results = zeros(4,size(cur_res,2));
    days_vec = days_cell{m};
    
    for d = 1:length(days_vec)
        for i = 1:4
            cur_m_results(i,d) = cur_res(3*(i-1)+1,d)/sum(cur_res(3*(i-1)+1:3*i,d));
        end
    end
    cur_m_results(3:4,1) = [nan;nan];
    subplot(1,length(trial_types_summary_cell),m)
    for i = 1:4
        if i < 3
            plot(days_vec(1),cur_m_results(i,1)',lines_vec(i),'LineWidth',2,'Color',ball_colour,'MarkerFaceColor',ball_colour,'MarkerSize',10)
            hold on
            plot(days_vec(2:end),cur_m_results(i,2:end)',lines_vec(i),'LineWidth',2,'Color',ball_colour,'MarkerFaceColor',ball_colour,'MarkerSize',10)
        else
            plot(days_vec(1),cur_m_results(i,1)',lines_vec(i),'LineWidth',2,'Color',bmi_colour,'MarkerFaceColor',bmi_colour,'MarkerSize',10)
            hold on
            plot(days_vec(2:end),cur_m_results(i,2:end)',lines_vec(i),'LineWidth',2,'Color',bmi_colour,'MarkerFaceColor',bmi_colour,'MarkerSize',10)  
        end
    end
    yline(0.5,'--','LineWidth',2);
    title("Mouse " + m)
    yticks([0, 0.25, 0.5, 0.75, 1])
    if m == 1
        ylabel("Fraction Correct")
        yticklabels([0, 0.25, 0.5, 0.75, 1])
    else
        yticklabels([])
    end
    % Some hardcoded stuff for labels etc.
    xlabel("Day")
    xlim([-0.5,4+0.5])
    ylim([0,1.1])
    xticks(0:4)
    xticklabels({'Train','1','2','3','4'})
    box off
    axis('square')
    if m==1
        legend('Ball Left', 'Ball Right', 'BMI Left', 'BMI Right')
    end
end

%% Now include summary results and statistical analysis
% Possibly specific days and or mice excluded...

% Bar plot showing mean across all animals, training, testing
% ball and testing BMI, as well as individual animals plotted, and linked

% figure
summary_results = zeros(length(trial_types_summary_cell),3);
tot_b = zeros(length(trial_types_summary_cell),2);
tot_n = zeros(length(trial_types_summary_cell),2);
for m = 1:length(trial_types_summary_cell)
    cur_res = trial_types_summary_cell{m};
    
    summary_results(m,1) = sum(cur_res([1,4],1))/sum(cur_res(1:6,1));
    summary_results(m,2) = sum(cur_res([1,4],2:end),'all')/sum(cur_res(1:6,2:end),'all');
    summary_results(m,3) = sum(cur_res([7,10],2:end),'all')/sum(cur_res(7:end,2:end),'all');
    
    tot_b(m,1) = sum(cur_res([7,10],2:end),'all');
    tot_b(m,2) = sum(cur_res(7:end,2:end),'all');
    tot_n(m,1) = sum(cur_res([1,4],2:end),'all');
    tot_n(m,2) = sum(cur_res(1:6,2:end),'all');
    
end

figure
% with connected points
bar(mean(summary_results),'FaceColor',[.7 .7 .7],'LineWidth',2)
xticklabels(["Train";"Ball";"BMI"])
yline(0.5,'--','LineWidth',2);
ylim([0,1.1])
yticks([0,0.25,0.5,0.75,1])
hold on

for m = 1:length(trial_types_summary_cell)
    plot(1:3,summary_results(m,:),'-o','Color',colour_vec(m,:),'MarkerFaceColor',colour_vec(m,:),'MarkerSize',10,'LineWidth',2)
end
box off
% axis('square')
ylabel("Fraction correct")
title("Summary of behavioral performance")