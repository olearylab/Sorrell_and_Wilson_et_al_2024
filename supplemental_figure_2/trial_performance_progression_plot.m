function [all_rho,all_p] = trial_performance_progression_plot(tbt_cell,smooth_n)
% 14/09/2024

% function for plotting BMI performance over time.

% Plot each session separately
% Also plot all sessions concatenated
% Plots best fit line.
% Calculate correlation. 
% Smooth with movmean smooth_n, also used for correlations.
% remove endpoints in movmean

% Numbered by BMI trial, does this make sense?

num_mice = size(tbt_cell,1);
num_days = size(tbt_cell,2);

% remove mouse 4 day 2
tbt_cell{4,2} = [];

%% load and get moving average of bmi trial performance
types_vec = [1,4,7,10]; % correct trials

m_d_performance = cell(num_mice,num_days);
m_performance = cell(num_mice,1);

num_m_ts = nan.*ones(num_mice,num_days);
for m = 1:num_mice
    cur_m_bmi_trials = [];
    for d = 1:num_days
        if ~isempty(tbt_cell{m,d})
            tbt_details = tbt_cell{m,d};
            
            bmi_trials = tbt_details(3,ismember(tbt_details(1,:),[3,4]));
            cur_m_bmi_trials = [cur_m_bmi_trials,bmi_trials];
            
            m_d_performance{m,d} = movmean(ismember(bmi_trials,types_vec([3,4])),smooth_n,"Endpoints","discard");
            num_m_ts(m,d) = length(bmi_trials);
        end
        
    end
    m_performance{m} = movmean(ismember(cur_m_bmi_trials,types_vec([3,4])),smooth_n,"Endpoints","discard");
end


% sample period
dt = 1/30;
% store correlations and probabilities
all_rho = cell(2,1);
all_p = cell(2,1);
%% Line plots, concatenated sessions
% circle_size = 100;

% lr combined
figure
m_all = [];
x_all = [];
for m = 1:num_mice
    cur_plot = m_performance{m};
    plot(cur_plot,'Color',[0.5,0.5,0.5],'LineWidth',2)
    hold on
    m_all = [m_all;cur_plot'];
    x_all = [x_all;(1:length(cur_plot))'];
end

mdl = fitlm(x_all,m_all);
xx = linspace(min(x_all),max(x_all),20)';
[ypred,yci] = predict(mdl,xx);
%plot(xx,yci(:,1),'--','LineWidth',2,'Color','k')
%plot(xx,yci(:,2),'--','LineWidth',2,'Color','k')
plot(xx,ypred,'-','LineWidth',2,'Color','k')
% h = fill([xx',fliplr(xx')],[yci(:,2)',fliplr(yci(:,1)')],'k','LineStyle','none');
% set(h,'facealpha',.1)

% or
% P = polyfit(x_all, m_all, 1);
% xx = linspace(min(x_all),max(x_all),20)';
% Pfit = polyval(P, xx);
% plot(xx,Pfit,'LineWidth',3,'Color','k')

[all_rho{1},all_p{1}] = corr(x_all,m_all,'rows','pairwise');

%xticks([10,15,20])
%yticks([10,15,20])
axis('square')
xlabel("Smoothed BMI trial number")
ylabel("Fraction correct")
title("Progression of BMI trial performance")
ylim([0,1])

%% Line plots, each session separately

% lr combined
titles = ["Session 1";"Session 2";"Session 3";"Session 4"];
figure
cur_rho = nan.*ones(num_days,1);
cur_p = nan.*ones(num_days,1);
for d = 1:num_days
    m_d_all = [];
    x_d_all = [];
    subplot(1,num_days,d)
    for m = 1:num_mice
        cur_plot = m_d_performance{m,d};
        plot(cur_plot,'Color',[0.5,0.5,0.5],'LineWidth',2)
        hold on
        m_d_all = [m_d_all;cur_plot'];
        x_d_all = [x_d_all;(1:length(cur_plot))'];
    end
    
    mdl = fitlm(x_d_all,m_d_all);
    xx = linspace(min(x_d_all),max(x_d_all),20)';
    [ypred,yci] = predict(mdl,xx);
    %plot(xx,yci(:,1),'--','LineWidth',2,'Color','k')
    %plot(xx,yci(:,2),'--','LineWidth',2,'Color','k')
    plot(xx,ypred,'-','LineWidth',2,'Color','k')
    % h = fill([xx',fliplr(xx')],[yci(:,2)',fliplr(yci(:,1)')],'k','LineStyle','none');
    % set(h,'facealpha',.1)

    % or
    % P = polyfit(x_d_all, m_d_all, 1);
    % xx = linspace(min(x_d_all),max(x_d_all),20)';
    % Pfit = polyval(P, xx);
    % plot(xx,Pfit,'LineWidth',3,'Color','k')
    
    %xticks([10,15,20])
    %yticks([10,15,20])
    axis('square')
    xlabel("Smoothed BMI trial Number")
    ylabel("Fraction correct")
    title(titles(d))
    ylim([0,1])
    [cur_rho(d),cur_p(d)] = corr(x_d_all,m_d_all,'rows','pairwise');
end

all_rho{2} = cur_rho;
all_p{2} = cur_p;

%% Hierarchical bootstrap

% [p_boot,all_centres,all_sems] = run_H_boot_ets(squeeze(mean_trial_lengths(:,:,1)), squeeze(mean_trial_lengths(:,:,2)),true);
% 
% h_boots.all_p_boot = p_boot;
% h_boots.all_centres = all_centres;
% h_boots.all_sems = all_sems;