function [all_rho,all_p] = trial_length_progression_plot(virmen_cell,tbt_cell,smooth_n)
% 08/09/2024

% function for plotting BMI trial lengths over time.

% Plot each session separately
% Also plot all sessions concatenated
% Plot one with left and right separate, one together.
% Plots best fit line.
% Calculate correlation. 
% Smooth by smooth_n, also used for correlations.

% Only correct trials, and numbered by correct BMI trials only. 

num_mice = size(virmen_cell,1);
num_days = size(virmen_cell,2);

%% load and get trial length results
types_vec = [1,4,7,10];

m_d_trial_lengths = cell(num_mice,num_days);
m_d_trial_lengths_lrsep = cell(num_mice,num_days,2);
m_trial_lengths = cell(num_mice,1);
m_trial_lengths_lrsep = cell(num_mice,2);

num_m_ts = nan.*ones(num_mice,2);
for m = 1:num_mice
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
            m_d_trial_lengths{m,d} = trial_lengths(ismember(tbt_details(3,:),types_vec([3,4])));
            m_trial_lengths{m} = [m_trial_lengths{m};trial_lengths(ismember(tbt_details(3,:),types_vec([3,4])))];
            for j = 1:2
                m_d_trial_lengths_lrsep{m,d,j} = trial_lengths(ismember(tbt_details(3,:),types_vec(3+j-1)));
                m_trial_lengths_lrsep{m,j} = [m_trial_lengths_lrsep{m,j};trial_lengths(ismember(tbt_details(3,:),types_vec(3+j-1)))];
            end
        end
        
    end
    for j = 1:2
        num_m_ts(m,j) = length(m_d_trial_lengths{m,j});
    end
end


% sample period
dt = 1/30;
% store correlations and probabilities
all_rho = cell(4,1);
all_p = cell(4,1);
%% Line plots, concatenated sessions
% circle_size = 100;

% lr combined
figure
m_all = [];
x_all = [];
for m = 1:num_mice
    cur_plot = movmean(m_trial_lengths{m}.*dt,smooth_n);
    plot(cur_plot,'Color',[0.5,0.5,0.5],'LineWidth',2)
    hold on
    m_all = [m_all;cur_plot];
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
xlabel("Number correct BMI trial")
ylabel("BMI trial lengths (s)")
title("Progression of correct BMI trial lengths")
ylim([0,40])


% lr separate
titles = ["Progression of correct left BMI trial lengths";"Progression of correct right BMI trial lengths"];
x_labels = ["Number correct BMI left trial";"Number correct BMI right trial"];
figure
cur_rho = nan.*ones(2,1);
cur_p = nan.*ones(2,1);
for j = 1:2
    subplot(1,2,j)
    m_all_j = [];
    x_all_j = [];
    for m = 1:num_mice
        cur_plot = movmean(m_trial_lengths_lrsep{m,j}.*dt,smooth_n);
        plot(cur_plot,'Color',[0.5,0.5,0.5],'LineWidth',2)
        hold on
        m_all_j = [m_all_j;cur_plot];
        x_all_j = [x_all_j;(1:length(cur_plot))'];
    end
    %xticks([10,15,20])
    %yticks([10,15,20])
    axis('square')
    xlabel(x_labels(j))
    ylabel("BMI trial lengths (s)")
    title(titles(j))
    ylim([0,40])
    
    mdl = fitlm(x_all_j,m_all_j);
    xx = linspace(min(x_all_j),max(x_all_j),20)';
    [ypred,yci] = predict(mdl,xx);
    %plot(xx,yci(:,1),'--','LineWidth',2,'Color','k')
    %plot(xx,yci(:,2),'--','LineWidth',2,'Color','k')
    plot(xx,ypred,'-','LineWidth',2,'Color','k')
    % h = fill([xx',fliplr(xx')],[yci(:,2)',fliplr(yci(:,1)')],'k','LineStyle','none');
    % set(h,'facealpha',.1)
    
    % or
    % P = polyfit(x_all_j, m_all_j, 1);
    % xx = linspace(min(x_all_j),max(x_all_j),20)';
    % Pfit = polyval(P, xx);
    % plot(xx,Pfit,'LineWidth',3,'Color','k')
    
    [cur_rho(j),cur_p(j)] = corr(x_all_j,m_all_j,'rows','pairwise');
end
all_rho{2} = cur_rho;
all_p{2} = cur_p;
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
        cur_plot = movmean(m_d_trial_lengths{m,d}.*dt,smooth_n);
        plot(cur_plot,'Color',[0.5,0.5,0.5],'LineWidth',2)
        hold on
        m_d_all = [m_d_all;cur_plot];
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
    xlabel("Number correct BMI trial")
    ylabel("BMI trial lengths (s)")
    title(titles(d))
    ylim([0,40])
    [cur_rho(d),cur_p(d)] = corr(x_d_all,m_d_all,'rows','pairwise');
end

all_rho{3} = cur_rho;
all_p{3} = cur_p;

% lr serparate
titles = ["Session 1 left trials";"Session 1 right trials"];
x_labels = ["Number correct BMI left trial";"Number correct BMI right trial"];
figure
cur_rho = nan.*ones(num_days,2);
cur_p = nan.*ones(num_days,2);
for d = 1:num_days
    for j = 1:2
        m_d_all_j = [];
        x_d_all_j = [];
        subplot(2,num_days,d+(j-1)*num_days)
        for m = 1:num_mice
            cur_plot = movmean(m_d_trial_lengths_lrsep{m,d,j}.*dt,smooth_n);
            plot(cur_plot,'Color',[0.5,0.5,0.5],'LineWidth',2)
            hold on
            m_d_all_j = [m_d_all_j;cur_plot];
            x_d_all_j = [x_d_all_j;(1:length(cur_plot))'];
        end
        
        mdl = fitlm(x_d_all_j,m_d_all_j);
        xx = linspace(min(x_d_all_j),max(x_d_all_j),20)';
        [ypred,yci] = predict(mdl,xx);
        %plot(xx,yci(:,1),'--','LineWidth',2,'Color','k')
        %plot(xx,yci(:,2),'--','LineWidth',2,'Color','k')
        plot(xx,ypred,'-','LineWidth',2,'Color','k')
        % h = fill([xx',fliplr(xx')],[yci(:,2)',fliplr(yci(:,1)')],'k','LineStyle','none');
        % set(h,'facealpha',.1)

        % or
        % P = polyfit(x_d_all_j, m_d_all_j, 1);
        % xx = linspace(min(x_d_all_j),max(x_d_all_j),20)';
        % Pfit = polyval(P, xx);
        % plot(xx,Pfit,'LineWidth',3,'Color','k')
        
        %xticks([10,15,20])
        %yticks([10,15,20])
        axis('square')
        xlabel(x_labels(j))
        ylabel("BMI trial lengths (s)")
        title(titles(j))
        ylim([0,40])
        [cur_rho(d,j),cur_p(d,j)] = corr(x_d_all_j,m_d_all_j,'rows','pairwise');
    end
end

all_rho{4} = cur_rho;
all_p{4} = cur_p;
%% Hierarchical bootstrap

% [p_boot,all_centres,all_sems] = run_H_boot_ets(squeeze(mean_trial_lengths(:,:,1)), squeeze(mean_trial_lengths(:,:,2)),true);
% 
% h_boots.all_p_boot = p_boot;
% h_boots.all_centres = all_centres;
% h_boots.all_sems = all_sems;