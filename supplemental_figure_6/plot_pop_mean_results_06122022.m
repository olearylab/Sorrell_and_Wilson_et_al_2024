function [h_boots] = plot_pop_mean_results_06122022(all_popboot_mean,all_overall_n_means,centres)
% 06/12/2022

% plot results for mean value comparisons across trial types.

num_mice = size(all_popboot_mean,1);
num_days = size(all_popboot_mean,2);
nbins = length(centres);
colour_vec = [[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0]];

%% Plot mean changes by bin

% Plot average across mice and sessions, with faint lines for average for
% individual mice.

% Average for both directions:
figure
% mean of differences for each mouse
m_means = nan.*ones(num_mice,length(centres));

md_means = nan.*ones(num_mice,num_days,length(centres));
for m = 1:num_mice
    cur_m = nan.*ones(num_days,length(centres));
    for d = 1:num_days
        if ~isempty(all_popboot_mean{m,d})
            cur_means = all_popboot_mean{m,d};
            av_lr = zeros(2,length(centres));
            for i = 1:2
                % Average across left and right
                av_lr(i,:) = mean(cur_means([(i-1)*2+1,i*2],:));
            end  
            % BMI - ball means
            cur_m(d,:) = av_lr(2,:)- av_lr(1,:);
          
            
        end
        
    end
    m_means(m,:) = squeeze(mean(cur_m,'omitnan'));
    md_means(m,1:sum(~isnan(cur_m(:,1))),:) = cur_m(~isnan(cur_m(:,1)),:);
end

yline(0,'--','LineWidth',2);
hold on

for m = 1:num_mice

    plot(centres,squeeze(m_means(m,:)),'-','LineWidth',1.5,'Color',[0.5,0.5,0.5])

end


plot(centres,squeeze(mean(m_means(:,:),'omitnan')),'LineWidth',3,'Color','k')
xlabel("Maze position (cm)")
xticks([0,200])
ylabel(["Mean normalized activity";"differece, (BMI - ball) (a.u.)"])
box off
axis('square')

%% Plot overall mean changes by neuron

n_diff_cell = cell(num_mice,1);
% calculate cumulative distribution for each mouse
num_points = 1000;
cum_line = zeros(num_mice,num_points);
m_vals = zeros(num_mice,num_points);

% prepare for h bootstrapping. proportion increasing
prop_ready = nan.*ones(num_mice,num_days);
for m = 1:num_mice
    cur_m_diffs = [];
    d_ind = 0;
    for d = 1:num_days
        if ~isempty(all_overall_n_means{m,d})
            d_ind = d_ind+1;
            cur_means = all_overall_n_means{m,d};

            cur_diffs = cur_means(2,:)-cur_means(1,:);
            
            cur_m_diffs = [cur_m_diffs,cur_diffs];
            
            prop_ready(m,d_ind) = sum(cur_diffs>0)/length(cur_diffs);
            
        end
    end
    n_diff_cell{m} = cur_m_diffs;
    cur_vals = linspace(min(cur_m_diffs),max(cur_m_diffs),num_points);
    m_vals(m,:) = cur_vals;
    for i = 1:num_points
        cum_line(m,i) = sum(cur_m_diffs<cur_vals(i))/length(cur_m_diffs);
    end
end

all_diffs = [];
for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(all_overall_n_means{m,d})
            cur_means = all_overall_n_means{m,d};

            cur_diffs = cur_means(2,:)-cur_means(1,:);
            
            all_diffs = [cur_m_diffs,cur_diffs];
            
        end
    end
end
all_vals = linspace(min(all_diffs),max(all_diffs),num_points);
all_cum_line = zeros(1,num_points);
% get cumulative distribution
for i = 1:num_points
    all_cum_line(i) = sum(all_diffs<all_vals(i))/length(all_diffs);
end

% Run H-bootstrapping
rng(1);
boot_samps = 1000;
num_trials = 4;
bootstats = get_bootstrapped_equalsamples(prop_ready,boot_samps,num_trials,'mean');
%Get mean and SEM of bootstrapped samples:
h_boots.all_sems = std(bootstats);
h_boots.all_centres = mean(bootstats);

figure

for m = 1:num_mice
    plot(m_vals(m,:),cum_line(m,:),'-','LineWidth',1.5,'Color',[0.5,0.5,0.5])
    hold on
end

plot(all_vals,all_cum_line,'-','LineWidth',3,'Color','k')
xline(0,'--','LineWidth',2);
yline(0.5,'--','LineWidth',2);
box off
axis('square')
xlabel(["Mean normalized activity";"difference, (BMI - ball) (a.u.)"])
ylabel(["Cumulative fraction"; "of neurons"])

