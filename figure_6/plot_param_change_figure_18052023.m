function [h_boots] = plot_param_change_figure_18052023(sig_results_all,sig_summary_cell,all_full_results_both,full_all_param_means,sig_check_all,bootstrap_means,centres,ex_neu,ex_dir,ex_md)
% 18/05/2023

colours_vec = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];
centres = centres*0.74;
CI_vals = [2.5,97.5];
t_types = [1,4,7,10];
%% Summary plots
figure

num_mice = size(all_full_results_both,1);
num_days = size(all_full_results_both,2);

% Want proportions of: gain, loss, both no change, pos, amp, width
all_proportions_both = nan.*ones(num_mice,num_days,2,6);

num_tuned = zeros(num_mice,num_days,4);

for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(all_full_results_both{m,d})
            full_results_both = all_full_results_both{m,d};
            cur_sig = sig_results_all{m,d};
            cur_sig_summary = sig_summary_cell{m,d};
            total_tuned = sum(cur_sig==1)>0;
            % Separately for left and right tuned neurons
            for i = 1:2
                cur_neurons = sum(cur_sig([i,i+2],:)==1)>0;
                cur_both = squeeze(cur_sig_summary(i,1,:)==1);
                all_proportions_both(m,d,i,1) = sum(squeeze(cur_sig_summary(i,2,:))==1)./sum(cur_neurons);
                all_proportions_both(m,d,i,2) = sum(squeeze(cur_sig_summary(i,3,:))==1)./sum(cur_neurons);
                % Change to sum(cur_neurons) instead of sum(cur_both) to get
                % proportion of tuned to left/right matching above.
                all_proportions_both(m,d,i,3) = squeeze(sum(sum(squeeze(full_results_both(cur_both,i,:))==0,2)==3))./sum(cur_neurons);
                all_proportions_both(m,d,i,4:6) = squeeze(sum(squeeze(full_results_both(cur_both,i,:))))./sum(cur_neurons);
            end

            for i = 1:4
                num_tuned(m,d,i) = sum(cur_sig(i,:));
            end
        end
    end
end

% average over left and right
all_proportions_both_dir = squeeze(mean(all_proportions_both,3,'omitnan'));

%% Violin plot version
figure
x_labels = ["Gain";"Loss";"No change";"Location";"Amplitdue";"Width"];

violin_cell = cell(1,6);
all_proportions_both_dir = permute(all_proportions_both_dir,[3,1,2]);
for i = 1:6
    violin_cell{i} = all_proportions_both_dir(i,:);
end

violin(violin_cell,'facecolor',[0.5,0.5,0.5],'medc',[],'edgecolor',[])
title('Tuning curve changes')
ylabel('Fraction of neurons')
xticklabels(x_labels)
xtickangle(45)
box off
% axis('square')

%% Plot example tuning curves
figure
sig_results = sig_results_all{ex_md(1),ex_md(2)};
cur_sig_check = sig_check_all{ex_md(1),ex_md(2)};
cur_neurons = sum(sig_results==1)>0;

bootstrap_means = bootstrap_means(:,:,:,cur_neurons);
sig_results = sig_results(:,cur_neurons);
cur_sig_check = cur_sig_check(cur_neurons);
zdim = size(sig_results,2);

z_binned_means = zeros(size(bootstrap_means,3),zdim,length(t_types));
  
for i = 1:length(t_types)
    z_binned_means(:,:,i) = mean(squeeze(bootstrap_means(i,:,:,:)),'omitnan');
end

% CIs, upper/lower x ntypes x nbins x zdim
CIs_all = zeros(2,length(t_types),size(bootstrap_means,3),zdim);
for i = 1:size(bootstrap_means,3)
    for j = 1:length(t_types)
        CIs_all(:,j,i,:) = prctile(squeeze(bootstrap_means(j,:,i,:)),CI_vals);
    end
end

% Plot example neurons
titles = ["Gain";"Loss";"No Change";"Location";"Amplitude";"Width"];
num_examples = size(ex_neu,2);
for i = 1:6

    for j = 1:num_examples
        subplot(2,6,(j-1)*6+i)

        cur_n = ex_neu(i,j);
        cur_dir = ex_dir(i,j); % 1 or 2
        for t = 1:2
            hold on
            h = fill([centres,fliplr(centres)],[squeeze(CIs_all(1,cur_dir + 2*(t-1),:,cur_n))',fliplr(squeeze(CIs_all(2,cur_dir + 2*(t-1),:,cur_n))')],colours_vec(t,:),'LineStyle','none');
            set(h,'facealpha',.3)
            % Means
            plot(centres,z_binned_means(:,cur_n,cur_dir + 2*(t-1)),'Linewidth',3,'Color',colours_vec(t,:))

        end
        xticks([0,200])
        yline(cur_sig_check(cur_n),'--','LineWidth',2); % singificance threshold
    end
end

for i = 1:6
    subplot(2,6,i)
    title(titles(i))
    xticks([])
end

subplot(2,6,1)
ylabel("\DeltaF/F")

subplot(2,6,7)
ylabel("\DeltaF/F")
xlabel("Maze position (cm)")

%% Plot example neurons with SDs - alternative version

% CIs, upper/lower x ntypes x nbins x zdim
lims_all = zeros(2,length(t_types),size(bootstrap_means,3),zdim);
for i = 1:length(t_types)       
    lims_all(1,i,:,:) = z_binned_means(:,:,i) - squeeze(std(squeeze(bootstrap_means(i,:,:,:))));
    lims_all(2,i,:,:) = z_binned_means(:,:,i) + squeeze(std(squeeze(bootstrap_means(i,:,:,:))));
end

figure
% Plot example neurons
titles = ["gain";"loss";"no change";"location";"amplitude";"width"];
num_examples = size(ex_neu,2);
for i = 1:6

    for j = 1:num_examples
        subplot(2,6,(j-1)*6+i)

        cur_n = ex_neu(i,j);
        cur_dir = ex_dir(i,j); % 1 or 2
        for t = 1:2
            hold on
            h = fill([centres,fliplr(centres)],[squeeze(lims_all(1,cur_dir + 2*(t-1),:,cur_n))',fliplr(squeeze(lims_all(2,cur_dir + 2*(t-1),:,cur_n))')],colours_vec(t,:),'LineStyle','none');
            set(h,'facealpha',.3)
            % Means
            plot(centres,z_binned_means(:,cur_n,cur_dir + 2*(t-1)),'Linewidth',3,'Color',colours_vec(t,:))
        end
        xticks([0,200])
        yline(cur_sig_check(cur_n),'--','LineWidth',2); % Significance threshold
    end
end

for i = 1:6
    subplot(2,6,i)
    title(titles(i))
end

subplot(2,6,1)
ylabel("\DeltaF/F")

subplot(2,6,7)
ylabel("\DeltaF/F")
xlabel("maze position (cm)")

%% Hierarchical bootstrap

rng(1);
boot_samps = 1000;
num_trials = 4;

all_centres = nan.*ones(size(all_proportions_both_dir,1),1);
all_sems = nan.*ones(size(all_proportions_both_dir,1),1);

% Move Nans to the end of rows
stats_data = nan.*ones(size(all_proportions_both_dir,1),num_mice,num_days);
for m = 1:num_mice

    for i = 1:size(all_proportions_both_dir,1)
        cur_data = squeeze(all_proportions_both_dir(i,m,:));

        stats_data(i,m,1:sum(~isnan(cur_data))) = cur_data(~isnan(cur_data));
    end

end

for i = 1:size(all_proportions_both_dir,1)
    bootstats = get_bootstrapped_equalsamples(squeeze(stats_data(i,:,:)),boot_samps,num_trials,'mean');
    %Get mean and SEM of bootstrapped samples:
    all_sems(i) = std(bootstats);
    all_centres(i) = mean(bootstats);
end

h_boots.all_centres = all_centres;
h_boots.all_sems = all_sems;