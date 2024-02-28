function [bootstats_center,bootstats_sem,h_boots] = plot_many_noise_corrs(full_pairs_means_cell)
% 05/07/2023

% Plot example session clustered neurons
% Plot summary data

% Example session noise correlation visualisation
ex_md = [4,3];
ex_pairs = full_pairs_means_cell{ex_md(1),ex_md(2)};
titles = ["Ball Left";"Ball Right";"BMI Left";"BMI Right"];
% Cluster data for visualisation. Cluster according to ball left trial
T = clusterdata(squeeze(ex_pairs(1,:,:)),1);
for i = 1:4
    % Sort clusters
    [a,inds] = sort(T);
    figure 
    imagesc(squeeze(ex_pairs(i,inds,inds)),[-1,1])
    colormap('redblue')
    xticks([])
    yticks([])
    title(titles(i));
    axis('square')
end

%% Summary

num_mice = size(full_pairs_means_cell,1);
num_days = size(full_pairs_means_cell,2);

full_means = nan.*ones(4,num_mice,num_days);

full_abs_means = nan.*ones(4,num_mice,num_days);

% Mean of positive correlations
full_pos_means = nan.*ones(4,num_mice,num_days);

full_sim = nan.*ones(4,4,num_mice,num_days);


for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(full_pairs_means_cell{m,d})
            cur_pairs = full_pairs_means_cell{m,d};
            for i = 1:4
                cur_block = squeeze(cur_pairs(i,:,:));
                % remove repeated correlations
                low_t = tril(nan.*ones(size(cur_block)));
                
                cur_block = cur_block + low_t;
                
                cur_block_ready = cur_block(:);
                full_means(i,m,d) = mean(cur_block_ready,'omitnan');
                full_abs_means(i,m,d) = mean(abs(cur_block_ready),'omitnan');
                full_pos_means(i,m,d) = mean(cur_block_ready(cur_block_ready>=0),'omitnan');
                for j = 1:4
                    cur_block_j = squeeze(cur_pairs(j,:,:));
                    % remove repeated correlations
                    low_t = tril(nan.*ones(size(cur_block_j)));
                
                    cur_block_j = cur_block_j + low_t;
                
                    cur_block_ready_j = cur_block_j(:);
                    % Calculate cosine similairty
                    full_sim(i,j,m,d) = sum(cur_block_ready(~isnan(cur_block_ready)).*cur_block_ready_j(~isnan(cur_block_ready_j)))/(norm(cur_block_ready(~isnan(cur_block_ready)))*norm(cur_block_ready_j(~isnan(cur_block_ready_j))));
                    
                end
            end
        end
    end
end

% Plot cosine similarities
num_sess = num_mice*num_days;
figure
plot_off = linspace(-0.4,0.4,num_sess);

i = 1;

lr_combined_sims = mean([squeeze(full_sim(1,3,:))';squeeze(full_sim(2,4,:))'],'omitnan');

scatter(plot_off+i.*ones(1,num_sess),lr_combined_sims,'filled','k')
hold on
plot([i+plot_off(1),i+plot_off(end)],[mean(lr_combined_sims,'omitnan'),mean(lr_combined_sims,'omitnan')],'k','LineWidth',2)
title(["Similarity of";"noise correlations"])
ylabel("Cosine similarity")
ylim([0,1])
xlim([0.5,1.5])
xticklabels([])

% Magnitude of positive noise correlations
% l/R combined
scatter_size = 50;

% BMI ball scatter - positive only 
lr_combined_pos_means = nan.*ones(2,num_mice,num_days);
for i = 1:2
    lr_combined_pos_means(i,:,:) = mean(full_pos_means([(i-1)*2+1,i*2],:,:));
end

figure
scatter(lr_combined_pos_means(1,:),lr_combined_pos_means(2,:),scatter_size,'filled','k')
hold on
plot([0.15,0.35],[0.15,0.35],'--','Color',[0.5,0.5,0.5],'LineWidth',2)
% title("Noise Correlations")
title("Mean positive noise correlations")
xlabel("Ball trials")
ylabel("BMI trials")
axis('square')

%% Hierarchical bootstrap

lr_combined_sims = (squeeze(full_sim(1,3,:,:))+squeeze(full_sim(2,4,:,:)))/2;

stats_data = nan.*ones(num_mice,num_days);
for m = 1:num_mice

    cur_data = lr_combined_sims(m,:);

    stats_data(m,1:sum(~isnan(cur_data))) = cur_data(~isnan(cur_data));
end

% Similarities
rng(1);
boot_samps = 1000;
num_trials = 4;
bootstats = get_bootstrapped_equalsamples(stats_data,boot_samps,num_trials,'mean');
%Get mean and SEM of bootstrapped samples:
bootstats_sem = std(bootstats);
bootstats_center = mean(bootstats);

% only positive correlations
rng(1)
[p_boot_pos,h_centre_pos,h_sem_pos] = run_H_boot_ets(squeeze(lr_combined_pos_means(1,:,:)), squeeze(lr_combined_pos_means(2,:,:)), true);
h_boots.p_boot_pos = p_boot_pos;
h_boots.centres_pos = h_centre_pos;
h_boots.sem_pos = h_sem_pos;