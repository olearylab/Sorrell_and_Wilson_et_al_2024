function [h_boots] = plot_offline_neuron_results_all_w_example_01122022(mice,data_days,example_m_d,reps,training_trials,ex_weights,ex_stat,ex_im)
% 01/12/2022

% Function for making paper ready figures showing comparison between pixel
% and neuron decoding.

% For plotting combined
plot_rows = 2;

% path to results files
res_path = "offline_results/";
res_path_n = "offline_neuron_results_CNN_16112022/";

line_width = 3;
circle_size = 100;
label_size = 50;

colours_vec = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];

% example is from mouse 4
ex_mouse = 4;

num_mice = length(mice);
n_days = zeros(num_mice,1);
for m = 1:num_mice
    n_days(m) = length(data_days{m});
end
num_days = max(n_days);

%% Plot

figure

% set figure size if desired
set(gcf,'position',[1,42,1440,508])

plot_nums = [1,7,3,4]; % ypos, VA, forward velocity, angular velocity
ex_offsets = [1.4804,1.4844];
Titles = ["Y Position", "View Angle", "Forward Velocity", "Angular Velocity"];

ytick_vec = [[0 1 2];[-2,0,2];[0 0.2 0.4];[-1, 0, 1]];

% Example decoding

results_struct = importdata(res_path + "res_struct_" + example_m_d + "_rep1_" + training_trials + ".mat");
t = (1:size(results_struct.xtest(results_struct.test_valid,:),1))/results_struct.model_params.fs;

% subtract offset from pitch and yaw, and convert into velocities
circum = 64;
V = 0.32;
alpha = -50/75*circum/V;
beta = 0.05*circum/V/2.5;
results_struct.xtest(:,3) = alpha*(results_struct.xtest(:,3)-ex_offsets(1));
results_struct.xprediction(:,3) = alpha*(results_struct.xprediction(:,3)-ex_offsets(1));
results_struct.xtest(:,4) = -beta*(results_struct.xtest(:,4)-ex_offsets(2));
results_struct.xprediction(:,4) = -beta*(results_struct.xprediction(:,4)-ex_offsets(2));

% convert from vu into m
vu_conv = 0.0074;
results_struct.xtest(:,3) = results_struct.xtest(:,3)*vu_conv;
results_struct.xprediction(:,3) = results_struct.xprediction(:,3)*vu_conv;
results_struct.xtest(:,1) = results_struct.xtest(:,1)*vu_conv;
results_struct.xprediction(:,1) = results_struct.xprediction(:,1)*vu_conv;

ylabel_vec = ["(m)";"(rad)";"(m/s)";"(rad/s)"];

for i = 1:4
    ax = subplot(plot_rows,4,i);
    
    max1 = max(results_struct.xtest(:,plot_nums(i)));
    max2 = max(results_struct.xprediction(:,plot_nums(i)));
    
    min1 = min(results_struct.xtest(:,plot_nums(i)));
    min2 = min(results_struct.xprediction(:,plot_nums(i)));
    
    baseline = min([min1,min2]);
    max_val = max([max1,max2]);
    
    if i ==2
        max_val = pi;
        baseline = -pi;
    end
    
    % Specificly to make slightly nicer for DW113 20210701 example
    if i ==4
        max_val = 1.5;
        baseline = -1.5;
    end
    
    hold on
    plot(t,results_struct.xtest(results_struct.test_valid,plot_nums(i)),'LineWidth',line_width,'Color',[0.5,0.5,0.6])

    plot(t,results_struct.xprediction(results_struct.test_valid,plot_nums(i)),'LineWidth',line_width,'Color',colours_vec(ex_mouse,:))
    title(Titles(i))
    xlabel("Time (s)")
    xlim([0,60])
    ylim([baseline,max_val]);
    
    xticks([0, 60])
    xticklabels({'0','60'})
    yticks(ytick_vec(i,:))
    ylabel(ylabel_vec(i))
    
    box off
    axis('square')
end

subplot(plot_rows,4,1)

% legend on view angle figure
subplot(plot_rows,4,2)
legend('True', 'Decoded')

subplot(plot_rows,4,4)
ylim([-1.5,1.5])

% Summary results

%% R^2
p_R2 = [];
n_R2 = [];

p_R2_all = [];
n_R2_all = [];

p_R2_all_mat = nan.*ones(num_mice,num_days,4);
n_R2_all_mat = nan.*ones(num_mice,num_days,4);

for m = 1:length(mice)
    data_d = data_days{m};
    m_vec = zeros(length(data_d),4);
    m_vec_n = zeros(length(data_d),4);
    for d = 1:length(data_d)
        md_vec = zeros(reps,4);
        md_vec_n = zeros(reps,4);
        for r = 1:reps
            results_struct = importdata(res_path + "res_struct_" + mice(m) + "_" + data_d(d) + "_rep"+ r + "_" + training_trials + ".mat");
            md_vec(r,:) = results_struct.all_res(2,plot_nums);
            results_struct = importdata(res_path_n + "res_struct_" + mice(m) + "_" + data_d(d) + "_rep"+ r + "_" + training_trials + ".mat");
            md_vec_n(r,:) = results_struct.all_res(2,plot_nums);
        end
        m_vec(d,:) = mean(md_vec);
        m_vec_n(d,:) = mean(md_vec_n);
    end
    % Mice averages
    p_R2 = [p_R2;mean(m_vec,1,'omitnan')];
    n_R2 = [n_R2;mean(m_vec_n,1,'omitnan')];
    
    % all_res
    p_R2_all = [p_R2_all;m_vec];
    n_R2_all = [n_R2_all;m_vec_n];
    
    p_R2_all_mat(m,1:size(m_vec,1),:) = m_vec;
    n_R2_all_mat(m,1:size(m_vec_n,1),:) = m_vec_n;
    
    for i = 1:4
        ax = subplot(plot_rows,4,8+i-4);
        hold on
%% 07032023: scatter cells against pixels
        scatter(m_vec_n(:,i),m_vec(:,i),circle_size,'k','filled')

        yticks([0 0.5 1])
        yticklabels({'0','0.5','1'})
        xticks([0 0.5 1])
        xticklabels({'0','0.5','1'})
        
        xlabel("Cells R^{2}")
    end
end
for i = 1:4
    subplot(plot_rows,4,8+i-4)
    plot([0,1],[0,1],'--','color','k')
    axis('square')
end

subplot(plot_rows,4,8+1-4)
ylabel("Pixels R^{2}")

%% Example weights
ex_w_num = 1;
figure
set(gcf,'position',[1,42,1440,254])
% Mean Downsampled Image 
ax = subplot(1,4,1);
mm = mean(ex_im(:));
ss = std(ex_im(:));
ll = mm - ss;
uu = mm + 5*ss;
imagesc(ex_im,[ll,uu])
axis('square')
box on

xticks([])
yticks([])
colormap(ax,'gray')
title("Mean Downsampled Image")
set(gca,'Position',[0.13,0.11,0.126394067899881,0.815])
hold on
% scalebar (1.2 microns per pixel)
plot([10,10+(100/4.8)],[15,15],'LineWidth',4,'Color','w')
   
% Transpose needed for mice 1-3
cur_weights = reshape(ex_weights(1:end-1,plot_nums(ex_w_num)),128,128)';

up_lim = max(cur_weights(:));
low_lim = min(cur_weights(:));

lim_val = max(abs([low_lim,up_lim]));

% Need to downsample the neuron masks
[ex_mask] = downsample_neuron_masks(ex_stat);

% Weights
ax1 = subplot(1,4,2);
imagesc(cur_weights,[-lim_val,lim_val])
axis('square')
box on

xticks([])
yticks([])
colormap(ax1,'redblue')
colorbar
title("Y Position Weights")

cur_n_weights = cur_weights;
cur_n_weights(ex_mask==0) = nan;

% Only weights within neurons
ax2 = subplot(1,4,3);
h = imagesc(cur_n_weights,[-lim_val,lim_val]);
set(h, 'AlphaData', ~isnan(cur_n_weights))
axis('square')
box on
xticks([])
yticks([])
colormap(ax2,'redblue')
set(gca,'Color',[0.8,0.8,0.8])
colorbar
title("Weights in Neurons")


cur_o_weights = cur_weights;
cur_o_weights(ex_mask==1) = nan;

% Only weights outside neurons
ax3 = subplot(1,4,4);
h = imagesc(cur_o_weights,[-lim_val,lim_val]);
set(h, 'AlphaData', ~isnan(cur_o_weights))
axis('square')
box on
xticks([])
yticks([])
colormap(ax3,'redblue')
set(gca,'Color',[0.8,0.8,0.8])
colorbar
title("Weights Outside Neurons")

%% Hierarcical bootstrapping
all_p_boot = nan.*ones(4,1);
all_centres = nan.*ones(4,2);
all_sems = nan.*ones(4,2);

for i = 1:4
    [all_p_boot(i),all_centres(i,:),all_sems(i,:)] = run_H_boot_ets(squeeze(p_R2_all_mat(:,:,i)), squeeze(n_R2_all_mat(:,:,i)),false);
end

h_boots.all_p_boot = all_p_boot;
h_boots.all_centres = all_centres;
h_boots.all_sems = all_sems;

% also compare position to forward velocity
[h_boots.p_boot1,h_boots.centres1,h_boots.sems1] = run_H_boot_ets(squeeze(p_R2_all_mat(:,:,1)), squeeze(p_R2_all_mat(:,:,3)),false);

% also compare view angle to angular velocity
[h_boots.p_boot2,h_boots.centres2,h_boots.sems2] = run_H_boot_ets(squeeze(p_R2_all_mat(:,:,2)), squeeze(p_R2_all_mat(:,:,4)),false);