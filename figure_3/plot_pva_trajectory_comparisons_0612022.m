function [] = plot_pva_trajectory_comparisons_0612022(virmen_data,tbt_details,nbins,linearise_x,boot_samps,sub_sample,ex_offset)
% 0612022

% function for plotting comparison of forward velocity and view angle between ball and
% bmi trials (i.e. decoder output on bmi trials, to control signal on ball
% trials.

x_vec = [13,7,16,7];
% unit of cm/s
pos_scale = 0.74;
ex_mouse = 4;

circum = 64;
V = 0.32;
alpha = -50/75*circum/V;
beta = 0.05*circum/V/2.5;

virmen_data(13,:) = (alpha*(virmen_data(13,:)-ex_offset))*pos_scale;
virmen_data(16,:) = (alpha*(virmen_data(16,:)-ex_offset))*pos_scale;

virmen_data(7,:) = wrapToPi(virmen_data(7,:));

%% 
% num trials x nbins x numx
[x_binned,centres] = bin_kin_data(virmen_data,x_vec,linearise_x,nbins);
centres = centres*pos_scale;

types_vec = [1,4,7,10];


%% Bootstrapped with CIs

color_vec = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];

% Bootstrapping of trials
[bootstrap_means] = calc_bootsrapped_means_subsample(x_binned,tbt_details,boot_samps,types_vec,sub_sample);
x_binned_means = zeros(size(bootstrap_means,3),length(x_vec),length(types_vec));
  
for i = 1:length(types_vec)
    x_binned_means(:,:,i) = mean(squeeze(bootstrap_means(i,:,:,:)),'omitnan');
end

CI_vals = [2.5,97.5];

% Get confidence intervals
% CIs, upper/lower x ntypes x nbins x zdim
CIs_all = zeros(2,length(types_vec),size(bootstrap_means,3),length(x_vec));
for i = 1:size(bootstrap_means,3)
    for j = 1:length(types_vec)
        CIs_all(:,j,i,:) = prctile(squeeze(bootstrap_means(j,:,i,:)),CI_vals);
    end
end

figure
for i = 1:2
    for j = 1:2        
        subplot(2,2,(i-1)*2+j)
        hold on
        h = fill([centres,fliplr(centres)],[squeeze(CIs_all(1,j,:,i))',fliplr(squeeze(CIs_all(2,j,:,i))')],color_vec(1,:),'EdgeColor','none');
        set(h,'facealpha',.3)
        
        plot(centres,squeeze(x_binned_means(:,i,j)),'Color',color_vec(1,:),'LineWidth',2)
        
        hold on
        h = fill([centres,fliplr(centres)],[squeeze(CIs_all(1,j+2,:,i+2))',fliplr(squeeze(CIs_all(2,j+2,:,i+2))')],color_vec(2,:),'EdgeColor','none');
        set(h,'facealpha',.3)
        
        plot(centres,squeeze(x_binned_means(:,i+2,j+2)),'Color',color_vec(2,:),'LineWidth',2)
        
        if i == 2
            ylim([-2,2]);
            yline(0,'--','LineWidth',2);
        else
            % ylim for forward velocity
            ylim([0,40]);
        end
        box off
    end
end
subplot(2,2,1)
ylabel(["Forward velocity"; "(cm/s)"]);
title("Left Trials")
% axis('square')

subplot(2,2,2)
title("Right trials")
% axis('square')

subplot(2,2,3)
ylabel(["View angle"; "(rad)"]);
xlabel("Linearized position (cm)")
% axis('square')

subplot(2,2,4)
xlabel("Linearized position (cm)")
% axis('square')

%% Bootstrapped with SEM (alternative to CIs)

% CIs, upper/lower x ntypes x nbins x zdim
lims_all = zeros(2,length(types_vec),size(bootstrap_means,3),length(x_vec));
for i = 1:length(types_vec)       
    lims_all(1,i,:,:) = x_binned_means(:,:,i) - squeeze(std(squeeze(bootstrap_means(i,:,:,:))));
    lims_all(2,i,:,:) = x_binned_means(:,:,i) + squeeze(std(squeeze(bootstrap_means(i,:,:,:))));
end

figure
for i = 1:2
    for j = 1:2        
        subplot(2,2,(i-1)*2+j)
        hold on
        h = fill([centres,fliplr(centres)],[squeeze(lims_all(1,j,:,i))',fliplr(squeeze(lims_all(2,j,:,i))')],color_vec(1,:),'EdgeColor','none');
        set(h,'facealpha',.3)
        
        plot(centres,squeeze(x_binned_means(:,i,j)),'Color',color_vec(1,:),'LineWidth',2)
        
        hold on
        h = fill([centres,fliplr(centres)],[squeeze(lims_all(1,j+2,:,i+2))',fliplr(squeeze(lims_all(2,j+2,:,i+2))')],color_vec(2,:),'EdgeColor','none');
        set(h,'facealpha',.3)
        
        plot(centres,squeeze(x_binned_means(:,i+2,j+2)),'Color',color_vec(2,:),'LineWidth',2)
        
        if i == 2
            ylim([-2,2]);
            yline(0,'--','LineWidth',2);
        else
            % insert ylim for pitch
            ylim([0,40]);
        end
        box off
    end
end
subplot(2,2,1)
ylabel(["forward velocity"; "(cm/s)"]);
title("left trials")
% axis('square')

subplot(2,2,2)
title("right trials")
% axis('square')

subplot(2,2,3)
ylabel(["view angle"; "(rad)"]);
xlabel("linearized position (cm)")
% axis('square')

subplot(2,2,4)
xlabel("linearized position (cm)")
% axis('square')