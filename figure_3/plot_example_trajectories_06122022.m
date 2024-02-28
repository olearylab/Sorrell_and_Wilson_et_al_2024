function [] = plot_example_trajectories_06122022(virmen_cell, tbt_cell, ex_md, offsets, linearise_x, nbins)
% 06/12/2022

% Plot example trajectories for forward velocity and angular velocity to compare between bmi and
% control trials.

x_vec = [13,15];
%% 
% 2 x 2
types_vec = [1,4,7,10];

virmen_data = virmen_cell{ex_md(1),ex_md(2)};
tbt_details = tbt_cell{ex_md(1),ex_md(2)};

% remove ITI
ITI = virmen_data(8,:);
cleaned_valid = clean_valid_data(ITI);
virmen_data = virmen_data(:,cleaned_valid);

% convert into velocity
circum = 64;
V = 0.32;
alpha = -50/75*circum/V;
beta = 0.05*circum/V/2.5;
virmen_data(13,:) = alpha*(virmen_data(13,:)-offsets(1));
virmen_data(15,:) = -beta*(virmen_data(15,:)-offsets(2));

% Bin data: trials x nbins x xdim
[x_binned,centres] = bin_kin_data(virmen_data,x_vec,linearise_x,nbins);

%% Plot
% Plot bootstrapped mean and confidence interval

plot_colours = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];
figure
pos_scale = .74;
centres = centres*pos_scale;

boot_samps = 100;
sub_sample = false;
[bootstrap_means] = calc_bootsrapped_means_subsample(x_binned,tbt_details,boot_samps,types_vec,sub_sample);
x_binned_means = zeros(size(bootstrap_means,3),length(x_vec),length(types_vec));
  
for i = 1:length(types_vec)
    x_binned_means(:,:,i) = mean(squeeze(bootstrap_means(i,:,:,:)),'omitnan');
end

CI_vals = [2.5,97.5];

% CIs, upper/lower x ntypes x nbins x zdim
CIs_all = zeros(2,length(types_vec),size(bootstrap_means,3),length(x_vec));
for i = 1:size(bootstrap_means,3)
    for j = 1:length(types_vec)
        CIs_all(:,j,i,:) = prctile(squeeze(bootstrap_means(j,:,i,:)),CI_vals);
    end
end

% figure
for i = 1:2
    for j = 1:2        
        subplot(2,2,(i-1)*2+j)
        hold on
        h = fill([centres,fliplr(centres)],[squeeze(CIs_all(1,j,:,i))',fliplr(squeeze(CIs_all(2,j,:,i))')],plot_colours(1,:),'LineStyle','none');
        set(h,'facealpha',.3)
        
        plot(centres,squeeze(x_binned_means(:,i,j)),'Color',plot_colours(1,:),'LineWidth',2)
        
        hold on
        h = fill([centres,fliplr(centres)],[squeeze(CIs_all(1,j+2,:,i))',fliplr(squeeze(CIs_all(2,j+2,:,i))')],plot_colours(2,:),'LineStyle','none');
        set(h,'facealpha',.3)
        
        plot(centres,squeeze(x_binned_means(:,i,j+2)),'Color',plot_colours(2,:),'LineWidth',2)
        
        if i == 2
            yline(0,'--','LineWidth',2);
        end
        box off
    end
end

ax11 = subplot(2,2,1);
title("Left trials")
ylabel(["Forward velocity"; "(cm/s)"])
% axis('square')

ax12 = subplot(2,2,2);
title("Right trials")
% axis('square')

ax13 = subplot(2,2,3);
ylabel(["Angular velocity"; "(rad/s)"])
xlabel("Linearized position (cm)")
% axis('square')

ax14 = subplot(2,2,4);
xlabel("Linearized position (cm)")
yticks([-0.5,0,0.5,1])
% axis('square')

linkaxes([ax11,ax12])
linkaxes([ax13,ax14])

%% Alternative plot: Bootstrapped with SEM

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
        h = fill([centres,fliplr(centres)],[squeeze(lims_all(1,j,:,i))',fliplr(squeeze(lims_all(2,j,:,i))')],plot_colours(1,:),'EdgeColor','none');
        set(h,'facealpha',.3)
        
        plot(centres,squeeze(x_binned_means(:,i,j)),'Color',plot_colours(1,:),'LineWidth',2)
        
        hold on
        h = fill([centres,fliplr(centres)],[squeeze(lims_all(1,j+2,:,i))',fliplr(squeeze(lims_all(2,j+2,:,i))')],plot_colours(2,:),'EdgeColor','none');
        set(h,'facealpha',.3)
        
        plot(centres,squeeze(x_binned_means(:,i,j+2)),'Color',plot_colours(2,:),'LineWidth',2)
        
        if i == 2
            yline(0,'--','LineWidth',2);
        end
        box off
    end
end
ax11 = subplot(2,2,1);
title("left trials")
ylabel(["forward velocity"; "(cm/s)"])
% axis('square')

ax12 = subplot(2,2,2);
title("right trials")
yticklabels([])
% axis('square')

ax13 = subplot(2,2,3);
ylabel(["angular velocity"; "(rad/s)"])
xlabel("linearized position (cm)")
% axis('square')

ax14 = subplot(2,2,4);
xlabel("linearized position (cm)")
yticks([-0.5,0,0.5,1])
yticklabels([])
% axis('square')

linkaxes([ax11,ax12])
linkaxes([ax13,ax14])