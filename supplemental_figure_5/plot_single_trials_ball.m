function [] = plot_single_trials_ball(virmen_data,tbt_details,nbins,linearise_x,ex_offset,num_trials)
% 08/08/2023

% Plot example single trials for ball trajectories

x_vec = [13,13,15,15];
types_vec = [1,4,7,10];
% Unit cm/s
pos_scale = 0.74;

% Convert to velocities
circum = 64;
V = 0.32;
alpha = -50/75*circum/V;
beta = 0.05*circum/V/2.5;

virmen_data(13,:) = (alpha*(virmen_data(13,:)-ex_offset(1)))*pos_scale;
virmen_data(16,:) = (alpha*(virmen_data(16,:)-ex_offset(1)))*pos_scale;

virmen_data(15,:) = -beta*(virmen_data(15,:)-ex_offset(2));

virmen_data(7,:) = wrapToPi(virmen_data(7,:));

% Remove ITI
ITI = virmen_data(8,:);
cleaned_valid = clean_valid_data(ITI);
virmen_data_clean = virmen_data(:,cleaned_valid);
trial_num = virmen_data_clean(12,:);

% num trials x nbins x numx
[x_binned,centres] = bin_kin_data(virmen_data,x_vec,linearise_x,nbins);
centres = centres*pos_scale;

% Find most similar and most different trajectories
all_diffs = nan.*ones(max(trial_num),max(trial_num),2);
for n = 1:max(trial_num)
    for t = 1:max(trial_num)
        all_diffs(n,t,1) = sqrt(sum((x_binned(n,:,1) - x_binned(t,:,1)).^2));
        all_diffs(n,t,2) = sqrt(sum((x_binned(n,:,3) - x_binned(t,:,3)).^2));
    end
end

% Select trials
% fv left, fv right, av left, av right
plot_trials_ball = nan.*ones(4,num_trials);
plot_trials_min = nan.*ones(4,num_trials);
plot_trials_max = nan.*ones(4,num_trials);
for i = 1:2
    cur_trials_1 = find(tbt_details(3,:)==types_vec(i));
    cur_trials_2 = find(tbt_details(3,:)==types_vec(i+2));
    rn1 = randi(length(cur_trials_1),num_trials,1);
    rn2 = rn1;
    
    plot_trials_ball(i,:) = cur_trials_1(rn1);
    plot_trials_ball(i+2,:) = cur_trials_1(rn2);
    
    [~,inds] = mink(all_diffs(plot_trials_ball(i,1),cur_trials_2,1),num_trials);
    plot_trials_min(i,:) = cur_trials_2(inds);
    [~,inds] = maxk(all_diffs(plot_trials_ball(i,1),cur_trials_2,1),num_trials);
    plot_trials_max(i,:) = cur_trials_2(inds);
    [~,inds] = mink(all_diffs(plot_trials_ball(i+2,1),cur_trials_2,2),num_trials);
    plot_trials_min(i+2,:) = cur_trials_2(inds);
    [~,inds] = maxk(all_diffs(plot_trials_ball(i+2,1),cur_trials_2,2),num_trials);
    plot_trials_max(i+2,:) = cur_trials_2(inds);
        
end
%% Individual trials - binned

color_vec = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];

% Plot very similar ball and bmi trials
figure
for i = 1:4  
    subplot(2,2,i)   
        
    plot(centres,squeeze(x_binned(plot_trials_ball(i,:),:,i)),'Color',color_vec(1,:),'LineWidth',2)
    hold on
    plot(centres,squeeze(x_binned(plot_trials_min(i,:),:,i)),'Color',color_vec(2,:),'LineWidth',2)

    if i > 2
        yline(0,'--','LineWidth',2);
    end
    box off
end
subplot(2,2,1)
ylabel(["Forward velocity"; "(cm/s)"]);
title("Left Trials")
% axis('square')

subplot(2,2,2)
title("Right trials")
% axis('square')

subplot(2,2,3)
ylabel(["Angular Velocity"; "(rad/s)"]);
xlabel("Linearized position (cm)")
% axis('square')

subplot(2,2,4)
xlabel("Linearized position (cm)")
% axis('square')

% plot very different ball and bmi trials
figure
for i = 1:4  
    subplot(2,2,i)   
        
    plot(centres,squeeze(x_binned(plot_trials_ball(i,:),:,i)),'Color',color_vec(1,:),'LineWidth',2)
    hold on
    plot(centres,squeeze(x_binned(plot_trials_max(i,:),:,i)),'Color',color_vec(2,:),'LineWidth',2)

    if i > 2
        yline(0,'--','LineWidth',2);
    end
    box off
end
subplot(2,2,1)
ylabel(["Forward velocity"; "(cm/s)"]);
title("Left Trials")
% axis('square')

subplot(2,2,2)
title("Right trials")
% axis('square')

subplot(2,2,3)
ylabel(["Angular Velocity"; "(rad/s)"]);
xlabel("Linearized position (cm)")
% axis('square')

subplot(2,2,4)
xlabel("Linearized position (cm)")
% axis('square')