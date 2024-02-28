function [] = plot_pva_example_trajectories(virmen_data,tbt_details,nbins,linearise_x,ex_offset,num_trials)
% 09/06/2023

% Function for plotting example trajectories of forward velocity and view
% angle against linearised position.

x_vec = [13,7,16,7];
% units of cm and cm/s
pos_scale = 0.74;

circum = 64;
V = 0.32;
alpha = -50/75*circum/V;
beta = 0.05*circum/V/2.5;

% Convert forward velocity and decoded forward velocity from measured ball
% volatges into cm/s
virmen_data(13,:) = (alpha*(virmen_data(13,:)-ex_offset(1)))*pos_scale;
virmen_data(16,:) = (alpha*(virmen_data(16,:)-ex_offset(1)))*pos_scale;

virmen_data(7,:) = wrapToPi(virmen_data(7,:));

% Remove ITI
ITI = virmen_data(8,:);
cleaned_valid = clean_valid_data(ITI);
virmen_data_clean = virmen_data(:,cleaned_valid);
trial_num = virmen_data_clean(12,:);
virmen_data_clean(6,:) = pos_scale.*(virmen_data_clean(6,:) + abs(virmen_data_clean(5,:)));

% Find selection of trials with ok view angle
good_trials = [];
for n = 1:max(trial_num)
    if max(abs(virmen_data_clean(7,trial_num==n)))<3
        good_trials = [good_trials,n];
    end
end

%% Individual trials - binned
% num trials x nbins x numx
[x_binned,centres] = bin_kin_data(virmen_data,x_vec,linearise_x,nbins);
centres = centres*pos_scale;

types_vec = [1,4,7,10];
color_vec = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];

figure
for i = 1:2
    for j = 1:2    
        cur_trials1 = find(tbt_details(3,:)==types_vec(j));
        cur_trials2 = find(tbt_details(3,:)==types_vec(j+2));
        % rn1 = randperm(length(cur_trials1),num_trials);
        % rn2 = randperm(length(cur_trials2),num_trials);
        cur_trials1 = cur_trials1(ismember(cur_trials1,good_trials));
        cur_trials2 = cur_trials2(ismember(cur_trials2,good_trials));
        rn1 = num_trials+1:2*num_trials;
        rn2 = num_trials+1:2*num_trials;
        subplot(2,2,(i-1)*2+j)
        hold on
        
        plot(centres,squeeze(x_binned(cur_trials1(rn1),:,i)),'Color',color_vec(1,:),'LineWidth',2)
 
        plot(centres,squeeze(x_binned(cur_trials2(rn2),:,i+2)),'Color',color_vec(2,:),'LineWidth',2)
        
        if i == 2
            yline(0,'--','LineWidth',2);
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