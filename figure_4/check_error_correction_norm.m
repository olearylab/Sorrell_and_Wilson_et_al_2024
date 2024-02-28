function [error_mat,x_binned] = check_error_correction_norm(virmen_data,tbt_details,mean_binned,std_binned,nbins,offsets,plot_res)
% 09/09/2023

% Function for calculating normalised view angle error
% Normalise by std of view angle in each bin
% Use mean view angle from training data (passed as mean_binned)

% Correct trial indicators in tbt_details
types_vec = [1,4,7,10];

% convert ball volatages to velocities
circum = 64;
V = 0.32;
alpha = -50/75*circum/V;
beta = 0.05*circum/V/2.5;
virmen_data(13,:) = alpha*(virmen_data(13,:)-offsets(1));
virmen_data(15,:) = -beta*(virmen_data(15,:)-offsets(2));

% Bin behavioural data
linearise_x = true;
x_vec = [6,7,15,17]; % [y position, view angle, angular velocity, decoded view angle]
virmen_data(7,:) = wrapToPi(virmen_data(7,:));
[x_binned,centres] = bin_kin_data(virmen_data,x_vec,linearise_x,nbins);

% Calculate error signal

error_mat = nan.*ones(size(x_binned,1),size(x_binned,2));
% Edited to also calculate for ball trials - legacy
mean_ind = [1,2,1,2];
% Calculate normalised heading deviations for trials of each trial type
for i = 1:4
    cur_trials = find(tbt_details(3,:)==types_vec(i));
    for n = 1:length(cur_trials)
        error_mat(cur_trials(n),:) = (squeeze(x_binned(cur_trials(n),:,2))' - squeeze(mean_binned(:,2,mean_ind(i))))./squeeze(std_binned(:,2,mean_ind(i)));
    end
end

% some optional plotting
if plot_res
    figure

    for i = 1:2
        subplot(1,2,i)
        cur_trials = find(tbt_details(3,:)==types_vec(2+i));
        cur_errors = error_mat(cur_trials,:);
        cur_yaw = squeeze(x_binned(cur_trials,:,3));
        ro = corr(cur_errors(:),cur_yaw(:));
        scatter(cur_errors(:),cur_yaw(:));
        title(ro)
    end
end