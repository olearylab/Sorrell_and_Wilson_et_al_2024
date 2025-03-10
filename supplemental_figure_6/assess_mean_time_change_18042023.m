function [overall_trial_means,kept_trials] = assess_mean_time_change_18042023(virmen_data,zdata,tbt_details,model_params,sub_sample,normalise_z,plot_res,av_num)
% 18/04/2023

% assess whether mean neural activity differences change throughout a
% session.

%%
types_vec = [1,4,7,10];

%% remove invalid data
ITI = virmen_data(8,:);
cleaned_valid = clean_valid_data(ITI);

virmen_data = virmen_data(:,cleaned_valid);
%% Preprocess neural data

if ndims(zdata) == 2
    model_params.spatial = false;
    create_reg = false;
else
    create_reg = true;
end
model_params.reg_images = false;

[zdatafilt, train_mean, Rfixed] = preprocess_images_05_2021(zdata,model_params,[],[],cleaned_valid,create_reg);
  

%% Subsample
trial_num = virmen_data(12,:);
if sub_sample
    [kept_trials] = subsample_trials(tbt_details,types_vec);
    kept_trials = sort(kept_trials(:));
else
    kept_trials = 1:max(trial_num);
end

kept_trial_num = ismember(trial_num,kept_trials);

zdatafilt = zdatafilt(kept_trial_num,:);
virmen_data = virmen_data(:,kept_trial_num);

%% Normalise (z-score) post filtering, before bootstrapping.
% Normalise filtered data with invalid data removed
if normalise_z 
    % Normalise according to ball trial activity
    [zdatafilt] = normalise_z_by_ball(zdatafilt,virmen_data,tbt_details);
end

%% Calculate each trials mean:

overall_trial_means = zeros(length(kept_trials),1);

trial_num = virmen_data(12,:);

for i = 1:length(overall_trial_means)
    cur_trial = kept_trials(i);

    overall_trial_means(i) = mean(zdatafilt(ismember(trial_num,cur_trial),:),'all','omitnan');
end

% Optional plotting
if plot_res
    kept_tbt = tbt_details(:,kept_trials);
    figure
    for t = 1:2
        cur_t = ismember(kept_tbt(3,:),types_vec([2*(t-1)+1,2*(t-1)+2]));
        plot(kept_trials(cur_t),overall_trial_means(cur_t),'LineWidth',2)
        hold on
    end
    
    figure
    % running average
    for t = 1:2
        cur_t = ismember(kept_tbt(3,:),types_vec([2*(t-1)+1,2*(t-1)+2]));
        cur_kept = kept_trials(cur_t);
        cur_means = movmean(overall_trial_means(cur_t),av_num);
        plot(cur_kept,cur_means,'LineWidth',2)
        hold on
    end
end