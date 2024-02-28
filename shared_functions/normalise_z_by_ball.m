function [z_norm] = normalise_z_by_ball(zfilt,virmen_data,tbt_details)
% 09/09/2023

% z score neural data, using mean and std calculated from ball trials only

% Assume neural data is already filtered.

% Can't just normalise with all trials as ratio of correct ball to bmi is
% different each session (therefore would need to balance). 

% Use only correct ball trials, and use only valid samples

types_vec = [1,4,7,10];

ITI = virmen_data(8,:);
cleaned_valid = clean_valid_data(ITI);

virmen_data = virmen_data(:,cleaned_valid);
zfilt_clean = zfilt(cleaned_valid,:);
trial_num = virmen_data(12,:);

cur_trials = find(ismember(tbt_details(3,:),types_vec([1,2])));
kept_trial_num = ismember(trial_num,cur_trials);

zfilt_clean = zfilt_clean(kept_trial_num,:);

% neuron specific means and stds
z_means = mean(zfilt_clean);
z_stds = std(zfilt_clean);

z_norm = (zfilt - z_means)./z_stds;
