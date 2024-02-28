function [bootstrap_means,sig_results,sig_check,centres] = significant_tuning_check_v3_24102022(zdata, virmen_data, tbt_details, num_shuffles, linearise_x, nbins, boot_samps, shuff_limit, xnum)
% 24/10/2022

% Function for determining cells with significant tuning.
% Using circular shuffle check.

% significance check is whether max value is higher than 99% of
% shuffles max.

% 24/10/2022 change: Nan out all BMI trial timepoints

zdim = size(zdata,2);
t_types = [1,4,7,10];


%% remove invalid data
ITI = virmen_data(8,:);
cleaned_valid = clean_valid_data(ITI);

virmen_data_clean = virmen_data(:,cleaned_valid);
%% Preprocess and bin data

initialise_params;
if ndims(zdata) == 2
    model_params.spatial = false;
    create_reg = false;
else
    create_reg = true;
end
model_params.reg_images = false;
disp("Pre-Processing")
[zdatafilt, train_mean, Rfixed] = preprocess_images_05_2021(zdata,model_params,[],[],[],create_reg);
disp("Pre-Processing Complete")   

% Remove invalid samples
zdatafilt_clean = zdatafilt(cleaned_valid,:);

num_samples = size(zdata,1);

%% Nan out bmi trials
zdatafilt(ismember(virmen_data(1,:),[3,4]),:) = nan;

%% Perform shuffles
% shuffle limit as an input
shuff_inds = randi([shuff_limit,(num_samples-shuff_limit)],num_shuffles,1);

shuff_types = [1,4];
shuff_z_binned_means_max = zeros(length(shuff_types),num_shuffles,zdim);
% Shuffle using full data (not cleaned data).
for i = 1:num_shuffles
    shuff_zdatafilt = circshift(zdatafilt,shuff_inds(i));
    
    % Bin into trials
    [z_binned, centres] = bin_neural_data_general_any(virmen_data_clean,shuff_zdatafilt(cleaned_valid,:),xnum,linearise_x,nbins);
    
    % Get mean across trials 
    z_binned_means = zeros(length(shuff_types),nbins,zdim);
    for j = 1:length(shuff_types)
        z_binned_means(j,:,:) = mean(z_binned(tbt_details(3,:)==shuff_types(j),:,:),'omitnan');
    end
    % get maximum across bins
    shuff_z_binned_means_max(:,i,:) = squeeze(max(z_binned_means,[],2,'omitnan'));
    
end
% Get maximum across left and right
shuff_means_max = squeeze(max(shuff_z_binned_means_max,[],1,'omitnan'));
%% Unshuffled
% z_binned is trials x nbins x zdim
[z_binned, centres] = bin_neural_data_general_any(virmen_data_clean,zdatafilt_clean,xnum,linearise_x,nbins);

% Use bootstrapping to calculate distributions of sample means for each
% trial type in each bin. 
% Number of resamples defined by boot_samps
% bootstrap_means is num_types x boot_samps x bins x neurons

[bootstrap_means] = calc_bootsrapped_means(z_binned,tbt_details,boot_samps,t_types);
    
%% Determine significant tuning
sig_results = zeros(length(t_types),zdim);

boot_means_mean = squeeze(mean(bootstrap_means,2,'omitnan'));
boot_means_mean_max = squeeze(max(boot_means_mean,[],2,'omitnan'));

% Get threshold as 99th percentile across shuffles
sig_check = squeeze(prctile(shuff_means_max,99,1));

% Check significance for each trial type.
for i = 1:4
    sig_results(i,:) = boot_means_mean_max(i,:) > sig_check;
end