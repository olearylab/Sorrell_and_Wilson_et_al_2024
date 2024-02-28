function [full_results,full_results_overlap,full_results_both,all_params,mean_curve_params] = assess_tuning_curve_param_changes_alt_20022023(virmen_data,zdata,tbt_details,linearise_x,nbins,boot_samps,CI_vals,xnum,plot_res,sig_check,sub_sample)
% 20/02/2023

% Assess changes in tuning curves by looking at changes in peak location,
% peak width, and peak amplitude.
% Includes option to subsample

% New version, assess peak width as number of bins above significance
% threshold.

% Significance if BMI mean outside CI for ball, and ball mean outside CI
% for BMI

t_types = [1,4,7,10];

%% remove invalid data
ITI = virmen_data(8,:);
cleaned_valid = clean_valid_data(ITI);

virmen_data = virmen_data(:,cleaned_valid);
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
[zdatafilt, train_mean, Rfixed] = preprocess_images_05_2021(zdata,model_params,[],[],cleaned_valid,create_reg);
disp("Pre-Processing Complete")   

% Bin neural data
% z_binned is trials x nbins x zdim
[z_binned, centres] = bin_neural_data_general_any(virmen_data,zdatafilt,xnum,linearise_x,nbins);

%% Bootstrapping
% Use bootstrapping to calculate distributions of sample means for each
% trial type in each bin. 
% Number of resamples defined by boot_samps
% bootstrap_means is num_types x boot_samps x bins x neurons

[bootstrap_means] = calc_bootsrapped_means_subsample(z_binned,tbt_details,boot_samps,t_types,sub_sample);

%% Calculate distances between mean curves
% Not used anymore

num_neurons = size(z_binned,3);
z_binned_means = zeros(size(z_binned,2),size(z_binned,3),length(t_types));
z_means_dists = zeros(size(z_binned,3),length(t_types),length(t_types));
  
% Mean tuning curves
for i = 1:length(t_types)
    z_binned_means(:,:,i) = mean(squeeze(bootstrap_means(i,:,:,:)),'omitnan');
end

% distances between mean tuning curves
for i = 1:length(t_types)
    for j = 1:length(t_types)
        z_means_dists(:,i,j) = sqrt(sum((z_binned_means(:,:,i)-z_binned_means(:,:,j)).^2));
    end
end

%% Determine significant changes in parameters

z_widths = zeros(length(t_types),boot_samps,size(z_binned,3));

% Define width as number of bins above significance threshold
% Define amplitude as difference between max and min (not just max value)

[max_z,z_peaks] = max(bootstrap_means,[],3);
max_z = squeeze(max_z);
z_peaks = squeeze(z_peaks);
[min_z] = min(bootstrap_means,[],3);
min_z = squeeze(min_z);
z_amps = max_z - min_z;

for i = 1:length(t_types)
    for b = 1:boot_samps
        for n = 1:size(z_binned,3)
            z_widths(i,b,n) = sum(squeeze(bootstrap_means(i,b,:,n))>(sig_check(n)));
        end
    end
end

% Store paramters for each bootstrap
all_params = zeros(size(z_peaks,1),size(z_peaks,2),size(z_peaks,3),3);
all_params(:,:,:,1) = z_peaks;
all_params(:,:,:,2) = z_amps;
all_params(:,:,:,3) = z_widths;

all_params_means = squeeze(mean(all_params,2,'omitnan'));

% Check if mean for BMI lies outside of confidence interval for control:
% [lower-upper, neurons, left-right, parameter]
% get CIs for ball trials
CIs_all = zeros(2,size(z_binned,3),2,3);
for i = 1:2
    for j = 1:3
        CIs_all(:,:,i,j) = prctile(squeeze(all_params(i,:,:,j)),CI_vals);
    end
end

full_results = zeros(size(z_binned,3),2,3);
for i = 1:2
    for j = 1:3
        % Check if BMI mean lies in CI for normal
        cur_means = squeeze(mean(squeeze(all_params(i+2,:,:,j))));
        cur_CIs = squeeze(CIs_all(:,:,i,j));
        full_results(:,i,j) = (cur_means < cur_CIs(1,:)) | (cur_means > cur_CIs(2,:));
    end
end

% Get CIs for BMI trials
CIs_all_b = zeros(2,size(z_binned,3),2,3);
for i = 1:2
    for j = 1:3
        CIs_all_b(:,:,i,j) = prctile(squeeze(all_params(i+2,:,:,j)),CI_vals);
    end
end

full_results_overlap = zeros(size(z_binned,3),2,3);
for i = 1:2
    for j = 1:3
        % Check if BMI CI lies outside CI for normal
        cur_CIs_b = squeeze(CIs_all_b(:,:,i,j));
        cur_CIs = squeeze(CIs_all(:,:,i,j));
        full_results_overlap(:,i,j) = (cur_CIs_b(2,:) < cur_CIs(1,:)) | (cur_CIs_b(1,:) > cur_CIs(2,:));
    end
end

% Check if ball mean outside CI for bmi
full_results_b = zeros(size(z_binned,3),2,3);

for i = 1:2
    for j = 1:3
        cur_means = squeeze(mean(squeeze(all_params(i,:,:,j))));
        cur_CIs = squeeze(CIs_all_b(:,:,i,j));
        full_results_b(:,i,j) = (cur_means < cur_CIs(1,:)) | (cur_means > cur_CIs(2,:));     
    end
end

% Check if both means outside other CI - this is what we used
full_results_both = full_results & full_results_b;

%% Get values for change in mean amplitude, width, and peak location

z_widths_mean = zeros(size(z_binned,3),length(t_types));

[max_z_means,z_peaks_mean] = max(z_binned_means);
max_z_means = squeeze(max_z_means);
z_peaks_mean = squeeze(z_peaks_mean);
[min_z_means] = min(z_binned_means);
min_z_means = squeeze(min_z_means);
z_amps_mean = max_z_means - min_z_means;

for i = 1:length(t_types)
    for n = 1:size(z_binned,3)
        z_widths_mean(n,i) = sum(squeeze(z_binned_means(:,n,i))>(sig_check(n)));
    end
end

mean_curve_params = zeros(size(z_peaks_mean,1),size(z_peaks_mean,2),3);
mean_curve_params(:,:,1) = z_peaks_mean;
mean_curve_params(:,:,2) = z_amps_mean;
mean_curve_params(:,:,3) = z_widths_mean;

mean_curve_params = permute(mean_curve_params,[2,1,3]);

if plot_res
    % optional plotting
    figure
    histogram(z_amps_mean(:,3)-z_amps_mean(:,1))
    figure
    histogram(z_peaks_mean(:,3)-z_peaks_mean(:,1))
    figure
    histogram(z_widths_mean(:,3)-z_widths_mean(:,1))
end

if plot_res
    figure
    for i = 1:2
        subplot(1,2,i)
        bar(squeeze(sum(full_results(:,i,:)))./size(z_binned,3))
        ylim([0,1])
    end

    figure
    for i = 1:2
        subplot(1,2,i)
        bar(squeeze(sum(full_results_overlap(:,i,:)))./size(z_binned,3))
        ylim([0,1])
    end
    
    figure
    for i = 1:2
        subplot(1,2,i)
        bar(squeeze(sum(full_results_both(:,i,:)))./size(z_binned,3))
        ylim([0,1])
    end
end
