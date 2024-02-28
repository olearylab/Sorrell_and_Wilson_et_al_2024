function [bootstrap_means_1,bootstrap_means_2,centres] = calc_bootstrap_from_raw_subsample_norm_halves(zdata, virmen_data, tbt_details, linearise_x, nbins,boot_samps, xnum,sub_sample,normalise_z)
% 03/03/2022

% Function for calculating bootstrapped mean tuning curves from raw neural
% data.
types_vec = [1,4,7,10];
max_types = max(tbt_details(1,:));
types_vec = types_vec(1:max_types);
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

%% Normalise (z-score) post filtering, before bootstrapping.
% Normalise filtered data with invalid data removed
if normalise_z 
    zdatafilt = zscore(zdatafilt);
end

% z_binned is trials x nbins x zdim
[z_binned, centres] = bin_neural_data_general_any(virmen_data,zdatafilt,xnum,linearise_x,nbins);
num_trials = size(z_binned,1);
%% Bootstrapping
% Use bootstrapping to calculate distributions of sample means for each
% trial type in each bin. 
% Number of resamples defined by boot_samps
% bootstrap_means is num_types x boot_samps x bins x neurons

[bootstrap_means_1,bootstrap_means_2] = calc_bootsrapped_means_subsample_halves(z_binned,tbt_details,boot_samps,types_vec,sub_sample);