function [population_bootstrap_mean,neuron_means,overall_n_means] = assess_mean_population_change_10092023(virmen_data,zdata,tbt_details,model_params,linearise_x,nbins,boot_samps,xnum,normalise_z,type_groups)
% 10/09/2023

% Assess mean change in population activity in bins along trial. Is there
% overall increase (or decrease) and where?
% Assess overall mean change (single number)

% normalising just with ball trials.

% Combine left and right together

%%
types_vec = [1,4,7,10];

%% remove invalid data
ITI = virmen_data(8,:);
cleaned_valid = clean_valid_data(ITI);

virmen_data = virmen_data(:,cleaned_valid);
%% Preprocess data

if ndims(zdata) == 2
    model_params.spatial = false;
    create_reg = false;
else
    create_reg = true;
end
model_params.reg_images = false;

[zdatafilt, train_mean, Rfixed] = preprocess_images_05_2021(zdata,model_params,[],[],cleaned_valid,create_reg);

%% Normalise (z-score) post filtering, before bootstrapping.
% Normalise filtered data with invalid data removed
if normalise_z 
    % normalise according to ball trials
    [zdatafilt] = normalise_z_by_ball(zdatafilt,virmen_data,tbt_details);
end

%% Calculate overall neuron means
% mean value for each neuron,
% averaged across whole population for each trial type.

overall_n_means = zeros(size(type_groups,1),size(zdatafilt,2));

trial_num = virmen_data(12,:);

for t = 1:size(type_groups,1)
    cur_trials = find(ismember(tbt_details(3,:),type_groups(t,:)));
    
    overall_n_means(t,:) = mean(zdatafilt(ismember(trial_num,cur_trials),:),'omitnan');
end

%% Bin neural data
% z_binned is trials x nbins x zdim
[z_binned, centres] = bin_neural_data_general_any(virmen_data,zdatafilt,xnum,linearise_x,nbins);

z_binned = z_binned(~isnan(z_binned(:,1,1)),:,:);
%% Bootstrapping
% Use bootstrapping to calculate distributions of sample means for each
% trial type in each bin. 
% Number of resamples defined by boot_samps
% bootstrap_means is num_types x boot_samps x bins x neurons

[bootstrap_means] = calc_bootsrapped_means_subsample(z_binned,tbt_details,boot_samps,types_vec,false);

mean_bootstrap = squeeze(mean(bootstrap_means,2,'omitnan'));

% population boostrap mean is num_types x bins
population_bootstrap_mean = squeeze(mean(mean_bootstrap,3,'omitnan'));

final_mean = mean(population_bootstrap_mean,2,'omitnan');

% Possibly compare individual neurons means
neuron_means = squeeze(mean(mean_bootstrap,2,'omitnan'));