function [full_corr,binned_corr_pairs,binned_corrs,full_pairs_means] = calc_noise_correlations_04072023(zdata,virmen_data,tbt_details, model_params, nbins, linearise_x, sub_sample, normalise_z)
% 11/10/2022

% Calculate noise correlations 
% Correlations between pairs of neurons within each position bin across trials
zdim = size(zdata,2);
types_vec = [1,4,7,10];
xnum = 6;
num_types = length(types_vec);
%% remove invalid data
ITI = virmen_data(8,:);
cleaned_valid = clean_valid_data(ITI);

virmen_data = virmen_data(:,cleaned_valid);
%% Preprocess and bin data

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
% option to subsample
if sub_sample
    [kept_trials] = subsample_trials(tbt_details,types_vec);
    kept_trials = kept_trials(:);
else
    kept_trials = 1:max(trial_num);
end

kept_trial_num = ismember(trial_num,kept_trials);

zdatafilt = zdatafilt(kept_trial_num,:);
virmen_data = virmen_data(:,kept_trial_num);

%% Normalise (z-score) post filtering, before bootstrapping.
% Normalise filtered data with invalid data removed
if normalise_z 
    % Normalise according to neural data on ball trials
    [zdatafilt] = normalise_z_by_ball(zdatafilt,virmen_data,tbt_details);
end

% Bin neural data
% z_binned is trials x nbins x zdim
[z_binned, centres] = bin_neural_data_general_any(virmen_data,zdatafilt,xnum,linearise_x,nbins);
z_binned = z_binned(~isnan(z_binned(:,1,1)),:,:);
tbt_details = tbt_details(:,sort(kept_trials));

% store indices for correct trials of each type
full_trials = cell(4,1);
types_vec = [1,4,7,10];
for i = 1:4
    full_trials{i} = find(ismember(tbt_details(3,:),types_vec(i)));
end

% fast method for computing correlations. Set up neural data for this
zbt = permute(z_binned,[1,3,2]);
zb_ext = zbt(:,:);

% extract all combinations of trial type correlations
full_corr = nan*ones(nbins,nbins,4);
num_pairs = (zdim*zdim - zdim)/2;
binned_corr_pairs = nan*ones(nbins,4,num_pairs);
full_pairs_means = nan.*ones(4,zdim,zdim);

for k = 1:4
    cur_coeff = corrcoef(zb_ext(full_trials{k},:),'rows','pairwise');
    % Only care about same bin correlations
    cur_pairs = nan.*ones(nbins,zdim,zdim);
    for i = 1:nbins

        cur_block = cur_coeff((i-1)*zdim + 1:i*zdim,(i-1)*zdim + 1:i*zdim);
        
        cur_pairs(i,:,:) = cur_block;
        
        % Nan out lower triangle which is are repeats and the diagonal
        low_t = tril(nan.*ones(zdim));
        cur_block = cur_block + low_t;
        cur_block_ready = cur_block(:);
        % Averaged across neuron pairs
        full_corr(i,i,k) = mean(cur_block_ready(~isnan(cur_block_ready)),'omitnan');
        binned_corr_pairs(i,k,:) = cur_block_ready(~isnan(cur_block_ready));
        
    end
    full_pairs_means(k,:,:) = mean(cur_pairs,'omitnan');
end

% Extract the correlations we care about
binned_corrs = zeros(4,nbins);
for i = 1:4
    binned_corrs(i,:) = diag(squeeze(full_corr(:,:,i)));
end