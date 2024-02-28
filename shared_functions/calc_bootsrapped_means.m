function [bootstrap_means] = calc_bootsrapped_means(z_binned,tbt_details,boot_samps,types_vec)
% 07/03/2022

% Function for calculating bootstrapped means of binned neural data.

% resample trials for each trial type. Whole trials, same for all
% neurons sampled, rather than each bin and/or each neuron independently.

num_types = length(types_vec); % NL, NR, BL, BR. Could change
bootstrap_means = zeros(num_types,boot_samps,size(z_binned,2),size(z_binned,3));

for n = 1:num_types
    cur_trials = find(ismember(tbt_details(3,:),types_vec(n)));
    num_samples = length(cur_trials);
    for b = 1:boot_samps
        cur_samples = datasample(cur_trials,num_samples);
        cur_resamp = z_binned(cur_samples,:,:);
        bootstrap_means(n,b,:,:) = mean(cur_resamp,'omitnan');
    end
end
    