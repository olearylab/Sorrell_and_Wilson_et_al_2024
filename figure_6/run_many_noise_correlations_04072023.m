function [full_pairs_means_cell] = run_many_noise_correlations_04072023(z_cell,virmen_cell,tbt_cell, model_params, nbins, linearise_x, sub_sample, normalise_z)
% 28/03/2023

% Calculate noise correlations for all sessions

num_mice = size(z_cell,1);
num_days = size(z_cell,2);

full_pairs_means_cell = cell(num_mice,num_days);

for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(z_cell{m,d})
            [full_corr,binned_corr_pairs,binned_corrs,full_pairs_means] = calc_noise_correlations_04072023(z_cell{m,d},virmen_cell{m,d},tbt_cell{m,d}, model_params, nbins, linearise_x, sub_sample, normalise_z);
            full_pairs_means_cell{m,d} = full_pairs_means;
        end
    end
end
