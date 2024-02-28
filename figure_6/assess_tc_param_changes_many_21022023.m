function [all_proportions,all_proportions_overlap,all_proportions_both,all_full_results,all_full_results_overlap,all_full_results_both,all_mean_curve_params,full_all_params] = assess_tc_param_changes_many_21022023(z_cell,virmen_cell,tbt_cell,linearise_x,nbins,boot_samps,CI_vals,xnum,sig_check_all,sub_sample)
% 21/02/2023

% Function for assessing changes of parameters of tuning curves for many
% mice and many days

num_mice = size(z_cell,1);
num_days = size(z_cell,2);

%%
all_proportions = zeros(num_mice,num_days,2,3);
all_proportions_overlap = zeros(num_mice,num_days,2,3);
all_proportions_both = zeros(num_mice,num_days,2,3);
all_full_results = cell(num_mice,num_days);
all_full_results_overlap = cell(num_mice,num_days);
all_full_results_both = cell(num_mice,num_days);

all_mean_curve_params = cell(num_mice,num_days);
full_all_params = cell(num_mice,num_days);

for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(z_cell{m,d})

            [full_results,full_results_overlap,full_results_both,all_params,mean_curve_params] = assess_tuning_curve_param_changes_alt_20022023(virmen_cell{m,d},z_cell{m,d},tbt_cell{m,d},linearise_x,nbins,boot_samps,CI_vals,xnum,false,sig_check_all{m,d},sub_sample);
            all_full_results{m,d} = full_results;
            all_full_results_overlap{m,d} = full_results_overlap;
            all_full_results_both{m,d} = full_results_both;
            all_proportions(m,d,:,:) = squeeze(sum(full_results))./size(full_results,1);
            all_proportions_overlap(m,d,:,:) = squeeze(sum(full_results_overlap))./size(full_results_overlap,1);
            all_proportions_both(m,d,:,:) = squeeze(sum(full_results_both))./size(full_results_both,1);
            
            % added saving of mean curve parameters
            all_mean_curve_params{m,d} = mean_curve_params;
            full_all_params{m,d} = all_params;
        end
    end
end