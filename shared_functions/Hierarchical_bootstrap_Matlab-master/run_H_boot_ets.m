function [p_boot,h_centre,h_sem] = run_H_boot_ets(data1, data2, reorder_nans)
% 23/08/2023

% Calculate bootstrap means and sem
% Then calculate stats test for a-b. Always select b to have larger mean,
% then add up proportion greater than or equal to 0.

% Set random seed for consistency
rng(1);

if reorder_nans
    num_mice = size(data1,1);
    num_days = size(data2,2);
    % In case nan's are not in correct place
    orig_data1 = data1;
    orig_data2 = data2;
    for m = 1:num_mice
        d_ind = 0;
        for d = 1:num_days
            if ~isnan(orig_data1(m,d)) 
                d_ind = d_ind+1;
                data1(m,d_ind) = orig_data1(m,d);
                data2(m,d_ind) = orig_data2(m,d);
            end
        end
        if d_ind<num_days
            data1(m,d_ind+1:num_days) = nan;
            data2(m,d_ind+1:num_days) = nan;
        end
    end
end

% Hard coded for now, I don't think these will be changed
boot_samps = 1000;
num_trials = 4;

% Calculates probability that data with lower mean is >= data with higher
% mean. Therefore low probability is significant for lower<higher. 
% Get means and sems
[p_b, bootstats, h_centre, h_sem] = get_bootstrap_results_equalsamples(data1,data2,boot_samps,num_trials,'mean');

if h_centre(1) > h_centre(2)
    data_test = data2 - data1;
else
    data_test = data1 - data2;
end

[bootstats] = get_bootstrapped_equalsamples(data_test,boot_samps,num_trials,'mean');

p_boot = sum(bootstats>=0)/length(bootstats);

