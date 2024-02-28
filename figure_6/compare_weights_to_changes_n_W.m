function [va_corrs,yaw_corrs] = compare_weights_to_changes_n_W(n_Wout_online_cell,n_bmi_weights_cell,all_overall_n_means)
% 19/09/2023

% calculate correlation between decoder weights and changes in neural activity

% Use weights already extracted from neurons for specific variables 
% view angle for online, corrective angular velocity for offline
% Run for all subsamples/crossvalidations and average

num_mice = size(n_bmi_weights_cell,1);
num_days = size(n_bmi_weights_cell,2);
nmodels = 5; % Crossvalidations
nsubs = 30; % subsamples

% originally had alternative calculations, removed now, only use first
% value
va_corrs = nan.*ones(3,num_mice,num_days);
yaw_corrs = nan.*ones(3,nmodels*nsubs,num_mice,num_days);

for m = 1:num_mice

    for d = 1:num_days
        if ~isempty(n_bmi_weights_cell{m,d}) 
            
            n_weights_mean  = n_Wout_online_cell{m,d};
            
            cur_b = n_bmi_weights_cell{m,d};
            
            % Check correlation with neural changes
            cur_means = all_overall_n_means{m,d};
            cur_diffs = cur_means(2,:)-cur_means(1,:);
            
            % calculate correlation between absolute weights and changes
            va_corrs(1,m,d) = corr(abs(cur_diffs'),abs(n_weights_mean));
            ind = 0;
            for n = 1:nmodels
                for s = 1:nsubs
                    ind = ind +1;
                    n_weights_mean_n = squeeze(cur_b(s,n,:));
                    % Check correlation with neural changes
                    yaw_corrs(1,ind,m,d) = corr(abs(cur_diffs'),abs(n_weights_mean_n));
                end
            end
                
            
        end
    end
end