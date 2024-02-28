function [z_correction_corrs_cell,z_correction_corrs_p_cell] = run_many_neuron_correction_corrs_spaced_bins(z_binned_cell,x_binned_cell,error_cell,tbt_cell,error_thresh,bin_space)
% 27/10/2023

% Calculate individual neuron correlations with corrective movements, with
% spacing between bins kept for correlation.

% 
num_mice = size(z_binned_cell,1);
num_days = size(z_binned_cell,2);

types_vec = [1,4,7,10];

z_correction_corrs_cell = cell(num_mice,num_days);
z_correction_corrs_p_cell = cell(num_mice,num_days);

%% 
for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(z_binned_cell{m,d})
            
            
            cur_tbt = tbt_cell{m,d};
            
            % Calculate only for BMI trials
            cur_trials = find(ismember(cur_tbt(3,:),types_vec([3,4])));
            
            cur_z = z_binned_cell{m,d}(cur_trials,:,:);
            zdim = size(cur_z,3);
            cur_x = squeeze(x_binned_cell{m,d}(cur_trials,:,3))';
            cur_e = error_cell{m,d}(cur_trials,:)';
            
            % Calculate only for bins where movements are heading
            % correcting and heading deviaitons are above threshold
            signs_opposite = sign(cur_e) ~= sign(cur_x);
            above_thresh = abs(cur_e)>error_thresh;
            
            % Create single vector for bins to keep
            kept_bins = above_thresh&signs_opposite;
            % Remove bins too close together in each trial
            for t = 1:size(kept_bins,2)
                
                cur_t = kept_bins(:,t);
                b_count = 0;
                for b = 1:length(cur_t)
                    % Biased towards the start
                    if b_count==0
                        cur_b = cur_t(b);
                        if cur_b
                            b_count = 1;
                        end
                    else
                        b_count = b_count+1;
                        kept_bins(b,t) = 0;
                        if b_count > bin_space
                            b_count = 0;
                        end
                    end
                    
                end
                
            end
            cur_x = cur_x(kept_bins);

            z_correction_corrs = nan.*ones(zdim,1);
            z_correction_corrs_p = nan.*ones(zdim,1);
            
            % Calculate correlations separately for each neuron
            for n = 1:zdim
                cur_zn = squeeze(cur_z(:,:,n))';
                cur_zn = cur_zn(kept_bins);
                [z_correction_corrs(n),z_correction_corrs_p(n)] = corr(cur_zn,abs(cur_x));
            end
            z_correction_corrs_cell{m,d} = z_correction_corrs;
            z_correction_corrs_p_cell{m,d} = z_correction_corrs_p;
        end
    end
end