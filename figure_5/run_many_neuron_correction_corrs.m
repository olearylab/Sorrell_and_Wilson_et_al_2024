function [z_correction_corrs_cell,z_correction_corrs_p_cell] = run_many_neuron_correction_corrs(z_binned_cell,x_binned_cell,error_cell,tbt_cell,error_thresh)
% 28/07/2023

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
            cur_x = squeeze(x_binned_cell{m,d}(cur_trials,:,3));
            cur_e = error_cell{m,d}(cur_trials,:);
            
            % Calculate only for bins where movements are heading
            % correcting and heading deviaitons are above threshold
            signs_opposite = sign(cur_e) ~= sign(cur_x);
            signs_opposite = signs_opposite(abs(cur_e)>error_thresh);
            cur_x = cur_x(abs(cur_e)>error_thresh);

            z_correction_corrs = nan.*ones(zdim,1);
            z_correction_corrs_p = nan.*ones(zdim,1);
            
            % Calculate correlations separately for each neuron
            for n = 1:zdim
                cur_zn = squeeze(cur_z(:,:,n));
                cur_zn = cur_zn(abs(cur_e)>error_thresh);
                [z_correction_corrs(n),z_correction_corrs_p(n)] = corr(cur_zn(signs_opposite),abs(cur_x(signs_opposite)));
            end
            z_correction_corrs_cell{m,d} = z_correction_corrs;
            z_correction_corrs_p_cell{m,d} = z_correction_corrs_p;
        end
    end
end