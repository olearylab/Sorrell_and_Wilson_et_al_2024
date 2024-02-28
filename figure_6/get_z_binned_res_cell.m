function [z_binned_res_cell] = get_z_binned_res_cell(z_cell,virmen_cell,tbt_cell,model_params,normalise_z,nbins)
% 12/09/2023

% Bin residual (BMI - ball) neural activity for each session.

% 
num_mice = size(z_cell,1);
num_days = size(z_cell,2);

z_binned_res_cell = cell(num_mice,num_days);

for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(z_cell{m,d})
            
            [z_binned] = calc_binned_residual_activity(z_cell{m,d},virmen_cell{m,d},tbt_cell{m,d},model_params,normalise_z,nbins);

            z_binned_res_cell{m,d} = z_binned;
                
        end
    end
end
