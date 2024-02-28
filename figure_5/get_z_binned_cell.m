function [z_binned_cell] = get_z_binned_cell(z_cell,virmen_cell,tbt_cell,model_params,normalise_z,nbins)
% 03/04/2023

% Create cell containing each sessions neural activity binned according to
% linearised position.
% 
num_mice = size(z_cell,1);
num_days = size(z_cell,2);

z_binned_cell = cell(num_mice,num_days);

for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(z_cell{m,d})
            
            [z_binned] = get_norm_z_binned(z_cell{m,d},virmen_cell{m,d},tbt_cell{m,d},model_params,normalise_z,nbins);

            z_binned_cell{m,d} = z_binned;
                
        end
    end
end
