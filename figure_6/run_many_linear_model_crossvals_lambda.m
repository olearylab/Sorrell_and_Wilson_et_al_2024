function [lambda_cell] = run_many_linear_model_crossvals_lambda(z_cell,virmen_cell,tbt_cell,extra_norm,lambda_vec,zscore_x,nbins,accuracy_type)
% 20/11/2023

num_mice = size(z_cell,1);
num_days = size(z_cell,2);

initialise_params;
model_params.spatial = false;
model_params.reg_images = false;

lambda_cell = cell(num_mice,num_days); 

for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(z_cell{m,d})
            disp("Running mouse " + m + " day " + d)
            [lambda_r2_mat,lambda_mat] = linear_model_crossvals_lambda(virmen_cell{m,d},z_cell{m,d},tbt_cell{m,d},model_params,extra_norm,lambda_vec,zscore_x,nbins,accuracy_type);
            lambda_cell{m,d} = lambda_mat;
        end
    end
end
  