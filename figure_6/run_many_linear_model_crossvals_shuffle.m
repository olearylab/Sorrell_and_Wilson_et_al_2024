function [all_r2_cell_shuff,trial_avs_cell_shuff,trial_av_r2_cell_shuff,binned_r2_cell_shuff,trial_av_cd_r2_cell_shuff,trial_av_rmse_cell_shuff] = run_many_linear_model_crossvals_shuffle(z_cell,virmen_cell,tbt_cell,n_shuff,extra_norm,lambda_cell,zscore_x,nbins,shuff_limit)
% 20/11/2023

num_mice = size(z_cell,1);
num_days = size(z_cell,2);

initialise_params;
model_params.spatial = false;
model_params.reg_images = false;

all_r2_cell_shuff = cell(num_mice,num_days); 
trial_avs_cell_shuff = cell(num_mice,num_days); 
trial_av_r2_cell_shuff = cell(num_mice,num_days); 
binned_r2_cell_shuff = cell(num_mice,num_days); 
trial_av_cd_r2_cell_shuff = cell(num_mice,num_days); 
trial_av_rmse_cell_shuff = cell(num_mice,num_days);

for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(z_cell{m,d})
            disp("Running mouse " + m + " day " + d)
            [all_r2_cell_shuff{m,d},trial_avs_cell_shuff{m,d},trial_av_r2_cell_shuff{m,d},binned_r2_cell_shuff{m,d},trial_av_cd_r2_cell_shuff{m,d},trial_av_rmse_cell_shuff{m,d}] = linear_model_crossvals_shuffle(virmen_cell{m,d},z_cell{m,d},tbt_cell{m,d},model_params,n_shuff,extra_norm,lambda_cell{m,d},zscore_x,nbins,shuff_limit);
        end
    end
end
   
% save trial_av_rmse_cell_shuff.mat trial_av_rmse_cell_shuff