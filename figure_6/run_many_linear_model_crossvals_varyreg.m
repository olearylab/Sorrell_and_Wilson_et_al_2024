function [all_r2_cell,trial_avs_cell,trial_av_r2_cell,binned_r2_cell,trial_av_cd_r2_cell,trial_av_rmse_cell] = run_many_linear_model_crossvals_varyreg(z_cell,virmen_cell,tbt_cell,num_subs,extra_norm,lambda_cell,zscore_x,nbins)
% 20/11/2023

num_mice = size(z_cell,1);
num_days = size(z_cell,2);

initialise_params;
model_params.spatial = false;
model_params.reg_images = false;

all_r2_cell = cell(num_mice,num_days); 
trial_avs_cell = cell(num_mice,num_days); 
trial_av_r2_cell = cell(num_mice,num_days); 
binned_r2_cell = cell(num_mice,num_days); 
trial_av_cd_r2_cell = cell(num_mice,num_days); 
trial_av_rmse_cell = cell(num_mice,num_days);

for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(z_cell{m,d})
            disp("Running mouse " + m + " day " + d)
            [all_r2_cell{m,d},trial_avs_cell{m,d},trial_av_r2_cell{m,d},binned_r2_cell{m,d},trial_av_cd_r2_cell{m,d},trial_av_rmse_cell{m,d}] = linear_model_crossvals_varyreg(virmen_cell{m,d},z_cell{m,d},tbt_cell{m,d},model_params,num_subs,extra_norm,lambda_cell{m,d},zscore_x,nbins);
        end
    end
end
    