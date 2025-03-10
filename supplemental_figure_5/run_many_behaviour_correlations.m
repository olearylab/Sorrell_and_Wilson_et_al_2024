function [behav_corr_cell,behav_std_corr_cell,behav_full_corr_cell] = run_many_behaviour_correlations(virmen_cell,tbt_cell,offsets)
% 25/07/2023

% Calculate beahviour correlations within and across trial types
sub_sample = false;

num_mice = size(virmen_cell,1);
num_days = size(virmen_cell,2);

behav_corr_cell = cell(num_mice,num_days,2);
behav_full_corr_cell = cell(num_mice,num_days,2);
behav_std_corr_cell = cell(num_mice,num_days,2);


for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(virmen_cell{m,d})
            
            % Run for forward velocity
            [behav_corr_cell{m,d,1},behav_std_corr_cell{m,d,1},behav_full_corr_cell{m,d,1}] = calc_behaviour_correlations(virmen_cell{m,d},tbt_cell{m,d},offsets(m,:),sub_sample,13);
            
            % Run for angular velocity
            [behav_corr_cell{m,d,2},behav_std_corr_cell{m,d,2},behav_full_corr_cell{m,d,2}] = calc_behaviour_correlations(virmen_cell{m,d},tbt_cell{m,d},offsets(m,:),sub_sample,15);            
            
        end
    end
end