function [all_popboot_mean,all_neuron_means,all_overall_n_means] = run_many_mean_population_change(z_cell,virmen_cell,tbt_cell,model_params,linearise_x,nbins,boot_samps,xnum,normalise_z,type_groups)
% 09/09/2022

% Check mean population actvitiy change for many mice and days

num_mice = size(z_cell,1);
num_days = size(z_cell,2);

all_popboot_mean = cell(num_mice,num_days);
all_neuron_means = cell(num_mice,num_days);
all_overall_n_means = cell(num_mice,num_days);

for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(z_cell{m,d})
            
            [all_popboot_mean{m,d},all_neuron_means{m,d},all_overall_n_means{m,d}] = assess_mean_population_change_10092023(virmen_cell{m,d},z_cell{m,d},tbt_cell{m,d},model_params,linearise_x,nbins,boot_samps,xnum,normalise_z,type_groups);
        end
    end
end