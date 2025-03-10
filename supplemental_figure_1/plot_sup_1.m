%% plot supplement figure 1
addpath("../shared_functions")

results_struct = importdata("../figure_1/offline_results/res_struct_DW81_20200902_rep1_80.mat");
xfull = importdata("sup_xtest.mat");
tbt_details = importdata("sup_tbt_details.mat");

offsets = [1.4953,1.4961;1.4953,1.4961;1.4953,1.4961;1.4804,1.4844;1.4804,1.4844];
plot_trial = 290;

% Example supplemental figure 1 a
plot_example_pos_dis(results_struct,xfull,plot_trial)

% Plot supplemental figure 1 b
plot_example_error_accumulation(results_struct,xfull,tbt_details,[1,2],7,4,offsets(1,2),plot_trial);