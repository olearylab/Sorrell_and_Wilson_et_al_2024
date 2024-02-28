%% Run analysis for figure 5

addpath("../shared_functions")
addpath("../figure_4")

% Set up workspace
virmen_cell = importdata("../saved_files/virmen_cell_5.mat");
tbt_cell = importdata("../saved_files/tbt_cell_5.mat");
z_cell = importdata("../saved_files/z_cell_CNN_5_v2_cut.mat");

initialise_params;
normalise_z = true;
nbins = 50;
% bin neural data
[z_binned_cell] = get_z_binned_cell(z_cell,virmen_cell,tbt_cell,model_params,normalise_z,nbins);

error_thresh = 1;
virmen_train_cell = importdata("../saved_files/virmen_train_cell_5.mat");
tbt_train_cell = importdata("../saved_files/tbt_train_cell_5.mat");
% bin heading deviations and behvaiour
[error_cell,x_binned_cell,ro_all,sign_opposite,h_boots] = run_many_check_BMI_correct_errors_norm(virmen_cell,tbt_cell,virmen_train_cell,tbt_train_cell,nbins,error_thresh,false);

% Run neuron specific correlations with corrective movements
bin_space = 5;
[z_correction_corrs_cell,z_correction_corrs_p_cell] = run_many_neuron_correction_corrs_spaced_bins(z_binned_cell,x_binned_cell,error_cell,tbt_cell,error_thresh,bin_space);

%%
% Run decoding of angular velocity (and other variables) during corrective
% movements. Example here for single session. Run for all sessions.
train_mean_cell = importdata("../saved_files/train_mean_cell.mat");
Rfixed_cell = importdata("../saved_files/Rfixed_cell.mat");
zfull = importdata("../saved_files/DW81_20200928_test.mat");

ex_md = [1,1];
offsets = [1.4953,1.4961;1.4953,1.4961;1.4953,1.4961;1.4804,1.4844;1.4804,1.4844];
num_subs = 30;
initialise_params;
model_params.loops = 20;

% Get mean heading in training ball trials for heading deviation
% calculation
[mean_binned,std_binned] = calculate_mean_ball_va(virmen_train_cell{ex_md(1)},tbt_train_cell{ex_md(1)},nbins,offsets(ex_md(1),:));
% Run cross-validated decoding for many subsamples
[all_res_mat,bmi_weights,num_samps] = decoding_check_corrective(virmen_cell{ex_md(1),ex_md(2)},zfull,tbt_cell{ex_md(1),ex_md(2)},model_params,train_mean_cell{ex_md(1)},Rfixed_cell{ex_md(1)},nbins,offsets(ex_md(1),:),mean_binned,std_binned,num_subs);



