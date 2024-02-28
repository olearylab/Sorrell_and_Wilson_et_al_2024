%% Run figure 6 analysis

addpath("../shared_functions")
addpath("../figure_4")

% Run neural svms
z_cell_CNN = importdata("../saved_files/z_cell_CNN_5_v2_cut.mat");
virmen_cell = importdata("../saved_files/virmen_cell_5.mat");
tbt_cell = importdata("../saved_files/tbt_cell_5.mat");
kept_types = [1,2,3,4];
nbins = 50;
balance_types = true;
keep_incorrect = false;
l_r_sep = true;
num_subs = 100;
[accuracies_cell] = run_many_subs_svm_bvn_neural_classifier(z_cell_CNN,virmen_cell,tbt_cell,nbins,kept_types,balance_types,keep_incorrect,l_r_sep,num_subs);

% Run on shuffled trial labels
n_shuffles = 100;
[accuracies_cell_shuffle] = run_many_svm_bvn_neural_classifier_shuffle(z_cell_CNN,virmen_cell,tbt_cell,nbins,kept_types,balance_types,keep_incorrect,l_r_sep,n_shuffles);
save accuracies_cell_neural_shuffle_new.mat accuracies_cell_shuffle;

% Calculate mean neural activity
initialise_params;
linearise_x = true;
boot_samps = 100;
xnum = 6;
normalise_z = true;
type_groups = [1,4;7,10];
[all_popboot_mean,all_neuron_means,all_overall_n_means] = run_many_mean_population_change(z_cell_CNN,virmen_cell,tbt_cell,model_params,true,nbins,boot_samps,xnum,normalise_z,type_groups);

% Calculate significant tuning
shuff_limit = 300;
num_shuffles = 1000;
[sig_results_all,sig_check_all] = run_many_sig_tuning_check(z_cell_CNN, virmen_cell, tbt_cell, num_shuffles, linearise_x, nbins, boot_samps, shuff_limit, xnum);
% Calculate whether neurons gain or lost peaks 
[sig_summary_cell] = calc_peak_gain_loss(sig_results_all);

% Calculate changes in tuning
CI_vals = [2.5,97.5];
sub_sample = false;
[all_proportions,all_proportions_overlap,all_proportions_both,all_full_results,all_full_results_overlap,all_full_results_both,all_mean_curve_params,full_all_params] = assess_tc_param_changes_many_21022023(z_cell_CNN,virmen_cell,tbt_cell,linearise_x,nbins,boot_samps,CI_vals,xnum,sig_check_all,sub_sample);

% Calculate binned residual activity
normalise_z = true;
[z_binned_res_cell] = get_z_binned_res_cell(z_cell_CNN,virmen_cell,tbt_cell,model_params,normalise_z,nbins);

error_thresh = 1;
virmen_train_cell = importdata("../saved_files/virmen_train_cell_5.mat");
tbt_train_cell = importdata("../saved_files/tbt_train_cell_5.mat");
% bin heading deviations and behvaiour
[error_cell,x_binned_cell,ro_all,sign_opposite,h_boots] = run_many_check_BMI_correct_errors_norm(virmen_cell,tbt_cell,virmen_train_cell,tbt_train_cell,nbins,error_thresh,false);

% Project pixel fluorescence through view angle decoder and calculate
% differentiated output
% Results for example sesssion given here, need to run for all sessions
zfull = importdata("../saved_files/DW81_20200928_test.mat");
Wout_online_cell = importdata("../saved_files/Wout_online_cell.mat");
Wout = Wout_online_cell{ex_md(1)};
Wout = Wout(:,2);
train_mean_cell = importdata("../saved_files/train_mean_cell.mat");
Rfixed_cell = importdata("../saved_files/Rfixed_cell.mat");
ex_md = [1,1];
normalise_z = false;
use_bias = false;
% Output is stored in x_binned
[z_binned,z_binned_res,x_binned,all_means] = calc_binned_residual_activity_from_samps_pix(zfull,virmen_cell{ex_md(1),ex_md(2)},tbt_cell{ex_md(1),ex_md(2)},Wout,model_params,normalise_z,nbins,train_mean_cell{ex_md(1)},Rfixed_cell{ex_md(1)},use_bias);

% Calculate noise correlations
sub_sample = false;
[full_pairs_means_cell] = run_many_noise_correlations_04072023(z_cell_CNN,virmen_cell,tbt_cell, model_params, nbins, linearise_x, sub_sample, normalise_z);

% Calculate correlations of decoder weights and mean activity changes
stat_cell = importdata("../saved_files/stat_cell_CNN_v2_cut.mat");
bmi_weights_cell = importdata("../saved_files/bmi_weights_cell_corrective_pixels_full.mat");
is_offline = false;
w_ind = 2;
% Extract weights within neurons for online view angle weights
% Commented out as this takes a long time
% [n_Wout_online_cell] = get_all_n_weights(Wout_online_cell,stat_cell,is_offline,w_ind);
w_ind = 4;
% Extract weights within neurons for online view angle weights
% Commented out as this takes a long time
% [n_bmi_weights_cell] = get_all_n_weights_subbed(bmi_weights_cell,stat_cell,w_ind);

% Load in the output from above to avoid running again
n_Wout_online_cell = importdata("n_W_online_cell_cut.mat");
n_bmi_weights_cell = importdata("n_W_cell_corrective_subs_cut.mat");

[va_corrs,yaw_corrs] = compare_weights_to_changes_n_W(n_Wout_online_cell,n_bmi_weights_cell,all_overall_n_means);
