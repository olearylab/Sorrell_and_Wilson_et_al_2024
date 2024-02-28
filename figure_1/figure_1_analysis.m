%% Figure 1 analysis
addpath("../shared_functions")
%% Pixel-based decoding
% For each session load:

% Images, should be (samples x 128 x 128) e.g.
zfull = importdata("DW81_20200902_final.mat");

% virmen data, e.g.
xfull = importdata("trimmed_virmen_DW81_20200902_test_revised.mat");

% trial type and outcome information, e.g.
tbt_details = importdata("DW81_20200902_trial_by_trial_details.mat");

% Then run
reps = 5;
training_trials = 80;
intialise_params;

[results_struct] = run_save_offline_single(xfull, zfull, tbt_details, model_params, reps, training_trials);

%% Neuron-based decoding

% For each session load:

% neuron data, should be (samples x neurons) e.g.
zfull = importdata("DW81_20200902_neurons.mat");

% virmen data, e.g.
xfull = importdata("trimmed_virmen_DW81_20200902_test_revised.mat");

% trial type and outcome information, e.g.
tbt_details = importdata("DW81_20200902_trial_by_trial_details.mat");

% This can be run for each sessions using:
reps = 5;
training_trials = 80;
intialise_params;
model_params.spatial = false;
model_params.reg_images = false;
model_params.lrate = 10^-2;

[results_struct_n] = run_save_offline_single(xfull, zfull, tbt_details, model_params, reps, training_trials);

