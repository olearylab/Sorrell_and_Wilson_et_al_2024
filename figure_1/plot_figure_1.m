%% FIGURE 1 PLOTTING
addpath("../shared_functions")
% Prepare workspace for plotting figure 1

mice = ["DW81";"DW83";"DW86";"DW113";"DW129"];
data_days{1} = ["20200902";"20200907";"20200924"];
data_days{2} = ["20200902";"20200914";"20200922"];
data_days{3} = ["20200902";"20200907";"20200914"];
data_days{4} = ["20210623";"20210628";"20210630";"20210701"];
data_days{5} = ["20211111"];

example_m_d = "DW113_20210701";
reps = 5;
training_trials = 80;

results_struct = importdata("offline_results/res_struct_DW81_20200902_rep1_80.mat");
ex_weights = results_struct.Wout;
stat_cell_CNN_offline = importdata("stat_cell_CNN_offline.mat");

ex_stat = stat_cell_CNN_offline{1};

ex_im = importdata("mean_training_im_DW81_20200902.mat");

% Plot figure 1 c-n
[h_boots] = plot_offline_neuron_results_all_w_example_01122022(mice,data_days,example_m_d,reps,training_trials,ex_weights,ex_stat,ex_im);