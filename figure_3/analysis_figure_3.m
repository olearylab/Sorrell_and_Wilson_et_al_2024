%% Figure 3 analysis

% Run Behavioural SVM analysis
virmen_cell = importdata("../saved_files/virmen_cell_5.mat");
tbt_cell = importdata("../saved_files/tbt_cell_5.mat");
nbins = 50;
kept_types = [1,2,3,4];
balance_types = true;
keep_incorrect = false;
l_r_sep = true;
num_subs = 100;

% Run for many subsamples of trials
[accuracies_cell] = run_many_subs_svm_bvn_behaviour_classifier(virmen_cell,tbt_cell,nbins,kept_types,balance_types,keep_incorrect,l_r_sep,num_subs);

n_shuffles = 100;
% Run for many shuffles
[accuracies_cell_shuff] = run_many_svm_bvn_behaviour_classifier_shuffle(virmen_cell,tbt_cell,nbins,kept_types,balance_types,keep_incorrect,l_r_sep,n_shuffles);
