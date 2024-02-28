function [boot_stats,sig_results,sig_results_train,sig_results_tot] = performance_stats(trial_types_summary_cell)
% 10/08/2023

num_days = 5;

summary_results = nan.*ones(length(trial_types_summary_cell),num_days,2);
tot_b = nan.*ones(length(trial_types_summary_cell),num_days,2);
tot_n = nan.*ones(length(trial_types_summary_cell),num_days,2);
for m = 1:length(trial_types_summary_cell)
    cur_res = trial_types_summary_cell{m};
    n_days = size(cur_res,2);
    for d = 1:n_days

        summary_results(m,d,1) = sum(cur_res([1,4],d),'all')/sum(cur_res(1:6,d),'all');
        summary_results(m,d,2) = sum(cur_res([7,10],d),'all')/sum(cur_res(7:end,d),'all');

        tot_b(m,d,1) = sum(cur_res([7,10],d),'all');
        tot_b(m,d,2) = sum(cur_res(7:end,d),'all');
        tot_n(m,d,1) = sum(cur_res([1,4],d),'all');
        tot_n(m,d,2) = sum(cur_res(1:6,d),'all');
    end
end

sig_results = nan.*ones(length(trial_types_summary_cell),num_days-1,2);
sig_results_train = nan.*ones(length(trial_types_summary_cell),1);
% summed significance check
tot_b_all = squeeze(sum(tot_b,2,'omitnan'));
tot_n_all = squeeze(sum(tot_n,2,'omitnan'));

sig_results_tot = nan.*ones(length(trial_types_summary_cell),2);

for m = 1:length(trial_types_summary_cell)
    for d = 1:num_days
        if d ==1
            sig_results_train(m) = myBinomTest(tot_n(m,d,1),tot_n(m,d,2),0.5,'one');
        else
            % Each day separate
            sig_results(m,d-1,1) = myBinomTest(tot_n(m,d,1),tot_n(m,d,2),0.5,'one');
            sig_results(m,d-1,2) = myBinomTest(tot_b(m,d,1),tot_b(m,d,2),0.5,'one');
        end
    end
    % Using total across all sessions for each mouse
    sig_results_tot(m,1) = myBinomTest(tot_n_all(m,1),tot_n_all(m,2),0.5,'one');
    sig_results_tot(m,2) = myBinomTest(tot_b_all(m,1),tot_b_all(m,2),0.5,'one');
end

%% Hierarchical bootstrap
% nans already is correct place
% traing only has one day.
rng(1)
boot_samps = 1000;
num_trials = 1;

bootstats = get_bootstrapped_equalsamples(squeeze(summary_results(:,1,1)),boot_samps,num_trials,'mean');

boot_stats.train_center = mean(bootstats);
boot_stats.train_std = std(bootstats);

num_trials = 4;
bootstats = get_bootstrapped_equalsamples(squeeze(summary_results(:,2:end,1)),boot_samps,num_trials,'mean');

boot_stats.ball_center = mean(bootstats);
boot_stats.ball_std = std(bootstats);

bootstats = get_bootstrapped_equalsamples(squeeze(summary_results(:,2:end,2)),boot_samps,num_trials,'mean');

boot_stats.bmi_center = mean(bootstats);
boot_stats.bmi_std = std(bootstats);