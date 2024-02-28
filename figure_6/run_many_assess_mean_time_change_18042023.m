function [] = run_many_assess_mean_time_change_18042023(virmen_cell,z_cell,tbt_cell,model_params,av_num)
% 18/04/2023

% Run mean assessment througout sessions for many sessions

num_mice = size(virmen_cell,1);
num_days = size(virmen_cell,2);

sub_sample = false;
normalise_z = true;
plot_res = false;

means_cell = cell(num_mice,num_days);
trials_cell = cell(num_mice,num_days);

types_vec = [1,4,7,10];

for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(z_cell{m,d})
            
            [overall_trial_means,kept_trials] = assess_mean_time_change_18042023(virmen_cell{m,d},z_cell{m,d},tbt_cell{m,d},model_params,sub_sample,normalise_z,plot_res,av_num);
            
            means_cell{m,d} = overall_trial_means;
            trials_cell{m,d} = kept_trials;
        end
    end
end

%% Plotting

colour_vec = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];

figure
all_means = cell(2,1);
all_x = cell(2,1);
for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(means_cell{m,d})
            tbt_details = tbt_cell{m,d};
            overall_trial_means = means_cell{m,d};
            kept_trials = trials_cell{m,d};
            
            kept_tbt = tbt_details(:,kept_trials);
            % running average
            for t = 1:2
                cur_t = ismember(kept_tbt(3,:),types_vec([2*(t-1)+1,2*(t-1)+2]));
                cur_kept = kept_trials(cur_t);
                cur_means = movmean(overall_trial_means(cur_t),av_num);
                all_means{t} = [all_means{t};cur_means];
                all_x{t} = [all_x{t},cur_kept];
                plot(cur_kept,cur_means,'LineWidth',1,'Color',colour_vec(t,:))
                hold on
            end
            
        end
    end
end
for t = 1:2
    P = polyfit(all_x{t}, all_means{t}, 1);
    Pfit = polyval(P, sort(all_x{t}));
    plot(sort(all_x{t}),Pfit,'LineWidth',3,'Color',colour_vec(t,:))
end
xlabel("Trial Number")
ylabel("Mean Normalised Activity (a.u.)")
title("Mean Activity Over Time")
legend("Ball","BMI")
box off
axis('square')