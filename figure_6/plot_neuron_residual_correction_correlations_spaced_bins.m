function [h_boot,ex_md,ro,ro_p,n_kept,stats_xc] = plot_neuron_residual_correction_correlations_spaced_bins(z_binned_res_cell,error_cell,tbt_cell,x_binned_cell,error_thresh,bin_space,num_shuff,shuff_limit)
% 27/10/2023

% Function for plotting correlation of residual population neural activity with error
% correcting movements. 

% Plotting for a single session. Stats calculated across all sessions with h-boot.

% Z should be normalised.

% Include circular shuffles for H-boot sig test. Can I do any shuffles?

% Discarding bins too close together (to hopefully reduce dependence
% between samples). 

num_mice = size(error_cell,1);
num_days = size(error_cell,2);
types_vec = [1,4,7,10];
colour_vec = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];

% store correlations
ro = NaN(num_mice,num_days);
ro_p = NaN(num_mice,num_days);
n_kept = NaN(num_mice,num_days);
ro_shuff = NaN(num_shuff,num_mice,num_days);
ro_circ_shuff = NaN(num_shuff,num_mice,num_days);

z_plot = cell(num_mice,num_days);
x_plot = cell(num_mice,num_days);

for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(error_cell{m,d})
            
            z_binned = z_binned_res_cell{m,d};
            error_mat = error_cell{m,d};
            tbt_details = tbt_cell{m,d};
            x_binned = x_binned_cell{m,d};
            
            % Just bmi trials
            cur_trials = ismember(tbt_details(3,:),types_vec([3,4]));
            cur_error = error_mat(cur_trials,:)';
            % Only use bins where heading deviations are above
            % threshold
            above_thresh = (abs(cur_error)>error_thresh);
            
            cur_z = z_binned(cur_trials,:,:);
            % Calculate population average
            pop_av = squeeze(mean(cur_z,3,'omitnan'))';
            cur_x = squeeze(x_binned(cur_trials,:,3))'; % get just ball angular velocity

            % Check for error correcting movements
            signs_opposite = (sign(cur_x) ~= sign(cur_error));
            
            % Create single vector for bins to keep
            kept_bins = above_thresh&signs_opposite;
            % Remove bins too close together in each trial
            for t = 1:size(kept_bins,2)
                
                cur_t = kept_bins(:,t);
                b_count = 0;
                for b = 1:length(cur_t)
                    % Biased towards the start
                    if b_count==0
                        cur_b = cur_t(b);
                        if cur_b
                            b_count = 1;
                        end
                    else
                        b_count = b_count+1;
                        kept_bins(b,t) = 0;
                        if b_count > bin_space
                            b_count = 0;
                        end
                    end
                    
                end
                
            end
            
            pop_av = pop_av(kept_bins);
            cur_x = abs(cur_x(kept_bins));
            
            % Store results
            [ro(m,d),ro_p(m,d)] = corr(cur_x,pop_av,'rows','pairwise');
            n_kept(m,d) = length(cur_x);
            z_plot{m,d} = pop_av;
            x_plot{m,d} = cur_x;
            
            % Shuffled correlations - circular
            num_samples = length(cur_x);
            shuff_inds = randi([shuff_limit,(num_samples-shuff_limit)],num_shuff,1);
            for n = 1:num_shuff
                
                shuff_z = circshift(pop_av,shuff_inds(n));
                
                ro_circ_shuff(n,m,d) = corr(cur_x,shuff_z,'rows','pairwise');
                
            end
            
            % Shuffled correlations - complete shuffle
            for n = 1:num_shuff
                shuffle_vec = randperm(length(pop_av));
                shuff_z = pop_av(shuffle_vec);
                
                ro_shuff(n,m,d) = corr(cur_x,shuff_z,'rows','pairwise');
                
            end
            
        end
    end
end

ex_md = [0,0];
[max_vec,max_inds] = max(ro);
[~,ex_md(2)] = max(max_vec);
ex_md(1) = max_inds(ex_md(2));

%% example z with error correcting movements
figure
x_all = x_plot{ex_md(1),ex_md(2)}.*(180/pi);
scatter(x_all,z_plot{ex_md(1),ex_md(2)},'o','MarkerEdgeColor',colour_vec(2,:),'MarkerFaceColor',colour_vec(2,:))
hold on

xlabel(["Ball angular velocity"; "magnitude during corrections (deg/s)"])
ylabel(["Mean BMI trial normalized"; "residual population activity"])
yline(0,'--','LineWidth',2);

mdl = fitlm(x_all,z_plot{ex_md(1),ex_md(2)});
xx = linspace(min(x_all),max(x_all),20)';
[ypred,yci] = predict(mdl,xx);
plot(xx,yci(:,1),'--','LineWidth',2,'Color','k')
plot(xx,yci(:,2),'--','LineWidth',2,'Color','k')
plot(xx,ypred,'-','LineWidth',2,'Color','k')
h = fill([xx',fliplr(xx')],[yci(:,2)',fliplr(yci(:,1)')],'k','LineStyle','none');
set(h,'facealpha',.1)

box off
axis('square')
title(["Residual population neural activity";"vs heading correcting movement"]);

%% H-boots

ro_shuff_ready = squeeze(mean(ro_shuff,'omitnan'));
ro_circ_shuff_ready = squeeze(mean(ro_circ_shuff,'omitnan'));

[h_boot.all_p_boot,h_boot.all_centres,h_boot.all_sems] = run_H_boot_ets(ro,ro_shuff_ready,true);
[h_boot.all_p_boot_circ,h_boot.all_centres_circ,h_boot.all_sems_circ] = run_H_boot_ets(ro,ro_circ_shuff_ready,true);

%% R^2 (variance explained)
mdl_R = NaN(num_mice,num_days);
mdl_Radj = NaN(num_mice,num_days);
R2 = NaN(num_mice,num_days);
var_ex = NaN(num_mice,num_days);
for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(error_cell{m,d})
            x_all = x_plot{m,d}.*(180/pi);
            z_all = z_plot{m,d};
            mdl = fitlm(x_all,z_all);
            mdl_R(m,d) = mdl.Rsquared.Ordinary;
            mdl_Radj(m,d) = mdl.Rsquared.Adjusted;
            [ypred,yci] = predict(mdl,x_all);
            R2(m,d) = 1 - sum((z_all - ypred).^2)/sum((z_all - mean(z_all)).^2);
            var_ex(m,d) = sum((ypred - mean(z_all)).^2)/sum((z_all - mean(z_all)).^2);
        end
    end
end
stats_xc.mdl_R = mdl_R;
stats_xc.mdl_Radj = mdl_Radj;
stats_xc.R2 = R2;
stats_xc.var_ex = var_ex;

[stats_xc.h_boot_p_boot1,stats_xc.h_boot_all_centres1,stats_xc.h_boot_all_sems1] = run_H_boot_ets(mdl_R,mdl_Radj,true);
[stats_xc.h_boot_p_boot2,stats_xc.h_boot_all_centres2,stats_xc.h_boot_all_sems2] = run_H_boot_ets(R2,var_ex,true);