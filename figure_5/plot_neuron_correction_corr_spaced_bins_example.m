function [n_ind,n_samps] = plot_neuron_correction_corr_spaced_bins_example(z_binned_cell,x_binned_cell,error_cell,tbt_cell,z_correction_corrs_cell,ex_md,error_thresh,bin_space)
% 27/10/2023

% plot scatter of correlation for a single neuron between activity and
% error correcting behaviour. - using new spaced bins method

% Only for bins with heading deviations above error_thresh, and where
% movements are error correcting. Only keep bins separated by bin_space

types_vec = [1,4,7,10];

z_binned = z_binned_cell{ex_md(1),ex_md(2)};
error_mat = error_cell{ex_md(1),ex_md(2)};
tbt_details = tbt_cell{ex_md(1),ex_md(2)};
z_correction_corrs = z_correction_corrs_cell{ex_md(1),ex_md(2)};
x_binned = x_binned_cell{ex_md(1),ex_md(2)};

% get neuron with maximum correlation
[n_corr,n_ind] = max(z_correction_corrs);

% get bmi trials
cur_trials = ismember(tbt_details(3,:),types_vec([3,4]));

% Get data for single neuron during BMI trials
zz = squeeze(z_binned(cur_trials,:,n_ind))';
xx = squeeze(x_binned(cur_trials,:,3))'.*(180/pi);
ee = error_mat(cur_trials,:)';
% Index of when movements were heading correcting
signs_opposite = sign(ee) ~= sign(xx);
% Index of when heading deviations were above threshold
above_thresh = abs(ee)>error_thresh;
% Combine both conditions
kept_bins = signs_opposite&above_thresh;

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

% reduce data to only above threshold
zz = zz(kept_bins);
xx = xx(kept_bins);
n_samps = length(xx);

colour_vec = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];

% Plot only for corrective movements
figure
scatter(abs(xx),zz,'filled','MarkerEdgeColor',colour_vec(2,:),'MarkerFaceColor',colour_vec(2,:))

hold on

mdl = fitlm(abs(xx),zz);
xx = linspace(min(abs(xx)),max(abs(xx)),20)';
[ypred,yci] = predict(mdl,xx);
plot(xx,yci(:,1),'--','LineWidth',2,'Color','k')
plot(xx,yci(:,2),'--','LineWidth',2,'Color','k')
plot(xx,ypred,'-','LineWidth',2,'Color','k')
h = fill([xx',fliplr(xx')],[yci(:,2)',fliplr(yci(:,1)')],'k','LineStyle','none');
set(h,'facealpha',.1)

box off
xlabel("View angle velocity magnitude (deg/s)")
ylabel(["Mean BMI trial";"normalised neural activity"])
yline(0,'--','LineWidth',2);
title(["Example neural activity";"vs heading correcting movements"])
axis('square')