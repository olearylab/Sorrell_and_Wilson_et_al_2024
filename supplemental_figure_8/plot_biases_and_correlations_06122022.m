function [RMSE_all,corr_vals,p_vals] = plot_biases_and_correlations_06122022(virmen_cell,tbt_cell,summary_cell,nbins,linearise_x,ex_offsets)
% 06/12/2022

% show various biases of mice and potential correlations with performance.
% For now just main mice and days (maybe change to full for complete
% analysis). might make things easier if 1d

% 1. main corrider bias
% 2. magnitude of turn bias
% 3. Decoder weight bias.
% 4. control trial view angle RMSE
% 5. training data RSME

num_mice = size(virmen_cell,1);
num_days = size(virmen_cell,2);

% subtract offset from pitch and yaw, and convert into velocities
circum = 64;
V = 0.32;
alpha = -50/75*circum/V;
beta = 0.05*circum/V/2.5;

% convert from vu into m
vu_conv = 0.0074;

% colour_vec = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880];[0.3010 0.7450 0.9330];[0.6350 0.0780 0.1840]];
corr_vals = nan.*ones(3,1);
p_vals = nan.*ones(3,1);
%% Plot straight section bias

% Plot spatially binned lines, then average over all or specific bin?
% Could just plot whole trial?

va_means = nan.*ones(4,nbins,num_mice,num_days);
% Maybe should also save confidence intervals/ stds...
types_vec = [1,4,7,10];
RMSE_all = nan.*ones(3,num_mice,num_days);
RMSE_pitch_all = nan.*ones(3,num_mice,num_days);
R_square_all = nan.*ones(3,num_mice,num_days);
R_square_pitch_all = nan.*ones(3,num_mice,num_days);
performance_sum = nan.*ones(3,num_mice,num_days);
for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(virmen_cell{m,d})
            virmen_data = virmen_cell{m,d};
            virmen_data(7,:) = wrapToPi(virmen_data(7,:));
            tbt_details = tbt_cell{m,d};
            % trials x bins x xdim
            [x_binned,centres] = bin_kin_data(virmen_data,7,linearise_x,nbins);
            
            for i = 1:4
                va_means(i,:,m,d) = mean(x_binned(tbt_details(3,:) == types_vec(i),:,:),1,'omitnan');
            end
            
            % Make sure not wrapped for RMSE calculation.
            virmen_data = virmen_cell{m,d};
            summary_data = summary_cell{m,d};
            
            virmen_data(:,13) = alpha*(virmen_data(:,13)-ex_offsets(m));
            virmen_data(:,16) = alpha*(virmen_data(:,16)-ex_offsets(m));
            virmen_data(:,13) = virmen_data(:,13)*vu_conv;
            virmen_data(:,16) = virmen_data(:,16)*vu_conv;
            
            for i = 1:2
                [RMSE, R_square, r_p, test_valid_new, xtest_new, xprediction_new] = calc_errors_R_online_results(virmen_data,tbt_details,i,2,[],[],[13,7]);
                RMSE_all(i,m,d) = RMSE(2);
                RMSE_pitch_all(i,m,d) = RMSE(1);
                R_square_all(i,m,d) = R_square(2);
                R_square_pitch_all(i,m,d) = R_square(1);
                performance_sum(i,m,d) = summary_data(types_vec(2+i))/sum(summary_data(types_vec(2+i):types_vec(2+i)+2));
            end
            [RMSE, R_square, r_p, test_valid_new, xtest_new, xprediction_new] = calc_errors_R_online_results(virmen_data,tbt_details,[1,2],2,[],[],[13,7]);
            RMSE_all(3,m,d) = RMSE(2);
            RMSE_pitch_all(3,m,d) = RMSE(1);
            R_square_all(3,m,d) = R_square(2);
            R_square_pitch_all(3,m,d) = R_square(1);
            performance_sum(3,m,d) = (summary_data(types_vec(3))+summary_data(types_vec(4)))/sum(summary_data(types_vec(3):end));
        end
    end
end

pos_scale = 0.74;
centres = centres*pos_scale;
cue_end = 200*pos_scale;
turn_point = 300*pos_scale;
%% Plot comparison of decoder and performance trends
% uncomment to replace RMSE plots with R^2 plots
% RMSE_all = R_square_all;
% RMSE_pitch_all = R_square_pitch_all;
%
% colour_vec = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880];[0.3010 0.7450 0.9330];[0.6350 0.0780 0.1840]];

colour_vec = [[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0]];
ex_colours = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980]];
% colour_vec(ex_mice,:) = ex_colours;

line_width = 2;
marker_size = 10;
% figure
% Set figure size/position
% set(gcf,'position',[1,459,918,234])
xx = 1:num_days;
for m = 1:num_mice
    idx = ~isnan(RMSE_pitch_all(3,m,:));
    
    figure(1)
    plot(xx(idx),squeeze(RMSE_pitch_all(3,m,idx)),'-o','LineWidth',line_width,'Color',colour_vec(m,:),'MarkerFaceColor',colour_vec(m,:),'MarkerSize',marker_size)
    hold on
    box off
    yticks([0.06,0.09,0.12])
    figure(2)
    plot(xx(idx),squeeze(RMSE_all(3,m,idx)),'-o','LineWidth',line_width,'Color',colour_vec(m,:),'MarkerFaceColor',colour_vec(m,:),'MarkerSize',marker_size);
    hold on
    box off
    %yticks([0.3,0.6,0.9])
    figure(3)
    plot(xx(idx),squeeze(performance_sum(3,m,idx)),'-o','LineWidth',line_width,'Color',colour_vec(m,:),'MarkerFaceColor',colour_vec(m,:),'MarkerSize',marker_size)
    hold on
    box off
    yticks([0.5,0.75,1])
end
% Calculate correlations using all data points and corresponding day (maybe
% not best way)
mice_mat = repmat(1:num_days,num_mice,1);
[corr_vals(1),p_vals(1)] = corr(mice_mat(:),RMSE_all(3,:)','Rows','complete');
% corr_rmse = corr_rmse(1,2);
[corr_vals(2),p_vals(2)] = corr(mice_mat(:),RMSE_pitch_all(3,:)','Rows','complete');
% corr_pitch_rmse = corr_pitch_rmse(1,2);
[corr_vals(3),p_vals(3)] = corr(mice_mat(:),performance_sum(3,:)','Rows','complete');
% corr_bmi = corr_bmi(1,2);

num_dp = 3;
% subplot(1,3,1)
figure(1)
% plot(mean(squeeze(RMSE_pitch_all(3,:,:)),'omitnan'),'-o','Color','k','LineWidth',2)
mdl = fitlm(mice_mat(:),RMSE_pitch_all(3,:));
xx = linspace(min(mice_mat(:)),max(mice_mat(:)),100)';
[ypred,yci] = predict(mdl,xx);
plot(xx,yci(:,1),'--','LineWidth',2,'Color','k')
plot(xx,yci(:,2),'--','LineWidth',2,'Color','k')
plot(xx,ypred,'-','LineWidth',2,'Color','k')
h = fill([xx',fliplr(xx')],[yci(:,2)',fliplr(yci(:,1)')],'k','LineStyle','none');
set(h,'facealpha',.1)
ylabel(["Forward Velocity";"RMSE (m/s)"])
xlabel('Day')
% xlim([0.75,4.25])
xticks([1,2,3,4])
% title("\rho = " + round(corr_pitch_rmse,num_dp))

% subplot(1,3,2)
figure(2)
% plot(mean(squeeze(RMSE_all(3,:,:)),'omitnan'),'-o','Color','k','LineWidth',2)
mdl = fitlm(mice_mat(:),RMSE_all(3,:));
xx = linspace(min(mice_mat(:)),max(mice_mat(:)),100)';
[ypred,yci] = predict(mdl,xx);
plot(xx,yci(:,1),'--','LineWidth',2,'Color','k')
plot(xx,yci(:,2),'--','LineWidth',2,'Color','k')
plot(xx,ypred,'-','LineWidth',2,'Color','k')
h = fill([xx',fliplr(xx')],[yci(:,2)',fliplr(yci(:,1)')],'k','LineStyle','none');
set(h,'facealpha',.1)
ylabel(["View Angle"; "RMSE (rad)"])
xlabel('Day')
% xlim([0.75,4.25])
xticks([1,2,3,4])
% title("\rho = " + round(corr_rmse,num_dp))

% subplot(1,3,3)
figure(3)
% plot(mean(squeeze(performance_sum(3,:,:)),'omitnan'),'-o','Color','k','LineWidth',2)
mdl = fitlm(mice_mat(:),performance_sum(3,:));
xx = linspace(min(mice_mat(:)),max(mice_mat(:)),100)';
[ypred,yci] = predict(mdl,xx);
plot(xx,yci(:,1),'--','LineWidth',2,'Color','k')
plot(xx,yci(:,2),'--','LineWidth',2,'Color','k')
plot(xx,ypred,'-','LineWidth',2,'Color','k')
h = fill([xx',fliplr(xx')],[yci(:,2)',fliplr(yci(:,1)')],'k','LineStyle','none');
set(h,'facealpha',.1)
ylabel(["Fraction Correct"; "(BMI Trials)"])
xlabel('Day')
% xlim([0.75,4.25])
xticks([1,2,3,4])
% title("\rho = " + round(corr_bmi,num_dp))