function [] = plot_example_error_accumulation(results_struct,xfull,tbt_details,t_types,va_ind,yaw_ind,yaw_offset,plot_trial)
% 09/08/2023

[RMSE, R_square, r_p, test_valid_new, xtest_new, xprediction_new] = calc_errors_R_05_2021(results_struct,xfull,tbt_details,t_types,va_ind,yaw_ind,yaw_offset);

% xprediction_new 9 is the view angle from decoded angular velocity

d_col = [0.4940 0.1840 0.5560];

figure

subplot(2,1,1)

yline(0,'--','LineWidth',2);
hold on

ITI = xfull(8,:);
cleaned_valid = clean_valid_data(ITI);
xfull = xfull(:,cleaned_valid);
xtest_new = xtest_new(cleaned_valid,:);
xprediction_new = xprediction_new(cleaned_valid,:);
trial_num = xfull(12,:);

circum = 64;
V = 0.32;
alpha = -50/75*circum/V;
beta = 0.05*circum/V/2.5;

% Convert yaw into angular velocity
xtest_new(:,yaw_ind) = -beta*(xtest_new(:,yaw_ind)-yaw_offset);

xprediction_new(:,yaw_ind) = -beta*(xprediction_new(:,yaw_ind)-yaw_offset);

t = (1:sum(trial_num==plot_trial))/30;

plot(t,xtest_new(trial_num==plot_trial,yaw_ind),'Color',[0.5,0.5,0.5],'LineWidth',2)
hold on 
plot(t,xprediction_new(trial_num==plot_trial,yaw_ind),'Color',d_col,'LineWidth',2)
title("Decoded View Angle Velocity")
xlabel("Time (s)")
ylabel(["View Angle Velocity";"(rad/s)"])

subplot(2,1,2)

yline(0,'--','LineWidth',2);
hold on

plot(t,xtest_new(trial_num==plot_trial,va_ind),'Color',[0.5,0.5,0.5],'LineWidth',2)
hold on 
plot(t,xprediction_new(trial_num==plot_trial,9),'Color',d_col,'LineWidth',2)
title("View Angle from Decoded View Angle Velocity")
xlabel("Time (s)")
ylabel(["View Angle";"(rad)"])

figure
plot(xtest_new(:,va_ind),'Color',[0.5,0.5,0.5],'LineWidth',2)
hold on 
plot(xprediction_new(:,9),'Color',d_col,'LineWidth',2)
title("View Angle from Decoded View Angle Velocity")
