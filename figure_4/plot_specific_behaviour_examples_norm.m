function [] = plot_specific_behaviour_examples_norm(virmen_data,tbt_details,nbins,offsets,mean_binned,std_binned,centres)
% 01/05/2023

% Plot example ball angular velocity, decoder angular velocity, and heading
% deviation

% label for correct trials of each type in tbt_details(3,:)
types_vec = [1,4,7,10];
% convert bin centres to cm
centres = centres*0.74;

% Calculate heading deviations and binned behavioural data
[error_mat,x_binned] = BMI_correct_errors_norm(virmen_data,tbt_details,mean_binned,std_binned,nbins,offsets);

% Look at bmi left trials
ex_x_binned = x_binned(tbt_details(3,:)==types_vec(3),:,:);
ex_error_mat = error_mat(tbt_details(3,:)==types_vec(3),:);

% specify trial to plot
ex_ind = 6;
% plot
figure
yyaxis left
colour_vec = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];
plot(centres,ex_error_mat(ex_ind,:),'Color','k','LineWidth',3)
ylabel("Normalised heading deviation (a.u.)")
ylim([-5,5])
hold on
yyaxis right
plot(centres,ex_x_binned(ex_ind,:,3),'-','Color',colour_vec(2,:),'LineWidth',3)
hold on
plot(centres,ex_x_binned(ex_ind,:,5),':','Color',colour_vec(2,:),'LineWidth',3)
ylabel("Angular velocity (rad/s)")

xlabel("Linearised position (cm)")
title("Heading correction example")
yline(0,'--','LineWidth',2);
ylim([-1.5,1.5])
box off
axis('square')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = colour_vec(2,:);
