function [] = plot_example_position_trajectories(virmen_data,tbt_details,num_trials)
% 08/06/2023

% Plot example positional trajectories

types_vec = [1,4,7,10];

to_cm = 0.74;
xpos = to_cm*virmen_data(5,:);
ypos = to_cm*virmen_data(6,:);
trial_num = virmen_data(12,:);
ITI = virmen_data(8,:);
cleaned_valid = clean_valid_data(ITI);

%% PLot just maze
figure
mazex = to_cm*[-0.1,0.1,0.1,1,1,15,15,-15,-15,-1,-1,-0.1,-0.1];
mazey = to_cm*[0,0,300,300,305,305,307.5,307.5,305,305,300,300,0];
plot(mazex,mazey,'Color',[0.5,0.5,0.5],'LineWidth',2)
box off
axis equal
% axis off
hold on
xlim([-125,125])
ylim([0,250])

scalex = [-40,-40,-20];
scaley = [20,0,0];
plot(scalex,scaley,'k','LineWidth',2)

box off
axis off

%% Many plots single direction - corridor
figure

colour_vec = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];
plot_types = [2,4];
for i = 1:2
    cur_trials = find(tbt_details(3,:) == types_vec(plot_types(i)));
    rn = randperm(length(cur_trials),num_trials);
    % rn = 1:num_trials;
    for n = 1:num_trials
        subplot(1,2*num_trials,(i-1)*num_trials + n)
        plot(mazex,mazey,'Color',[0.5,0.5,0.5],'LineWidth',2)
        hold on
        plot(xpos(trial_num==cur_trials(rn(n))&cleaned_valid),ypos(trial_num==cur_trials(rn(n))&cleaned_valid),'-','Color',colour_vec(i,:),'LineWidth',2)
        xlim(to_cm*[-0.11,0.11])
        ylim(to_cm*[0,299])

        box off
        axis off
    end
end

scalex = [-0.05,-0.05,0.05];
scaley = [10,0,0];
plot(scalex,scaley,'k','LineWidth',2)

%% Single plot single direction - turn
figure

plot(mazex,mazey,'Color',[0.5,0.5,0.5],'LineWidth',2)
hold on

for i = 1:2
    cur_trials = find(tbt_details(3,:) == types_vec(plot_types(i)));
    rn = randperm(length(cur_trials),num_trials);
    for n = 1:num_trials
        plot(xpos(trial_num==cur_trials(rn(n))&cleaned_valid),ypos(trial_num==cur_trials(rn(n))&cleaned_valid),'-','Color',colour_vec(i,:),'LineWidth',2)
    end
end

scalex = [-9,-9,-7];
scaley = [224,222,222];
plot(scalex,scaley,'k','LineWidth',2)

xlim(to_cm*[-15,15])
ylim(to_cm*[299,309])

box off
axis off