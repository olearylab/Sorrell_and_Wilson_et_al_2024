function [] = plot_perf_across_days(summary_cell)
% 09/06/2023
num_mice = size(summary_cell,1);
num_days = size(summary_cell,2);

types_vec = [1,4,7,10];

line_width = 2;
marker_size = 10;

performance_sum = nan.*ones(3,num_mice,num_days);
performance_sum_ball = nan.*ones(3,num_mice,num_days);

for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(summary_cell{m,d})
            summary_data = summary_cell{m,d};
            for i = 1:2
                performance_sum(i,m,d) = summary_data(types_vec(2+i))/sum(summary_data(types_vec(2+i):types_vec(2+i)+2));
                performance_sum_ball(i,m,d) = summary_data(types_vec(i))/sum(summary_data(types_vec(i):types_vec(i)+2));
            end
            performance_sum(3,m,d) = (summary_data(types_vec(3))+summary_data(types_vec(4)))/sum(summary_data(types_vec(3):end));
            performance_sum_ball(3,m,d) = (summary_data(types_vec(1))+summary_data(types_vec(2)))/sum(summary_data(1:6));
        end
    end
end

% Plot BMI trial performance across days
figure
xx = 1:num_days;
for m = 1:num_mice
    idx = ~isnan(performance_sum(1,m,:));
    plot(xx(idx),squeeze(performance_sum(3,m,idx)),'-o','LineWidth',line_width,'Color','k','MarkerFaceColor','k','MarkerSize',marker_size)
    hold on
end
box off
yticks([0,0.25,0.5,0.75,1])
ylim([0,1])
yline(0.5,'--','LineWidth',2);
xticks([1,2,3,4])
xlim([0.75,4.25])
title(["BMI Trial Performance Across Days"])
xlabel("Day")
ylabel("Fraction Correct")

% Same but for Ball trials
figure
xx = 1:num_days;
for m = 1:num_mice
    idx = ~isnan(performance_sum_ball(1,m,:));
    plot(xx(idx),squeeze(performance_sum_ball(3,m,idx)),'-o','LineWidth',line_width,'Color','k','MarkerFaceColor','k','MarkerSize',marker_size)
    hold on
end
box off
yticks([0,0.25,0.5,0.75,1])
ylim([0,1])
yline(0.5,'--','LineWidth',2);
xticks([1,2,3,4])
xlim([0.75,4.25])
title(["Ball Trial Performance Across Days"])
xlabel("Day")
ylabel("Fraction Correct")