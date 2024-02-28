% Script for making supplement figure 2 legend

lines_vec = ["-o";"--^";"-o";"--^"];
colour_vec = [[0.4940 0.1840 0.5560];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880]];
for i = 1:4
    
    plot(1:10,1:10,lines_vec(i),'LineWidth',2,'Color',colour_vec(i,:),'MarkerFaceColor',colour_vec(i,:),'MarkerSize',10)
    hold on
    
end

legend("Ball Left","Ball Right","BMI Left","BMI Right")

figure
lines_vec = ["-";"--";"-";"--"];
colour_vec = [[0.4940 0.1840 0.5560];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880]];
for i = 1:4
    
    plot(1:10,1:10,lines_vec(i),'LineWidth',2,'Color',colour_vec(i,:),'MarkerFaceColor',colour_vec(i,:),'MarkerSize',10)
    hold on
    
end

legend("Ball Left","Ball Right","BMI Left","BMI Right")