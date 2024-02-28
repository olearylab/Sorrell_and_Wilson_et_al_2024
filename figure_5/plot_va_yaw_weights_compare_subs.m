function [dot_prods,dot_prods_yaw,h_boots] = plot_va_yaw_weights_compare_subs(Wout_online_cell,bmi_weights_cell)
% 29/06/2023

% Plot example weights, and results of dot products of weights.
num_mice = size(bmi_weights_cell,1);
num_days = size(bmi_weights_cell,2);

ex_md = [5,1];
yaw_ind = 4;
nmodels = 5;

% Plot image of weights
ex_w = squeeze(Wout_online_cell{ex_md(1)}(:,2));
exbmi_w = squeeze(bmi_weights_cell{ex_md(1),ex_md(2)}(30,1,:,yaw_ind));

figure

cur_weights = reshape(ex_w(1:end-1),128,128);

up_lim = max(cur_weights(:));
low_lim = min(cur_weights(:));

lim_val = max(abs([low_lim,up_lim]));

imagesc(cur_weights,[-lim_val,lim_val])
axis('square')
box on
xticks([])
yticks([])
colormap('redblue')
colorbar
title("View Angle Weights")

figure
% Transpose needed for mice 1-3. Not transpose for 4-5
% Invert as beta is negative
cur_weights = -1.0*reshape(exbmi_w(1:end-1),128,128);

up_lim = max(cur_weights(:));
low_lim = min(cur_weights(:));

lim_val = max(abs([low_lim,up_lim]));


imagesc(cur_weights,[-lim_val,lim_val])
axis('square')
box on
xticks([])
yticks([])
colormap('redblue')
colorbar
title("Ball Angular Velocity Weights")

%% Calculate dot product and plot
nsubs = size(bmi_weights_cell{1,1},1);

dot_prods = nan.*ones(nmodels*nsubs,num_mice,num_days);
num_shuff = 100;
dot_prods_yaw = nan.*ones(nmodels*nsubs,nmodels*nsubs,num_mice,num_days);

for m = 1:num_mice
    Wout  = squeeze(Wout_online_cell{m}(1:end-1,2));
    % Add in transpose for m1-3 due to way images were saved so online
    % weights are transposed compared to offline weights
    if m < 4
        ww = reshape(Wout,128,128)';
        Wout = ww(:);
    end
    for d = 1:num_days
        if ~isempty(bmi_weights_cell{m,d})
            cur_b = bmi_weights_cell{m,d};
            cur_b = permute(cur_b,[3,4,1,2]);
            cur_b = cur_b(:,:,:);
            for n = 1:nmodels*nsubs
                % -1.0 for offline corrective angular velocity weights
                % since there is an inversion from ball voltages to angular
                % velocities
                cur_n = -1.0*squeeze(cur_b(1:end-1,yaw_ind,n));
                dot_prods(n,m,d) = sum((cur_n.*Wout))/(norm(cur_n)*norm(Wout));
                for i = 1:nmodels*nsubs
                    if i ~= n 
                        cur_ni = -1.0*squeeze(cur_b(1:end-1,yaw_ind,i));
                        dot_prods_yaw(n,i,m,d) = sum((cur_n.*cur_ni))/(norm(cur_n)*norm(cur_ni));
                    end
                end
            end
                
            
        end
    end
end

%% Plot
num_sess = num_mice*num_days;
% offset for plotting
plot_off = linspace(-0.4,0.4,num_sess);

figure
mean_dot_prods = squeeze(mean(dot_prods,'omitnan'));
mean_yaw_prod = squeeze(mean(dot_prods_yaw,[1,2],'omitnan'));
scatter(plot_off+ones(1,num_sess),mean_dot_prods(:)','filled','k')
hold on
plot([plot_off(1)+1,plot_off(end)+1],[mean(mean_dot_prods(:),'omitnan'),mean(mean_dot_prods(:),'omitnan')],'k','LineWidth',2)

scatter(plot_off+2.*ones(1,num_sess),squeeze(mean_yaw_prod(:))','filled','k')
hold on
plot([plot_off(1)+2,plot_off(end)+2],[mean(squeeze(mean_yaw_prod(:)),'omitnan'),mean(squeeze(mean_yaw_prod(:)),'omitnan')],'k','LineWidth',2)

title(["View Angle and Ball Angular Velocity"; "Weight Comparison"])
ylabel("Cosine Similarity")
yline(0,'--','LineWidth',2);
xticks([1,2])
% xlabel("Weights")
xticklabels([])
axis('square')

% Get mean and sem with H-bootstrapping
dots_ready = nan.*ones(2,num_mice,num_days);
for m = 1:num_mice
    cur_dots = mean_dot_prods(m,:);
    cur_dots_yaw = mean_yaw_prod(m,:);
    dots_ready(1,m,1:sum(~isnan(cur_dots))) = cur_dots(~isnan(cur_dots));
    dots_ready(2,m,1:sum(~isnan(cur_dots_yaw))) = cur_dots_yaw(~isnan(cur_dots_yaw));
end

[all_p_boot,all_centres,all_sems] = run_H_boot_ets(squeeze(dots_ready(1,:,:)), squeeze(dots_ready(2,:,:)),false);

h_boots.all_p_boot = all_p_boot;
h_boots.all_centres = all_centres;
h_boots.all_sems = all_sems;
