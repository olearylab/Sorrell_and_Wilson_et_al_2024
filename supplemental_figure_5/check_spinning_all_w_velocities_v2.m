function [h_boots] = check_spinning_all_w_velocities_v2(virmen_cell,tbt_cell,offsets,ex_days,ex_trials)
% 30/09/2024

% Function for checking the number of bmi trials where the mice spin the
% ball significantly (|theta|>2pi when integrated), and plotting examples.

% Separate figure for each mouse. Plot binned instead of time

% offsets contains both pitch and yaw offsets.

num_mice = size(virmen_cell,1);
num_days = size(virmen_cell,2);
types_vec = [1,4,7,10];
spin_cell = cell(num_mice,num_days);
spin_cell_dir = cell(num_mice,num_days,2);
spin_sum = nan.*ones(num_mice,num_days);
spin_sum_fraction = nan.*ones(num_mice,num_days);

color_vec_1 = [[0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880]];
color_vec_2 = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];

pos_scale = 0.74;
circum = 64;
V = 0.32;
alpha = -50/75*circum/V;
beta = 0.05*circum/V/2.5;

% get edges for binning
nbins = 50;
virmen_data = virmen_cell{1,1};
ITI = virmen_data(8,:);
test_valid = clean_valid_data(ITI);
virmen_data_clean = virmen_data(:,test_valid);
[bb,edges] = discretize(virmen_data_clean(6,:),nbins);
centres = (edges(2:end)+edges(1:end-1))/2;
centres = centres*pos_scale;
%% Find trials and plot examples
for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(virmen_cell{m,d})
            
            virmen_data = virmen_cell{m,d};
            yaw_vec = virmen_data(15,:);
            pitch_vec = virmen_data(13,:);
            va_vec = virmen_data(7,:);
            lin_pos = virmen_data(6,:) + abs(virmen_data(5,:));
            ITI = virmen_data(8,:);
            test_valid = clean_valid_data(ITI);
            
            [va] = calc_va_alt_05_2021(yaw_vec,va_vec,offsets(m,2),test_valid);
            
            va = va(test_valid);
            va_vec = va_vec(test_valid);
            
            lin_pos = lin_pos(test_valid);
            
            pitch_vec = pitch_vec(test_valid);
            yaw_vec = yaw_vec(test_valid);
            % Convert into true forward and angular velocity
            pitch_vec = (alpha*(pitch_vec-offsets(m,1)))*pos_scale;
            yaw_vec = -beta*(yaw_vec-offsets(m,2));
            
            trial_num = virmen_data(12,test_valid);
            
            tbt_details = tbt_cell{m,d};
            
            % determine largest view angle on each trial
            max_vas = nan.*ones(size(tbt_details,2),1);
            for i = 1:size(tbt_details,2)
                max_vas(i) = max(abs(va(trial_num==i)));
            end
            % get correct bmi trial numbers
            bmi_trials = find(ismember(tbt_details(3,:),types_vec([3,4])));
            dir_vec = tbt_details(1,ismember(tbt_details(3,:),types_vec([3,4])));
            % get trial numbers for correct bmi trials with spins
            spin_list = bmi_trials(max_vas(bmi_trials)>(2*pi));
            dir_list = dir_vec(max_vas(bmi_trials)>(2*pi));
            spin_cell{m,d} = spin_list;
            spin_sum(m,d) = length(spin_list);
            spin_sum_fraction(m,d) = length(spin_list)/length(bmi_trials);
            % left and right separate
            for i = 1:2
                bmi_trials = find(ismember(tbt_details(3,:),types_vec([2+i])));
                spin_list_dir = bmi_trials(max_vas(bmi_trials)>(2*pi));
                spin_cell_dir{m,d,i} = spin_list_dir;
            end
            
            % get mean ball trial velocities
            virmen_data(13,:) = alpha*(virmen_data(13,:)-offsets(m,1))*pos_scale;
            virmen_data(15,:) = -beta*(virmen_data(15,:)-offsets(m,2));
            
            % Bin data: trials x nbins x xdim
            x_vec = [13,15];
            linearise_x=true;
            [x_binned,cc] = bin_kin_data(virmen_data,x_vec,linearise_x,nbins);
            
            boot_samps = 100;
            sub_sample = false;
            [bootstrap_means] = calc_bootsrapped_means_subsample(x_binned,tbt_details,boot_samps,types_vec,sub_sample);
            x_binned_means = zeros(size(bootstrap_means,3),length(x_vec),length(types_vec));

            for i = 1:length(types_vec)
                x_binned_means(:,:,i) = mean(squeeze(bootstrap_means(i,:,:,:)),'omitnan');
            end

            CI_vals = [2.5,97.5];

            % CIs, upper/lower x ntypes x nbins x zdim
            CIs_all = zeros(2,length(types_vec),size(bootstrap_means,3),length(x_vec));
            for i = 1:size(bootstrap_means,3)
                for j = 1:length(types_vec)
                    CIs_all(:,j,i,:) = prctile(squeeze(bootstrap_means(j,:,i,:)),CI_vals);
                end
            end
            % plot examples for example day
            if d==ex_days(m)
                figure

                cur_spin = spin_list(ex_trials(m)); 
                [binned] = discretize(lin_pos(trial_num==cur_spin),edges);
                
%                 % get correct ball trials for direction of example
%                 ball_trials = find(ismember(tbt_details(3,:),types_vec(dir_list(ex_trials(m)))));
%                 % get correct ball trials with sensible trajectories
%                 good_ball_list = ball_trials(max_vas(ball_trials)<(pi));
                
                % Plot integrated view angle and view angle
                subplot(1,3,1)
                cur_trial = va(trial_num==cur_spin);
                [binned_va] = bin_single_trial(cur_trial,binned,nbins);
                % t = (1:length(cur_trial))./30;
                plot(centres,wrapToPi(binned_va),':','Color',color_vec_1(1,:),'LineWidth',2)
                hold on
                cur_trial = va_vec(trial_num==cur_spin);
                [binned_va_vec] = bin_single_trial(cur_trial,binned,nbins);
                plot(centres,wrapToPi(binned_va_vec),'Color',color_vec_1(2,:),'LineWidth',2)
                yline(0,'--','LineWidth',2);
                box off
                legend("Treadmill","Decoder")

%                 cur_ball = good_ball_list(ex_trials(m));
%                 [ball_binned] = discretize(lin_pos(trial_num==cur_ball),edges);
                 subplot(1,3,2)
%                 cur_trial = pitch_vec(trial_num==cur_ball);
%                 [binned_pitch] = bin_single_trial(cur_trial,ball_binned,nbins);
                % t = (1:length(cur_trial))./30;
%                 plot(centres,binned_pitch,'Color',[0.5,0.5,0.6],'LineWidth',2)
%                 hold on
                
                hold on
                h = fill([centres,fliplr(centres)],[squeeze(CIs_all(1,dir_list(ex_trials(m)),:,1))',fliplr(squeeze(CIs_all(2,dir_list(ex_trials(m)),:,1))')],color_vec_2(1,:),'LineStyle','none');
                set(h,'facealpha',.3)
        
                plot(centres,squeeze(x_binned_means(:,1,dir_list(ex_trials(m)))),'Color',color_vec_2(1,:),'LineWidth',2)
                
                cur_trial = pitch_vec(trial_num==cur_spin);
                [binned_pitch] = bin_single_trial(cur_trial,binned,nbins);
                % t = (1:length(cur_trial))./30;
                plot(centres,binned_pitch,'Color',color_vec_2(2,:),'LineWidth',2)
                yline(0,'--','LineWidth',2);
                box off
                legend("Mean ball trial","BMI trial")
                
                subplot(1,3,3)
%                 cur_trial = yaw_vec(trial_num==cur_ball);
%                 [binned_yaw] = bin_single_trial(cur_trial,ball_binned,nbins);
%                 % t = (1:length(cur_trial))./30;
%                 plot(centres,binned_yaw,'Color',[0.5,0.5,0.6],'LineWidth',2)
%                 hold on
                hold on
                h = fill([centres,fliplr(centres)],[squeeze(CIs_all(1,dir_list(ex_trials(m)),:,2))',fliplr(squeeze(CIs_all(2,dir_list(ex_trials(m)),:,2))')],color_vec_2(1,:),'LineStyle','none');
                set(h,'facealpha',.3)

                plot(centres,squeeze(x_binned_means(:,2,dir_list(ex_trials(m)))),'Color',color_vec_2(1,:),'LineWidth',2)

                cur_trial = yaw_vec(trial_num==cur_spin);
                [binned_yaw] = bin_single_trial(cur_trial,binned,nbins);
                % t = (1:length(cur_trial))./30;
                plot(centres,binned_yaw,'Color',color_vec_2(2,:),'LineWidth',2)
                yline(0,'--','LineWidth',2);
                box off

                ax1 = subplot(1,3,1);
                ylabel("(rad)");
                title("Heading direction")
                xlabel("Maze position (cm)")
                ylim([-3.5,3.5])
                % axis('square')

                ax2 = subplot(1,3,2);
                ylabel("(cm/s)");
                title("Ball forward velocity")
                xlabel("Maze position (cm)")
                % axis('square')
                
                ax3 = subplot(1,3,3);
                ylabel("(rad/s)");
                title("Ball angular velocity")
                xlabel("Maze position (cm)")
                % axis('square')

                % linkaxes([ax1,ax2,ax3],'x')
            end
        end
    end
end


%% Plot summary

% Scatter plot with mean line
num_sess = num_mice*num_days;
plot_off = linspace(-0.4,0.4,num_sess);

plot_res = spin_sum_fraction;

figure
scatter(plot_off+1.*ones(1,num_sess),squeeze(plot_res(:))','filled','k')
hold on
plot([1+plot_off(1),1+plot_off(end)],[mean(squeeze(plot_res(:)),'omitnan'),mean(squeeze(plot_res(:)),'omitnan')],'k','LineWidth',2)
title(["Frequency of large ball rotations";"during correct BMI trials"])
ylabel("Fraction with rotations >2\pi")
xticks([])
% xticklabels()
xlabel("Session")
xlim([0.5,1.5])
ylim([0,1])
axis('square')

%% Hierarchical bootstrapping
% What kind of statistical test do we run for this.

rng(1);
[all_p_boot,all_centres,all_sems] = run_H_boot_ets(plot_res, plot_res,true);

h_boots.all_p_boot = all_p_boot;
h_boots.all_centres = all_centres;
h_boots.all_sems = all_sems;