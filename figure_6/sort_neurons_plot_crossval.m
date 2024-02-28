function [all_sorted_means,all_sorted_means_joined,I_store,cur_neurons_set] = sort_neurons_plot_crossval(zdata,virmen_data,tbt_details,sig_results,linearise_x,nbins,boot_samps,xnum,sub_sample,normalise_z,plot_res,kept_types)
% 07/12/2022

% Make visualisations of neurons sorted according to ball or bmi trials

% Use bootstrapped mean tuning curves. normalise between 0 and 1

% Provide as input just CNN neurons. 
% subsample (optional) and normalise (z-score, optional)

%% Bootstrap
% num types x bootsamps x nbins x neurons
[bootstrap_means_1,bootstrap_means_2, centres] = calc_bootstrap_from_raw_subsample_norm_halves(zdata, virmen_data, tbt_details, linearise_x, nbins,boot_samps, xnum,sub_sample,normalise_z);

bootstrap_means_1 = bootstrap_means_1(kept_types,:,:,:);
bootstrap_means_2 = bootstrap_means_2(kept_types,:,:,:);
% Get means
mean_bootstrap_means_1 = squeeze(mean(bootstrap_means_1,2,'omitnan'));
mean_bootstrap_means_2 = squeeze(mean(bootstrap_means_2,2,'omitnan'));
%% (Possibly) keep only significanty tuned neurons
% Option to keep only significantly tuned neurons
if ~isempty(sig_results)
    cur_neurons = sum(sig_results==1)>0;
    
    cur_neurons_l = sum(sig_results([1,3],:)==1)>0;
    cur_neurons_r = sum(sig_results([2,4],:)==1)>0;
    cur_neurons_set = [cur_neurons;cur_neurons_l;cur_neurons_r];
else
    cur_neurons = ones(1,size(mean_bootstrap_means_1,3))==1;
    cur_neurons_set = ones(1,size(mean_bootstrap_means_1,3))==1;
end

% For normalising overall. Group 2 is for plotting
[max_val] = squeeze(max(mean_bootstrap_means_2,[],[1,2]))';
[min_val] = squeeze(min(mean_bootstrap_means_2,[],[1,2]))';

%% sort according to each trial type
num_types = size(mean_bootstrap_means_1,1);
all_sorted_means = zeros(num_types,num_types,size(mean_bootstrap_means_1,3),nbins);
I_store = zeros(num_types,size(mean_bootstrap_means_1,3));
for i = 1:size(mean_bootstrap_means_1,1)
    % Sort using mean_bootstrap_means_1
    [sorted_means,I] = sort_means_04_2022(squeeze(mean_bootstrap_means_1(i,:,:)),false,[]);
    I_store(i,:) = I;
    for j = 1:size(mean_bootstrap_means_2,1)
        % Plot using mean_bootstrap_means_2
        % Normalise
        cur_bootstrap_means = (squeeze(mean_bootstrap_means_2(j,:,:)) - min_val)./(max_val - min_val);
        [sorted_means,II] = sort_means_04_2022(cur_bootstrap_means,false,I);
        all_sorted_means(i,j,:,:) = sorted_means;
        all_sorted_means_joined((i-1)*size(mean_bootstrap_means_2,3)+1:i*size(mean_bootstrap_means_2,3),(j-1)*nbins+1:j*nbins) = sorted_means;
    end
end

%% Plot

% Rows are normalised 0 to 1

centres = centres*0.74;

if plot_res
    
    figure
    I_vec = [];
    for k = 1:size(I_store,1)
        I_vec = [I_vec,I_store(k,:)];
    end
    imagesc(all_sorted_means_joined(cur_neurons(I_vec),:),[0,1])
    x100 = (100/(centres(end)+centres(1)))*(nbins)+0.5;
    x200 = (200/(centres(end)+centres(1)))*(nbins)+0.5;
    xticks([0.5,x100,x200])
    xticklabels([0,100,200])
    yticks([size(all_sorted_means_joined,1)-(150-1),size(all_sorted_means_joined,1)-(100-1),size(all_sorted_means_joined,1)-(50-1)])
    yticklabels([150,100,50])
    box off
    axis('square')
end