function [z_binned, centres] = bin_neural_data_general_any(virmen_data,zdatafilt,xnum,linearise_x,nbins)
% 02/12/2021
% Function for binning neural data according to any variable
% assumes input is preprocessed
% for neurons or pixels
% Only bins valid data

ITI = virmen_data(8,:);
cleaned_valid = clean_valid_data(ITI);
zdim = size(zdatafilt,2);
trial_num_full = virmen_data(12,:);
% make sure ITI assigned to previous trial
trial_num_full = reassign_ITI_trial_num(cleaned_valid,trial_num_full);
trial_num_clean = trial_num_full(cleaned_valid);

% if linearising maze, add modulus of x position to y position
if linearise_x
    virmen_data(6,:) = virmen_data(6,:) + abs(virmen_data(5,:));
end

% wrap VA to pi
virmen_data(7,:) = wrapToPi(virmen_data(7,:));

% remove invalid samples
zdatafilt_clean = zdatafilt(cleaned_valid,:);
virmen_data_clean = virmen_data(:,cleaned_valid);
% Bin data according to xnum variable
z_binned = zeros(max(trial_num_clean),nbins,zdim);
[binned,edges] = discretize(virmen_data_clean(xnum,:),nbins);
for i = 1:max(trial_num_clean)
    [binned] = discretize(virmen_data_clean(xnum,trial_num_clean==i),edges);
    ztrial = zdatafilt_clean(trial_num_clean == i,:);
    for j = 1:nbins
        if sum(binned==j) == 1
            z_binned(i,j,:) = ztrial(binned==j,:);
        else
            z_binned(i,j,:) = mean(ztrial(binned==j,:),'omitnan');
        end
    end
end

% Get bin centres
centres = (edges(2:end)+edges(1:end-1))/2;