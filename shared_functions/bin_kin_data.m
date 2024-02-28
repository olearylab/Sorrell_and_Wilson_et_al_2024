function [x_binned,centres] = bin_kin_data(virmen_data,x_vec,linearise_x,nbins)
% 17/11/2021
xdim = length(x_vec);
ITI = virmen_data(8,:);
cleaned_valid = clean_valid_data(ITI);
% zdim = size(zdata,2);
trial_num_full = virmen_data(12,:);
% make sure ITI assigned to previous trial
trial_num_full = reassign_ITI_trial_num(cleaned_valid,trial_num_full);
trial_num_clean = trial_num_full(cleaned_valid);
num_trials = max(trial_num_clean);

% if linearising maze, add modulus of x position to y position
if linearise_x
    virmen_data(6,:) = virmen_data(6,:) + abs(virmen_data(5,:));
end

% bin valid data according to position
virmen_data_clean = virmen_data(:,cleaned_valid);
x_binned = zeros(num_trials,nbins,xdim);
[binned,edges] = discretize(virmen_data_clean(6,:),nbins);
for i = 1:num_trials
    [binned] = discretize(virmen_data_clean(6,trial_num_clean==i),edges);
    xtrial = virmen_data_clean(x_vec,trial_num_clean == i)';
    for j = 1:nbins
        if sum(binned==j) == 1
            x_binned(i,j,:) = xtrial(binned==j,:);
        else
            x_binned(i,j,:) = mean(xtrial(binned==j,:));
        end
    end
end

centres = (edges(2:end)+edges(1:end-1))/2;