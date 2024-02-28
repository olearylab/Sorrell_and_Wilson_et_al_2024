function [full_corr,std_corr,full_corr_cell] = calc_behaviour_correlations(virmen_data,tbt_details,offsets,sub_sample,x_num)
% 25/07/2023

% Calculate correlations between behaviour on ball and bmi trials.

types_vec = [1,4,7,10];
linearise_x = true;
nbins = 50;

% remove ITI
ITI = virmen_data(8,:);
cleaned_valid = clean_valid_data(ITI);
virmen_data = virmen_data(:,cleaned_valid);

% convert into velocity
circum = 64;
V = 0.32;
alpha = -50/75*circum/V;
beta = 0.05*circum/V/2.5;
virmen_data(13,:) = alpha*(virmen_data(13,:)-offsets(1));
virmen_data(15,:) = -beta*(virmen_data(15,:)-offsets(2));

virmen_data(16,:) = alpha*(virmen_data(16,:)-offsets(1));

virmen_data(7,:) = wrapToPi(virmen_data(7,:));

%% Subsample
trial_num = virmen_data(12,:);
if sub_sample
    [kept_trials] = subsample_trials(tbt_details,types_vec);
    kept_trials = kept_trials(:);
else
    kept_trials = 1:max(trial_num);
end

kept_trial_num = ismember(trial_num,kept_trials);

virmen_data = virmen_data(:,kept_trial_num);

% Bin data: trials x nbins x xdim
[x_binned,centres] = bin_kin_data(virmen_data,x_num,linearise_x,nbins);

tbt_details = tbt_details(:,sort(kept_trials));
num_trials = length(kept_trials);

full_trials = cell(4,1);
types_vec = [1,4,7,10];
for i = 1:4
    full_trials{i} = find(ismember(tbt_details(3,:),types_vec(i)));
end

% need to make bins first variable. Also only one variable in x_num.
cur_coeff = corrcoef(x_binned','rows','pairwise');

full_corr = nan*ones(4,4);

full_corr_cell = cell(4,4);

std_corr = nan*ones(4,4);

% Extract correlations within and across trial types
for i = 1:4
    for j = 1:4
        cur_block = cur_coeff(full_trials{i},full_trials{j});
        if i == j
            % Remove diagonals of 1s, and double up of correlations
            low_t = tril(nan.*ones(length(full_trials{i})));
            cur_block = cur_block + low_t;
        end
        cur_block_ready = cur_block(:);
        full_corr(i,j) = mean(cur_block_ready(~isnan(cur_block_ready)),'omitnan');
        full_corr_cell{i,j} = cur_block_ready(~isnan(cur_block_ready));
        std_corr(i,j) = std(cur_block_ready(~isnan(cur_block_ready)));
    end
end