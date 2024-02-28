function [ss,centres] = prep_data_for_bvn_lstm(virmen_data,tbt_details,nbins,keep_incorrect,kept_types,balance_types)
% 06/05/2022

% Prepare data for input into LSTM function.

types_vec = [1,4,7,10];
virmen_data(7,:) = wrapToPi(virmen_data(7,:));
% bin kinematic data
[x_binned,centres] = bin_kin_data(virmen_data,1:size(virmen_data,1),true,nbins);
% Option for keeping incorrect trials
if ~keep_incorrect
    x_binned = x_binned(tbt_details(2,:)==1,:,:);
else
    x_binned = x_binned(tbt_details(2,:)~=0,:,:);
end
x_binned_perm = permute(x_binned,[1,3,2]);
ss.binnedVirmenData2 = x_binned_perm;



trial_type = tbt_details(1,:);

if ~keep_incorrect
    trial_type = trial_type(tbt_details(2,:)==1);
    ss.correctVec = tbt_details(2,tbt_details(2,:)==1);
    tbt_details_cut = tbt_details(:,tbt_details(2,:)==1);
else
    trial_type = trial_type(tbt_details(2,:)~=0);
    ss.correctVec = tbt_details(2,tbt_details(2,:)~=0);
    tbt_details_cut = tbt_details(:,tbt_details(2,:)~=0);
end

% Store trial type for kept trials
ss.trialType = trial_type(ismember(trial_type,kept_types));
% Store whether trials are correct
ss.correctVec = ss.correctVec(ismember(trial_type,kept_types));
% Stored binned virmen data
ss.binnedVirmenData2 = ss.binnedVirmenData2(ismember(trial_type,kept_types),:,:);
% trim the trial details to only kept trial types
tbt_details_cut = tbt_details_cut(:,ismember(trial_type,kept_types));

% subsample to balance the number of each trial type
if balance_types
    if ~keep_incorrect
        
        sep_trial_nums = zeros(length(kept_types),1);
        for i = 1:length(kept_types)
            sep_trial_nums(i) = sum(tbt_details_cut(3,:) == types_vec(kept_types(i)));
        end
        kept_num = min(sep_trial_nums(sep_trial_nums~=0));

        all_kept = [];
        for i = 1:length(kept_types)
            cur_trials = find(tbt_details_cut(3,:)==types_vec(kept_types(i)));
            if ~isempty(cur_trials)
                cur_kept = datasample(cur_trials,kept_num,'Replace',false);
                all_kept = [all_kept,cur_kept];
            end
        end
    else

        sep_trial_nums = zeros(length(kept_types),1);
        for i = 1:length(kept_types)
            sep_trial_nums(i) = sum(ss.trialType == kept_types(i));
        end
        kept_num = min(sep_trial_nums(sep_trial_nums~=0));

        all_kept = [];
        for i = 1:length(kept_types)
            cur_trials = find(ss.trialType==kept_types(i));
            if ~isempty(cur_trials)
                cur_kept = datasample(cur_trials,kept_num,'Replace',false);
                all_kept = [all_kept,cur_kept];
            end
        end  
        
    end
    
    
    ss.binnedVirmenData2 = ss.binnedVirmenData2(all_kept,:,:);
    ss.trialType = ss.trialType(all_kept);
    ss.correctVec = ss.correctVec(all_kept);
    
end

ss.origTrialType = ss.trialType;
% Set ball trials to type 1, and bmi trials to type 2
ss.trialType(ss.trialType==2) = 1;
ss.trialType(ss.trialType==3) = 2;
ss.trialType(ss.trialType==4) = 2;

ss.correct_only = ~keep_incorrect;