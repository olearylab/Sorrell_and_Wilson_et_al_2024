function [ITI_trial_num] = reassign_ITI_trial_num(cleaned_valid,trial_num)

% function for reassigining trial numbers such that the ITI isn't split
% across 2 trial numbers but assigned entirely to the previous trial.
ITI_trial_num = trial_num;
for i = 2:length(trial_num)
    if ~cleaned_valid(i)
        ITI_trial_num(i) = ITI_trial_num(i-1);
    end
end