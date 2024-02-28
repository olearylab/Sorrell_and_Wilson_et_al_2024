function [kept_trials] = subsample_trials(tbt_details,types_vec)
% 05/09/2022

% Subsample by number of trials to have equal number of each trial type.

num_types = length(types_vec);
%%
num_trials = zeros(num_types,1);
for i = 1:num_types
    num_trials(i) = sum(tbt_details(3,:)==types_vec(i));
end

kept_num = min(num_trials);

kept_trials = zeros(num_types,kept_num);

for i = 1:num_types
    cur_trials = find(tbt_details(3,:)==types_vec(i));
    kept_trials(i,:) = datasample(cur_trials,kept_num,'Replace',false);
end

