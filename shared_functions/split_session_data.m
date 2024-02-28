function [ztrain, ztest, xtrain, xtest] = split_session_data(zdata, virmen_data, tbt_details, t_types, training_trials, balanced, by_test)

% create a training set and test set from a single session.

% zdata is whole session neural data
% virmen_data is whole session virmen data
% tbt_details has details of trial types and outcomes
% t_types tells us what trial types we are after.
% training trials is total number of trials used for training. If balanced
% this must be even
% balanced determines whether we require an equal number of left and right
% trials for training. All trials are actually kept (not discarded here)
% but will be balanced in DW_Online_training. This ensures there are enough
% of each to give the total desired training trials.
% by_test determines if we are setting number of training or testing trials

% Train and test sets will be split in time. Therefore if balancing left
% and right trials, there will likely be some trials lost (that essentially
% remain in the training set but are removed in training, rather than being
% put into the test set). This more accurately matches what would happen in
% experiments (i.e. a training session and a testing session).

% Added consideration of removal of poor trials (including whatever
% the current criterion for this is), otherwise may be less and unbalanced
% in the end.

types_vec = [1,4,7,10]; % correct trial indicators for types 1,2,3,4
types_vec = types_vec(t_types);
trial_num = virmen_data(12,:);

% find trial numbers of correct trials we want to use
cur_trials = find(ismember(tbt_details(3,:),types_vec));

% remove poor trials: trials where magnitude of view angle increases
% above 2pi. 
poor_trials = [];
for i = 1:length(cur_trials)
    if max(abs(virmen_data(7,trial_num==cur_trials(i)))) > 2*pi
        poor_trials = [poor_trials,cur_trials(i)];
    end
end

% If setting number of test rather than training trials
if by_test
    cur_trials = find(ismember(tbt_details(3,:),types_vec));
    cur_trials = cur_trials(~ismember(cur_trials,poor_trials));
    training_end = cur_trials(max(trial_num) - training_trials);
    train_length = find(trial_num==training_end,1,'last');
    
% if balancing left and right trials in training set
elseif balanced
    assert(~mod(training_trials,2),'Odd number of training trials')
    
    cur_l_trials = find(ismember(tbt_details(3,:),types_vec(1:2:end)));
    cur_r_trials = find(ismember(tbt_details(3,:),types_vec(2:2:end)));
    cur_l_trials = cur_l_trials(~ismember(cur_l_trials,poor_trials));
    cur_r_trials = cur_r_trials(~ismember(cur_r_trials,poor_trials));
    
    l_end = cur_l_trials(training_trials/2);
    r_end = cur_r_trials(training_trials/2);
    training_end = max([l_end, r_end]);
    train_length = find(trial_num==training_end,1,'last');
else
    % if not balancing left and right trials
    cur_trials = find(ismember(tbt_details(3,:),types_vec));
    cur_trials = cur_trials(~ismember(cur_trials,poor_trials));
    training_end = cur_trials(training_trials);
    train_length = find(trial_num==training_end,1,'last');
    
end

% Split into training and testing sets
ztrain = zdata(1:train_length,:,:);
xtrain = virmen_data(:,1:train_length);
ztest = zdata(train_length+ 1:end,:,:);
xtest = virmen_data(:,train_length + 1:end);


    
    