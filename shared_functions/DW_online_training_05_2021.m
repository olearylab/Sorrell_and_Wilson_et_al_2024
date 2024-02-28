function [Wout, xpred, origztrainfilt, train_mean, Rfixed, train_valid, model_params] = DW_online_training_05_2021(ztrain, virmen_data, model_params)
% Outer training function for cleaning training data and downsampling
% images (if need be).

% zimages: samplesxheightxwidth (possibly downsampled)
% if neurons, zimages: samplesxneurons
% virmen_data: samplesxvariables (full virmen data matrix)

%%

% trim data if necessary
assert(size(virmen_data,2)>=size(ztrain,1),'Check detected frame triggers')
orig_length = size(ztrain,1);
ztrain = ztrain(~isnan(virmen_data(1,1:orig_length)),:,:);
virmen_data = virmen_data(:,~isnan(virmen_data(1,1:orig_length)));

%%
% downsample images
% still assumes downsampling from 512x512 to 128x128
if ~model_params.downsampled
    ztrain_new = zeros(size(ztrain,1),128,128);
    for i = 1:size(ztrain,1)
        zimage = squeeze(ztrain(i,:,:));
        ztrain_new(i,:,:) = imresize(zimage,1/4);
    end
    ztrain = ztrain_new;
    clear ztrain_new
end

% For calculating mean training image for registration
if  ndims(ztrain) == 3
    train_vec = double(ztrain(:,:));
    train_mean = mean(train_vec);
    train_mean = reshape(train_mean,128,128);
    Rfixed = imref2d(size(train_mean));
    clear train_vec
else
    train_mean = 0;
    Rfixed = 0;
end

% Ensure all pixel values are positive for dff filtering
if model_params.dff
    ztrain = ztrain + model_params.dff_offset;
end

% Obtain valid data indicator vector
train_invalid_sample = virmen_data(8,:);
% Obtain true valid vector
[train_valid] = clean_valid_data(train_invalid_sample);

%% Determine correct trials.
% This is overly complicated due to missed reward indicators in binned data
trial_num = virmen_data(12,:);
trial_type = virmen_data(1,:);
trial_correct = virmen_data(9,:);
ITI = virmen_data(8,:);
trial_by_trial_type = zeros(max(trial_num),1);
trial_by_trial_correct = zeros(max(trial_num),1);
% timeout trials stay as 0
for j = 1:max(trial_num)
    cur_trial_type = trial_type(trial_num == j);
    cur_trial_xpos = virmen_data(5,trial_num==j);
    trial_by_trial_type(j) = cur_trial_type(end);
    if sum(trial_correct(trial_num == j) ~= 0) % reward present
        trial_by_trial_correct(j) = 1; 
    elseif ismember(-1,ITI(trial_num == j)) % trial completed
        if (((trial_by_trial_type(j) == 1) || (trial_by_trial_type(j) == 3)) && cur_trial_xpos(end) < 0) % animal in correct direction for left trials
            trial_by_trial_correct(j) = 1;
        elseif (((trial_by_trial_type(j) == 2) || (trial_by_trial_type(j) == 4)) && cur_trial_xpos(end) > 0) % animal in correct direction for right trials
            trial_by_trial_correct(j) = 1;
        else
        trial_by_trial_correct(j) = -1; %incorrect trials
        end
    end
end

% combine trial type and outcome
trial_by_trial_details = zeros(3,length(trial_by_trial_type));
    
trial_by_trial_details(1,:) = trial_by_trial_type;
trial_by_trial_details(2,:) = trial_by_trial_correct;

trial_by_trial_details(3,trial_by_trial_correct == 1) = (trial_by_trial_type(trial_by_trial_correct == 1)-1)*3+1;
trial_by_trial_details(3,trial_by_trial_correct == -1) = (trial_by_trial_type(trial_by_trial_correct == -1)-1)*3+2;
trial_by_trial_details(3,trial_by_trial_correct == 0) = (trial_by_trial_type(trial_by_trial_correct == 0)-1)*3+3;

%%
% keep only correct trials
cur_trials = find(ismember(trial_by_trial_details(3,:),model_params.types_vec));
kept_trials = ismember(trial_num,cur_trials);

train_valid = train_valid & kept_trials;

% remove poor trials: trials where magnitude of view angle increases
% above 2pi.

poor_trials = [];
if model_params.remove_poor
    for i = 1:length(cur_trials)
        if max(abs(virmen_data(7,trial_num==cur_trials(i)))) > 2*pi
            train_valid(trial_num==cur_trials(i)) = 0;
            poor_trials = [poor_trials,cur_trials(i)];
        end
    end
end
% Balance left and right trials. 
% Discarding any excess trials from the start of the session
if model_params.balance_lr
    cur_l_trials = find(ismember(trial_by_trial_details(3,:),model_params.types_vec(1:2:end)));
    cur_r_trials = find(ismember(trial_by_trial_details(3,:),model_params.types_vec(2:2:end)));
    cur_l_trials = cur_l_trials(~ismember(cur_l_trials,poor_trials));
    cur_r_trials = cur_r_trials(~ismember(cur_r_trials,poor_trials));
    kept_num = min([length(cur_l_trials),length(cur_r_trials)]);
    cur_l_trials = cur_l_trials(length(cur_l_trials)-kept_num+1:end);
    cur_r_trials = cur_r_trials(length(cur_r_trials)-kept_num+1:end);
    all_cur_trials = [cur_l_trials;cur_r_trials]; % Should be all_cur_trials = [cur_l_trials,cur_r_trials]; but makes no difference.

    kept_trials = ismember(trial_num,all_cur_trials);

    train_valid = train_valid & kept_trials;
    
end

% wrap to pi view angle
virmen_data(7,:) = wrapToPi(virmen_data(7,:));

% extract only decoded variables
xtrain = virmen_data(model_params.xnums,:)';

ztrain = double(ztrain);

% Ensure functions with persistent variables are reset
clear dff_filt
clear Predict_and_update
clear Predict_and_update_ZA_LMS_2022

% Optional registering to training mean.
if model_params.train_reg
    for i = 1:size(ztrain,1)
        tformEstimate= imregcorr(squeeze(ztrain(i,:,:)),train_mean,'translation');
        movingReg = imwarp(squeeze(ztrain(i,:,:)),tformEstimate,'OutputView',Rfixed);
        ztrain(i,:,:) = movingReg;
    end
end

%% Training
% perform decoder training
[Wout, xpred, origztrainfilt, model_params] = LMS_training_Mod_New_05_2021(ztrain, xtrain, model_params, train_valid);
