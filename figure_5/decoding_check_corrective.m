function [all_res_mat,bmi_weights,num_samps] = decoding_check_corrective(virmen_data,zdata,tbt_details,model_params,train_mean,Rfixed,nbins,offsets,mean_binned,std_binned,num_subs)
% 08/09/2023

% run cross-validated decoding on samples only when behaviour is error
% correcting. 
% Then also calculate error separately for large and small error, where
% error is normalised. Large is defined as great than 1 (this is 1 std).

% Only corrective movements, and only care about BMI trials
nModels = 5; % number of cv partitions
types_vec = [1,4,7,10];
kept_types = [3,4];

%% Preproccess neural data
[zfilt, train_mean, Rfixed] = preprocess_images_05_2021(zdata,model_params,train_mean,Rfixed,[],false);
model_params.dff = false;
model_params.spatial = false;
model_params.reg_images = false;

% Store trained weights
bmi_weights = nan.*ones(num_subs,nModels,size(zfilt,2)+1,length(model_params.xnums));

%% Get normalised errors and correcting vec
% remove invalid data from here on
ITI = virmen_data(8,:);
cleaned_valid = clean_valid_data(ITI);
virmen_data = virmen_data(:,cleaned_valid);
zfilt = zfilt(cleaned_valid,:);
trial_num = virmen_data(12,:);

% Calculate heading deviations and whether movements were heading
% correcting
[error_mat,correcting_mat] = check_error_correction_norm_samples(virmen_data,tbt_details,nbins,offsets,mean_binned,std_binned);

orig_virmen_data = virmen_data;
orig_zfilt = zfilt;
orig_tbt_details = tbt_details;
orig_error_mat = error_mat;
orig_correcting_mat = correcting_mat;
orig_trial_num = trial_num;

% Store all decoding accuracies
all_res_mat = nan.*ones(num_subs,3,3,length(model_params.xnums));

max_vas = nan.*ones(size(orig_tbt_details,2),1);
for i = 1:size(orig_tbt_details,2)
    max_vas(i) = max(abs(orig_virmen_data(7,orig_trial_num==i)));
end
% artificially set poor trials to 13 (ignored in training and testing).
orig_tbt_details(3,max_vas>(2*pi)) = 13; 

% Get lowest number of trials across trial types for subsampling
sep_trial_nums = zeros(length(kept_types),1);
for i = 1:length(kept_types)
    sep_trial_nums(i) = sum(orig_tbt_details(3,:) == types_vec(kept_types(i)));
end
kept_num = min(sep_trial_nums(sep_trial_nums~=0));

for n = 1:num_subs
    %% Balance trial types (subsample)
    % Only choose trials that are correct and "good"

    all_kept = [];
    for i = 1:length(kept_types)
        cur_trials = find(orig_tbt_details(3,:)==types_vec(kept_types(i)));
        if ~isempty(cur_trials)
            cur_kept = datasample(cur_trials,kept_num,'Replace',false);
            all_kept = [all_kept,cur_kept];
        end
    end

    virmen_data = orig_virmen_data(:,ismember(orig_trial_num,all_kept));
    zfilt = orig_zfilt(ismember(orig_trial_num,all_kept),:);
    sorted_kept = sort(all_kept);
    tbt_details = orig_tbt_details(:,sorted_kept);

    error_mat = orig_error_mat(ismember(orig_trial_num,all_kept));
    correcting_mat = orig_correcting_mat(ismember(orig_trial_num,all_kept));

    %% Keep only samples where movements were corrective


    virmen_data = virmen_data(:,correcting_mat);
    zfilt = zfilt(correcting_mat,:);
    error_mat = error_mat(correcting_mat);
    trial_num = virmen_data(12,:);

    %% Train on BMI trials
    bmi_trials = [find(tbt_details(3,:) == types_vec(3)),find(tbt_details(3,:) == types_vec(4))];
    bmi_trial_nums = sorted_kept(bmi_trials);
    % Choosing to repeat one split, rather than use 2, to balance left and
    % right
    c = cvpartition(length(bmi_trials)/2,'KFold',nModels);
    testSplitID = NaN(size(tbt_details,2));

    %% Train cross-validated decoders
    Xprediction = NaN(size(virmen_data,2),length(model_params.xnums));
    for cvSplitNumber=1:nModels
        trainIdx = [c.training(cvSplitNumber);c.training(cvSplitNumber)];
        testIdx = [c.test(cvSplitNumber);c.test(cvSplitNumber)];

        XTrain = virmen_data(:,ismember(trial_num,bmi_trial_nums(trainIdx)));
        ZTrain = zfilt(ismember(trial_num,bmi_trial_nums(trainIdx)),:);

        [Wout, xpred, origztrainfilt, train_mean, Rfixed, train_valid, model_params] = DW_online_training_06_2023(ZTrain, XTrain, model_params);

        bmi_weights(n,cvSplitNumber,:,:) = Wout;

        XTest = virmen_data(:,ismember(trial_num,bmi_trial_nums(testIdx)));
        ZTest = zfilt(ismember(trial_num,bmi_trial_nums(testIdx)),:);

        xprediction = zeros(size(XTest,2),length(model_params.xnums));

        % Test decoder
        for i = 1:size(ZTest,1)
        [xpred,zreg] = Decoder_Online_Mod_05_2021(squeeze(ZTest(i,:)), Wout, model_params, train_mean,Rfixed);
        xprediction(i,:) = xpred;
        end

        % Store only the accuracy on the testing trials for each
        % cross-validation
        Xprediction(ismember(trial_num,bmi_trial_nums(testIdx)),:) = xprediction;

        % Not used
        testSplitID(testIdx)=cvSplitNumber;
    end

    %% Calculate decoding accuracies
    test_ITI = virmen_data(8,:);
    [test_valid] = clean_valid_data(test_ITI);

    xtrain = virmen_data';
    xtrain = xtrain(:,model_params.xnums);
    xfull = virmen_data;
    xtest = virmen_data';
    xtest = xtest(:,model_params.xnums);

    results_struct.xtest = xtest;
    results_struct.xtrain = xtrain;
    results_struct.xprediction = Xprediction;
    results_struct.test_valid = test_valid;
    results_struct.Wout = Wout;
    results_struct.train_mean = train_mean;
    results_struct.Rfixed = Rfixed;
    results_struct.model_params = model_params;

    % Use orig_tbt_details as the trial numbers haven't changed
    [RMSE, R_square, r_p, test_valid_new, xtest_new, xprediction_new] = calc_errors_R_05_2021(results_struct,virmen_data,orig_tbt_details,[3,4],7,[],[]);
    all_res_mat(n,3,:,:) = [RMSE;R_square;r_p];
    num_samps = NaN(3,1);
    num_samps(3) = sum(results_struct.test_valid);

    %% low error samples
    results_struct.test_valid = test_valid&(abs(error_mat)<=1);
    [RMSE, R_square, r_p, test_valid_new, xtest_new, xprediction_new] = calc_errors_R_05_2021(results_struct,virmen_data,orig_tbt_details,[3,4],7,[],[]);
    all_res_mat(n,1,:,:) = [RMSE;R_square;r_p];

    num_samps(1) = sum(results_struct.test_valid);
    %% high error samples
    results_struct.test_valid = test_valid&(abs(error_mat)>1);
    [RMSE, R_square, r_p, test_valid_new, xtest_new, xprediction_new] = calc_errors_R_05_2021(results_struct,virmen_data,orig_tbt_details,[3,4],7,[],[]);
    all_res_mat(n,2,:,:) = [RMSE;R_square;r_p];

    num_samps(2) = sum(results_struct.test_valid);

end