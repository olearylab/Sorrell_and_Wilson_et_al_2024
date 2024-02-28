function [trial_av_r2,lambda_mat] = linear_model_crossvals_lambda(virmen_data,zdata,tbt_details,model_params,extra_norm,lambda_vec,zscore_x,nbins,accuracy_type)
% 22/11/2023

% fit linear regression models of activity as alternative to GLM. Get R^2
% zdata should be neurons. 
% This function for computing optimal lambda

% Run for each direction separately
rng(1);
nModels = 5; % number of cv partitions
types_vec = [1,4,7,10];
kept_types = [1,2,3,4];
trial_num = virmen_data(12,:);
zdim = size(zdata,2);
num_lam = length(lambda_vec);

%% Preproccess neural data
[zfilt, train_mean, Rfixed] = preprocess_images_05_2021(zdata,model_params,[],[],[],false);
model_params.dff = false;
model_params.spatial = false;
model_params.reg_images = false; 

%% Balance trial types (subsample)
% Only choose trials that are correct and "good"
max_vas = nan.*ones(size(tbt_details,2),1);
for i = 1:size(tbt_details,2)
    % NB: Should change to only look at valid samples!!! Need to change
    % everywhere.
    max_vas(i) = max(abs(virmen_data(7,trial_num==i)));
end
tbt_details(3,max_vas>(2*pi)) = 13; % artificially set poor trials to 13.

sep_trial_nums = zeros(length(kept_types),1);
for i = 1:length(kept_types)
    sep_trial_nums(i) = sum(tbt_details(3,:) == types_vec(kept_types(i)));
end
kept_num = min(sep_trial_nums(sep_trial_nums~=0));

%% remove invalid data
ITI = virmen_data(8,:);
cleaned_valid = clean_valid_data(ITI);

orig_virmen_data = virmen_data(:,cleaned_valid);
orig_zfilt = zfilt(cleaned_valid,:);
orig_tbt_details = tbt_details;

orig_trial_num = orig_virmen_data(12,:);

%% Prepare behavioural features
% For now will use pitch, yaw and 10 bumps of linearised position
line_pos = orig_virmen_data(6,:) + abs(orig_virmen_data(5,:));
num_bumps = 10;
bump_means = linspace(min(line_pos),max(line_pos),num_bumps);
bump_std = (bump_means(2)-bump_means(1))/2;

all_bumps = NaN(num_bumps,length(line_pos));

for i = 1:num_bumps
    
    all_bumps(i,:) = normpdf(line_pos,bump_means(i),bump_std);

end

v_length = size(orig_virmen_data,1);
orig_virmen_data = [orig_virmen_data;all_bumps];
%ind_vec = [13,15,v_length+1:size(orig_virmen_data,1)];
ind_vec = [13:15];
%ind_vec = [v_length+1:size(orig_virmen_data,1)];

% Get edges for binning data later
[binned,edges] = discretize(line_pos,nbins);

% for convenience save all in one.
lambda_r2_mat = nan.*ones(num_lam,2,2,zdim);
lambda_binned_r2 = nan.*ones(num_lam,2,2,zdim);

% Compute trial average r2s
trial_avs = NaN(2,2,nbins,zdim);
trial_avs_in = NaN(2,2,nbins,zdim,num_lam);
trial_av_r2 = NaN(num_lam,2,2,zdim);
trial_av_cd_r2 = NaN(num_lam,2,2,zdim);
trial_av_rmse = NaN(num_lam,2,2,zdim);

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

trial_num = virmen_data(12,:);

% Run separately for each direction
for i = 1:2
    %% Train on ball trials
    ball_trials = find(tbt_details(3,:) == types_vec(i));
    c = cvpartition(length(ball_trials),'KFold',nModels);
    % No need to split in 2
    % c = cvpartition(length(ball_trials)/2,'KFold',nModels);
    testSplitID = NaN(size(tbt_details,2));

    % Get original trial numbers
    ball_trial_nums = sorted_kept(ball_trials);

    bmi_trials = find(tbt_details(3,:) == types_vec(i+2));
    bmi_trial_nums = sorted_kept(bmi_trials);

    if extra_norm
        zfilt(ismember(trial_num,ball_trial_nums),:) = zscore(zfilt(ismember(trial_num,ball_trial_nums),:));
        zfilt(ismember(trial_num,bmi_trial_nums),:) = zscore(zfilt(ismember(trial_num,bmi_trial_nums),:));
    else
        zfilt(ismember(trial_num,[ball_trial_nums,bmi_trial_nums]),:) = zscore(zfilt(ismember(trial_num,[ball_trial_nums,bmi_trial_nums]),:));
    end

    %% zscore behaviour as well
    if zscore_x
        for j = 1:length(ind_vec)
            virmen_data(ind_vec(j),ismember(trial_num,[ball_trial_nums,bmi_trial_nums])) = zscore(virmen_data(ind_vec(j),ismember(trial_num,[ball_trial_nums,bmi_trial_nums])));
        end
    end

    %% Train network
    zballtest_balltrain = NaN(size(zfilt,1),size(zfilt,2),num_lam);
    zbmitest_bmitrain = NaN(size(zfilt,1),size(zfilt,2),num_lam);
    for cvSplitNumber=1:nModels
        trainIdx = c.training(cvSplitNumber);
        testIdx = c.test(cvSplitNumber);

        XTrain = virmen_data(ind_vec,ismember(trial_num,ball_trial_nums(trainIdx)))';
        ZTrain = zfilt(ismember(trial_num,ball_trial_nums(trainIdx)),:);

        % add bias term
        XTrain = [XTrain,ones(size(XTrain,1),1)];
        
        XTest = virmen_data(ind_vec,ismember(trial_num,ball_trial_nums(testIdx)))';
        ZTest = zfilt(ismember(trial_num,ball_trial_nums(testIdx)),:);
        
        % add bias term
        XTest = [XTest,ones(size(XTest,1),1)];

        for l = 1:num_lam

            if lambda_vec(l)==0
                W = XTrain\ZTrain;
            else
                W = inv(XTrain'*XTrain + lambda_vec(l).*eye(size(XTrain,2)))*XTrain'*ZTrain;
            end
            
            zballtest_balltrain(ismember(trial_num,ball_trial_nums(testIdx)),:,l) = XTest*W;
            
        end
        % Not used
        testSplitID(testIdx)=cvSplitNumber;
    end

    % Calculate R^2
%     for l = 1:num_lam
%         lambda_r2_mat(l,i,1,:) = 1 - sum((zfilt(ismember(trial_num,ball_trial_nums),:) - zballtest_balltrain(ismember(trial_num,ball_trial_nums),:,l)).^2)./sum((zfilt(ismember(trial_num,ball_trial_nums),:) - mean(zfilt(ismember(trial_num,ball_trial_nums),:))).^2);
%     end
    % squared pearson method
    for l = 1:num_lam
        lambda_r2_mat(l,i,1,:) = diag(corr(zfilt(ismember(trial_num,ball_trial_nums),:),zballtest_balltrain(ismember(trial_num,ball_trial_nums),:,l))).^2;
    end

    %% Train on BMI trials
    c = cvpartition(length(bmi_trials),'KFold',nModels);
    % No need to split in 2
    % c = cvpartition(length(bmi_trials)/2,'KFold',nModels);
    testSplitID = NaN(size(tbt_details,2));

    % for convenience save all in one.
    % 1 = train ball/ test ball. 2 = train ball/ test bmi. 3 = train bmi/ test
    % bmi. 4 = train bmi/ test ball

    %% Train network
    for cvSplitNumber=1:nModels
        trainIdx = c.training(cvSplitNumber);
        testIdx = c.test(cvSplitNumber);

        XTrain = virmen_data(ind_vec,ismember(trial_num,bmi_trial_nums(trainIdx)))';
        ZTrain = zfilt(ismember(trial_num,bmi_trial_nums(trainIdx)),:);

        % add bias term
        XTrain = [XTrain,ones(size(XTrain,1),1)];
        
        XTest = virmen_data(ind_vec,ismember(trial_num,bmi_trial_nums(testIdx)))';
        ZTest = zfilt(ismember(trial_num,bmi_trial_nums(testIdx)),:);

        % add bias term
        XTest = [XTest,ones(size(XTest,1),1)];

        for l = 1:num_lam
            if lambda_vec(l) == 0
                W = XTrain\ZTrain;
            else
                W = inv(XTrain'*XTrain + lambda_vec(l).*eye(size(XTrain,2)))*XTrain'*ZTrain;
            end
            
            zbmitest_bmitrain(ismember(trial_num,bmi_trial_nums(testIdx)),:,l) = XTest*W;
            
        end

        % Not used
        testSplitID(testIdx)=cvSplitNumber;
    end
    
    % Calculate R^2
%     for l = 1:num_lam
%         lambda_r2_mat(l,i,2,:) = 1 - sum((zfilt(ismember(trial_num,bmi_trial_nums),:) - zbmitest_bmitrain(ismember(trial_num,bmi_trial_nums),:,l)).^2)./sum((zfilt(ismember(trial_num,bmi_trial_nums),:) - mean(zfilt(ismember(trial_num,bmi_trial_nums),:))).^2);
%     end
    % Squared Pearson method
    for l = 1:num_lam
        lambda_r2_mat(l,i,2,:) = diag(corr(zfilt(ismember(trial_num,bmi_trial_nums),:),zbmitest_bmitrain(ismember(trial_num,bmi_trial_nums),:,l))).^2;
    end
    
    % ball test
    z_binned = zeros(length(ball_trial_nums),nbins,zdim);
    z_binned_in = zeros(length(ball_trial_nums),nbins,zdim,num_lam);
    cur_lin_pos = virmen_data(6,:) + abs(virmen_data(5,:));
    for j = 1:length(ball_trial_nums)
        [binned] = discretize(cur_lin_pos(trial_num==ball_trial_nums(j)),edges);
        ztrial = zfilt(trial_num==ball_trial_nums(j),:);
        ztrial_in = zballtest_balltrain(trial_num==ball_trial_nums(j),:,:);
        for jj = 1:nbins
            if sum(binned==jj) == 1
                z_binned(j,jj,:) = ztrial(binned==jj,:);
                z_binned_in(j,jj,:,:) = ztrial_in(binned==jj,:,:);
            else
                z_binned(j,jj,:) = mean(ztrial(binned==jj,:),'omitnan');
                z_binned_in(j,jj,:,:) = mean(ztrial_in(binned==jj,:,:),'omitnan');
            end
        end
    end
    
    % Calculate binned R2 over all trials 
    z_b = permute(z_binned,[3,1,2]);
    z_b_in = permute(z_binned_in,[3,4,1,2]);
    for l = 1:num_lam
        lambda_binned_r2(l,i,1,:) = diag(corr(z_b(:,:)',squeeze(z_b_in(:,l,:))')).^2;
    end

    trial_avs(i,1,:,:) = mean(z_binned,'omitnan');
    trial_avs_in(i,1,:,:,:) = mean(z_binned_in,'omitnan');

    % bmi test
    z_binned = zeros(length(bmi_trial_nums),nbins,zdim);
    z_binned_in = zeros(length(bmi_trial_nums),nbins,zdim,num_lam);
    cur_lin_pos = virmen_data(6,:) + abs(virmen_data(5,:));
    for j = 1:length(bmi_trial_nums)
        [binned] = discretize(cur_lin_pos(trial_num==bmi_trial_nums(j)),edges);
        ztrial = zfilt(trial_num==bmi_trial_nums(j),:);
        ztrial_in = zbmitest_bmitrain(trial_num==bmi_trial_nums(j),:,:);
        for jj = 1:nbins
            if sum(binned==jj) == 1
                z_binned(j,jj,:) = ztrial(binned==jj,:);
                z_binned_in(j,jj,:,:) = ztrial_in(binned==jj,:,:);
            else
                z_binned(j,jj,:) = mean(ztrial(binned==jj,:),'omitnan');
                z_binned_in(j,jj,:,:) = mean(ztrial_in(binned==jj,:,:),'omitnan');
            end
        end
    end
    
    z_b = permute(z_binned,[3,1,2]);
    z_b_in = permute(z_binned_in,[3,4,1,2]);
    for l = 1:num_lam
        lambda_binned_r2(l,i,2,:) = diag(corr(z_b(:,:)',squeeze(z_b_in(:,l,:))')).^2;
    end
    
    trial_avs(i,2,:,:) = mean(z_binned,'omitnan');
    trial_avs_in(i,2,:,:,:) = mean(z_binned_in,'omitnan');
    
    for l = 1:num_lam
        trial_av_r2(l,i,1,:) = diag(corr(squeeze(trial_avs(i,1,:,:)),squeeze(trial_avs_in(i,1,:,:,l)))).^2;
        trial_av_r2(l,i,2,:) = diag(corr(squeeze(trial_avs(i,2,:,:)),squeeze(trial_avs_in(i,2,:,:,l)))).^2;
        trial_av_cd_r2(l,i,1,:) = 1 - sum((squeeze(trial_avs(i,1,:,:)) - squeeze(trial_avs_in(i,1,:,:,l))).^2)./sum((squeeze(trial_avs(i,1,:,:)) - mean(squeeze(trial_avs(i,1,:,:)))).^2);
        trial_av_cd_r2(l,i,2,:) = 1 - sum((squeeze(trial_avs(i,2,:,:)) - squeeze(trial_avs_in(i,2,:,:,l))).^2)./sum((squeeze(trial_avs(i,2,:,:)) - mean(squeeze(trial_avs(i,2,:,:)))).^2);  
        trial_av_rmse(l,i,1,:) = sqrt(mean((squeeze(trial_avs(i,1,:,:)) - squeeze(trial_avs_in(i,1,:,:,l))).^2));
        trial_av_rmse(l,i,2,:) = sqrt(mean((squeeze(trial_avs(i,2,:,:)) - squeeze(trial_avs_in(i,2,:,:,l))).^2));   
    end

end

% all r2
if accuracy_type=="r2_corrsquare"
    [~,lambda_inds] = max(lambda_r2_mat);
% trial averaged corrsquare
elseif accuracy_type=="tav_r2_corrsquare"
    [~,lambda_inds] = max(trial_av_r2);
% binned corrsquare
elseif accuracy_type=="bin_r2_corrsquare"
    [~,lambda_inds] = max(lambda_binned_r2);
% trial averaged coefficient of determination
elseif accuracy_type=="tav_r2_cd"
    [~,lambda_inds] = max(trial_av_cd_r2);
% trial averaged RMSE
elseif accuracy_type=="tav_rmse"
    [~,lambda_inds] = max(trial_av_rmse);
end
lambda_mat = lambda_vec(squeeze(lambda_inds));