function [all_r2_mat,trial_avs_all,trial_av_r2,binned_r2,trial_av_cd_r2,trial_av_rmse] = linear_model_crossvals_shuffle(virmen_data,zdata,tbt_details,model_params,n_shuff,extra_norm,lambda_mat,zscore_x,nbins,shuff_limit)
% 23/11/2023

% fit linear regression models of activity as alternative to GLM. Get R^2
% zdata should be neurons

% Run for each direction separately
rng(1);
nModels = 5; % number of cv partitions
types_vec = [1,4,7,10];
kept_types = [1,2,3,4];
trial_num = virmen_data(12,:);
zdim = size(zdata,2);

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
% 1 = train ball/ test ball. 2 = train ball/ test bmi. 3 = train bmi/ test
% bmi. 4 = train bmi/ test ball
all_r2_mat = nan.*ones(n_shuff,2,4,zdim);
binned_r2 = nan.*ones(n_shuff,2,4,zdim);
% 1 = ball trial true, 2 = train ball/ test ball. 3 = train ball/ test bmi. 
% 4 = bmi trial true, 5 = train bmi/ test bmi. 6 = train bmi/ test ball
trial_avs = NaN(n_shuff,2,2,nbins,zdim);
trial_avs_in = NaN(n_shuff,2,2,nbins,zdim);
trial_avs_cross = NaN(n_shuff,2,2,nbins,zdim,nModels);
trial_av_r2 = NaN(n_shuff,2,4,zdim);
trial_av_cd_r2 = NaN(n_shuff,2,4,zdim);
trial_av_rmse = NaN(n_shuff,2,4,zdim);

for n = 1:n_shuff

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

        % shuffle behavioural data
        num_samples = sum(ismember(trial_num,ball_trial_nums));
        shuff_ind = randi([shuff_limit,(num_samples-shuff_limit)]);
        virmen_data(ind_vec,ismember(trial_num,ball_trial_nums)) = circshift(virmen_data(ind_vec,ismember(trial_num,ball_trial_nums)),shuff_ind,2);
        num_samples = sum(ismember(trial_num,bmi_trial_nums));
        shuff_ind = randi([shuff_limit,(num_samples-shuff_limit)]);
        virmen_data(ind_vec,ismember(trial_num,bmi_trial_nums)) = circshift(virmen_data(ind_vec,ismember(trial_num,bmi_trial_nums)),shuff_ind,2);       
        
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
        zballtest_balltrain = NaN(size(zfilt));
        zbmitest_balltrain = NaN(size(zfilt,1),size(zfilt,2),nModels);
        zbmitest_bmitrain = NaN(size(zfilt));
        zballtest_bmitrain = NaN(size(zfilt,1),size(zfilt,2),nModels);
        for cvSplitNumber=1:nModels
            trainIdx = c.training(cvSplitNumber);
            testIdx = c.test(cvSplitNumber);

            XTrain = virmen_data(ind_vec,ismember(trial_num,ball_trial_nums(trainIdx)))';
            ZTrain = zfilt(ismember(trial_num,ball_trial_nums(trainIdx)),:);
            
            % add bias term
            XTrain = [XTrain,ones(size(XTrain,1),1)];

            % run for neurons regularisation separately
            W = NaN(length(ind_vec)+1,zdim);
            for neu = 1:zdim
                if lambda_mat(i,1,neu) == 0
                    W(:,neu) = XTrain\ZTrain(:,neu);
                else
                    W(:,neu) = inv(XTrain'*XTrain + lambda_mat(i,1,neu).*eye(size(XTrain,2)))*XTrain'*ZTrain(:,neu);
                end
            end
              

            XTest = virmen_data(ind_vec,ismember(trial_num,ball_trial_nums(testIdx)))';
            ZTest = zfilt(ismember(trial_num,ball_trial_nums(testIdx)),:);
            
            % add bias term
            XTest = [XTest,ones(size(XTest,1),1)];

            zballtest_balltrain(ismember(trial_num,ball_trial_nums(testIdx)),:) = XTest*W;

            bmi_XTest = virmen_data(ind_vec,ismember(trial_num,bmi_trial_nums))';
            bmi_ZTest = zfilt(ismember(trial_num,bmi_trial_nums),:);
            
            % add bias term
            bmi_XTest = [bmi_XTest,ones(size(bmi_XTest,1),1)];

            zbmitest_balltrain(ismember(trial_num,bmi_trial_nums),:,cvSplitNumber) = bmi_XTest*W;

            % Not used
            testSplitID(testIdx)=cvSplitNumber;
        end

        % Calculate R^2
        % Explained variance method would only work for the training
        % dataset?
        % all_r2_mat(n,i,1,:) = sum((zballtest_balltrain(ismember(trial_num,ball_trial_nums),:) - mean(zfilt(ismember(trial_num,ball_trial_nums),:))).^2)./sum((zfilt(ismember(trial_num,ball_trial_nums),:) - mean(zfilt(ismember(trial_num,ball_trial_nums),:))).^2);
        
        % coefficient of determination method
%         all_r2_mat(n,i,1,:) = 1 - sum((zfilt(ismember(trial_num,ball_trial_nums),:) - zballtest_balltrain(ismember(trial_num,ball_trial_nums),:)).^2)./sum((zfilt(ismember(trial_num,ball_trial_nums),:) - mean(zfilt(ismember(trial_num,ball_trial_nums),:))).^2);
%         temp_r2 = NaN(nModels,zdim);
%         for nn = 1:nModels
%             temp_r2(nn,:) = 1 - sum((zfilt(ismember(trial_num,bmi_trial_nums),:) - zbmitest_balltrain(ismember(trial_num,bmi_trial_nums),:,nn)).^2)./sum((zfilt(ismember(trial_num,bmi_trial_nums),:) - mean(zfilt(ismember(trial_num,bmi_trial_nums),:))).^2);
%         end
%         all_r2_mat(n,i,2,:) = mean(temp_r2);

        % Squared pearson method
        all_r2_mat(n,i,1,:) = diag(corr(zfilt(ismember(trial_num,ball_trial_nums),:),zballtest_balltrain(ismember(trial_num,ball_trial_nums),:))).^2;
        temp_r2 = NaN(nModels,zdim);
        for nn = 1:nModels
            temp_r2(nn,:) = diag(corr(zfilt(ismember(trial_num,bmi_trial_nums),:),zbmitest_balltrain(ismember(trial_num,bmi_trial_nums),:,nn))).^2;
        end
        all_r2_mat(n,i,2,:) = mean(temp_r2);

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
            
            % run for neurons regularisation separately
            W = NaN(length(ind_vec)+1,zdim);
            for neu = 1:zdim
                if lambda_mat(i,2,neu) == 0
                    W(:,neu) = XTrain\ZTrain(:,neu);
                else
                    W(:,neu) = inv(XTrain'*XTrain + lambda_mat(i,2,neu).*eye(size(XTrain,2)))*XTrain'*ZTrain(:,neu);
                end
            end

            XTest = virmen_data(ind_vec,ismember(trial_num,bmi_trial_nums(testIdx)))';
            ZTest = zfilt(ismember(trial_num,bmi_trial_nums(testIdx)),:);
            
            % add bias term
            XTest = [XTest,ones(size(XTest,1),1)];

            zbmitest_bmitrain(ismember(trial_num,bmi_trial_nums(testIdx)),:) = XTest*W;

            ball_XTest = virmen_data(ind_vec,ismember(trial_num,ball_trial_nums))';
            ball_ZTest = zfilt(ismember(trial_num,ball_trial_nums),:);
            
            % add bias term
            ball_XTest = [ball_XTest,ones(size(ball_XTest,1),1)];

            zballtest_bmitrain(ismember(trial_num,ball_trial_nums),:,cvSplitNumber) = ball_XTest*W;

            % Not used
            testSplitID(testIdx)=cvSplitNumber;
        end

        % Calculate R^2
%         all_r2_mat(n,i,3,:) = 1 - sum((zfilt(ismember(trial_num,bmi_trial_nums),:) - zbmitest_bmitrain(ismember(trial_num,bmi_trial_nums),:)).^2)./sum((zfilt(ismember(trial_num,bmi_trial_nums),:) - mean(zfilt(ismember(trial_num,bmi_trial_nums),:))).^2);
%         temp_r2 = NaN(nModels,zdim);
%         for nn = 1:nModels
%             temp_r2(nn,:) = 1 - sum((zfilt(ismember(trial_num,ball_trial_nums),:) - zballtest_bmitrain(ismember(trial_num,ball_trial_nums),:,nn)).^2)./sum((zfilt(ismember(trial_num,ball_trial_nums),:) - mean(zfilt(ismember(trial_num,ball_trial_nums),:))).^2);
%         end
%         all_r2_mat(n,i,4,:) = mean(temp_r2);  

        % Squared pearson method
        all_r2_mat(n,i,3,:) = diag(corr(zfilt(ismember(trial_num,bmi_trial_nums),:),zbmitest_bmitrain(ismember(trial_num,bmi_trial_nums),:))).^2;
        temp_r2 = NaN(nModels,zdim);
        for nn = 1:nModels
            temp_r2(nn,:) = diag(corr(zfilt(ismember(trial_num,ball_trial_nums),:),zballtest_bmitrain(ismember(trial_num,ball_trial_nums),:,nn))).^2;
        end
        all_r2_mat(n,i,4,:) = mean(temp_r2);

        
        % Calculate trial_averages
        z_binned = zeros(length(ball_trial_nums),nbins,zdim);
        z_binned_in = zeros(length(ball_trial_nums),nbins,zdim);
        z_binned_cross = zeros(length(ball_trial_nums),nbins,zdim,nModels);
        cur_lin_pos = virmen_data(6,:) + abs(virmen_data(5,:));
        for j = 1:length(ball_trial_nums)
            [binned] = discretize(cur_lin_pos(trial_num==ball_trial_nums(j)),edges);
            ztrial = zfilt(trial_num==ball_trial_nums(j),:);
            ztrial_in = zballtest_balltrain(trial_num==ball_trial_nums(j),:);
            ztrial_cross = zballtest_bmitrain(trial_num==ball_trial_nums(j),:,:);
            for jj = 1:nbins
                if sum(binned==jj) == 1
                    z_binned(j,jj,:) = ztrial(binned==jj,:);
                    z_binned_in(j,jj,:) = ztrial_in(binned==jj,:);
                    z_binned_cross(j,jj,:,:) = ztrial_cross(binned==jj,:,:);
                else
                    z_binned(j,jj,:) = mean(ztrial(binned==jj,:),'omitnan');
                    z_binned_in(j,jj,:) = mean(ztrial_in(binned==jj,:),'omitnan');
                    z_binned_cross(j,jj,:,:) = mean(ztrial_cross(binned==jj,:,:),'omitnan');
                end
            end
        end
        
        % Calculate binned R2 over all trials 
        z_b = permute(z_binned,[3,1,2]);
        z_b_in = permute(z_binned_in,[3,1,2]);
        z_b_cross = permute(z_binned_cross,[3,4,1,2]);
        binned_r2(n,i,1,:) = diag(corr(z_b(:,:)',z_b_in(:,:)')).^2;
        temp_r2 = NaN(nModels,zdim);
        for nn = 1:nModels
            temp_r2(nn,:) = diag(corr(z_b(:,:)',squeeze(z_b_cross(:,nn,:))')).^2;
        end
        binned_r2(n,i,2,:) = mean(temp_r2);
        
        trial_avs(n,i,1,:,:) = mean(z_binned,'omitnan');
        trial_avs_in(n,i,1,:,:) = mean(z_binned_in,'omitnan');
        trial_avs_cross(n,i,1,:,:,:) = mean(z_binned_cross,'omitnan');
        
        % bmi test
        z_binned = zeros(length(bmi_trial_nums),nbins,zdim);
        z_binned_in = zeros(length(bmi_trial_nums),nbins,zdim);
        z_binned_cross = zeros(length(bmi_trial_nums),nbins,zdim,nModels);
        cur_lin_pos = virmen_data(6,:) + abs(virmen_data(5,:));
        for j = 1:length(bmi_trial_nums)
            [binned] = discretize(cur_lin_pos(trial_num==bmi_trial_nums(j)),edges);
            ztrial = zfilt(trial_num==bmi_trial_nums(j),:);
            ztrial_in = zbmitest_bmitrain(trial_num==bmi_trial_nums(j),:);
            ztrial_cross = zbmitest_balltrain(trial_num==bmi_trial_nums(j),:,:);
            for jj = 1:nbins
                if sum(binned==jj) == 1
                    z_binned(j,jj,:) = ztrial(binned==jj,:);
                    z_binned_in(j,jj,:) = ztrial_in(binned==jj,:);
                    z_binned_cross(j,jj,:,:) = ztrial_cross(binned==jj,:,:);
                else
                    z_binned(j,jj,:) = mean(ztrial(binned==jj,:),'omitnan');
                    z_binned_in(j,jj,:) = mean(ztrial_in(binned==jj,:),'omitnan');
                    z_binned_cross(j,jj,:,:) = mean(ztrial_cross(binned==jj,:,:),'omitnan');
                end
            end
        end
        
        % Calculate binned R2 over all trials 
        z_b = permute(z_binned,[3,1,2]);
        z_b_in = permute(z_binned_in,[3,1,2]);
        z_b_cross = permute(z_binned_cross,[3,4,1,2]);
        binned_r2(n,i,3,:) = diag(corr(z_b(:,:)',z_b_in(:,:)')).^2;
        temp_r2 = NaN(nModels,zdim);
        for nn = 1:nModels
            temp_r2(nn,:) = diag(corr(z_b(:,:)',squeeze(z_b_cross(:,nn,:))')).^2;
        end
        binned_r2(n,i,4,:) = mean(temp_r2);
        
        trial_avs(n,i,2,:,:) = mean(z_binned,'omitnan');
        trial_avs_in(n,i,2,:,:) = mean(z_binned_in,'omitnan');
        trial_avs_cross(n,i,2,:,:,:) = mean(z_binned_cross,'omitnan');
        
        % squared pearson r2 method
        % Changed order: 1 ball/ball 2 bmi/ball 3 bmi/bmi 4 ball/bmi
        trial_av_r2(n,i,1,:) = diag(corr(squeeze(trial_avs(n,i,1,:,:)),squeeze(trial_avs_in(n,i,1,:,:)))).^2;
        temp_r2 = NaN(nModels,zdim);
        for nn = 1:nModels
            temp_r2(nn,:) = diag(corr(squeeze(trial_avs(n,i,1,:,:)),squeeze(trial_avs_cross(n,i,1,:,:,nn)))).^2;
        end
        trial_av_r2(n,i,2,:) = mean(temp_r2);
        trial_av_r2(n,i,3,:) = diag(corr(squeeze(trial_avs(n,i,2,:,:)),squeeze(trial_avs_in(n,i,2,:,:)))).^2;
        temp_r2 = NaN(nModels,zdim);
        for nn = 1:nModels
            temp_r2(nn,:) = diag(corr(squeeze(trial_avs(n,i,2,:,:)),squeeze(trial_avs_cross(n,i,2,:,:,nn)))).^2;
        end
        trial_av_r2(n,i,4,:) = mean(temp_r2);
        
        % coefficient of determination
        % Changed order: 1 ball/ball 2 bmi/ball 3 bmi/bmi 4 ball/bmi
        trial_av_cd_r2(n,i,1,:) = 1 - sum((squeeze(trial_avs(n,i,1,:,:)) - squeeze(trial_avs_in(n,i,1,:,:))).^2)./sum((squeeze(trial_avs(n,i,1,:,:)) - mean(squeeze(trial_avs(n,i,1,:,:)))).^2);
        temp_r2 = NaN(nModels,zdim);
        for nn = 1:nModels
            temp_r2(nn,:) = 1 - sum((squeeze(trial_avs(n,i,1,:,:)) - squeeze(trial_avs_cross(n,i,1,:,:,nn))).^2)./sum((squeeze(trial_avs(n,i,1,:,:)) - mean(squeeze(trial_avs(n,i,1,:,:)))).^2);
        end
        trial_av_cd_r2(n,i,2,:) = mean(temp_r2);
        trial_av_cd_r2(n,i,3,:) = 1 - sum((squeeze(trial_avs(n,i,2,:,:)) - squeeze(trial_avs_in(n,i,2,:,:))).^2)./sum((squeeze(trial_avs(n,i,2,:,:)) - mean(squeeze(trial_avs(n,i,2,:,:)))).^2);
        temp_r2 = NaN(nModels,zdim);
        for nn = 1:nModels
            temp_r2(nn,:) = 1 - sum((squeeze(trial_avs(n,i,2,:,:)) - squeeze(trial_avs_cross(n,i,2,:,:,nn))).^2)./sum((squeeze(trial_avs(n,i,2,:,:)) - mean(squeeze(trial_avs(n,i,2,:,:)))).^2);
        end
        trial_av_cd_r2(n,i,4,:) = mean(temp_r2);
        
        % rmse
        % Changed order: 1 ball/ball 2 bmi/ball 3 bmi/bmi 4 ball/bmi
        trial_av_rmse(n,i,1,:) = sqrt(mean((squeeze(trial_avs(n,i,1,:,:)) - squeeze(trial_avs_in(n,i,1,:,:))).^2));
        temp_r2 = NaN(nModels,zdim);
        for nn = 1:nModels
            temp_r2(nn,:) = sqrt(mean((squeeze(trial_avs(n,i,1,:,:)) - squeeze(trial_avs_cross(n,i,1,:,:,nn))).^2));
        end
        trial_av_rmse(n,i,2,:) = mean(temp_r2);
        trial_av_rmse(n,i,3,:) = sqrt(mean((squeeze(trial_avs(n,i,2,:,:)) - squeeze(trial_avs_in(n,i,2,:,:))).^2));
        temp_r2 = NaN(nModels,zdim);
        for nn = 1:nModels
            temp_r2(nn,:) = sqrt(mean((squeeze(trial_avs(n,i,2,:,:)) - squeeze(trial_avs_cross(n,i,2,:,:,nn))).^2));
        end
        trial_av_rmse(n,i,4,:) = mean(temp_r2);
        
    end 
end
trial_avs_all.trial_avs = trial_avs;
trial_avs_all.trial_avs_in = trial_avs_in;
trial_avs_all.trial_avs_cross = trial_avs_cross;