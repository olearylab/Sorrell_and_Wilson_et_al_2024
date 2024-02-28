function [results_struct] = DW_train_test_offline_func_05_2021(ztrain, ztest, xtrain, xtest, model_params, tbt_details, t_types, va_ind, yaw_ind)

% Function for training and testing decoders offline
% requires train and testing data already separated
%%

% option to decode sin and cos theta, then combine back to theta at the end
if model_params.va_2d
    sin_va = sin(xtrain(7,:));
    cos_va = cos(xtrain(7,:));
    xtrain = [xtrain; sin_va; cos_va; xtrain(7,:)];
    model_params.xnums = [model_params.xnums,size(xtrain,1)-2,size(xtrain,1)-1];
end

tic
disp('training start')
% Run decoder training
[Wout, xpred, origztrainfilt, train_mean, Rfixed, train_valid, model_params] = DW_online_training_05_2021(ztrain, xtrain, model_params);
toc;
disp('training complete')
% calculate decoder output on training data
xtrainprediction = origztrainfilt*Wout;

% initalise testing decoder output matrix
xprediction = zeros(size(xtest,2),length(model_params.xnums));

% For storing registered test images
ztest_reg = zeros(size(ztest));

tic;
disp('testing start')
% run decoding on testing data
for i = 1:size(ztest,1)
[xpred,zreg] = Decoder_Online_Mod_05_2021(squeeze(ztest(i,:,:)), Wout, model_params, train_mean,Rfixed);
xprediction(i,:) = xpred;
ztest_reg(i,:,:) = zreg;
end
toc;
disp('testing complete')

% If using sin and cos of view angle
if model_params.va_2d
    xprediction = [xprediction, atan2(xprediction(:,end-1),xprediction(:,end))];
    sin_va = sin(xtest(7,:));
    cos_va = cos(xtest(7,:));
    xtest = [xtest; sin_va; cos_va; xtest(7,:)];
    model_params.xnums = [model_params.xnums,size(xtest,1)];
    
end

% get vector indicating valid samples in testing 
test_ITI = xtest(8,:);
[test_valid] = clean_valid_data(test_ITI);

% create reduced matrices for saving
xtrain = xtrain';
xtrain = xtrain(:,model_params.xnums);
xfull = xtest;
xtest = xtest';
xtest = xtest(:,model_params.xnums);

% store results in a struct
results_struct.xtest = xtest;
results_struct.xtrain = xtrain;
results_struct.xprediction = xprediction;
results_struct.xtrainprediction = xtrainprediction;
results_struct.train_valid = train_valid;
results_struct.test_valid = test_valid;
results_struct.Wout = Wout;
results_struct.train_mean = train_mean;
results_struct.Rfixed = Rfixed;
results_struct.model_params = model_params;

yaw_offset = 1.492; % Should be set according to the mouse (not actually needed)
% calculate decoder accuracy in terms of RMSE, R^2 and pearson correlation.
[RMSE, R_square, r_p, test_valid_new, xtest_new, xprediction_new] = calc_errors_R_05_2021(results_struct,xfull,tbt_details,t_types,va_ind,yaw_ind,yaw_offset);
all_res = [RMSE;R_square;r_p];
results_struct.all_res = all_res;