function [results_struct] = run_save_offline_single(xfull, zfull, tbt_details, model_params, reps, training_trials)
% 06/10/2023

% Split into training and testing data sets
[ztrain, ztest, xtrain, xtest] = split_session_data(zfull, xfull, tbt_details, [1,2], training_trials, true, false);
        
% Preprocess training data
[ztrainfilt, train_mean, Rfixed] = preprocess_images_05_2021(ztrain,model_params,[],[],[],true);
  
% Preprocess testing data
[ztestfilt, train_mean, Rfixed] = preprocess_images_05_2021(ztest,model_params,train_mean,Rfixed,[],false);
  
% Change parameters since data is now filtered
model_params.dff = false;
model_params.reg_images = false;
model_params.spatial = false; 
        
for i = 1:reps

    % run decoder training and testing
    [results_struct] = DW_train_test_offline_func_05_2021(ztrainfilt, ztestfilt, xtrain, xtest, model_params, tbt_details, [1,2], 7, 4);

end
        
  