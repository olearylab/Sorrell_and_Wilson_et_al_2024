function [Wout, xpred, origztrainfilt, train_mean, Rfixed, train_valid, model_params] = DW_online_training_06_2023(ztrain, virmen_data, model_params)
% 26/06/2023

% Edited from 05_2021 to deal with training datasets that don't have all
% trials included. Assumes only correct trials already 

% Outer training function for cleaning training data and downsampling
% images (if need be).

% zimages: samplesxheightxwidth (possibly downsampled)
% if neurons, zimages: samplesxneurons
% virmen_data: samplesxvariables (full virmen data matrix)

%%

% trim data
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
%tic
%disp('training start')
[Wout, xpred, origztrainfilt, model_params] = LMS_training_Mod_New_05_2021(ztrain, xtrain, model_params, train_valid);
%toc;
%disp('training complete')