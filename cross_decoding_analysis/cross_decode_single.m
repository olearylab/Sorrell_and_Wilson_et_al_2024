function [xprediction] = cross_decode_single(W,ztest,model_params,train_mean,Rfixed)
% 01/10/2024

% for each subsample and cross-validation calculate the decoder output
n_subs = size(W,1);
n_crossval = size(W,2);
n_samps = size(ztest,1);
yaw_ind = 4;
%% Preproccess neural data
tic
disp('Preprocessing start')
[zfilt, train_mean, Rfixed] = preprocess_images_05_2021(ztest,model_params,train_mean,Rfixed,[],false);
disp('Preprocessing stop')
toc
model_params.dff = false;
model_params.spatial = false;
model_params.reg_images = false;
% Test decoder

xprediction = nan.*ones(n_subs,n_crossval,n_samps);
zfilt = [zfilt,ones(size(zfilt,1),1)];

tic
disp('decoding start')
for s = 1:n_subs
    for c = 1:n_crossval
        Wout = squeeze(W(s,c,:,yaw_ind));
        % takes too long
        % [xpred,zreg] = Decoder_Online_Mod_05_2021(squeeze(zfilt(i,:)), Wout, model_params, train_mean,Rfixed);

        xprediction(s,c,:) = zfilt*Wout;
    end
end
disp('decoding stop')
toc