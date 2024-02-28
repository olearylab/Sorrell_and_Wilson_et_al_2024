function [ztrainfilt, train_mean, Rfixed] = preprocess_images_05_2021(ztrain,model_params,train_mean,Rfixed,train_valid,create_reg)

% Get paramters
fs = model_params.fs;
afast = model_params.afast;
aslow = model_params.aslow;

blursmall = model_params.blursmall;
blurlarge = model_params.blurlarge;

% reset df/f filter
clear dff_filt
% inialise output matrix depeding on whether input is neurons or pixels
if ndims(ztrain) == 3
    ztrainfilt = zeros(size(ztrain,1),size(ztrain,2)^2);
else
    ztrainfilt = zeros(size(ztrain,1),size(ztrain,2));
end
ztrain = double(ztrain);

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

% create training mean for image registration
if create_reg
    train_vec = double(ztrain(:,:));
    train_mean = mean(train_vec);
    train_mean = reshape(train_mean,128,128);
    Rfixed = imref2d(size(train_mean));
    clear train_vec
end


%%
% burn in DF/F filter
for i = 1:size(ztrain,1)

    zvar = ztrain((end+1)-i,:);


    % DF/F filtering
    if model_params.dff
        zvar = zvar + model_params.dff_offset;
        [ztrainfilt(i,:)] = dff_filt(zvar, afast, aslow);
    else
        ztrainfilt(i,:) = zvar;
    end

end
%%
% preprocessing
for i = 1:size(ztrain,1)

    zvar = ztrain(i,:);

    
    % If registering images
    if model_params.reg_images
        test_image = reshape(zvar',128,128);
        tformEstimate= imregcorr(test_image,train_mean,'translation');
        movingReg = imwarp(test_image,tformEstimate,'OutputView',Rfixed);
        zreg = movingReg;
        zvar = movingReg(:)';
    else
        zreg = zvar;
    end

    % DF/F filtering
    if model_params.dff
        zvar = zvar + model_params.dff_offset;
        [ztrainfilt(i,:)] = dff_filt(zvar, afast, aslow);
    else
        ztrainfilt(i,:) = zvar;
    end

    % Spatial bandpass filter of data
    if model_params.spatial
        [ztrainfilt(i,:)] = spatial_bandpass(ztrainfilt(i,:),blursmall,blurlarge);
    end
end
% option to remove invalid samples here or not
if ~isempty(train_valid)
    ztrainfilt = ztrainfilt(train_valid,:);
end

% option to zscore data if desired
if model_params.zscore_z
    zdim = size(ztrainfilt,2);
    zmeans = zeros(zdim,1);
    z_vars = zeros(zdim,1);
    for n = 1:zdim
        zmeans(n) = mean(ztrainfilt(:,n));
        z_vars(n) = var(ztrainfilt(:,n));
        ztrainfilt(:,n) = (ztrainfilt(:,n) - zmeans(n))/sqrt(z_vars(n));
    end
    model_params.zmeans = zmeans;
    model_params.z_vars = z_vars;
end