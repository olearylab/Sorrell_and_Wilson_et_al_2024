function [xpred,zreg] = Decoder_Online_Mod_05_2021(zimage, Win, model_params, train_mean,Rfixed)
%downsampling of image
if ~model_params.downsampled
    zimage = imresize(zimage,1/4);
end

zvar = double(zimage);

%%
% If registering images
if model_params.reg_images
    tformEstimate= imregcorr(zvar,train_mean,'translation');
    movingReg = imwarp(zvar,tformEstimate,'OutputView',Rfixed);
    zreg = movingReg;
    zvar = movingReg(:)';
else
    zreg = zimage;
    zvar = zvar(:)';
end
%%

%DF/F filtering
if model_params.dff
    zvar = zvar + model_params.dff_offset;
    [zvar] = dff_filt(zvar, model_params.afast, model_params.aslow);
end

%Spatial bandpass filtering
if model_params.spatial
    [zvar] = spatial_bandpass(zvar,model_params.blursmall,model_params.blurlarge);
end
zvar = [zvar,1];

% optional zscoring
if model_params.zscore_z
    for n = 1:(length(zvar)-1)
        zvar(n) = (zvar(n) - model_params.zmeans(n))/sqrt(model_params.z_vars(n));
    end
end
% Predict 
[xpred]  = decoder_step(zvar, Win);

% recover true values if zscored
if model_params.zscore_x
    for n = 1:length(xpred)
        xpred(n) = xpred(n)*sqrt(model_params.x_vars(n)) + model_params.xmeans(n);
    end
end