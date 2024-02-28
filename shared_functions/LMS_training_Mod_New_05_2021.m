function [Wout, xpred, ztrainfilt, model_params] = LMS_training_Mod_New_05_2021(zimages, xvars, model_params, train_valid)
% Inner training function that pre-processes training data then performs
% LMS training

% zimages: samplesxheightxwidth (downsampled)
% if neurons, zimages: samplesxneurons
% xvars: samplesxvariables

% legacy parameter always true for LMS training
update = true;

train_length = size(xvars,1);

% allows for use with neurons or images
if ndims(zimages) == 3
    ztrainfilt = zeros(size(zimages,1),size(zimages,2)^2+1);
    % Initialise weight matrix
    if length(model_params.Win)==1
        W = zeros(size(zimages,2)^2+1,size(xvars,2));
        W(end,:) = mean(xvars,1);
    else
        W = model_params.Win;
    end
else
    ztrainfilt = zeros(size(zimages,1),size(zimages,2)+1);
    if length(model_params.Win)==1
        W = zeros(size(zimages,2)+1,size(xvars,2));
        W(end,:) = mean(xvars,1);
    else
        W = model_params.Win;
    end
end

% burn in dff filter (in reverse).
if model_params.dff
    for i = 1:size(xvars,1)
        % turn into array
        zvar = zimages((end+1)-i,:);

        % DF/F filtering
        [ztrainfilt(i,1:end-1)] = dff_filt(zvar, model_params.afast, model_params.aslow);

    end
end

% Now actually filter

for i = 1:size(xvars,1)
    % turn into array
    zvar = zimages(i,:);

    % DF/F filtering
    if model_params.dff
        [ztrainfilt(i,1:end-1)] = dff_filt(zvar, model_params.afast, model_params.aslow);
    else
        ztrainfilt(i,1:end-1) = zvar;
    end

    % Spatial bandpass filter of data
    if model_params.spatial
        [ztrainfilt(i,1:end-1)] = spatial_bandpass(ztrainfilt(i,1:end-1),model_params.blursmall,model_params.blurlarge);
    end
    ztrainfilt(i,end) = 1;
end

% Optional z-scoring
if model_params.zscore_z
    zdim = size(ztrainfilt,2)-1;
    zmeans = zeros(zdim,1);
    z_vars = zeros(zdim,1);
    for n = 1:zdim
        zmeans(n) = mean(ztrainfilt(train_valid,n));
        z_vars(n) = var(ztrainfilt(train_valid,n));
        ztrainfilt(:,n) = (ztrainfilt(:,n) - zmeans(n))/sqrt(z_vars(n));
    end
    model_params.zmeans = zmeans;
    model_params.z_vars = z_vars;
end
    
if model_params.zscore_x
    xdim = size(xvars,2);
    xmeans = zeros(xdim,1);
    x_vars = zeros(xdim,1);
    for n = 1:xdim
        xmeans(n) = mean(xvars(train_valid,n));
        x_vars(n) = var(xvars(train_valid,n));
        xvars(:,n) = (xvars(:,n) - xmeans(n))/sqrt(x_vars(n));
    end
    model_params.xmeans = xmeans;
    model_params.x_vars = x_vars;
    W(end,:) = mean(xvars,1);
end

% Optional offset subtraction - for pitch and yaw
if model_params.subtract_offsets && sum(ismember(model_params.xnums,13))>0
    xvars(:,model_params.xnums==13) = xvars(:,model_params.xnums==13) - model_params.offsets(1);
    if sum(ismember(model_params.xnums,15))>0 && ~model_params.balance_hist
        xvars(:,model_params.xnums==15) = xvars(:,model_params.xnums==15) - model_params.offsets(2);
    end
end

for i = 1:model_params.loops
    % randomly shuffle training data
    shuffle = randperm(train_length);

    for j = 1:train_length
        
        % LMS
        [xpred, Wout]  = Predict_and_update(ztrainfilt(shuffle(j),:), xvars(shuffle(j),:), W, model_params.lrate, train_valid(shuffle(j)), update);

    end

end