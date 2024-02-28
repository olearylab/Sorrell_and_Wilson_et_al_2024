function [accuracies_mat] = svm_neural_classifier(s)
% 02/08/2023

% Updated SVM classification.

% bin specific

tempNeuralData = s.binnedNeuralData2;

% Z-score each neuron separately
for n = 1:size(tempNeuralData,2)
    tempNeuralData(:,n,:) = zscore(tempNeuralData(:,n,:),0,'all');
end

nbins = size(tempNeuralData,3);

c_type = 'svm';

% set up partitions
% train only on correct trials
s.controlIdx = find(s.trialType<3);
if s.correct_only
    s.controlIdx = s.controlIdx(s.correctVec(s.controlIdx)==1);
end
s.controlBinnedData = tempNeuralData(s.controlIdx,:,:);
s.controlLabels = s.trialType(s.controlIdx);

s.controlLabelSequence = zeros(size(s.controlBinnedData,1),size(s.controlBinnedData,3));

for controlTrialNumber = 1:size(s.controlBinnedData,1);
    s.controlLabelSequence(controlTrialNumber,:) = categorical(repmat(s.controlLabels(controlTrialNumber),[1, size(tempNeuralData,3)]),[1 2]);
end

nModels = 5; % number of cv partitions

c = cvpartition(size(s.controlLabelSequence,1)/2,'KFold',nModels);
testSplitID = NaN(size(s.trialType));

many_Mdls = cell(nModels,nbins);

%% Train svms
for cvSplitNumber=1:nModels
trainIdx = [c.training(cvSplitNumber);c.training(cvSplitNumber)];
testIdx = [c.test(cvSplitNumber);c.test(cvSplitNumber)];
   
XTrain = s.controlBinnedData(trainIdx,:,:);
YTrain = s.controlLabelSequence(trainIdx,:);

XTest = s.controlBinnedData(testIdx,:,:);
YTest = s.controlLabelSequence(testIdx,:);

XTrain = permute(XTrain,[1,3,2]);
XTest = permute(XTest,[1,3,2]);

tidx = find(testIdx==1);
for i = 1:nbins
    many_Mdls{cvSplitNumber,i} = fitclinear(squeeze(XTrain(:,i,:)),YTrain(:,i),'Learner',c_type);
end

testSplitID(s.controlIdx(testIdx))=cvSplitNumber;
end

predict_vec = zeros(length(s.trialType),nbins);
accuracies_mat = zeros(4,nbins);
correct_mat = zeros(length(s.trialType),nbins);

testTrialTypes = s.trialType;
testTrialTypes(testTrialTypes>2) = testTrialTypes(testTrialTypes>2) - 2;

% Test svms on test data
for trialNumber= 1:length(s.trialType);
if ~isnan(testSplitID(trialNumber));
    for i = 1:nbins
        bin_test = predict(many_Mdls{testSplitID(trialNumber),i},squeeze(tempNeuralData(trialNumber,:,i)));
        predict_vec(trialNumber,i) = bin_test;
        correct_mat(trialNumber,i) = predict_vec(trialNumber,i)==testTrialTypes(trialNumber);
    end
else
    temp1 = [];
    for i = 1:nbins
        for modelNumber = 1:nModels
            bin_test = predict(many_Mdls{modelNumber,i},squeeze(tempNeuralData(trialNumber,:,i)));
            temp1 = [temp1;bin_test==testTrialTypes(trialNumber)];
        end
        correct_mat(trialNumber,i) = squeeze(nanmean(temp1,1));
    end
end
end;

% Calculate classification accuracies
for i = 1:4
    accuracies_mat(i,:) = sum(correct_mat(s.trialType==i & s.correctVec==1,:),1)./sum(s.trialType==i & s.correctVec==1);
end
