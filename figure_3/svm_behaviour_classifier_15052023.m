function [accuracies_mat] = svm_behaviour_classifier_15052023(s)
% 15/05/2023

% Updated bin specific SVM classification of behaviour.

% zscore input data
tempVirmenData = s.binnedVirmenData2(:,[13,15],:);
tempVirmenData(:,1,:) = zscore(tempVirmenData(:,1,:),0,'all');
tempVirmenData(:,2,:) = zscore(tempVirmenData(:,2,:),0,'all');

nbins = size(tempVirmenData,3);

c_type = 'svm';

% set up partitions
% Train only on correct trials.
s.controlIdx = find(s.trialType<3);
if s.correct_only
    s.controlIdx = s.controlIdx(s.correctVec(s.controlIdx)==1);
end
s.controlBinnedData = tempVirmenData(s.controlIdx,:,:);
s.controlLabels = s.trialType(s.controlIdx);

s.controlLabelSequence = zeros(size(s.controlBinnedData,1),size(s.controlBinnedData,3));

for controlTrialNumber = 1:size(s.controlBinnedData,1);
    s.controlLabelSequence(controlTrialNumber,:) = categorical(repmat(s.controlLabels(controlTrialNumber),[1, size(tempVirmenData,3)]),[1 2]);
end

nModels = 5; % number of cv partitions

c = cvpartition(size(s.controlLabelSequence,1)/2,'KFold',nModels);
testSplitID = NaN(size(s.trialType));

many_Mdls = cell(nModels,nbins);

%% Train network
% Run for 5-fold cross-validation
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
    % In testing now
    % many_Mdls{cvSplitNumber,i} = fitclinear(squeeze(XTrain(:,i,:)),YTrain(:,i),'Learner',c_type,'Lambda',0);
    many_Mdls{cvSplitNumber,i} = fitclinear(squeeze(XTrain(:,i,:)),YTrain(:,i),'Learner',c_type);
    % many_Mdls{cvSplitNumber,i} = fitcsvm(squeeze(XTrain(:,i,:)),YTrain(:,i)); %,'KernelFunction','rbf');
end

testSplitID(s.controlIdx(testIdx))=cvSplitNumber;
end

predict_vec = zeros(length(s.trialType),nbins);
accuracies_mat = zeros(4,nbins);
correct_mat = zeros(length(s.trialType),nbins);

testTrialTypes = s.trialType;
testTrialTypes(testTrialTypes>2) = testTrialTypes(testTrialTypes>2) - 2;

% Test SVMs on test data for each cross-validation
for trialNumber= 1:length(s.trialType);
if ~isnan(testSplitID(trialNumber));
    for i = 1:nbins
        bin_test = predict(many_Mdls{testSplitID(trialNumber),i},squeeze(s.controlBinnedData(trialNumber,:,i)));
        predict_vec(trialNumber,i) = bin_test;
        correct_mat(trialNumber,i) = predict_vec(trialNumber,i)==testTrialTypes(trialNumber);
    end
else
    temp1 = [];
    for i = 1:nbins
        for modelNumber = 1:nModels
            bin_test = predict(many_Mdls{modelNumber,i},squeeze(tempVirmenData(trialNumber,:,i)));
            temp1 = [temp1;bin_test==testTrialTypes(trialNumber)];
        end
        correct_mat(trialNumber,i) = squeeze(nanmean(temp1,1));
    end
end
end;


for i = 1:4
    accuracies_mat(i,:) = sum(correct_mat(s.trialType==i & s.correctVec==1,:),1)./sum(s.trialType==i & s.correctVec==1);
end