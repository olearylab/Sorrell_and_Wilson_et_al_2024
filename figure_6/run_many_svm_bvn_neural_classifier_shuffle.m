function [accuracies_cell] = run_many_svm_bvn_neural_classifier_shuffle(z_cell,virmen_cell,tbt_cell,nbins,kept_types,balance_types,keep_incorrect,l_r_sep,n_shuffles)
% 02/08/2023
%%
rng(1);

num_mice = size(virmen_cell,1);
num_days = size(virmen_cell,2);

accuracies_cell = cell(num_mice,num_days);

initialise_params;

model_params.spatial = false;
create_reg = false;

model_params.reg_images = false;

for m = 1:num_mice
    disp("Mouse " + m)
    for d = 1:num_days
        if ~isempty(virmen_cell{m,d})
            disp("Day " + d)
            cur_acc = nan.*ones(n_shuffles,4,nbins);
            
            % Preprocess neural data
            [zdatafilt, train_mean, Rfixed] = preprocess_images_05_2021(z_cell{m,d},model_params,[],[],[],create_reg);
            
            for i = 1:n_shuffles
                % Prepare structure for svms
                [ss,centres] = prep_data_for_neural_bvn_shuffle(zdatafilt,virmen_cell{m,d},tbt_cell{m,d},nbins,keep_incorrect,kept_types,balance_types);

                % Additonal option to run separately on left and right trials
                if l_r_sep
                    s = ss;
                    kept_l = ismember(ss.origTrialType,[1,3]);
                    s.binnedNeuralData2 = ss.binnedNeuralData2(kept_l,:,:);
                    s.trialType = ss.trialType(kept_l);
                    s.correctVec = ss.correctVec(kept_l);
                    
                    % Run svms on left trials
                    [accuracies_l] = svm_neural_classifier(s);

                    s = ss;
                    kept_r = ismember(ss.origTrialType,[2,4]);
                    s.binnedNeuralData2 = ss.binnedNeuralData2(kept_r,:,:);
                    s.trialType = ss.trialType(kept_r);
                    s.correctVec = ss.correctVec(kept_r);
                    % Run svms on right trials
                    [accuracies_r] = svm_neural_classifier(s);
                    % average of left and right accuracies.
                    accuracies_mat = (accuracies_l + accuracies_r)/2;
                else
                    % run svms on left and right trials together
                    [accuracies_mat] = svm_neural_classifier(s);

                end
                
                cur_acc(i,:,:) = accuracies_mat;
                
            end

            accuracies_cell{m,d} = cur_acc;
        end
    end
end