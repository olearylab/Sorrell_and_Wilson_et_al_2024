function [accuracies_cell] = run_many_subs_svm_bvn_neural_classifier(z_cell,virmen_cell,tbt_cell,nbins,kept_types,balance_types,keep_incorrect,l_r_sep,num_subs)
% 02/08/2023
%%
rng(1);

num_mice = size(virmen_cell,1);
num_days = size(virmen_cell,2);

accuracies_cell = cell(num_mice,num_days);

for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(virmen_cell{m,d})
            accuracies_temp = NaN(num_subs,length(kept_types),nbins);
            disp("Running mouse " + m + " day " + d)
            
            zdata = z_cell{m,d};
            % Preprocess once
            initialise_params;
            if ndims(zdata) == 2
                model_params.spatial = false;
                create_reg = false;
            else
                create_reg = true;
            end
            model_params.reg_images = false;

            [zdatafilt, train_mean, Rfixed] = preprocess_images_05_2021(zdata,model_params,[],[],[],create_reg); 
            
            % Run for many subsamples of trials
            for n = 1:num_subs
                % prepare structure for running svms
                [ss,centres] = prep_data_for_neural_bvn_lstm_processed(zdatafilt,virmen_cell{m,d},tbt_cell{m,d},nbins,keep_incorrect,kept_types,balance_types);

                % Additonal option to run separately on left and right trials
                if l_r_sep
                    s = ss;
                    kept_l = ismember(ss.origTrialType,[1,3]);
                    s.binnedNeuralData2 = ss.binnedNeuralData2(kept_l,:,:);
                    s.trialType = ss.trialType(kept_l);
                    s.correctVec = ss.correctVec(kept_l);
                    
                    % Run svm classification on left trials
                    [accuracies_l] = svm_neural_classifier(s);

                    s = ss;
                    kept_r = ismember(ss.origTrialType,[2,4]);
                    s.binnedNeuralData2 = ss.binnedNeuralData2(kept_r,:,:);
                    s.trialType = ss.trialType(kept_r);
                    s.correctVec = ss.correctVec(kept_r);
                    % run svm classification on right trials
                    [accuracies_r] = svm_neural_classifier(s);
                    % average of left and right accuracies.
                    accuracies_mat = (accuracies_l + accuracies_r)/2;
                else
                    % Run svm classification on left and right trials
                    % together
                    [accuracies_mat] = svm_neural_classifier(s);

                end

                accuracies_temp(n,:,:) = accuracies_mat;
            end
            % Only keep average of all subsamples.
            accuracies_cell{m,d} = squeeze(mean(accuracies_temp));
        end
    end
end