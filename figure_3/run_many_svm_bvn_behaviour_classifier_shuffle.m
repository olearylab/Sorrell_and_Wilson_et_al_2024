function [accuracies_cell] = run_many_svm_bvn_behaviour_classifier_shuffle(virmen_cell,tbt_cell,nbins,kept_types,balance_types,keep_incorrect,l_r_sep,n_shuffles)

%%
% Set random seed
rng(1);
% Run shuffled svm classification

num_mice = size(virmen_cell,1);
num_days = size(virmen_cell,2);

accuracies_cell = cell(num_mice,num_days);

for m = 1:num_mice
    disp("Mouse "+ m)
    for d = 1:num_days
        if ~isempty(virmen_cell{m,d})
            disp("Day "+ d)
            cur_acc = nan.*ones(n_shuffles,4,nbins);
            % Run for many different shuffles
            for i = 1:n_shuffles
                [ss,centres] = prep_data_for_bvn_lstm_shuffled(virmen_cell{m,d},tbt_cell{m,d},nbins,keep_incorrect,kept_types,balance_types);

                % Additonal option to run separately on left and right trials
                if l_r_sep
                    s = ss;
                    kept_l = ismember(ss.origTrialType,[1,3]);
                    s.binnedVirmenData2 = ss.binnedVirmenData2(kept_l,:,:);
                    s.trialType = ss.trialType(kept_l);
                    s.correctVec = ss.correctVec(kept_l);
                    % Run svms on left trials only
                    [accuracies_l] = svm_behaviour_classifier_15052023(s);

                    s = ss;
                    kept_r = ismember(ss.origTrialType,[2,4]);
                    s.binnedVirmenData2 = ss.binnedVirmenData2(kept_r,:,:);
                    s.trialType = ss.trialType(kept_r);
                    s.correctVec = ss.correctVec(kept_r);
                    % Run svms on right trials only
                    [accuracies_r] = svm_behaviour_classifier_15052023(s);
                    % average of left and right accuracies.
                    accuracies_mat = (accuracies_l + accuracies_r)/2;
                else
                    % Run svms on left and right trials together
                    [accuracies_mat] = svm_behaviour_classifier_15052023(s);

                end
                cur_acc(i,:,:) = accuracies_mat;
                
            end
            accuracies_cell{m,d} = cur_acc;
        end
    end
end