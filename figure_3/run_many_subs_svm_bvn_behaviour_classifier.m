function [accuracies_cell] = run_many_subs_svm_bvn_behaviour_classifier(virmen_cell,tbt_cell,nbins,kept_types,balance_types,keep_incorrect,l_r_sep,num_subs)

%%
% Set random seed
rng(1);

num_mice = size(virmen_cell,1);
num_days = size(virmen_cell,2);

accuracies_cell = cell(num_mice,num_days);

for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(virmen_cell{m,d})
            disp("Running mouse " + m + " day " + d)
            accuracies_temp = NaN(num_subs,length(kept_types),nbins);
            % run for many subsamples of trials
            for n = 1:num_subs
                % set up data structure for running svms 
                [ss,centres] = prep_data_for_bvn_lstm(virmen_cell{m,d},tbt_cell{m,d},nbins,keep_incorrect,kept_types,balance_types);

                % Additonal option to run separately on left and right trials
                if l_r_sep
                    s = ss;
                    kept_l = ismember(ss.origTrialType,[1,3]);
                    s.binnedVirmenData2 = ss.binnedVirmenData2(kept_l,:,:);
                    s.trialType = ss.trialType(kept_l);
                    s.correctVec = ss.correctVec(kept_l);
                    % run svms on left trials
                    [accuracies_l] = svm_behaviour_classifier_15052023(s);

                    s = ss;
                    kept_r = ismember(ss.origTrialType,[2,4]);
                    s.binnedVirmenData2 = ss.binnedVirmenData2(kept_r,:,:);
                    s.trialType = ss.trialType(kept_r);
                    s.correctVec = ss.correctVec(kept_r);
                    % run svm on left trials
                    [accuracies_r] = svm_behaviour_classifier_15052023(s);
                    % average of left and right accuracies.
                    accuracies_mat = (accuracies_l + accuracies_r)/2;
                else
                    % run on both left and right trials together
                    [accuracies_mat] = svm_behaviour_classifier_15052023(s);

                end
                accuracies_temp(n,:,:) = accuracies_mat;
                
            end
            % Only keep average of all subsamples.
            accuracies_cell{m,d} = squeeze(mean(accuracies_temp));
        end
    end
end