function [all_res_mat] = calc_va_cross_decoding_all(virmen_cell,tbt_cell,virmen_train_cell,tbt_train_cell,nbins,offsets)
% 21/09/2024

% decode corrective angular velocity using online view angle decoder.

num_mice = size(virmen_cell,1);
num_days = size(virmen_cell,2);

% Store all decoding accuracies
all_res_mat = nan.*ones(num_mice,num_days,3,2);

for m = 1:num_mice
    [mean_binned,std_binned] = calculate_mean_ball_va(virmen_train_cell{m},tbt_train_cell{m},nbins,offsets(m,:));
    for d = 1:num_days
        if ~isempty(virmen_cell{m,d})
            virmen_data = virmen_cell{m,d};
            tbt_details = tbt_cell{m,d};
            
            %% Get normalised errors and correcting vec
            % remove invalid data from here on
            ITI = virmen_data(8,:);
            cleaned_valid = clean_valid_data(ITI);
            virmen_data = virmen_data(:,cleaned_valid);
            trial_num = virmen_data(12,:);

            % Calculate heading deviations and whether movements were heading
            % correcting
            [error_mat,correcting_mat] = check_error_correction_norm_samples(virmen_data,tbt_details,nbins,offsets,mean_binned,std_binned);

            %% Keep only samples where movements were corrective

            virmen_data = virmen_data(:,correcting_mat);

            
            %% Calculate decoding accuracies
            test_ITI = virmen_data(8,:);
            [test_valid] = clean_valid_data(test_ITI);
            test_nums = [17,15];
            prediction_nums = [17,17];

            xtrain = virmen_data';
            xtrain = xtrain(:,7);
            xtest = virmen_data';
            xtest = xtest(:,test_nums);
            
            Xprediction = virmen_data';
            Xprediction = Xprediction(:,prediction_nums);

            results_struct.xtest = xtest;
            results_struct.xtrain = xtrain;
            results_struct.xprediction = Xprediction;
            results_struct.test_valid = test_valid;
%            results_struct.Wout = Wout;
%            results_struct.train_mean = train_mean;
%            results_struct.Rfixed = Rfixed;
%            results_struct.model_params = model_params;

            % Use orig_tbt_details as the trial numbers haven't changed
            [RMSE, R_square, r_p, test_valid_new, xtest_new, xprediction_new] = calc_errors_R_05_2021(results_struct,virmen_data,tbt_details,[3,4],[],[],[]);
            all_res_mat(m,d,:,:) = [RMSE;R_square;r_p];
        end
    end
end