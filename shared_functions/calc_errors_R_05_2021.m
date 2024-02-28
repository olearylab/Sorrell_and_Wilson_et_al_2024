function [RMSE, R_square, r_p, test_valid_new, xtest_new, xprediction_new] = calc_errors_R_05_2021(results_struct,xfull,tbt_details,t_types,va_ind,yaw_ind,yaw_offset)
xtest = results_struct.xtest;
xprediction = results_struct.xprediction;
% For getting quick decoding accuracy

if ~isempty(xfull) && ~isempty(tbt_details)
    % correct trials only
    types_vec = [1,4,7,10];
    types_vec = types_vec(t_types);
    
    % keep only correct trials
    trial_num = xfull(12,:);
    cur_trials = find(ismember(tbt_details(3,:),types_vec));
    kept_trials = ismember(trial_num,cur_trials);
    
    test_valid = results_struct.test_valid & kept_trials;
    
    % remove poor trials: trials where magnitude of view angle increases
    % above 2pi.
    for i = 1:length(cur_trials)
        if max(abs(xfull(7,trial_num==cur_trials(i)))) > 2*pi
            test_valid(trial_num==cur_trials(i)) = 0;
        end
    end
    
end


if ~isempty(va_ind)
    xtest(:,va_ind) = wrapToPi(xtest(:,va_ind));
    if ~isempty(yaw_ind)
        % Include yaw to VA check
        ball_va = calc_va_alt_05_2021(xtest(:,yaw_ind),xtest(:,va_ind),yaw_offset,test_valid); 
        decode_yaw_va = calc_va_alt_05_2021(xprediction(:,yaw_ind),xtest(:,va_ind),yaw_offset,test_valid); 

        ball_va = wrapToPi(ball_va);
        decode_yaw_va = wrapToPi(decode_yaw_va);

        % make final value a check that va predicted from yaw matches true va
        xtest = [xtest,ball_va,xtest(:,va_ind)];
        xprediction = [xprediction,decode_yaw_va,ball_va];
    end
end

num_x = size(xtest,2);

% calculate decoder accuracy
RMSE = sqrt(mean((xtest(test_valid,:) - xprediction(test_valid,:)).^2,1));
r_p = zeros(1,num_x);
R_square = zeros(1,num_x);
for i = 1:num_x
    cur_cor = corrcoef(xtest(test_valid,i),xprediction(test_valid,i));
    r_p(i) = cur_cor(1,2);
    R_square(i) = 1 - sum((xtest(test_valid,i) - xprediction(test_valid,i)).^2)/sum((xtest(test_valid,i) - mean(xtest(test_valid,i))).^2);
end

test_valid_new = test_valid;
xtest_new = xtest;
xprediction_new = xprediction;