function [rel_results, va_rel] = basic_va_control_sim_23062023(virmen_data, tbt_details, sim_model,alpha)
% 08/05/2024

% set alpha = 0 (purely sensory) or 1 (purely step ahead).

% Run va control simulation using separate motor and sensory neurons

% Simplify to remove noise and other complicating factors.
sim_model.useful_neu = 10;
sim_model.alpha = alpha;
sim_model.neurons = 20;
sim_model.repeat_neu = 1;

types_vec = [1,4];
%% Set up simulation

num_motor = floor(sim_model.alpha*sim_model.useful_neu);

motor_neu = randperm(sim_model.useful_neu,num_motor);

% set constants
silent_prob = sim_model.silent_prob;
scale_std = sim_model.scale_std;
width_std = sim_model.width_std;
centre_std = sim_model.centre_std;
va_std_scale = sim_model.va_std_scale;
% train_length = sim_model.train_length;
train_trials = sim_model.train_trials;
ITI = virmen_data(8,:);
cleaned_valid = clean_valid_data(ITI);
virmen_data = virmen_data(:,cleaned_valid);
trial_num = virmen_data(12,:);
% wrap va to pi

test_valid = virmen_data(8,:);

cur_trials = find(ismember(tbt_details(3,:),types_vec));

new_tnum = trial_num;
% Keep only trials that don't rotate past pi and face the right direction.
ii = 1;
kept_trials = [];
for i = 1:length(cur_trials)
    if max(abs(virmen_data(7,trial_num==cur_trials(i)))) < pi
        test_valid(trial_num==cur_trials(i)) = 1;
        new_tnum(trial_num==cur_trials(i)) = ii;
        kept_trials = [kept_trials,i];
        ii = ii+1;
    end
end

tbt_details = tbt_details(:,kept_trials);
virmen_data(12,:) = new_tnum;
virmen_data = virmen_data(:,test_valid==1);
virmen_data(7,:) = wrapToPi(virmen_data(7,:));
trial_num = virmen_data(12,:);

% set up training and testing
if sim_model.use_all
    train_length = 0;
    train_data = virmen_data(7,:);
    test_data = virmen_data(7,:);

    test_trial_num = trial_num;
else
    train_length = find(trial_num==train_trials+1,1)-1;
    train_data = virmen_data(7,1:train_length);
    test_data = virmen_data(7,train_length+1:end);

    test_trial_num = trial_num(train_length+1:end);
end

fs = 30;

dt = 1/fs;
t = linspace(0,(length(train_data)-1)*dt,length(train_data));
t_test = linspace(0,(length(test_data)-1)*dt,length(test_data));

p = train_data;

% Encoding
% Each neurons activity is a gaussian bump function, with some noise added.
% Assume a constant noise for both groups.
% For non-variable group, this is it.
% For variable group, introduce also variability in peak centre, peak
% height, and a silencing probability
% Standard deviation set to 3 times distance between means
% random noise drawn from gaussian with 0 mean, variance set to ??
neurons = sim_model.neurons;
useful_neu = sim_model.useful_neu;
repeat_neu = sim_model.repeat_neu;
va_means = [];
for i = 1:repeat_neu
    va_means = [va_means,linspace(-pi,pi,useful_neu/repeat_neu)];
end
va_std = va_std_scale*(va_means(2)-va_means(1));

% NO noise
noise_std = sim_model.noise_std;

% Training neural data
n_rel_train = zeros(length(t),neurons);
n_unrel_train = zeros(length(t),neurons);
% Testing neural data
n_rel_test = zeros(length(t_test),neurons);
n_unrel_test = zeros(length(t_test),neurons);

% Create training and testing neural vectors. Using same positional
% trajectories but recreating neural vectores for each.
% mag_scale = normrnd(ones(1,length(va_means)),scale_std);
% width_rand = normrnd(va_std*ones(1,length(va_means)),width_std);
% centre_rand = normrnd(va_means,centre_std);
% silent = binornd(1,silent_prob,[1,useful_neu]);

% mag_scale_test = normrnd(ones(1,length(va_means)),scale_std);
% width_rand_test = normrnd(va_std*ones(1,length(va_means)),width_std);
% centre_rand_test = normrnd(va_means,centre_std);
% silent_test = binornd(1,silent_prob,[1,useful_neu]);
trial_prev = trial_num(1);
for i = 2:length(train_data)
%     if trial_num(i-1) ~= trial_prev
%         mag_scale = normrnd(ones(1,length(va_means)),scale_std);
%         width_rand = normrnd(va_std*ones(1,length(va_means)),width_std);
%         centre_rand = normrnd(va_means,centre_std);
%         silent = binornd(1,silent_prob,[1,useful_neu]);
%         
%         mag_scale_test = normrnd(ones(1,length(va_means)),scale_std);
%         width_rand_test = normrnd(va_std*ones(1,length(va_means)),width_std);
%         centre_rand_test = normrnd(va_means,centre_std);
%         silent_test = binornd(1,silent_prob,[1,useful_neu]);
%     end
    
    % sensory
    % n_rel_train(i-1,1:length(va_means)) = normpdf((p(i-1)-va_means)/va_std) + normrnd(0,noise_std,[1,useful_neu]);
    n_rel_train(i-1,1:length(va_means)) = normpdf((p(i-1)-va_means)/va_std);
    % motor
    % n_rel_train(i-1,motor_neu) = normpdf((p(i)-va_means(motor_neu))/va_std) + normrnd(0,noise_std,[1,num_motor]);
    n_rel_train(i-1,motor_neu) = normpdf((p(i)-va_means(motor_neu))/va_std);
    
    % sensory
    % n_unrel_train(i-1,1:length(va_means)) = (1-silent).*mag_scale.*normpdf((p(i-1)-centre_rand)./width_rand) + normrnd(0,noise_std,[1,useful_neu]);
    
    % motor
    % n_unrel_train(i-1,motor_neu) = (1-silent(motor_neu)).*mag_scale(motor_neu).*normpdf((p(i)-centre_rand(motor_neu))./width_rand(motor_neu)) + normrnd(0,noise_std,[1,num_motor]);
    
%     if useful_neu<neurons
%         n_rel_train(i,length(va_means)+1:end) = normrnd(0,noise_std,[1,neurons-useful_neu]);
%         n_unrel_train(i,length(va_means)+1:end) = normrnd(0,noise_std,[1,neurons-useful_neu]);
%     end
    % n_unrel_test(i,:) = (1-silent_test).*mag_scale_test.*normpdf((p(i)-centre_rand_test)./width_rand_test) + normrnd(0,noise_std,[1,neurons]);
    trial_prev = trial_num(i-1);
end

%% repeat for testing - open loop
% mag_scale = normrnd(ones(1,length(va_means)),scale_std);
% width_rand = normrnd(va_std*ones(1,length(va_means)),width_std);
% centre_rand = normrnd(va_means,centre_std);
% silent = binornd(1,silent_prob,[1,useful_neu]);
% 
% mag_scale_test = normrnd(ones(1,length(va_means)),scale_std);
% width_rand_test = normrnd(va_std*ones(1,length(va_means)),width_std);
% centre_rand_test = normrnd(va_means,centre_std);
% silent_test = binornd(1,silent_prob,[1,useful_neu]);
% Is this ok?
trial_prev = test_trial_num(1);
for i = 2:length(test_data)
%     if trial_num(train_length+i-1) ~= trial_prev
%         mag_scale = normrnd(ones(1,length(va_means)),scale_std);
%         width_rand = normrnd(va_std*ones(1,length(va_means)),width_std);
%         centre_rand = normrnd(va_means,centre_std);
%         silent = binornd(1,silent_prob,[1,useful_neu]);
%         
%         mag_scale_test = normrnd(ones(1,length(va_means)),scale_std);
%         width_rand_test = normrnd(va_std*ones(1,length(va_means)),width_std);
%         centre_rand_test = normrnd(va_means,centre_std);
%         silent_test = binornd(1,silent_prob,[1,useful_neu]);
%     end
    
    % sensory
    % n_rel_test(i-1,1:length(va_means)) = normpdf((test_data(i-1)-va_means)/va_std) + normrnd(0,noise_std,[1,useful_neu]);
    n_rel_test(i-1,1:length(va_means)) = normpdf((test_data(i-1)-va_means)/va_std);
    % motor
    % n_rel_test(i-1,motor_neu) = normpdf((test_data(i)-va_means(motor_neu))/va_std) + normrnd(0,noise_std,[1,num_motor]);
    n_rel_test(i-1,motor_neu) = normpdf((test_data(i)-va_means(motor_neu))/va_std);
    % sensory
    % n_unrel_test(i-1,1:length(va_means)) = (1-silent_test).*mag_scale_test.*normpdf((test_data(i-1)-centre_rand_test)./width_rand_test) + normrnd(0,noise_std,[1,useful_neu]);
    
    % motor
    % n_unrel_test(i-1,motor_neu) = (1-silent_test(motor_neu)).*mag_scale_test(motor_neu).*normpdf((test_data(i)-centre_rand_test(motor_neu))./width_rand_test(motor_neu)) + normrnd(0,noise_std,[1,num_motor]);
    
%     if useful_neu<neurons
%         n_rel_test(i,length(va_means)+1:end) = normrnd(0,noise_std,[1,neurons-useful_neu]);
%         n_unrel_test(i,length(va_means)+1:end) = normrnd(0,noise_std,[1,neurons-useful_neu]);
%     end
    trial_prev = test_trial_num(i-1);
end

%% Train decoders. 
% NB: I probably need some kind of neural netowrk to really
% test my hypotheses (as I am attempting to model brain, so linreg is
% probably not enough).

% If want to use inbuilt RLS - using SGD
% W = zeros(neurons+1,2);
% mdl = fitrlinear(n_rel_train,p,'OptimizeHyperparameters','auto');
% W(:,1) = [mdl.Beta;mdl.Bias];
% mdl = fitrlinear(n_unrel_train,p,'OptimizeHyperparameters','auto');
% W(:,2) = [mdl.Beta;mdl.Bias];

% Instead use OLS - might fail with large number of neurons?
W = zeros(neurons+1,2);
W(:,1) = [n_rel_train(1:end-1,:),ones(size(n_rel_train,1)-1,1)]\p(1:end-1)';
% W(:,2) = [n_unrel_train(1:end-1,:),ones(size(n_unrel_train,1)-1,1)]\p(1:end-1)';

rel_train_decode = [n_rel_train,ones(size(n_rel_train,1),1)]*W(:,1);
% unrel_train_decode = [n_unrel_train,ones(size(n_unrel_train,1),1)]*W(:,2);

% Test decoders
rel_test_decode = [n_rel_test,ones(size(n_rel_test,1),1)]*W(:,1);
% unrel_test_decode = [n_unrel_test,ones(size(n_unrel_test,1),1)]*W(:,2);
% rel_test_decode_unrel = [n_unrel_test,ones(size(n_unrel_test,1),1)]*W(:,1);
% unrel_test_decode_rel = [n_rel_test,ones(size(n_rel_test,1),1)]*W(:,2);

%% closed loop testing

% entirely sensory has no control signal
% simulate va control where va is purely sensory
% simulate va control when it is predictive? - what is control signal?
% Where I want my next view angle to be I guess? Should work on converting
% this into a control signal e.g. via PID control. 
% simulate some combination of these - alpha
va_rel = zeros(length(test_data),1);
% va_unrel = zeros(length(test_data),1);
alpha = sim_model.alpha;
lp = sim_model.lp;

% mag_scale = normrnd(ones(1,length(va_means)),scale_std);
% width_rand = normrnd(va_std*ones(1,length(va_means)),width_std);
% centre_rand = normrnd(va_means,centre_std);
% silent = binornd(1,silent_prob,[1,useful_neu]);
% 
% mag_scale_test = normrnd(ones(1,length(va_means)),scale_std);
% width_rand_test = normrnd(va_std*ones(1,length(va_means)),width_std);
% centre_rand_test = normrnd(va_means,centre_std);
% silent_test = binornd(1,silent_prob,[1,useful_neu]);
% Is this ok?
trial_prev = test_trial_num(1);
for i = 2:length(test_data)
    if trial_num(train_length + i-1) ~= trial_prev
%         mag_scale = normrnd(ones(1,length(va_means)),scale_std);
%         width_rand = normrnd(va_std*ones(1,length(va_means)),width_std);
%         centre_rand = normrnd(va_means,centre_std);
%         silent = binornd(1,silent_prob,[1,useful_neu]);
%         
%         mag_scale_test = normrnd(ones(1,length(va_means)),scale_std);
%         width_rand_test = normrnd(va_std*ones(1,length(va_means)),width_std);
%         centre_rand_test = normrnd(va_means,centre_std);
%         silent_test = binornd(1,silent_prob,[1,useful_neu]);
        
        % reset view angle
        va_rel(i-1) = 0;
        va_unrel(i-1) = 0;
    end
   
    % sensory
    % n_rel_test(i-1,1:length(va_means)) = normpdf((va_rel(i-1)-va_means)/va_std) + normrnd(0,noise_std,[1,useful_neu]);
    n_rel_test(i-1,1:length(va_means)) = normpdf((va_rel(i-1)-va_means)/va_std);
    % motor
    % n_rel_test(i-1,motor_neu) = normpdf((test_data(i)-va_means(motor_neu))/va_std) + normrnd(0,noise_std,[1,num_motor]);
    n_rel_test(i-1,motor_neu) = normpdf((test_data(i)-va_means(motor_neu))/va_std);
    
    % sensory
    % n_unrel_test(i-1,1:length(va_means)) = (1-silent_test).*mag_scale_test.*normpdf((va_unrel(i-1)-centre_rand_test)./width_rand_test) + normrnd(0,noise_std,[1,useful_neu]);
    
    % motor
    % n_unrel_test(i-1,motor_neu) = (1-silent_test(motor_neu)).*mag_scale_test(motor_neu).*normpdf((test_data(i)-centre_rand_test(motor_neu))./width_rand_test(motor_neu)) + normrnd(0,noise_std,[1,num_motor]);
    
%     if useful_neu<neurons
%         n_rel_test(i,length(va_means)+1:end) = normrnd(0,noise_std,[1,neurons-useful_neu]);
%         n_unrel_test(i,length(va_means)+1:end) = normrnd(0,noise_std,[1,neurons-useful_neu]);
%     end
    trial_prev = test_trial_num(i-1);
    
    % va_rel(i) = lp*[n_rel_test(i-1,:),1]*W(:,1) + (1-lp)*va_rel(i-1);
    va_rel(i) = [n_rel_test(i-1,:),1]*W(:,1);
    % va_unrel(i) = lp*[n_unrel_test(i-1,:),1]*W(:,2) + (1-lp)*va_unrel(i-1);
end    
% 

%% Calculate results

num_left = sum(tbt_details(3,test_trial_num(1):max(test_trial_num)) == types_vec(1));
num_right = sum(tbt_details(3,test_trial_num(1):max(test_trial_num)) == types_vec(2));
rel_cor = [0,0];
for i = test_trial_num(1):max(test_trial_num)
    if tbt_details(3,i) == types_vec(1)
        cur_rel = va_rel(test_trial_num==i);
        if cur_rel(end) > 0
            rel_cor(1) = rel_cor(1) + 1;
        end
%         cur_unrel = va_unrel(test_trial_num==i);
%         if cur_unrel(end) > 0
%             unrel_cor(1) = unrel_cor(1) + 1;
%         end
    elseif tbt_details(3,i) == types_vec(2)
        cur_rel = va_rel(test_trial_num==i);
        if cur_rel(end) < 0
            rel_cor(2) = rel_cor(2) + 1;
        end
%         cur_unrel = va_unrel(test_trial_num==i);
%         if cur_unrel(end) < 0
%             unrel_cor(2) = unrel_cor(2) + 1;
%         end
    end
end     

% Calculate proportion "correct" for each direction
rel_results = sum(rel_cor)/sum([num_left,num_right]);
% unrel_results = sum(unrel_cor)/sum([num_left,num_right]);

%% Plot

% Tuning Curves
figure
p_space = linspace(-pi,pi,100);
for n = 1:sim_model.useful_neu
    plot(p_space,normpdf(p_space-va_means(n)/va_std),'LineWidth',2);
    hold on
end
title("View Angle Tuning Curves")
ylabel("Neural Activity")
xlabel("View Angle (rad)")

% Example open loop and closed loop decoding
colour_vec = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];
num_examples = 3;
figure
l_trials = find(tbt_details(3,:)==types_vec(1));
r_trials = find(tbt_details(3,:)==types_vec(2));

subplot(2,2,1)
hold on
title("Left Trials")
ylabel("View Angle (rad)")
yline(0,'--','LineWidth',2,'Color','k');

subplot(2,2,2)
hold on
title("Right Trials")
yline(0,'--','LineWidth',2,'Color','k');

subplot(2,2,3)
hold on
ylabel("View Angle (rad)")
yline(0,'--','LineWidth',2,'Color','k');
xlabel("Sample")

subplot(2,2,4)
hold on
ylabel("View Angle (rad)")
yline(0,'--','LineWidth',2,'Color','k');
xlabel("Sample")

for i = 1:num_examples
    subplot(2,2,1)
    plot(test_data(test_trial_num==l_trials(i)),'LineWidth',2,'Color',[0.5,0.5,0.5])
    plot(rel_test_decode(test_trial_num==l_trials(i)),'--','LineWidth',2,'Color',colour_vec(1,:))
    
    subplot(2,2,2)
    plot(test_data(test_trial_num==r_trials(i)),'LineWidth',2,'Color',[0.5,0.5,0.5])
    plot(rel_test_decode(test_trial_num==r_trials(i)),'--','LineWidth',2,'Color',colour_vec(1,:))
    
    subplot(2,2,3)
    plot(test_data(test_trial_num==l_trials(i)),'LineWidth',2,'Color',[0.5,0.5,0.5])
    plot(va_rel(test_trial_num==l_trials(i)),'LineWidth',2,'Color',colour_vec(2,:))
    
    subplot(2,2,4)
    plot(test_data(test_trial_num==r_trials(i)),'LineWidth',2,'Color',[0.5,0.5,0.5])
    plot(va_rel(test_trial_num==r_trials(i)),'LineWidth',2,'Color',colour_vec(2,:))
end

% Example trajectories
figure
subplot(2,2,1)
hold on
title("Left Trial")
ylabel("View Angle (rad)")
xlabel("Sample")
plot(test_data(test_trial_num==l_trials(1)),'LineWidth',2,'Color',[0.5,0.5,0.5])
yline(0,'--','LineWidth',2,'Color','k');

subplot(2,2,2)
hold on
title("Right Trial")
xlabel("Sample")
plot(test_data(test_trial_num==r_trials(1)),'LineWidth',2,'Color',[0.5,0.5,0.5])
yline(0,'--','LineWidth',2,'Color','k');