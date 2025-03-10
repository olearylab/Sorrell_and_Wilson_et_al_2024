% set model_parameters

% set alpha paramter, proportion of signal that is motor
sim_model.alpha = 0.1;

% set parameters
sim_model.va_std_scale = 3;
sim_model.silent_prob = 0.2;
sim_model.scale_std = 0.2;
sim_model.width_std = 0.2;
sim_model.centre_std = 0.2;

% training lenth
% or include all data for train and test
% sim_model.train_length = 15000;
sim_model.train_trials = 30;
sim_model.use_all = true;

% set neuron constants
sim_model.neurons = 100;
sim_model.useful_neu = 100;
sim_model.repeat_neu = 2;

% set noise std
sim_model.noise_std = 0.1;

% include low pass filter
sim_model.lp = 0.2;