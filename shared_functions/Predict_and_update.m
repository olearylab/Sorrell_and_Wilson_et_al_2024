function [xpred, Wout]  = Predict_and_update(zin, xin, W0, lrate, valid, update)
% Run LMS training update step
    persistent Win
    if isempty(Win)
        Win = W0;
    end
    xpred = zin*Win;
    if valid && update
        error = xin - xpred;
        dW = lrate*zin'*error;
        Wout = Win + dW;
        Win = Wout;
    else
        Wout = Win;
    end