function [xpred]  = decoder_step(zin, W0)
% linear decoding for single sample
    xpred = zin*W0;
    end