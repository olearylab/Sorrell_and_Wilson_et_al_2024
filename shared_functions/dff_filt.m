function [zfilt] = dff_filt(zdata, afast, aslow)
% Perform online pseudo-DF/F filtering
    persistent first zfprev zsprev
    if isempty(first)
        first = false;
        zfprev = zdata;
        zsprev = zdata;
    end
    zfast = afast*zdata + (1-afast)*zfprev;
    zslow = aslow*zdata + (1-aslow)*zsprev;
    
    zfilt = (zfast - zslow)./zslow;
    zfilt(isnan(zfilt))=0;
    
    zfprev = zfast;
    zsprev = zslow;
end