function [va] = calc_va_alt_05_2021(yaw_vec,va_vec,offset,test_valid) 
% Function for calculating view angle from ball yaw velocity.
dt = 1/30;
circum = 64;
V = 0.32;
beta = 0.05*circum/V/2.5;

yaw_vec = yaw_vec-offset;
va = zeros(size(va_vec));
va(1) = va_vec(1);
for j = 2:length(yaw_vec)
    if (~test_valid(j-1))&&(test_valid(j))
        va(j) = va_vec(j);
    else
        va(j) = va(j-1) + yaw_vec(j-1)*-beta*dt;
    end
end