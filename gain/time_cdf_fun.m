function [cdf,nx] = time_cdf_fun(path)
%TIME_CDF_FUN  Cumulative distribution of travel time for each joystick gain.
%   [cdf, nx] = time_cdf_fun(path) pools the sessions in path with con_gain
%   and returns the CDF of movement duration (0-4 s grid nx) for gains 1, 1.5
%   and 2 in cdf.g1, cdf.g15, cdf.g2.

m = con_gain(path);
%% Gain 1
nx = linspace(0,4,100);
[ny_g1,~] =  hist(m.split.g1.t_travel,nx);
cdf.g1 = cumsum(ny_g1(:))./sum(ny_g1);

%% Gain 1.5
[ny_g15,~] =  hist(m.split.g15.t_travel,nx);
cdf.g15 = cumsum(ny_g15(:))./sum(ny_g15);

%% Gain 2
[ny_g2,~] =  hist(m.split.g2.t_travel,nx);
cdf.g2 = cumsum(ny_g2(:))./sum(ny_g2);

end