function [cdf,nx] = time_cdf_fun(path)

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