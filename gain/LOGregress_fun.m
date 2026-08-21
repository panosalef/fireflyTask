function   [b_v,bci_v] = LOGregress_fun(path)
%LOGREGRESS_FUN  Log-log regression of travel time on target distance and speed.
%   [b, bci] = LOGregress_fun(path) pools all gain conditions (con_gain) and
%   regresses log(travel time) on log(target distance) and log(mean speed).
%   Returns the coefficients and their 95% confidence intervals.
m = con_gain(path);

%% Regress Velocity all

Yv = [m.split.g1.t_travel' ; m.split.g15.t_travel' ; m.split.g2.t_travel'];

Xv = [[m.split.g1.targ.r' ; m.split.g15.targ.r' ; m.split.g2.targ.r'] ... 
    [m.split.g1.meanspeed'; m.split.g15.meanspeed'; m.split.g2.meanspeed']];

[b_v,bci_v] = regress(log(abs(Yv)),log(abs(Xv)));

end






