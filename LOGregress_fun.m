function   [b_v,bci_v] = LOGregress_fun(path)
m = con_gain(path);

%% Regress Velocity all

Yv = [m.split.g1.t_travel' ; m.split.g15.t_travel' ; m.split.g2.t_travel'];

Xv = [[m.split.g1.targ.r' ; m.split.g15.targ.r' ; m.split.g2.targ.r'] ... 
    [m.split.g1.meanspeed'; m.split.g15.meanspeed'; m.split.g2.meanspeed']];

[b_v,bci_v] = regress(log(abs(Yv)),log(abs(Xv)));

end






