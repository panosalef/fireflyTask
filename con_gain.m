function [m] = con_gain(genpath)
d=dir(fullfile(genpath,'*mat'));
fnames={d.name};

%% Initialize
m.all.targ.r = [];          m.all.targ.theta = [];
m.all.resp.r = [];          m.all.resp.theta = [];
m.all.t_stop = [];           m.all.t_move = [];
m.inds.good = [];           m.inds.sesno = [];
m.inds.g1 = [];             m.inds.g15 = [];       m.inds.g2 = [];
m.inds.rew = [];            m.all.meanspeed = [];

for i=1:numel(fnames)
    load(fullfile(genpath,fnames{i}));
    
    t_move = []; t_stop = []; rew = []; meanspeed = [];
    for ntrls = 1: numel(behaviours.trials)
        
        t_move(ntrls) = behaviours.trials(ntrls).events.t_move;
        t_stop(ntrls) = behaviours.trials(ntrls).events.t_stop;
        rew(ntrls) = behaviours.trials(ntrls).logical.reward;
        meanspeed(ntrls) = mean(behaviours.trials(ntrls).continuous.v(behaviours.trials(ntrls).continuous.ts>max(0,behaviours.trials(ntrls).events.t_move) &...
            behaviours.trials(ntrls).continuous.ts<behaviours.trials(ntrls).events.t_end ));
      
    end
    

    m.inds.good = [ m.inds.good behaviours.stats.trialtype.all.trlindx ];
    m.all.meanspeed = [m.all.meanspeed meanspeed];
    m.all.t_move = [m.all.t_move t_move];
    m.all.t_stop = [m.all.t_stop t_stop];
    m.all.targ.r = [ m.all.targ.r behaviours.stats.pos_final.r_targ ];
    m.all.targ.theta = [ m.all.targ.theta behaviours.stats.pos_final.theta_targ ];
    m.all.resp.r = [ m.all.resp.r behaviours.stats.pos_final.r_monk ];
    m.all.resp.theta = [ m.all.resp.theta behaviours.stats.pos_final.theta_monk ];
    m.inds.rew = [m.inds.rew rew];
    m.inds.g1 = [m.inds.g1 behaviours.stats.trialtype.controlgain(1).trlindx];
    m.inds.g15 = [m.inds.g15 behaviours.stats.trialtype.controlgain(2).trlindx];
    m.inds.g2 = [m.inds.g2 behaviours.stats.trialtype.controlgain(3).trlindx];
    m.inds.sesno = [m.inds.sesno repmat(i,[1 numel(behaviours.stats.pos_final.r_targ)])];
    
end
m.all.t_travel = m.all.t_stop - max(0,m.all.t_move);

%% Cartesian kolpa
m.all.targ.x = m.all.targ.r .* sind(m.all.targ.theta);
m.all.targ.y = m.all.targ.r .* cosd(m.all.targ.theta);
m.all.resp.x = m.all.resp.r .* sind(m.all.resp.theta);
m.all.resp.y = m.all.resp.r .* cosd(m.all.resp.theta);
m.all.error = sqrt((m.all.targ.x-m.all.resp.x).^2 + (m.all.targ.y-m.all.resp.y).^2);


%% %%%%% %%
%% SPLIT %%
%% %%%%% %%

%% Densities Gain 1
m.split.g1.targ.r = m.all.targ.r(m.inds.good & m.inds.g1);
m.split.g1.targ.theta = m.all.targ.theta(m.inds.good & m.inds.g1);
m.split.g1.resp.r = m.all.resp.r(m.inds.good & m.inds.g1);
m.split.g1.resp.theta = m.all.resp.theta(m.inds.good & m.inds.g1);
m.split.g1.meanspeed = m.all.meanspeed(m.inds.good & m.inds.g1);

%% Densities Gain 1.5
m.split.g15.targ.r = m.all.targ.r(m.inds.good & m.inds.g15);
m.split.g15.targ.theta = m.all.targ.theta(m.inds.good & m.inds.g15);
m.split.g15.resp.r = m.all.resp.r(m.inds.good & m.inds.g15);
m.split.g15.resp.theta = m.all.resp.theta(m.inds.good & m.inds.g15);
m.split.g15.meanspeed = m.all.meanspeed(m.inds.good & m.inds.g15);

%% Densities Gain 2
m.split.g2.targ.r = m.all.targ.r(m.inds.good & m.inds.g2);
m.split.g2.targ.theta = m.all.targ.theta(m.inds.good & m.inds.g2);

m.split.g2.resp.r = m.all.resp.r(m.inds.good & m.inds.g2);
m.split.g2.resp.theta = m.all.resp.theta(m.inds.good & m.inds.g2);
m.split.g2.meanspeed = m.all.meanspeed(m.inds.good & m.inds.g2);

%% Session Numbers
m.split.g1.sesno = m.inds.sesno(m.inds.good & m.inds.g1);
m.split.g15.sesno = m.inds.sesno(m.inds.good & m.inds.g15);
m.split.g2.sesno = m.inds.sesno(m.inds.good & m.inds.g2);

%% travel time
m.split.g1.t_travel = m.all.t_travel(m.inds.good & m.inds.g1);
m.split.g15.t_travel = m.all.t_travel(m.inds.good & m.inds.g15);
m.split.g2.t_travel = m.all.t_travel(m.inds.good & m.inds.g2);

%% Densities Gain 1
m.rew.g1.targ.r = m.all.targ.r(m.inds.good & m.inds.g1 & m.inds.rew);
m.rew.g1.targ.theta = m.all.targ.theta(m.inds.good & m.inds.g1 & m.inds.rew);
m.rew.g1.resp.r = m.all.resp.r(m.inds.good & m.inds.g1 & m.inds.rew);
m.rew.g1.resp.theta = m.all.resp.theta(m.inds.good & m.inds.g1 & m.inds.rew);

%% Densities Gain 1.5
m.rew.g15.targ.r = m.all.targ.r(m.inds.good & m.inds.g15 & m.inds.rew);
m.rew.g15.targ.theta = m.all.targ.theta(m.inds.good & m.inds.g15 & m.inds.rew);
m.rew.g15.resp.r = m.all.resp.r(m.inds.good & m.inds.g15 & m.inds.rew);
m.rew.g15.resp.theta = m.all.resp.theta(m.inds.good & m.inds.g15 & m.inds.rew);


%% Densities Gain 2
m.rew.g2.targ.r = m.all.targ.r(m.inds.good & m.inds.g2 & m.inds.rew);
m.rew.g2.targ.theta = m.all.targ.theta(m.inds.good & m.inds.g2 & m.inds.rew);

m.rew.g2.resp.r = m.all.resp.r(m.inds.good & m.inds.g2 & m.inds.rew);
m.rew.g2.resp.theta = m.all.resp.theta(m.inds.good & m.inds.g2 & m.inds.rew);

%% travel time
m.rew.g1.t_travel = m.all.t_travel(m.inds.good & m.inds.g1 & m.inds.rew);
m.rew.g15.t_travel = m.all.t_travel(m.inds.good & m.inds.g15 & m.inds.rew);
m.rew.g2.t_travel = m.all.t_travel(m.inds.good & m.inds.g2 & m.inds.rew);

%% travel time
m.rew.g1.meanspeed = m.all.meanspeed(m.inds.good & m.inds.g1 & m.inds.rew);
m.rew.g15.meanspeed = m.all.meanspeed(m.inds.good & m.inds.g15 & m.inds.rew);
m.rew.g2.meanspeed = m.all.meanspeed(m.inds.good & m.inds.g2 & m.inds.rew);



%% Keep only UNrewarded trials %%

%% Densities Gain 1
m.no_rew.g1.targ.r = m.all.targ.r(m.inds.good & m.inds.g1 & ~m.inds.rew);
m.no_rew.g1.targ.theta = m.all.targ.theta(m.inds.good & m.inds.g1 & ~m.inds.rew);

m.no_rew.g1.resp.r = m.all.resp.r(m.inds.good & m.inds.g1 & ~m.inds.rew);
m.no_rew.g1.resp.theta = m.all.resp.theta(m.inds.good & m.inds.g1 & ~m.inds.rew);


%% Densities Gain 1.5
m.no_rew.g15.targ.r = m.all.targ.r(m.inds.good & m.inds.g15 & ~m.inds.rew);
m.no_rew.g15.targ.theta = m.all.targ.theta(m.inds.good & m.inds.g15 & ~m.inds.rew);
m.no_rew.g15.resp.r = m.all.resp.r(m.inds.good & m.inds.g15 & ~m.inds.rew);
m.no_rew.g15.resp.theta = m.all.resp.theta(m.inds.good & m.inds.g15 & ~m.inds.rew);


%% Densities Gain 2
m.no_rew.g2.targ.r = m.all.targ.r(m.inds.good & m.inds.g2 & ~m.inds.rew);
m.no_rew.g2.targ.theta = m.all.targ.theta(m.inds.good & m.inds.g2 & ~m.inds.rew);
m.no_rew.g2.resp.r = m.all.resp.r(m.inds.good & m.inds.g2 & ~m.inds.rew);
m.no_rew.g2.resp.theta = m.all.resp.theta(m.inds.good & m.inds.g2 & ~m.inds.rew);

%% travel time
m.no_rew.g1.t_travel = m.all.t_travel(m.inds.good & m.inds.g1 & ~m.inds.rew);
m.no_rew.g15.t_travel = m.all.t_travel(m.inds.good & m.inds.g15 & ~m.inds.rew);
m.no_rew.g2.t_travel = m.all.t_travel(m.inds.good & m.inds.g2 & ~m.inds.rew);

%% travel time
m.no_rew.g1.meanspeed = m.all.meanspeed(m.inds.good & m.inds.g1 & ~m.inds.rew);
m.no_rew.g15.meanspeed = m.all.meanspeed(m.inds.good & m.inds.g15 & ~m.inds.rew);
m.no_rew.g2.meanspeed = m.all.meanspeed(m.inds.good & m.inds.g2 & ~m.inds.rew);


%% Convert to Cartesian



end