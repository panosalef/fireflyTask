function [m,m2] = concat_basisF2(genpath)
d=dir(fullfile(genpath,'*mat'));
fnames={d.name};

m.ptb_trials = struct('continuous',([]),'events',([]),'logical',([]),'prs',([]));

m.h_index = [];
m.targ.no_ptb.r = []; m.targ.ptb.r = []; m.targ.no_ptb.theta = []; m.targ.ptb.theta = [];
m.resp.no_ptb.r = []; m.resp.ptb.r = []; m.resp.no_ptb.theta = []; m.resp.ptb.theta = [];
m.theta_ptb = [];           % theta at the beginning of the perturbation
m.rew = [];
for i=1:numel(fnames)
    
    load(fullfile(genpath,fnames{i}));
    for j=1:numel(behaviours.trials);reward(j)= behaviours.trials(j).logical.reward;end
    theta_ptb = add_theta_ptb(behaviours);
    h_index = add_hindex1(behaviours);
    good_indx = behaviours.stats.trialtype.all.trlindx;
    ptb_indx = behaviours.stats.trialtype.ptb(2).trlindx;
    
    behaviours = add_ptbn2(behaviours);
    behaviours = add_continuous(behaviours);
    

    
    ptb_trials = behaviours.trials(good_indx & ptb_indx);
    behaviours.trials = behaviours.trials(good_indx);
    
    m.ptb_trials((end+1):(end+numel(ptb_trials))) = ptb_trials;
    m.targ.ptb.r = [m.targ.ptb.r behaviours.stats.pos_final.r_targ(good_indx & ptb_indx)];
    m.targ.ptb.theta = [m.targ.ptb.theta behaviours.stats.pos_final.theta_targ(good_indx & ptb_indx)];
    m.resp.ptb.r = [m.resp.ptb.r behaviours.stats.pos_final.r_monk(good_indx & ptb_indx)];
    m.resp.ptb.theta = [m.resp.ptb.theta behaviours.stats.pos_final.theta_monk(good_indx & ptb_indx)];
    
    
    m.theta_ptb = [m.theta_ptb theta_ptb];
    m.h_index = [m.h_index h_index];
    m.rew = [m.rew reward(good_indx & ptb_indx)];
end

m.ptb_trials(1) = [] ;
con_ptb = [m.ptb_trials.continuous];
prs_ptb = [m.ptb_trials.prs];
%%
m_add = con_ptbprs(genpath);
ms = split_trials(m_add);
[V_out,W_out]=gen_phantom(ms.np.all.targ.r,ms.np.all.targ.theta,ms.np.all.ts.v,ms.np.all.ts.w,...
    ms.p.all.targ.r,ms.p.all.targ.theta,m_add.t.t_ptbn);
%%
m2.ptb.continuous.v = {con_ptb.v_ptb};
m2.ptb.continuous.w = {con_ptb.w_ptb};
m2.ptb.targ.r = m.targ.ptb.r;
m2.ptb.targ.theta = m.targ.ptb.theta;
m2.ptb.resp.r = m.resp.ptb.r;
m2.ptb.resp.theta = m.resp.ptb.theta;

m2.ptb.prs.A_v = [prs_ptb.ptb_linear];
m2.ptb.prs.A_w = [prs_ptb.ptb_angular];
m2.ptb.prs.theta_ptb = [m.theta_ptb];

%no ptb
m2.no_ptb.continuous.v = V_out;
m2.no_ptb.continuous.w = W_out;
m2.ptb.rew = m.rew;
end