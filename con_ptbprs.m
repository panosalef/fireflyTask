% concatenate all sessions for spacial map prosessing

function  m = con_ptbprs(genpath)
d=dir(fullfile(genpath,'*mat'));
fnames={d.name};

%% Initialize
m.targ.r = [];          m.targ.theta = [];
m.resp.r = [];          m.resp.theta = [];
m.inds.good = [];       m.inds.ptb = [];        m.inds.sesno = [];      m.inds.forw = [];
m.ptb.w = [];           m.ptb.v = [];      m.ptb.h_index2 = [];
m.ts.v = {};            m.ts.w = {};
m.ptb.dx = [];          m.ptb.dy = [];
m.t.stop = [];          m.t.t_ptbn = []; m.t.t_ptbo = [];
m.theta_ptb = [];           % theta at the beginning of the perturbation
m.inds.rew = [];
%% Concatenate
for i=1:numel(fnames)
    load(fullfile(genpath,fnames{i}));
    behaviours = struct(behaviours);
    
    [ptb_dx,ptb_dy]=add_ptbDisp(behaviours);
    [v_ts,w_ts] = get_speedts(behaviours);
    h_index2 = add_hindex2(behaviours);   h_index2 = cellfun(@mean,h_index2);
    for j = 1:numel(behaviours.trials), t_ptbo(j) = behaviours.trials(j).events.t_ptb; end
    
    theta_ptb = add_theta_ptb(behaviours);% theta at the beginning of ptb
    
    events = [behaviours.trials.events];    t_stop = [events.t_stop];
    params = [ behaviours.trials.prs ];     v = [params.ptb_linear ];      w = [params.ptb_angular];
    m.inds.rew = [m.inds.rew behaviours.stats.trialtype.reward(2).trlindx];
    
    %% Get t_ptbn
    [behs]=add_ptbn2(behaviours) ; beh = behs.trials(behaviours.stats.trialtype.ptb(2).trlindx & behaviours.stats.trialtype.all.trlindx);
    
    event_ptbn = [beh.events] ; t_ptbn = [event_ptbn.t_ptbn];
    prs=[behaviours.trials.prs]; ptb_linear=[prs.ptb_linear]; forw = ptb_linear > 0;
    
    m.inds.ptb = [ m.inds.ptb behaviours.stats.trialtype.ptb(2).trlindx ];
    m.inds.good = [ m.inds.good behaviours.stats.trialtype.all.trlindx ];
    m.inds.forw = [m.inds.forw forw];
    
    m.targ.r = [ m.targ.r behaviours.stats.pos_final.r_targ ];
    m.targ.theta = [ m.targ.theta behaviours.stats.pos_final.theta_targ ];
    m.resp.r = [ m.resp.r behaviours.stats.pos_final.r_monk ];
    m.resp.theta = [ m.resp.theta behaviours.stats.pos_final.theta_monk ];
    
    m.ts.v = [m.ts.v v_ts];            m.ts.w = [m.ts.w w_ts];
    m.t.t_ptbn = [m.t.t_ptbn t_ptbn];
    m.t.t_ptbo = [m.t.t_ptbo t_ptbo(behaviours.stats.trialtype.ptb(2).trlindx & behaviours.stats.trialtype.all.trlindx)];
    m.t.stop = [ m.t.stop t_stop ];   m.inds.sesno = [m.inds.sesno repmat(i,[1 numel(behaviours.stats.pos_final.r_targ)])];
    m.ptb.w = [ m.ptb.w w ];           m.ptb.v = [ m.ptb.v v ];         m.ptb.h_index2 = [ m.ptb.h_index2 h_index2 ];
    m.ptb.dx = [m.ptb.dx ptb_dx'];          m.ptb.dy = [m.ptb.dy ptb_dy'];
    m.theta_ptb = [m.theta_ptb theta_ptb];
end

%% Calculate cartesian Coordinates
m.targ.x = sind(m.targ.theta).*m.targ.r;      m.targ.y = cosd(m.targ.theta).*m.targ.r;
m.resp.x = sind(m.resp.theta).*m.resp.r;      m.resp.y = cosd(m.resp.theta).*m.resp.r;

%% Calculate error vector
x_err = m.resp.x - m.targ.x;
y_err = m.resp.y - m.targ.y;

m.err.x = x_err;
m.err.y = y_err;
m.err.abs = sqrt((x_err).^2 + (y_err).^2);

end