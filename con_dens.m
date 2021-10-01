function [m]=con_dens(genpath)

d=dir(fullfile(genpath,'*mat'));
fnames={d.name};

%% Initialize
m.targ.r = [];          m.targ.theta = [];
m.resp.r = [];          m.resp.theta = [];
m.inds.good = [];       m.inds.sesno = [];
m.prs.dens_all = [];

for i=1:numel(fnames)
    
    load(fullfile(genpath,fnames{i}));
    d = []; for j=1:numel(behaviours.trials), d(j) = behaviours.trials(j).prs.floordensity; end
    m.inds.good = [ m.inds.good behaviours.stats.trialtype.all.trlindx ];
    m.targ.r = [ m.targ.r behaviours.stats.pos_final.r_targ ];
    m.targ.theta = [ m.targ.theta behaviours.stats.pos_final.theta_targ ];
    m.resp.r = [ m.resp.r behaviours.stats.pos_final.r_monk ];
    m.resp.theta = [ m.resp.theta behaviours.stats.pos_final.theta_monk ];
    m.prs.dens_all = [m.prs.dens_all d];
    m.inds.sesno = [m.inds.sesno repmat(i,[1 numel(behaviours.stats.pos_final.r_targ)])];
    
end

m.prs.dens_uniq = sort(unique(m.prs.dens_all));

for g = 1:numel(m.prs.dens_uniq)
    m.inds.dens(g,:) = m.prs.dens_all == m.prs.dens_uniq(g);
end

%% Densities 0.0001
m.keep.min.targ.r = m.targ.r(m.inds.good & m.prs.dens_all == .0001);
m.keep.min.targ.theta = m.targ.theta(m.inds.good & m.prs.dens_all == .0001);

m.keep.min.targ.x = m.keep.min.targ.r.*sind(m.keep.min.targ.theta);
m.keep.min.targ.y = m.keep.min.targ.r.*cosd(m.keep.min.targ.theta);


m.keep.min.resp.r = m.resp.r(m.inds.good & m.prs.dens_all == .0001);
m.keep.min.resp.theta = m.resp.theta(m.inds.good & m.prs.dens_all == .0001);

m.keep.min.resp.x = m.keep.min.resp.r.*sind(m.keep.min.resp.theta);
m.keep.min.resp.y = m.keep.min.resp.r.*cosd(m.keep.min.resp.theta);

m.keep.min.error = sqrt((m.keep.min.resp.y - m.keep.min.targ.y).^2 + (m.keep.min.resp.x - m.keep.min.targ.x).^2);

%% Densities 0.005
m.keep.max.targ.r = m.targ.r(m.inds.good & m.prs.dens_all == .005 );
m.keep.max.targ.theta = m.targ.theta(m.inds.good & m.prs.dens_all == .005);

m.keep.max.targ.x = m.keep.max.targ.r.*sind(m.keep.max.targ.theta);
m.keep.max.targ.y = m.keep.max.targ.r.*cosd(m.keep.max.targ.theta);


m.keep.max.resp.r = m.resp.r(m.inds.good & m.prs.dens_all == .005);
m.keep.max.resp.theta = m.resp.theta(m.inds.good & m.prs.dens_all == .005);

m.keep.max.resp.x = m.keep.max.resp.r.*sind(m.keep.max.resp.theta);
m.keep.max.resp.y = m.keep.max.resp.r.*cosd(m.keep.max.resp.theta);

m.keep.max.error = sqrt((m.keep.max.resp.y - m.keep.max.targ.y).^2 + (m.keep.max.resp.x - m.keep.max.targ.x).^2);

%%
m.keep.min.sesno = m.inds.sesno(m.inds.good & m.prs.dens_all == .0001) ;
m.keep.max.sesno = m.inds.sesno(m.inds.good & m.prs.dens_all == .005) ;


end
