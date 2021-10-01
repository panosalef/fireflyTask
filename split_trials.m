%Goes with con_ptbprs() function;

function m_split = split_trials(m)   
%%%%%%%%%%%%%%%%
%% Ptb Trials %%
%%%%%%%%%%%%%%%%

%% Target
m_split.p.all.targ.r = m.targ.r(m.inds.good & m.inds.ptb);
m_split.p.all.targ.theta = m.targ.theta(m.inds.good & m.inds.ptb);

m_split.p.all.targ.x = m_split.p.all.targ.r .* sind(m_split.p.all.targ.theta);
m_split.p.all.targ.y = m_split.p.all.targ.r .* cosd(m_split.p.all.targ.theta);

%% Responce
m_split.p.all.resp.r = m.resp.r(m.inds.good & m.inds.ptb);
m_split.p.all.resp.theta = m.resp.theta(m.inds.good & m.inds.ptb);

m_split.p.all.resp.x = m_split.p.all.resp.r .* sind(m_split.p.all.resp.theta);
m_split.p.all.resp.y = m_split.p.all.resp.r .* cosd(m_split.p.all.resp.theta);

%% Time Series
m_split.p.all.ts.v = m.ts.v(m.inds.good & m.inds.ptb);
m_split.p.all.ts.w = m.ts.w(m.inds.good & m.inds.ptb); 

%% H_index Slpit
m_split.prs.Hi = m.ptb.h_index2;
Hi_indx = m.ptb.h_index2 > 0;
m_split.prs.Hi_indx = Hi_indx;
m_split.prs.pdx = m.ptb.dx;
m_split.prs.pdy = m.ptb.dy;
%% target

%Helpfull
m_split.p.h.targ.r = m_split.p.all.targ.r(Hi_indx);
m_split.p.h.targ.theta = m_split.p.all.targ.theta(Hi_indx);

m_split.p.h.targ.x = m_split.p.h.targ.r .*sind(m_split.p.h.targ.theta);
m_split.p.h.targ.y = m_split.p.h.targ.r .*cosd(m_split.p.h.targ.theta);

%NotHelpfull
m_split.p.nh.targ.r = m_split.p.all.targ.r(~Hi_indx);
m_split.p.nh.targ.theta = m_split.p.all.targ.theta(~Hi_indx);

m_split.p.nh.targ.x = m_split.p.nh.targ.r .*sind(m_split.p.nh.targ.theta);
m_split.p.nh.targ.y = m_split.p.nh.targ.r .*cosd(m_split.p.nh.targ.theta);

%% responce

%Helpfull
m_split.p.h.resp.r = m_split.p.all.resp.r(Hi_indx);
m_split.p.h.resp.theta = m_split.p.all.resp.theta(Hi_indx);

m_split.p.h.resp.x = m_split.p.h.resp.r .*sind(m_split.p.h.resp.theta);
m_split.p.h.resp.y= m_split.p.h.resp.r .*cosd(m_split.p.h.resp.theta);

%NotHelpfull
m_split.p.nh.resp.r = m_split.p.all.resp.r(~Hi_indx);
m_split.p.nh.resp.theta = m_split.p.all.resp.theta(~Hi_indx);

m_split.p.nh.resp.x = m_split.p.nh.resp.r .*sind(m_split.p.nh.resp.theta);
m_split.p.nh.resp.y= m_split.p.nh.resp.r .*cosd(m_split.p.nh.resp.theta);


%%%%%%%%%%%%%%%%%%%
%% NO Ptb Trials %%
%%%%%%%%%%%%%%%%%%%

%% Target
m_split.np.all.targ.r = m.targ.r(m.inds.good & ~m.inds.ptb);
m_split.np.all.targ.theta = m.targ.theta(m.inds.good & ~m.inds.ptb);

m_split.np.all.targ.x = m_split.np.all.targ.r .* sind(m_split.np.all.targ.theta);
m_split.np.all.targ.y = m_split.np.all.targ.r .* cosd(m_split.np.all.targ.theta);

%% Responce
m_split.np.all.resp.r = m.resp.r(m.inds.good & ~m.inds.ptb);
m_split.np.all.resp.theta = m.resp.theta(m.inds.good & ~m.inds.ptb);

m_split.np.all.resp.x = m_split.np.all.resp.r .* sind(m_split.np.all.resp.theta);
m_split.np.all.resp.y = m_split.np.all.resp.r .* cosd(m_split.np.all.resp.theta);

%% Time Series
m_split.np.all.ts.v = m.ts.v(m.inds.good & ~m.inds.ptb);
m_split.np.all.ts.w = m.ts.w(m.inds.good & ~m.inds.ptb);


%% session number
m_split.prs.sesno.all = m.inds.sesno(logical(m.inds.good));
m_split.prs.sesno.p = m.inds.sesno(m.inds.good & m.inds.ptb);
m_split.prs.sesno.np = m.inds.sesno(m.inds.good & ~m.inds.ptb);
%% Other Params
m_split.prs.t_ptbn = m.t.t_ptbn;
m_split.prs.w = m.ptb.w(m.inds.good & m.inds.ptb);
m_split.prs.v = m.ptb.v(m.inds.good & m.inds.ptb);

end