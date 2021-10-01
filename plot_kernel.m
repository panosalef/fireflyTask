[~,m] = concat_basisF2('/Users/panos/Documents/MATLAB/quarantine/Firefly-ptb/m51');

f_idx = m.ptb.prs.A_v > 0 ;
con_idx = sign((m.ptb.targ.theta) .* m.ptb.prs.A_w) > 0 ;

%% Forward
continuous_f = m.ptb.continuous.v(f_idx);
Yresp_f  = [continuous_f{:}]; Yresp_f = Yresp_f(:);
ph_continuous_f = m.no_ptb.continuous.v(f_idx);
ph_Yresp_f  = [ph_continuous_f{:}]; ph_Yresp_f = ph_Yresp_f(:);
Amplitudes_f = m.ptb.prs.A_v(f_idx);

[kernel.f,ts] = get_kernel_Base(Yresp_f,Amplitudes_f,ph_Yresp_f);



%% Backward
continuous_b = m.ptb.continuous.v(~f_idx);
Yresp_b  = [continuous_b{:}]; Yresp_b = Yresp_b(:);
ph_continuous_b = m.no_ptb.continuous.v(~f_idx);
ph_Yresp_b  = [ph_continuous_b{:}]; ph_Yresp_b = ph_Yresp_b(:);
Amplitudes_b = m.ptb.prs.A_v(~f_idx);

[kernel.b,~] = get_kernel_Base(Yresp_b,Amplitudes_b,ph_Yresp_b);


%% Congruent
continuous_c = m.ptb.continuous.w(con_idx);
Yresp_c  = [continuous_c{:}]; Yresp_c = Yresp_c(:);
ph_continuous_c = m.no_ptb.continuous.w(con_idx);
ph_Yresp_c  = [ph_continuous_c{:}]; ph_Yresp_c = ph_Yresp_c(:);
Amplitudes_c = m.ptb.prs.A_w(con_idx);

[kernel.c,~] = get_kernel_Base(Yresp_c,Amplitudes_c,ph_Yresp_c);


%% Incongruent
continuous_ic = m.ptb.continuous.w(~con_idx);
Yresp_ic  = [continuous_ic{:}]; Yresp_ic = Yresp_ic(:);
ph_continuous_ic = m.no_ptb.continuous.w(~con_idx);
ph_Yresp_ic  = [ph_continuous_ic{:}]; ph_Yresp_ic = ph_Yresp_ic(:);
Amplitudes_ic = m.ptb.prs.A_w(~con_idx);

[kernel.ic,~] = get_kernel_Base(Yresp_ic,Amplitudes_ic,ph_Yresp_ic);

%% Plot
figure,
subplot(1,2,1),plot(ts,kernel.b,ts,kernel.f)
subplot(1,2,2),plot(ts,kernel.ic,ts,kernel.c)

