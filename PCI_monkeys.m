clear
monks = [44 51 53];

%% Monkeys
for i=1:numel(monks)
    [m] = con_ptbprs(strcat('/Users/panos/Documents/MATLAB/quarantine/Firefly-ptb/m',num2str(monks(i))));
    ms = split_trials(m);
    ses = unique(ms.prs.sesno.all);
    
    for s = 1:numel(ses)
        
        [ms.sp.all.resp,ms.sp.all.targ] = gen_simptb_final3(ms.np.all.targ.r(ms.prs.sesno.np == s),ms.np.all.targ.theta(ms.prs.sesno.np == s),...
            ms.np.all.ts.v(ms.prs.sesno.np == s), ms.np.all.ts.w(ms.prs.sesno.np == s),...
            ms.p.all.targ.r(ms.prs.sesno.p == s),ms.p.all.targ.theta(ms.prs.sesno.p == s),...
            ms.prs.t_ptbn(ms.prs.sesno.p == s),ms.prs.v(ms.prs.sesno.p == s),ms.prs.w(ms.prs.sesno.p == s));

        %% ROC Analysis
        [~,true.p,shuffled.p] = ComputeROCFirefly([ms.p.all.targ.r' ms.p.all.targ.theta'.*pi/180],[ms.p.all.resp.r' ms.p.all.resp.theta'.*pi/180],1500,10);
        [~,true.sp,shuffled.sp] = ComputeROCFirefly([ms.sp.all.targ.r' ms.sp.all.targ.theta'.*pi/180],[ms.sp.all.resp.r' ms.sp.all.resp.theta'.*pi/180],1500,10);
        [~,true.np,shuffled.np] = ComputeROCFirefly([ms.np.all.targ.r' ms.np.all.targ.theta'.*pi/180],[ms.np.all.resp.r' ms.np.all.resp.theta'.*pi/180],1500,10);
        
        %% Calculate AUC
        
        AUC.p(s,i) = trapz(shuffled.p,true.p);
        AUC.sp(s,i) = trapz(shuffled.sp,true.sp);
        AUC.np(s,i) = trapz(shuffled.np,true.np);
        
    end  
end
AUC.p(AUC.p == 0) = NaN;
AUC.sp(AUC.sp == 0) = NaN;

PCI_m = (AUC.p - AUC.sp)./(AUC.np - AUC.sp);
PCI_m(PCI_m < 0) = NaN;

PCI_q_std = nanstd(PCI_m(:,1));
PCI_b_mean = ceil(nanmean(PCI_m(:,2)*100))./100;    PCI_b_std = nanstd(PCI_m(:,2));
PCI_s_mean = ceil(nanmean(PCI_m(:,3)*100))./100;    PCI_s_std = nanstd(PCI_m(:,3));

%% Plot

figure,hold on
bar(PCI_q_mean)
err_nf = errorbar(1,PCI_q_mean,PCI_q_std);
err_nf.Color = [0 0 0];
err_nf.LineStyle = 'none';
err_nf.CapSize = 0;
err_nf.LineWidth = 1.25;
ylim([0 1.2])
title('mQ')

figure,hold on
bar(PCI_b_mean)
err_nf = errorbar(1,PCI_b_mean,PCI_s_std);
err_nf.Color = [0 0 0];
err_nf.LineStyle = 'none';
err_nf.CapSize = 0;
err_nf.LineWidth = 1.25;
ylim([0 1.2])
title('mB')

figure,hold on
bar(PCI_s_mean)
err_nf = errorbar(1,PCI_s_mean,PCI_s_std);
err_nf.Color = [0 0 0];
err_nf.LineStyle = 'none';
err_nf.CapSize = 0;
err_nf.LineWidth = 1.25;
ylim([0 1.2])
title('mS')


