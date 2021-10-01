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
        [~,true.p,shuffled.p] = ComputeROCFirefly([ms.p.all.targ.r' ms.p.all.targ.theta'.*pi/180],[ms.p.all.resp.r' ms.p.all.resp.theta'.*pi/180],1200,10);
        [~,true.sp,shuffled.sp] = ComputeROCFirefly([ms.sp.all.targ.r' ms.sp.all.targ.theta'.*pi/180],[ms.sp.all.resp.r' ms.sp.all.resp.theta'.*pi/180],1200,10);
        [~,true.np,shuffled.np] = ComputeROCFirefly([ms.np.all.targ.r' ms.np.all.targ.theta'.*pi/180],[ms.np.all.resp.r' ms.np.all.resp.theta'.*pi/180],1200,10);
        
        %% Calculate AUC
        
        AUC.p(s,i) = trapz(shuffled.p,true.p);
        AUC.sp(s,i) = trapz(shuffled.sp,true.sp);
        AUC.np(s,i) = trapz(shuffled.np,true.np);
        
    end
    
end
AUC.p(AUC.p == 0) = NaN;
AUC.sp(AUC.sp == 0) = NaN;

%% Humans


%with feedback
[Htarg_f,Hresp_f] = prep_data3_f('/Users/panos/Documents/MATLAB/quarantine/Firefly-ptb/Humans/Data/data_ptbwithfb.mat');

load '/Users/panos/Documents/MATLAB/quarantine/Firefly-ptb/Humans/Data/data_ptbwithfb.mat'
for i = 1:numel(subjects)
    
    for b = 1:20
        p_shuf = randi(numel(Htarg_f.p.r{i}),[1 numel(Htarg_f.p.r{i})]);
        sp_shuf = randi(numel(Htarg_f.sp.r{i}),[1 numel(Htarg_f.sp.r{i})]);
        np_shuf = randi(numel(Htarg_f.np.r{i}),[1 numel(Htarg_f.np.r{i})]);
        
        % ROC Analysis
        [~,Htrue_f.p(i,:),Hshuffled_f.p(i,:)] = ComputeROCFirefly([Htarg_f.p.r{i}(p_shuf)' Htarg_f.p.theta{i}(p_shuf)'.*pi/180],...
            [Hresp_f.p.r{i}(p_shuf)' Hresp_f.p.theta{i}(p_shuf)'.*pi/180],2000,10);
        [~,Htrue_f.sp(i,:),Hshuffled_f.sp(i,:)] = ComputeROCFirefly([Htarg_f.sp.r{i}(sp_shuf)' Htarg_f.sp.theta{i}(sp_shuf)'.*pi/180],...
            [Hresp_f.sp.r{i}(sp_shuf)' Hresp_f.sp.theta{i}(sp_shuf)'.*pi/180],2000,10);
        
        [~,Htrue_f.np(i,:),Hshuffled_f.np(i,:)] = ComputeROCFirefly([Htarg_f.np.r{i}(np_shuf)' Htarg_f.np.theta{i}(np_shuf)'.*pi/180],...
            [Hresp_f.np.r{i}(np_shuf)' Hresp_f.np.theta{i}(np_shuf)'.*pi/180],2000,10);
        
        
        % Calculate AUC
        HAUC_f.p(b,i) = trapz(Hshuffled_f.p(i,:),Htrue_f.p(i,:));
        HAUC_f.sp(b,i) = trapz(Hshuffled_f.sp(i,:),Htrue_f.sp(i,:));
        HAUC_f.np(b,i) = trapz(Hshuffled_f.np(i,:),Htrue_f.np(i,:));
        
        
    end
    
end




%without feedback
[Htarg_nf,Hresp_nf] = prep_data3_nf('/Users/panos/Documents/MATLAB/quarantine/Firefly-ptb/Humans/Data/data_ptbwithoutfb.mat');

load '/Users/panos/Documents/MATLAB/quarantine/Firefly-ptb/Humans/Data/data_ptbwithoutfb.mat'
for i = 1:numel(Htarg_nf.p.r)
    
    for b = 1:20
        p_shuf = randi(numel(Htarg_nf.p.r{i}),[1 numel(Htarg_nf.p.r{i})]);
        sp_shuf = randi(numel(Htarg_nf.sp.r{i}),[1 numel(Htarg_nf.sp.r{i})]);
        np_shuf = randi(numel(Htarg_nf.np.r{i}),[1 numel(Htarg_nf.np.r{i})]);
        
        % ROC Analysis
        [~,Htrue_nf.p(i,:),Hshuffled_nf.p(i,:)] = ComputeROCFirefly([Htarg_nf.p.r{i}(p_shuf)' Htarg_nf.p.theta{i}(p_shuf)'.*pi/180],...
            [Hresp_nf.p.r{i}(p_shuf)' Hresp_nf.p.theta{i}(p_shuf)'.*pi/180],2000,10);
        
        [~,Htrue_nf.sp(i,:),Hshuffled_nf.sp(i,:)] = ComputeROCFirefly([Htarg_nf.sp.r{i}(sp_shuf)' Htarg_nf.sp.theta{i}(sp_shuf)'.*pi/180],...
            [Hresp_nf.sp.r{i}(sp_shuf)' Hresp_nf.sp.theta{i}(sp_shuf)'.*pi/180],2000,10);
        
        [~,Htrue_nf.np(i,:),Hshuffled_nf.np(i,:)] = ComputeROCFirefly([Htarg_nf.np.r{i}(np_shuf)' Htarg_nf.np.theta{i}(np_shuf)'.*pi/180],...
            [Hresp_nf.np.r{i}(np_shuf)' Hresp_nf.np.theta{i}(np_shuf)'.*pi/180],2000,10);
        
        % Calculate AUC
        HAUC_nf.p(b,i) = trapz(Hshuffled_nf.p(i,:),Htrue_nf.p(i,:));
        HAUC_nf.sp(b,i) = trapz(Hshuffled_nf.sp(i,:),Htrue_nf.sp(i,:));
        HAUC_nf.np(b,i) = trapz(Hshuffled_nf.np(i,:),Htrue_nf.np(i,:));
        
    end
    
end


%% Plotting

colr = brewermap(4,'PiYG');

scatter(nanmean(AUC.sp(:,1)),nanmean(AUC.p(:,1)),50,'^','MarkerEdgeColor',colr(4,:),'MarkerFaceColor',colr(4,:)),hold on       %m44
scatter(nanmean(AUC.sp(:,2)),nanmean(AUC.p(:,2)),50,'s','MarkerEdgeColor',colr(4,:),'MarkerFaceColor',colr(4,:)),hold on   %m51
scatter(nanmean(AUC.sp(:,3)),nanmean(AUC.p(:,3)),50,'d','MarkerEdgeColor',colr(4,:),'MarkerFaceColor',colr(4,:)),hold on   %m53

scatter(mean(HAUC_nf.sp),mean(HAUC_nf.p),50,'MarkerEdgeColor',colr(1,:)),hold on       %humans without feedback
scatter(mean(HAUC_f.sp),mean(HAUC_f.p),50,'MarkerEdgeColor',colr(1,:),'MarkerFaceColor',colr(1,:)),hold on       %humans with feedback

plot(0:.1:1,0:.1:1,'--k')
axis([.5 .9 .5 .9])
xlabel('Simulated'),ylabel('True')

% Errorbars

% Humans
%noFeedback
hold on,
err_nf = errorbar(mean(HAUC_nf.sp),mean(HAUC_nf.p),std(HAUC_nf.p)./sqrt(20),std(HAUC_nf.p)./sqrt(20),std(HAUC_nf.sp)./sqrt(20),std(HAUC_nf.sp)./sqrt(20));
err_nf.Color = [0 0 0];
err_nf.LineStyle = 'none';
err_nf.CapSize = 0;
err_nf.LineWidth = 1.25;

%feedback
hold on,
err_f = errorbar(mean(HAUC_f.sp),mean(HAUC_f.p),std(HAUC_f.p)./sqrt(20),std(HAUC_f.p)./sqrt(20),std(HAUC_f.sp)./sqrt(20),std(HAUC_f.sp)./sqrt(20));
err_f.Color = [0 0 0];
err_f.LineStyle = 'none';
err_f.CapSize = 0;
err_f.LineWidth = 1.25;

% Monkeys

err_44 = errorbar(nanmean(AUC.sp(:,1)),nanmean(AUC.p(:,1)),...
    nanstd(AUC.p(:,1))./sqrt(10),nanstd(AUC.p(:,1))./sqrt(10),nanstd(AUC.sp(:,1))./sqrt(10),nanstd(AUC.sp(:,1))./sqrt(10)); %m44
err_44.Color = [0 0 0];
err_44.LineStyle = 'none';
err_44.CapSize = 0;
err_44.LineWidth = 1.25;

err_51 = errorbar(nanmean(AUC.sp(:,2)),nanmean(AUC.p(:,2)),...
    nanstd(AUC.p(:,2))./sqrt(10),nanstd(AUC.p(:,2))./sqrt(10),nanstd(AUC.sp(:,2))./sqrt(10),nanstd(AUC.sp(:,2))./sqrt(10)); %m51
err_51.Color = [0 0 0];
err_51.LineStyle = 'none';
err_51.CapSize = 0;
err_51.LineWidth = 1.25;

err_53 = errorbar(nanmean(AUC.sp(:,3)),nanmean(AUC.p(:,3)),...
    nanstd(AUC.p(:,3))./sqrt(9),nanstd(AUC.p(:,3))./sqrt(9),nanstd(AUC.sp(:,3))./sqrt(9),nanstd(AUC.sp(:,3))./sqrt(9)); %m53
err_53.Color = [0 0 0];
err_53.LineStyle = 'none';
err_53.CapSize = 0;
err_53.LineWidth = 1.25;



