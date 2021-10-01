[cdf44,nx] = time_cdf_fun('C:\Users\pa77\Documents\Github\Firefly-gain\m44');
[cdf51,~] = time_cdf_fun('C:\Users\pa77\Documents\Github\Firefly-gain\m51');
[cdf53,~] = time_cdf_fun('C:\Users\pa77\Documents\Github\Firefly-gain\m53');

cdf_g1 = [cdf44.g1'; cdf51.g1'; cdf53.g1'];
cdf_g15 = [cdf44.g15'; cdf51.g15'; cdf53.g15'];
cdf_g2 = [cdf44.g2'; cdf51.g2'; cdf53.g2'];


%% Plot
colr = brewermap(3,'Set1');

plot(nx,mean(cdf_g1),'color',colr(1,:),'linewidth',1.3), hold on
plot(nx,mean(cdf_g15),'color',colr(3,:),'linewidth',1.3), hold on
plot(nx,mean(cdf_g2),'color',colr(2,:),'linewidth',1.3), hold on

xlabel('Movement Duration')
ylabel('Cumulative freq')