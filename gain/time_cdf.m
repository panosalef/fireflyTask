%% time_cdf.m  -  travel-time CDFs per joystick gain, averaged over animals
% Expects one folder per animal under DATA_ROOT (m44, m51, m53), each holding
% the session .mat files. Set FIREFLY_DATA_ROOT or edit DATA_ROOT below.
DATA_ROOT = getenv('FIREFLY_DATA_ROOT'); if isempty(DATA_ROOT), DATA_ROOT = fullfile(pwd,'data','gain'); end

[cdf44,nx] = time_cdf_fun(fullfile(DATA_ROOT,'m44'));
[cdf51,~] = time_cdf_fun(fullfile(DATA_ROOT,'m51'));
[cdf53,~] = time_cdf_fun(fullfile(DATA_ROOT,'m53'));

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