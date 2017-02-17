% Log validation

close all
clear all
saveit=0;

filename={'TestLOG'}
load([cd,'/Result/',filename{1}])       % Load filename

dudz2=(diff(usave(:,end))/deltaz).^2;

ustar00 = sqrt(nutTot(1).*sqrt(dudz2(1)));
z0 = H/exp(boundValNu*kar/ustar00);
u00=ustar00/kar.*log((-flipud(faceZ)+H+z0)/z0)

z0 = H/exp(boundValNu*kar/ustar00)

figure
plot(usave(:,end),nodeZ,'o','linewidth',1)
grid on
hold on
plot(u00,faceZ,'linewidth',1)
hold off
legend('SICT','Log')
ylabel('z')
xlabel('u [m/s]')
title(['Simulated time [s]:', num2str(max(tsave))])

if saveit
    set(gcf,'color','white')
    export_fig('LogValidation','-pdf')
end