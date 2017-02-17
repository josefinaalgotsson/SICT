%Validation according to mixed layer depth 2-layer stratification
%Finds mld from model simulation. Uses KPA expression to find dhdt and the
%mld for KPA.
clc
close all
clear all
saveit=0;

movmeansiz=20;
filename={'Test'};            % Filenames to evaluate
legendentry={};                     % Empty cell array to save legententries
legendentry2={};
legendentry3={};


%For mld
f1=figure;
ax1=gca;
xlabel('Time [hrs]','FontSize',16)
ylabel('h [m]','FontSize',16)
title('Mixed layer depth, h','Fontsize',20)
axis ij
set(gcf,'color','white')
grid on
hold on
    
f2=figure;
ax2=gca;
xlabel('Time [hrs]','FontSize',16)
ylabel('h_{m} : h_{KPA}','FontSize',16)
title('Model to experimental ratio h_{m} : h_{KPA}','FontSize',20)
set(gcf,'color','white')
grid on
hold on

f3=figure;
ax3=gca;
xlabel('Time [hrs]','FontSize',16)
ylabel('h'' [m/s]','FontSize',16)
title('Entrainment velocity h''','FontSize',20)
set(gcf,'color','white')
grid on
hold on 

f4=figure;
ax4=gca;
xlabel('Time [hrs]','FontSize',16)
ylabel('h''_m : h''_{KPA}','FontSize',16)
title('Model to experimental ratio ratio h''_m : h''_{KP}','FontSize',20)
set(gcf,'color','white')
%ylim([0 2])
grid on
hold on

rho0=1000;                      % Reference density
t0=10;                          % Reference temperature
s0=35;                          % Reference salinity
alfa=0.17e-4;                   % Epansion coefficient
beta=7.6e-4;                    % Expansion coefficient

for n=1:length(filename)
    load([cd,'/Result/',filename{n}])       % Load filename
    mld=[];                                 % Empty vector for model mixed layer depth
    for k=1:size(ksave,2)                   % For all timesteps
        for m=1:size(ksave,1)               % For all depths
            if ksave(m,k)<1e-6              % If element reaches criteria for mixed layer
                mld(k,1)=deltaz*m-deltaz/2; % Calculate depth
                break
            end
            % If we have reached the bottom of the profile and no new mld
            % has been found put in nan
            if m==size(ksave,1) && ~isequal(length(mld),k)
                mld(k)=nan;
            end
        end
    end
    ustarA=sqrt(abs(boundValNu));                        % Analytical/ hard set ustar
    legendentry{end+1}=filename{n};             % Update legendentry with filename
    legendentry2{end+1}=filename{n};
    legendentry3{end+1}=filename{n};
    
    %Figure3: dh/dt 
    %dh/dt for model
    dhdt=diff(mld)./diff(tsave');               % dh dt of model
    dhdtmov=movmean(dhdt,movmeansiz);                      % dh/dt with moving mean
    newtime=(tsave(2:end)+tsave(1:end-1))/2';   % New timestamp
    % dh/dt for KPA
    rhosave=rho0.*(1-alfa.*(Tsave-t0)+beta.*(Ssave-s0)); % Calculate density from modeled T and S
    %t=min(tsave):1:max(tsave);                          % Set up a time vector
    
    Rv=0.6;                                             % Constant
    h0=find(diff(rhosave(:,1))==max(diff(rhosave(:,1))))*deltaz; %find breakpoint in initiation density
    drho0=rhosave(1,1)-rhosave(end,1);                  % Density difference at start
    rho=rhosave(end,1);                                 % Density in bottom layer
    g=9.81;                                             % Gravitational acceleration
    B=abs(g*h0*drho0/rho);                              % Buoyancy
    dhdt_KPA=ustarA^2.*(Rv/B)^(1/2)';                   % dhdt according to Kantha Phillips Azhad
    
    p2=plot(ax3,newtime./(60*60),dhdt);
    plot(ax3,newtime./(60*60),dhdt_KPA.*ones(size(newtime)))
    movcolor2=get(p2,'color');
    plot(ax3,newtime./(60*60),dhdtmov,'*-','color',movcolor2)
    
    % Figure 1: Mixed layer depths
    %Create mld curve from KPA
    mld_KPA=dhdt_KPA.*tsave'+h0;                 % Create line with slope dhdt
    plot(ax1,tsave./(60*60),mld,'-o')% Plot MLD of model
    plot(ax1,tsave./(60*60),mld_KPA)
    plot(ax1,tsave./(60*60),mld_KPA-5.5,'--r')
    
    % Figure 2: Quote between mixed layer depths
    quote_mld=mld./(mld_KPA-5);
    plot(ax2,tsave./(60*60),quote_mld);
    
    % Figure 4: Quote between dh/dts
    quote_dhdt=dhdt/dhdt_KPA;
    quote_dhdtmov=dhdtmov/dhdt_KPA;
    p=plot(ax4,newtime./(60*60),quote_dhdt);
    movcolor=get(p,'color');
    plot(ax4,newtime./(60*60),quote_dhdtmov,'*-','color',movcolor);
    
    legendentry{end+1}='KPA';             % Update legendentry with KPA
    legendentry2{end+1}='KPA';
    legendentry2{end+1}='moving mean';
    legendentry3{end+1}='moving mean';
   
end

%Figure 1 mld
legend(ax1,legendentry)

%Figure 2 mld quote
ylim(ax2,[0 2])
legend(ax2,filename)
plot(ax2,[tsave(1)./(60*60) 7],[1 1],'k--')

%Figure 3 dhdt
legend(ax3,legendentry2)

%Figure 4 dhdt quote
legend(ax4,legendentry3)
plot(ax4,[tsave(1)./(60*60) 7],[1 1],'k--')

ylim(ax1,[0 20])
ylim(ax2,[0 2])
ylim(ax3,[0 6e-3])
ylim(ax4,[0 2])
if saveit
timestamp=datestr(date,'yymmdd');
export_fig(f1,[cd,'/Result/Figures/',char(filename),'valMLD',timestamp],'-pdf')
export_fig(f2,[cd,'/Result/Figures/',char(filename),'valquote1',timestamp],'-pdf')
export_fig(f3,[cd,'/Result/Figures/',char(filename),'valdhdt',timestamp],'-pdf')
export_fig(f4,[cd,'/Result/Figures/',char(filename),'valdhdtquote2',timestamp],'-pdf')
savefig(f1,[cd,'/Result/Figures/',char(filename),'valMLD',timestamp])
savefig(f2,[cd,'/Result/Figures/',char(filename),'valquote1',timestamp])
savefig(f3,[cd,'/Result/Figures/',char(filename),'valdhdt',timestamp])
savefig(f4,[cd,'/Result/Figures/',char(filename),'valdhdtquote2',timestamp])
end


si=tsave(2)-tsave(1);
movmeantimespan=(si*movmeansiz);

disp(['Timespan for moving mean is ',num2str(movmeantimespan/60),' minutes'])