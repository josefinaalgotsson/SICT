% Author: MSc Josefina Algotsson, University of Gothenburg
% This software is licensed under the GNU General Public License v.3
%Validation according to mixed layer depth 2-layer stratification
%Finds mld from model simulation. Uses KPA expression to find dhdt and the
%mld for KPA.
clc
close all
saveit=0;


movmeansiz=20;
filename={'Test'};                  % Filenames to evaluate
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
ylabel('h_{m} : h_{KP}','FontSize',16)
title('Model to experimental ratio h_{m} : h_{KP}','FontSize',20)
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
ylabel('h''_m : h''_{KP}','FontSize',16)
title('Model to experimental ratio ratio h''_m : h''_{KP}','FontSize',20)
set(gcf,'color','white')
%ylim([0 2])
grid on
hold on


for n=1:length(filename)
    load([cd,'/Result/',filename{n}])       % Load filename
    mld=[];                                 % Empty vector for model mixed layer depth
    ustarA=sqrt(abs(boundValNu));                        % Analytical/ hard set ustar
    rho=rho0.*(1-alfa.*(Tsave(:,1)-T0)+beta.*(Ssave(:,1)-S0)); 
    % Density gradient
    drhodz = -diff(rho)/deltaz;    
    % Buoyancy frequency squared 
    N2=-g/rho0*drhodz;

    NA=mean(sqrt(abs(N2)));                 % Analytical/hard set N

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
    legendentry{end+1}=filename{n};             % Update legendentry with filename
    legendentry2{end+1}=filename{n};
    legendentry3{end+1}=filename{n};
    % Figure 1: Mixed layer depths
    %Create mld curve from KPA
    mld_KP=1.047.*ustarA.*sqrt((tsave./NA))';
    plot(ax1,tsave./(60*60),mld,'-o')% Plot MLD of model
    plot(ax1,tsave./(60*60),mld_KP)
     
    %Figure3: dh/dt 
    %dh/dt for model
    dhdt=diff(mld)./diff(tsave');               % dh dt of model
    dhdtmov=movmean(dhdt,movmeansiz);                    % dh/dt with moving mean
    newtime=(tsave(2:end)+tsave(1:end-1))/2';   % New timestamp
    % dh/dt for KPA
    dhdt_KP=1.047.*ustarA.*0.5*newtime'.^(-0.5)./sqrt(NA);   
    p=plot(ax3,newtime./(60*60),dhdt);
    plot(ax3,newtime./(60*60),dhdt_KP);
    movcolor=get(p,'color');
    plot(ax3,newtime./(60*60),dhdtmov,'*-','color',movcolor);

    % Figure 2: Quote between mixed layer depths
    quote_mld=mld./(mld_KP);
    plot(ax2,tsave./(60*60),quote_mld);
    
    % Figure 4: Quote between dh/dts
    quote_dhdt=dhdt./dhdt_KP;
    quote_dhdtmov=dhdtmov./dhdt_KP;
    p2=plot(ax4,newtime./(60*60),quote_dhdt);
    movcolor2=get(p2,'color');
    plot(ax4,newtime./(60*60),quote_dhdtmov,'*-','color',movcolor2);
    
    legendentry{end+1}='KP';             % Update legendentry with KPA
    legendentry2{end+1}='KP';
    legendentry2{end+1}='moving mean';
    legendentry3{end+1}='moving mean';
end

%Figure 1 mld
legend(ax1,legendentry)

%Figure 2 mld quote
plot(ax2,[tsave(1)./(60*60) 7],[1 1],'k--')
ylim(ax2,[0 2])
legend(ax2,filename)

%Figure 3 dhdt
legend(ax3,legendentry2)

%Figure 4 dhdt quote
plot(ax4,[tsave(1)./(60*60) 7],[1 1],'k--')
legend(ax4,legendentry3)

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