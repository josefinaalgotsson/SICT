function updateplots
global F_iu FiuLine
global F_ik FikLine
global Eps  EpsLine
global nutP nutLine
global Ps PsLine
global Pb PbLine
global momFlux momFluxLine
global rho pLine
global reset
global xlabel1 ylabel1
global xlabel2 ylabel2
global xlabel3
global xlabel4 ylabel4
global xlabel5 ylabel5
global xlabel6 ylabel6
global xlabel8 ylabel8
global xlabel9 ylabel9
global xlim1 ylim1
global xlim2 ylim2
global xlim3 ylim3
global xlim4 ylim4
global xlim5 ylim5
global xlim6 ylim6
global xlim8 ylim8
global xlim9 ylim9
global hdls nodeZ
global lockx locky
global faceZ
global u k
%----------------------PLOTTING--------------------------------------------
%---------------------concentration-------------------

if ~reset
    set(FiuLine,'Xdata',u);

    set(FikLine,'Xdata',k);
    
    set(PsLine,'Xdata',Ps);
 
    set(EpsLine,'Xdata',Eps);

    set(nutLine,'Xdata',nutP);

    set(momFluxLine,'Xdata',momFlux);

    set(PbLine,'Xdata',Pb);

    set(pLine,'Xdata',rho);

elseif reset==1 | reset==2
     
    FiuLine=plot(hdls.axis1,u,nodeZ,'-ob',...
        'markersize',4);
    set(hdls.axis1.XLabel,'string',xlabel1);
    set(hdls.axis1.YLabel,'string',ylabel1);
    
    set(hdls.axis1,'Ygrid','on','YMinorgrid','on');
    set(hdls.axis1,'Xgrid','on');
    
    FikLine=plot(hdls.axis3,k,faceZ,'-o');
    set(hdls.axis3.XLabel,'string',xlabel3);
    set(hdls.axis3,'Xgrid','on');
    set(hdls.axis3,'Ygrid','on','YMinorgrid','on');
    %set(hdls.axis3,'YtickLabel',{});
    
    PsLine=plot(hdls.axis2,Ps,faceZ,'-o');
    set(hdls.axis2.XLabel,'string',xlabel2);
    set(hdls.axis2,'Xgrid','on');
    set(hdls.axis2,'Ygrid','on','YMinorgrid','on');
    %set(hdls.axis2,'YtickLabel',{});
    set(hdls.axis2.YLabel,'string',ylabel2);
    
    EpsLine=plot(hdls.axis4,Eps,faceZ,'-o');
    set(hdls.axis4.XLabel,'string',xlabel4);
    set(hdls.axis4,'Xgrid','on');
    set(hdls.axis4,'Ygrid','on','YMinorgrid','on');
    %set(hdls.axis4,'YtickLabel',{});
    set(hdls.axis4.YLabel,'string',ylabel4);
    
    nutLine=plot(hdls.axis5,nutP,nodeZ(2:end-1),'-o');
    set(hdls.axis5.XLabel,'string',xlabel5);
    set(hdls.axis5,'Xgrid','on');
    set(hdls.axis5,'Ygrid','on','YMinorgrid','on');
    %set(hdls.axis5,'YtickLabel',{});
    set(hdls.axis5.YLabel,'string',ylabel5);
    
    
    momFluxLine=plot(hdls.axis6,momFlux,faceZ,'-o');
    set(hdls.axis6.XLabel,'string',xlabel6);
    set(hdls.axis6,'Xgrid','on');
    set(hdls.axis6,'Ygrid','on','YMinorgrid','on');
    %set(hdls.axis6,'YtickLabel',{});
    set(hdls.axis6.YLabel,'string',ylabel6);
    
    PbLine=plot(hdls.axis8,Pb,faceZ,'-o');
    set(hdls.axis8.XLabel,'string',xlabel8);
    set(hdls.axis8,'Xgrid','on');
    set(hdls.axis8,'Ygrid','on','YMinorgrid','on');
    %set(hdls.axis8,'YtickLabel',{});
    set(hdls.axis8.YLabel,'string',ylabel8);
    
    pLine=plot(hdls.axis9,rho,nodeZ,'-o');
    set(hdls.axis9.XLabel,'string',xlabel9);
    set(hdls.axis9,'Xgrid','on');
    set(hdls.axis9,'Ygrid','on','YMinorgrid','on');
    %set(hdls.axis9,'YtickLabel',{});
    set(hdls.axis9.YLabel,'string',ylabel9);

    

        xlim1=get(hdls.axis1,'XLim');
        ylim1=get(hdls.axis1,'YLim');
        xlim2=get(hdls.axis2,'XLim');
        ylim2=get(hdls.axis2,'YLim');
        xlim3=get(hdls.axis3,'XLim');
        ylim3=get(hdls.axis3,'YLim');
        xlim4=get(hdls.axis4,'XLim');
        ylim4=get(hdls.axis4,'YLim');
        xlim5=get(hdls.axis5,'XLim');
        ylim5=get(hdls.axis5,'YLim');
        xlim6=get(hdls.axis6,'XLim');
        ylim6=get(hdls.axis6,'YLim');
        xlim8=get(hdls.axis8,'XLim');
        ylim8=get(hdls.axis8,'YLim');
        xlim9=get(hdls.axis9,'XLim');
        ylim9=get(hdls.axis9,'YLim');
        end
        


%-----------------------------------------------------
