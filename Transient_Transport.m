function varargout = Transient_Transport(varargin)
% TRANSIENT_TRANSPORT MATLAB code for Transient_Transport.fig
%      TRANSIENT_TRANSPORT, by itself, creates a new TRANSIENT_TRANSPORT or raises the existing
%      singleton*.
%
%      H = TRANSIENT_TRANSPORT returns the handle to a new TRANSIENT_TRANSPORT or the handle to
%      the existing singleton*.
%
%      TRANSIENT_TRANSPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRANSIENT_TRANSPORT.M with the given input arguments.
%
%      TRANSIENT_TRANSPORT('Property','Value',...) creates a new TRANSIENT_TRANSPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Transient_Transport_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Transient_Transport_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Transient_Transport

% Last Modified by GUIDE v2.5 17-Feb-2017 14:55:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Transient_Transport_OpeningFcn, ...
    'gui_OutputFcn',  @Transient_Transport_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Transient_Transport is made visible.
function Transient_Transport_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Transient_Transport (see VARARGIN)

% Choose default command line output for Transient_Transport
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
% Define global variables
global  xlabel1  ylabel1
global  xlabel2  ylabel2
global  xlabel3  ylabel3
global  xlabel4  ylabel4
global  xlabel5  ylabel5
global  xlabel6  ylabel6
global  xlabel8  ylabel8
global  xlabel9  ylabel9
global  hs hdls
% Save all handles to global variable
hdls=handles;
% Velocity u plot
set(hdls.axis1.XLabel,'string','Velocity u [m/s]');
set(hdls.axis1.YLabel,'string','M.A.B');
set(hdls.axis1,'Ygrid','on');
set(hdls.axis1,'Xgrid','on');
% Production of t.k.e through shear plot
set(hdls.axis2.XLabel,'string',{'P_s','[m^2/s^3]'});
set(hdls.axis2.YLabel,'string','');
set(hdls.axis2,'Ygrid','on');
set(hdls.axis2,'Xgrid','on');
set(hdls.axis2,'yticklabels',{});
% Turbulent Kinetic Energy (t.k.e) plot
set(hdls.axis3.YLabel,'string','');
set(hdls.axis3.XLabel,'String',{'k','[m^2/s^2]'});
set(hdls.axis3,'Ygrid','on');
set(hdls.axis3,'Xgrid','on');
set(hdls.axis3,'yticklabels',{});
% Dissipation Epsilon plot
set(hdls.axis4.XLabel,'string',{char(949),'[m^2/s^3]'}');
set(hdls.axis4.YLabel,'string','');
set(hdls.axis4,'Ygrid','on');
set(hdls.axis4,'Xgrid','on');
set(hdls.axis4,'yticklabels',{});
% Turbulent viscosity nutP plot
set(hdls.axis5.XLabel,'string',{'\nu_t','[m^2/s]'});
set(hdls.axis5.YLabel,'string','');
set(hdls.axis5,'Ygrid','on');
set(hdls.axis5,'Xgrid','on');
set(hdls.axis5,'yticklabels',{});
% Momentum flux Tau plot
set(hdls.axis6.XLabel,'string',{'\tau','[m^2/s^2]'});
set(hdls.axis6.YLabel,'string','');
set(hdls.axis6,'Ygrid','on');
set(hdls.axis6,'Xgrid','on');
set(hdls.axis6,'yticklabels',{});
% Production of Turbulent Kinetic Energy, Pb, plot
set(hdls.axis8.XLabel,'string',{'P_b','[m^2/s^3]'});
set(hdls.axis8.YLabel,'string','');
set(hdls.axis8,'Ygrid','on');
set(hdls.axis8,'Xgrid','on');
set(hdls.axis8,'yticklabels',{});
% Density F_ip plot
set(hdls.axis9.XLabel,'string',{'\rho','[kg/m^3]'});
set(hdls.axis9.YLabel,'string','');
set(hdls.axis9,'Ygrid','on');
set(hdls.axis9,'Xgrid','on');
set(hdls.axis9,'yticklabels',{});

% Get and save x -and ylabels for plots for later
xlabel1=get(hdls.axis1.XLabel,'string');
ylabel1=get(hdls.axis1.YLabel,'string');
xlabel2=get(hdls.axis2.XLabel,'string');
ylabel2=get(hdls.axis2.YLabel,'string');
xlabel3=get(hdls.axis3.XLabel,'string');
ylabel3=get(hdls.axis3.YLabel,'string');
xlabel4=get(hdls.axis4.XLabel,'string');
ylabel4=get(hdls.axis4.YLabel,'string');
xlabel5=get(hdls.axis5.XLabel,'string');
ylabel5=get(hdls.axis5.YLabel,'string');
xlabel6=get(hdls.axis6.XLabel,'string');
ylabel6=get(hdls.axis6.YLabel,'string');
xlabel8=get(hdls.axis8.XLabel,'string');
ylabel8=get(hdls.axis8.YLabel,'string');
xlabel9=get(hdls.axis9.XLabel,'string');
ylabel9=get(hdls.axis9.YLabel,'string');
% Set up invisible figure
hs=addcomponents;
% Get default inputs from figure
getinputs %Get user inputs
% Signal that user has not yet set initiation
reset=0;
% Evaluate stratification
popStrat_Callback

% --- Executes on button press in advancedButton.
function advancedButton_Callback(hObject, eventdata, handles)
% hObject    handle to advancedButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Handle to external UI
global hs
% Make figure visible
hs.fig.Visible = 'on';

% --- Executes on button press in browseOutputFile.
function browseOutputFile_Callback(hObject, eventdata, handles)
% hObject    handle to browseOutputFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hdls j filename
% Open a UI popup window to select a filename
[filename] = uigetfile;
if ~filename==0
    set(hdls.editOutputFile,'string',filename);  % Writes filename text in UI
end
% Scroll to the right in output folder text box
j=findjobj(hdls.editOutputFolder);
j.setCaretPosition(length(get(hdls.editOutputFolder,'string')))
set(j,'FocusLostCallback',{@scrollright,hdls.editOutputFolder})

% --- Executes on button press in browseOutputFolder.
function browseOutputFolder_Callback(hObject, eventdata, handles)
% hObject    handle to browseOutputFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hdls pathname
pathname = uigetdir;
% Write output folder text in UI
if ~ pathname == 0
    set(hdls.editOutputFolder,'string',pathname);
elseif pathname==0
    set(hdls.editOutputFolder,'string','')
end
% Scroll to the right in output folder text box
jk=findjobj(hdls.editOutputFolder);
jk.setCaretPosition(length(get(hdls.editOutputFolder,'string')));
set(jk,'FocusLostCallback',{@scrollright,hdls.editOutputFolder});

% --- Executes on selection change in editNorth.
function editNorth_Callback(hObject, eventdata, handles)
% hObject    handle to editNorth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns editNorth contents as cell array
%        contents{get(hObject,'Value')} returns selected item from editNorth
% All handles
global hdls
% Initial- and boundary conditions
global bcNu bcSu
% Get inputs from User Interface
getinputs
% --------------- Evaluate BC and change unit text -----------------------
% If Dirichlet BC on north boundary
if bcNu==1
    % Change text for units
    set(hdls.editNorthUnit,'string','[m/s]')
    % If Neuman BC on north boundary
elseif bcNu==2
    % Change text for units
    set(hdls.editNorthUnit,'string','[m^2/s^2]')
end
% If Dirichlet BC on south boundary
if bcSu==1
    % Change text for units
    set(hdls.editSouthUnit,'string','[m/s]')
    % If Neuman BC on south boundary
elseif bcSu==2
    % Change text for units
    set(hdls.editSouthUnit,'string','[m^2/s^2]')
end
% ------------------------------------------------------------------------

% --- Executes on selection change in editSouth.
function editSouth_Callback(hObject, eventdata, handles)
% hObject    handle to editSouth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns editSouth contents as cell array
%        contents{get(hObject,'Value')} returns selected item from editSouth
editNorth_Callback; % Run editNorth_Callback


function hs = addcomponents
global hs
% Add components to advanced settings window
% Save handles in a struct

%Create figure and keep invisible
hs.fig = figure('Visible','off',...
    'Resize','on',...
    'Tag','fig','position',[360 278 300 220]);
%Create text and textboxes
hs.eostext= uicontrol('Style','text',...
    'Position',[30 200 100 20],...
    'String','Equation of state','fontweight','bold');

hs.roughnestext= uicontrol('Style','text',...
    'Position',[180 200 100 20],...
    'String','Roughness lengths','fontweight','bold');

hs.Miscstext= uicontrol('Style','text',...
    'Position',[180 120 100 20],...
    'String','Miscellaneous','fontweight','bold');

hs.alpha = uicontrol(hs.fig,'Position',[50 180 70 20],...
    'Style','edit',...
    'Tag','editAlpha','string','0.17e-4');

hs.alphatext= uicontrol('Style','text',...
    'Position',[0 180 50 20],...
    'String','Alpha');

hs.beta = uicontrol(hs.fig,'Position',[50 150 70 20],...
    'Style','edit',...
    'Tag','editBeta','string','7.6e-4');

hs.betatext= uicontrol('Style','text',...
    'Position',[0 150 50 20],...
    'String','Beta');

hs.rho0 = uicontrol(hs.fig,'Position',[50 120 70 20],...
    'Style','edit',...
    'Tag','editRho0','string','1000');

hs.rho0text= uicontrol('Style','text',...
    'Position',[0 120 50 20],...
    'String','Rho_0');

hs.T0 = uicontrol(hs.fig,'Position',[50 90 70 20],...
    'Style','edit',...
    'Tag','editT0','string','10');

hs.T0text= uicontrol('Style','text',...
    'Position',[0 90 50 20],...
    'String','T0');

hs.S0 = uicontrol(hs.fig,'Position',[50 60 70 20],...
    'Style','edit',...
    'Tag','editS0','string','35');

hs.S0text= uicontrol('Style','text',...
    'Position',[0 60 50 20],...
    'String','S0');

hs.z0s = uicontrol(hs.fig,'Position',[200 180 70 20],...
    'Style','edit',...
    'Tag','editz0s','string','0.05');

hs.z0stext= uicontrol('Style','text',...
    'Position',[150 180 50 20],...
    'String','Surface');

hs.z0b = uicontrol(hs.fig,'Position',[200 150 70 20],...
    'Style','edit',...
    'Tag','editz0b','string','0.05');

hs.z0btext= uicontrol('Style','text',...
    'Position',[150 150 50 20],...
    'String','Bottom');

hs.g = uicontrol(hs.fig,'Position',[200 90 70 20],...
    'Style','edit',...
    'Tag','editg','string','9.81');

hs.gtext= uicontrol('Style','text',...
    'Position',[150 90 50 20],...
    'String','g');

hs.btn = uicontrol(hs.fig,'Position',[100 30 70 20],...
    'String','Save & Close',...
    'Tag','SaveConstButton',...
    'Callback',@saveconst);

% --- Executes on selection change in popStrat.
function popStrat_Callback(hObject, eventdata, handles)
% hObject    handle to popStrat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hdls stratChoice N2
getinputs

% Linear stratification chosen
if stratChoice==1
    set(hdls.breakText,'string','Calculated N^2 [1/s^2]')
    set(hdls.editBreakDep,'string',num2str(mean(N2)))
    % 2-layer stratification chosen
elseif stratChoice==2
    %set(hdls.editBreakDep,'string','5')
    set(hdls.breakText,'string','Top layer thickness [m]:')
end

% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global reset
% Geometry of grid
global deltaz lg l nodeZ lg2 faceZ lb H
% Marching
global t deltat tmax iteratedtime  iterCounter tend
global ksave Tsave Ssave  usave tsave m  si Pbsave Pssave Epssave nutPsave
global momFluxsave
% Initial- and boundary conditions
global boundValSu bcSu
global boundValNu bcNu
global boundValNk bcNk
global boundValSk bcSk
global boundValST bcST
global boundValNT bcNT
global boundValSS bcSS
global boundValNS bcNS
% Modeled/transport equations
global F_iu u
global F_ik k
global T S
% Source terms
global Ps Eps Pb
% Model constants and stability functions
global  Cmu Cmu2 cb  momFlux Cmu0
% Stratification and equation os state
global alfa beta T0 S0 eosChoice rho rho0
% Miscellaneous
global N2 nutTot nut2Tot nutP nu nup  g Rt drhodz  vcase kar kmin

%Measure time it takes to reach tmax
tic;
% To reset plots
reset=0;
% Get inputs from UI
getinputs

% Reset the iterated time
iteratedtime=0;
% While next timestep won't make iterated time exceed tmax
while iteratedtime+deltat<=tmax
    % Counter for number of iterations
    iterCounter=iterCounter+1;
    % Update timevector
    timetracker
    % --------------------- Length scales --------------------------------
    % Density gradient
    drhodz = -diff(rho)/deltaz;
    % Buoyancy frequency squared
    N2=-g/rho0*drhodz;
    % Byoyancy length scale squared
    lb2=cb.^2*k./N2;
    % Limiting lb to the size of domain
    % (Length scale can not be larger than domain)
    lb2=min(lb2, max(abs(nodeZ)));
    % Buoyancy length scale with sign
    lb = sign(lb2).*sqrt(abs(lb2));
    % Total length scale, neutral
    l = lg;
    % Total length scale, neutral and squared
    l2 = l.^2;                      % Mixing length, total, neutral, squared.
    % Total length scale, stable and squared
    if ~vcase==0
        l(lb2>0) = sqrt(1./(1./lg2(lb>0) + 1./lb2(lb>0))); % Mixing length, total, stable, squared
        l(lb2<0) = sqrt(lg2(lb2<0) - l2(lb2<0).*lg2(lb2<0)./lb2(lb2<0));
    end
    % --------------------------------------------------------------------
    % -------------------- Stability functions ---------------------------
    % Calculate dissipation Epsilon
    Eps=Cmu0.^3.*k.^(3/2)./l;
    % Calculate Turbulent Richardson number
    Rt=k.^2.*N2./(Eps.^2);
    %Stability functions
    Cmu=(Cmu0+0.108*Rt)./(1+0.308.*Rt+0.00837.*Rt.^2);
    Cmu2=Cmu0./(1+0.277.*Rt);
    % --------------------------------------------------------------------
    % ------------ Turbulent viscosity and tirbulent diffusivity ---------
    % Calculate turbulent viscosity nu_t
    nut=Cmu.*sqrt(k).*l;
    % Calculate turbulent diffusivity nu_t^'
    nut2=Cmu2.*sqrt(k).*l;
    % Total turbulent viscosity
    nutTot=nut+nu;
    % Total turbulent diffusivity
    nut2Tot=nut2+nup;
    % Total turbulent viscosity in nodes
    nutP=interp1(faceZ,nutTot,nodeZ(2:end-1));
    % --------------------------------------------------------------------
    % -------------------- Source terms ----------------------------------
    % Velocity shear du/dz in faces
    dudz2=(diff(u)/deltaz).^2;
    % Production, P_s, of turbulent kinetic energy from shear
    Ps=nutTot.*dudz2;
    % Production, Pb, of turbulent kinetic energy from buoyancy
    Pb=nut2Tot.*N2;
    % Calculate momentum flux
    momFlux=nutTot.*-diff(u)./deltaz;
    %---------------------------------------------------------------------
    %%%%%%%%%%%%%%%
    if vcase==0
        ustar00 = sqrt(nutTot(1).*sqrt(dudz2(1)));
        z0 = H/exp(boundValNu*kar/ustar00);
        %l = kar.*min(flipud(faceZ)+z0,faceZ-faceZ(end)+z0);
    end
    %%%%%%%%%%%%%%%
    
    % Mask for positive Pb
    pmask=Pb>0;
    % Mask for negative Pb
    nmask=Pb<0;
    
    %---------------------Discretization coefficients----------------------
    % Common coefficient
    a_P0=1;
    %----------------U--------------------------------------------------
    A_iu=deltat.*nutTot(1:end-1)./(deltaz.^2);
    C_iu=deltat*nutTot(2:end)./(deltaz.^2);
    B_iu=a_P0+A_iu+C_iu;
    F_iu=u(2:end-1);
    % Dirichlet boundary condition, north
    if bcNu==1
        F_iu(1)=F_iu(1)+A_iu(1)*boundValNu;
    end
    if vcase>0
    A_iu(1)=0;
    end
    % Neuman boundary condition, north
    if bcNu==2
        B_iu(1)=a_P0+A_iu(1)+C_iu(1);
        F_iu(1,1)=F_iu(1,1)-boundValNu.*deltat./deltaz;
    end
    % Dirichlet boundary condition, south
    if bcSu==1
        F_iu(end)=F_iu(end)+C_iu(end)*boundValSu;
    end
    
    if vcase==0
        F_iu(end)=u(end);
    end
    
    C_iu(end)=0;
    % Neuman boundary condition, south
    if bcSu==2
        B_iu(end)=a_P0+A_iu(end)+C_iu(end);
        F_iu(end,1)=F_iu(end,1)+boundValSu.*deltat./deltaz;
    end
    % Calculate new timestep
    u(2:end-1,1)=tridiag(B_iu,-A_iu,-C_iu,F_iu);
    % Explicit Neuman boundary condition north
    if bcNu==2
        u(1)=u(2)-boundValNu.*deltaz./nutTot(1);
    end
    % Explicit Neuman boundary condition south
    if bcSu==2
        u(end)=u(end-1)+boundValSu.*deltaz./nutTot(end);
    end
    
    %---------------------------------------------------------------------
    %-------------k-----------------------------------------------------
    if bcNu==2
        ustar2=boundValNu;
    elseif vcase==0
        ustar2=momFlux(1);
    end
    
    boundValNk=Cmu0^(-2)*ustar2;
    A_ik=deltat.*nutP(1:end-1)./deltaz^2;
    C_ik=deltat.*nutP(2:end)./deltaz.^2;
    B_ik=a_P0+A_ik+C_ik+deltat.*Eps(2:end-1)./...
        k(2:end-1)+deltat.*Pb(2:end-1).*pmask(2:end-1)./k(2:end-1);
    F_ik=k(2:end-1)+Ps(2:end-1).*deltat+deltat.*Pb(2:end-1).*nmask(2:end-1);
    % Dirichlet boundary condition, north
    if vcase>0
        F_ik(1)=F_ik(1)+A_ik(1)*boundValNk;
    end

    if vcase==0
        B_ik(1)=B_ik(1)-A_ik(1);
        B_ik(end)=B_ik(end)-C_ik(end);
        F_ik(end)=F_ik(end)+C_ik(end)*k(end);
    end
    
    % Neuman boundary condition, north
    if bcNk==2
        B_ik(1)=B_ik(1)-A_ik(1);
        F_ik(1,1)=F_ik(1,1)-boundValNk.*deltat./deltaz;
    end
    if vcase>0
        A_ik(1)=0;
    end
    % Dirichlet boundary condition, south
    if bcSk==1
        F_ik(end)=F_ik(end)+C_ik(end)*boundValSk;
    end
    if vcase>0
        % Neuman boundary condition, south
        if bcSk==2
            B_ik(end)=B_ik(end)-C_ik(end);
            F_ik(end,1)=F_ik(end,1)+boundValSk.*deltat./deltaz;
        end
        if vcase>0
            C_ik(end)=0;
        end
    end
    % Calculate new timestep
    k(2:end-1,1)=tridiag(B_ik,-A_ik,-C_ik,F_ik);
    if vcase==0
        k(1)=k(2);
        k(end)=nutTot(end).*sqrt(dudz2(end));              %Bottom boundary value for k
        % Neuman boundary condition, south
        k=max(k,kmin);
        if bcSk==2
            %F_ik(end,1)=F_ik(end,1)+boundValSk.*deltat./deltaz;
        end
    end
    
    % Explicit Neuman boundary condition north
    if bcNk==2
        k(1)=k(2)-boundValNk.*deltaz./nutTot(1);
    end
    % Explicit Neuman boundary condition south
    if vcase>0
    if bcSk==2
        k(end)=k(end-1)+boundValSk.*deltaz./nutTot(end);
    end
    end
    % Correct negative k if necessary
    if vcase>0
    k=max(k,1e-20);
    end
    
    %---------------------------------------------------------------------
    %----------------Temp--------------------------------------------------
    A_iT=deltat.*nut2Tot(1:end-1)./(deltaz.^2);
    C_iT=deltat*nut2Tot(2:end)./(deltaz.^2);
    B_iT=a_P0+A_iT+C_iT ;
    F_iT=T(2:end-1);
    % Dirichlet boundary condition, north
    if bcNT==1
        F_iT(1)=F_iT(1)+A_iT(1)*boundValNT;
    end
    A_iT(1)=0;
    % Neuman boundary condition, north
    if bcNT==2
        B_iT(1)=a_P0+A_iT(1)+C_iT(1);
        F_iT(1,1)=F_iT(1,1)-boundValNT.*deltat./deltaz;
    end
    % Dirichlet boundary condition, south
    if bcST==1
        F_iT(end)=F_iT(end)+C_iT(end)*boundValST;
    end
    C_iT(end)=0;
    % Neuman boundary condition, south
    if bcST==2
        B_iT(end)=a_P0+A_iT(end)+C_iT(end);
        F_iT(end,1)=F_iT(end,1)+boundValST.*deltat./deltaz;
    end
    % Calculate new timestep
    T(2:end-1,1)=tridiag(B_iT,-A_iT,-C_iT,F_iT);
    % Explicit Neuman boundary condition north
    if bcNT==2
        T(1)=T(2)-boundValNT.*deltaz./nut2Tot(1);
    end
    % Explicit Neuman boundary condition south
    if bcST==2
        T(end)=T(end-1)+boundValST.*deltaz./nut2Tot(end);
    end
    
    %---------------------------------------------------------------------
    %----------------Salinity--------------------------------------------------
    A_iS=deltat.*nut2Tot(1:end-1)./(deltaz.^2);
    C_iS=deltat*nut2Tot(2:end)./(deltaz.^2);
    B_iS=a_P0+A_iS+C_iS ;
    F_iS=S(2:end-1);
    % Dirichlet boundary condition, north
    if bcNS==1
        F_iS(1)=F_iS(1)+A_iS(1)*boundValNS;
    end
    A_iS(1)=0;
    % Neuman boundary condition, north
    if bcNS==2
        B_iS(1)=a_P0+A_iS(1)+C_iS(1);
        F_iS(1,1)=F_iS(1,1)-boundValNS.*deltat./deltaz;
    end
    % Dirichlet boundary condition, south
    if bcSS==1
        F_iS(end)=F_iS(end)+C_iS(end)*boundValSS;
    end
    C_iS(end)=0;
    % Neuman boundary condition, south
    if bcSS==2
        B_iS(end)=a_P0+A_iS(end)+C_iS(end);
        F_iS(end,1)=F_iS(end,1)+boundValSS.*deltat./deltaz;
    end
    % Calculate new timestep
    S(2:end-1,1)=tridiag(B_iS,-A_iS,-C_iS,F_iS);
    % Explicit Neuman boundary condition north
    if bcNS==2
        S(1)=S(2)-boundValNS.*deltaz./nut2Tot(1);
    end
    % Explicit Neuman boundary condition south
    if bcSS==2
        S(end)=S(end-1)+boundValSS.*deltaz./nut2Tot(end);
    end
    %---------------------------------------------------------------------
    % Calculate density
    % Use linear equation of state if chosen by user
    if strcmp(eosChoice,'buttonlinear')
        rho=rho0.*(1-alfa.*(T-T0)+beta.*(S-S0));
        % Use sw_dens if chosen by user
    elseif strcmp(eosChoice,'buttonswdens')
        rho=sw_dens(S,T,0);
    end
    % Save timestep if appropriate
    if t/si==round(t/si)
        m=m+1;
        % Velocity
        usave(:,m)=u;
        % Time
        tsave(m)=t;
        % Turbulent kinetic energy
        ksave(:,m)=k;
        % Temperature
        Tsave(:,m)=T;
        % Salinity
        Ssave(:,m)=S;
        % Buoyancy production
        Pbsave(:,m)=Pb;
        % Shear production
        Pssave(:,m)=Ps;
        % Dissipation Epsilon
        Epssave(:,m)=Eps;
        % Turbulent viscosity
        nutPsave(:,m)=nutP;
        % Momentum flux
        momFluxsave(:,m)=momFlux;
    end
end

% Update plots in user interface
updateplots
tend=toc;
% Update time counter in user interface
showtrackers

%Executes on button press in SaveConstButton
function saveconst(hObject,event)
global alfa beta rho0 T0 S0 hs z0s z0b
%Collect values from textboxes in advanced settings figure window
%Coefficient of thermal expansion
alfa=str2num(hs.alpha.String);
%Coefficient of saline contraction
beta=str2num(hs.beta.String);
%Reference density
rho0=str2num(hs.rho0.String);
%Reference temprature
T0=str2num(hs.T0.String);
%Reference salinity
S0=str2num(hs.S0.String);
%Roughness length surface
z0s=str2num(hs.z0s.String);
%Roughness length bottom
z0b=str2num(hs.z0b.String);
%Make figure invisible
hs.fig.Visible = 'off';

% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Handle to all handles
global hdls
% Marching
global deltat
global usave tsave ksave Ssave Tsave Pbsave Pssave Epssave nutPsave
global momFluxsave
% Geometry of grid
global deltaz nodeZ faceZ z0s z0b H breakDep
% Initial- and boundary conditions
global boundValSu bcSu
global boundValNu bcNu
global boundValSk bcSk
global boundValNk bcNk
global boundValST bcST
global boundValNT bcNT
global boundValSS bcSS
global boundValNS bcNS
% Stratification and equation of state
global eosChoice stratChoice rho0 T0 S0 alfa beta
% Model constants and stability functions
global cb Cmu0 kar
% Miscellaneous
global g
global nu nup nutP nutTot nut2Tot
%Gets filename from UI
selectedFilename=get(handles.editOutputFile,'string');
%Gets directory from UI
selectedFolder=get(handles.editOutputFolder,'string');
% If user didn't select an output folder choose working directory as output
% folder
if isempty(selectedFolder)
    selectedFolder = pwd;
end
if isempty(selectedFilename)
    % Datenum of current date and time
    timestamp=clock;
    timestampstr=datestr(timestamp,'yyyymmHHMM')
    selectedFilename = ['Output',timestampstr];
end
%Save variables to file
save([selectedFolder '/' selectedFilename],...
    'deltat',...
    'usave', 'tsave', 'ksave', 'Ssave', 'Tsave', 'Pbsave', 'Pssave', 'Epssave', 'nutPsave',...
    'momFluxsave',...
    'deltaz', 'nodeZ', 'faceZ', 'z0s', 'z0b', 'H', 'breakDep',...
    'boundValSu', 'bcSu',...
    'boundValNu', 'bcNu',...
    'boundValSk', 'bcSk',...
    'boundValNk', 'bcNk',...
    'boundValST', 'bcST',...
    'boundValNT', 'bcNT',...
    'boundValSS', 'bcSS',...
    'boundValNS', 'bcNS',...
    'eosChoice', 'stratChoice', 'rho0', 'T0', 'S0', 'alfa', 'beta',...
    'cb', 'Cmu0', 'kar',...
    'g','nu','nup','nutP','nutTot','nut2Tot');
% Print message for user
msgbox(['Output saved to file ', selectedFilename,'.mat in ',selectedFolder])

% --- Executes on button press in setCreate.
function setCreate_Callback(hObject, eventdata, handles)
% hObject    handle to setCreate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hdls
% Modelled/Transport equations
global k T S u
% Geometry of grid and length scale
global nodeZ deltaz faceZ Nnodes lb2 lb l lg lg2
% Source terms
global Pb Ps Eps
% Model constants and stability functions
global Cmu0 Cmu2 Cmu cb
% Initial- and boundary conditions
global T1 T2 S1 S2 inu boundValNk boundValNu bcNu
global nutTot nutP nut2Tot  nu nup
% Equation of state and stratification
global alfa beta  T0 S0 rho0 rho breakDep eosChoice stratChoice
% Matrices for saving output
global Tsave Ssave usave tsave ksave Pbsave Pssave Epssave nutPsave
global momFluxsave
% Miscellaneous variables
global g drhodz  Rt N2 momFlux kmin t reset m dudz2 vcase

% Counter to save u and t
m=1;
% To reset plots
reset=1;
%Matrices for saving modelled/transport equations values
usave=[];
tsave=[];
ksave=[];
Tsave=[];
Ssave=[];
Pbsave=[];
Pssave=[];
Epssave=[];
nutPsave=[];
momFluxsave=[];

% Get user inputs
getinputs
% Set up geometry of grid spanning domain
setupgrid
% Set velocity
u=ones(Nnodes,1).*inu;
T=[];
S=[];
if bcNu==1
    u(1)=boundValNu;
end
% Minimum turbulent kinetic energy k
if vcase==0
kmin=1e-10;
else
    kmin=1e-6;
end
% Minimum turbulent kinetic energy in faces
k=zeros(Nnodes-1,1)+kmin;

% ----------------Set up initial stratification----------------------------
% 2-layer stratification chosen
if stratChoice==2
    % Where to place interface between the layers
    breakNode=round(breakDep/deltaz);
    if breakNode>Nnodes
        disp('Break too deep. 2-layer stratification was not set. Choose a smaller break depth.');
        return
    end
    % Set up initial temperature and salinity distribution
    T(1:breakNode,1)=T1;
    T(breakNode+1:Nnodes,1)=T2;
    S(1:breakNode,1)=S1;
    S(breakNode+1:Nnodes,1)=S2;
    % Linear stratification chosen
elseif stratChoice==1
    % Set up initial temperature and salinity distribution
    T=linspace(T1,T2,Nnodes)';
    S=linspace(S1,S2,Nnodes)';
end
% -------------------------------------------------------------------------
% ---------------------Calculate density-----------------------------------
% Linear equation of state chosen
if strcmp(eosChoice,'buttonlinear')
    rho=rho0.*(1-alfa.*(T-T0)+beta.*(S-S0));
    % sw_dens chosen
elseif strcmp(eosChoice,'buttonswdens')
    rho=sw_dens(S,T,0);
end
% -------------------------------------------------------------------------

% --------------------- Length scales -------------------------------------
% Density gradient
drhodz = -diff(rho)/deltaz;
% Buoyancy frequency squared
N2=-g/rho0*drhodz;
% Geometric length scale squared
lg2=lg.^2;
% Byoyancy length scale squared
lb2=cb.^2*k./N2;
% Limiting lb to size of domain (Length scale can not be larger than domain)
lb2=min(lb2, max(abs(nodeZ)));
% Buoyancy length scale with sign
lb = sign(lb2).*sqrt(abs(lb2));
% Total length scale, neutral.
l = lg;
% Total length scale, neutral and squared
l2 = l.^2;
% Total length scale, stable and squared
l(lb2>0) = sqrt(1./(1./lg2(lb>0) + 1./lb2(lb>0)));
l(lb2<=0) = sqrt(lg2(lb2<=0) - l2(lb2<=0).*lg2(lb2<=0)./lb2(lb2<=0));
% ------------------------------------------------------------------------

% -------------------- Stability functions -------------------------------
% Calculate dissipation Epsilon
Eps=Cmu0.^3.*k.^(3/2)./l;
% Calculate Turbulent Richardson number
Rt=k.^2.*N2./(Eps.^2);
%Stability functions
Cmu=(Cmu0+0.108*Rt)./(1+0.308.*Rt+0.00837.*Rt.^2);
Cmu2=Cmu0./(1+0.277.*Rt);
% ------------------------------------------------------------------------

% ------------ Turbulent viscosity and tirbulent diffusivity -------------
% Molecular viscocity
nu = 1e-6;
% Molecular diffusivity
nup = 1e-7;
% Calculate turbulent viscosity nu_t
nut=Cmu.*sqrt(k).*l;
% Calculate turbulent diffusivity nu_t^'
nut2=Cmu2.*sqrt(k).*l;
% Total turbulent viscosity
nutTot=nut+nu;
% Total turbulent diffusivity
nut2Tot=nut2+nup;
% Total turbulent viscosity in nodes
nutP=interp1(faceZ,nutTot,nodeZ(2:end-1));
% ------------------------------------------------------------------------
% -------------------- Source terms --------------------------------------
% Velocity shear du/dz in faces
dudz2=(diff(u)/deltaz).^2;
% Production, P_s, of turbulent kinetic energy from shear
Ps=nutTot.*dudz2;
% Production, Pb, of turbulent kinetic energy from buoyancy
Pb=-nut2Tot.*N2;
% Calculate momentum flux
momFlux=nutTot.*-diff(u)./deltaz;
%-------------------------------------------------------------------------

% Friction velocity cubed
%ustar3 = sign(boundValNk)*abs(boundValNu)^(3/2);
% Reset trackers
resettrackers
% Show trackers in UI
showtrackers
% Update plots in UI
updateplots
% Update text on reset button in UI


% Save initial conditions for velocity, turbulent kinetic energy,
% Temperature and Salinity as well as production terms and dissipation.
% Velocity
usave(:,m)=u;
% Time
tsave(m)=t;
% Turbulent kinetic energy
ksave(:,m)=k;
% Temperature
Tsave(:,m)=T;
% Salinity
Ssave(:,m)=S;
% Buoyancy production
Pbsave(:,m)=Pb;
% Shear production
Pssave(:,m)=Ps;
% Dissipation Epsilon
Epssave(:,m)=Eps;
% Turbulent viscosity
nutPsave(:,m)=nutP;
% Momentum flux
momFluxsave(:,m)=momFlux;
% Evaluate stratification
popStrat_Callback


% --- Outputs from this function are returned to the command line.
function varargout = Transient_Transport_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
% ------------------------------------------------------------------------



% --- Executes on button press in RebootButton.
function RebootButton_Callback(hObject, eventdata, handles)
% hObject    handle to RebootButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
runSICT
