function getinputs
global hdls hs
% Geometry of grid
global deltaz H z0b z0s 
%Model constants and stability functions
global Cmu0 
% Initial- and boundary conditions
global T1 T2 S1 S2 
global bcSu bcNu    
global boundValSu boundValNu
global bcSk bcNk
global boundValSk boundValNk
global bcNT bcST
global boundValST boundValNT
global bcSS bcNS
global boundValSS boundValNS
global inu
% Equation of state and stratification
global alfa beta T0 S0 rho0 breakDep eosChoice stratChoice 
% Marching
global deltat tmax si
% Miscellaneous
global g cb ustar2 kar vcase

%----------------- GET INPUTS FROM USER INTERFACE -------------------------
% ------- Get inputted boundary condition for velocity u -----------------
% Boundary condition kind on south side
bcSu=get(hdls.editSouth,'value');
% Value (of u or of ustar) on south side
boundValNu=str2num(get(hdls.editNorthValue,'string'));
% Boundary condition kind on north side
bcNu=get(hdls.editNorth,'value');
% Value (of u or of ustar) on south side
boundValSu=str2num(get(hdls.editSouthValue,'string'));
% ------------------------------------------------------------------------
% ----------- Get grid dimension and marching ----------------------------
% Get the inputed grid size (height of each cell)
deltaz=str2num(get(hdls.editGridSize,'String')); 
% Get the inputted height H of the domain
H=str2num(get(hdls.editDomHeight,'String')); 
% Get timestep
deltat=str2num(get(hdls.editTimeStep,'string'));     
% Get t-max (time to run the model)
tmax=str2num(get(hdls.editTmax,'String')); 
% Get save interval
si=str2num(get(hdls.editSaveInterval,'string'));  
% Get roughness length on south side
z0s=str2num(hs.z0s.String);  
% Get roughness length on south side
z0b=str2num(hs.z0b.String);    
% ------------------------------------------------------------------------
% ------------  Get inputs about stratification --------------------------
% Get choice about equation of state (either use linear or sw_dens)
eosChoice=get(hdls.eosGrouping.SelectedObject,'Tag');   
% Get choice about stratification 1) Linear temp and sal 2) 2-layer temp and sal 
stratChoice=get(hdls.popStrat,'value');
% Get depth at which 2-layer interface is to be put
breakDep=str2num(get(hdls.editBreakDep,'String')); 
% Reference density for linear equation of state
rho0=str2num(hs.rho0.String);
%Reference temperature for linear equation of state
T0=str2num(hs.T0.String);
%Reference salinity for linear equation of state
S0=str2num(hs.S0.String);
% Coefficient of thermal expansion for linear equation of state
alfa=str2num(hs.alpha.String);
% Coefficient of saline contraction for linear equation of state
beta=str2num(hs.beta.String);
% ------------------------------------------------------------------------
% ------ Get initial conditions for temperature and salinity -------------
T1=str2num(get(hdls.editT1,'string'));
T2=str2num(get(hdls.editT2,'string'));
S1=str2num(get(hdls.editS1,'string'));
S2=str2num(get(hdls.editS2,'string'));
% ------------------------------------------------------------------------
% Get gravitational acceleration
g=str2num(hs.g.String);
%--------------------------------------------------------------------------

%------------------------------SET VALUES----------------------------------
% --------- Set boundary conditions for turbulent kinetic energy k -------
% Boundary condition kind on north side
bcNk=1;  
% Boundary condition kind on south side
bcSk=2;         %BC kind for k
% Value of kflux on south side
boundValSk=0; 
%Friction velocity
if bcNu==2
    ustar2 = boundValNu;
elseif bcNu==1
    ustar2=boundValNu;
end

%ustar3 = sign(ustar3)*abs(ustar3)^(3/2);
%Stability function
Cmu0=0.556;
% Value of k on north side
boundValNk=Cmu0^(-2)*ustar2;
% ------------------------------------------------------------------------
% ---------------- Set boundary conditions for temperature T -------------
% Boundary condition kind on north side
bcNT=2;
% Value of T on north side
boundValNT=0;
% Boundary condition kind on south side
bcST=2;
% Value of T on south side
boundValST=0;
% ------------------------------------------------------------------------
% ---------------- Set boundary conditions for salinity S ----------------
% Boundary condition kind on north side
bcNS=2;
% Value of S on north side
boundValNS=0;
% Boundary condition kind on south side
bcSS=2;
% Value of S on south side
boundValSS=0;
% ------------------------------------------------------------------------
% -------------------- Set miscellaneous ---------------------------------
%Set initial velocity 
inu=0;          
% Model constant
cb=0.35;
% Von karmans constant   
kar=0.4;                            

if T1==T2 && S1==S2 && bcNu==1
    vcase=0;
elseif stratChoice==1
    vcase=1;
elseif stratChoice==2
    vcase=2;
end
%--------------------------------------------------------------------------
