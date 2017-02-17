function timetracker
global iteratedtime deltat 
global t 
global tvector iterCounter

% The updated elapsed modeled time                            
iteratedtime=iteratedtime+deltat;       
% Calc. accumulated time for clock display
t=t+deltat;   
% Put time into timevector for saving
tvector(iterCounter+1)=t;
