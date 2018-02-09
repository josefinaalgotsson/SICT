% Author: MSc Josefina Algotsson, University of Gothenburg
% This software is licenced under the GNU General Public License v.3
function resettracker
global t iterCounter
global tvector 
global reset tend

if reset
% Time
t=0;  
% Counter for number of interations
iterCounter=0;
% Vector for saving time
tvector=t;
% Real time whn model is finished running
tend=[];
end
