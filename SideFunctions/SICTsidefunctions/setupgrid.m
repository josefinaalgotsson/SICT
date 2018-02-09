% Author: MSc Josefina Algotsson, University of Gothenburg
% This software is licenced under the GNU General Public License v.3
function setupgrid
global H                  
global Nnodes deltaz                
global faceZ nodeZ                                     
global lg z0b z0s                   
global kar vcase
% ---------------------- SET UP GRID -------------------------------------
% Calculate number of nodes in grid
Nnodes=(H/deltaz)+1; 
% Face positions
faceZ=[H-deltaz/2:-deltaz:deltaz/2]'; 
% Node positions
nodeZ=linspace(H,0,Nnodes)';  
% Geometric length scale
lg=kar.*min(faceZ+z0b,H-faceZ+z0s);   
if vcase==0
    lg=kar*(faceZ-faceZ(end)+z0b);
end
%..........................................................................
