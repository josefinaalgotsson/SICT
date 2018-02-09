% Author: MSc Josefina Algotsson, University of Gothenburg
% This software is licenced under the GNU General Public License v.3
function scrollright(j,ev,hObject)
j.setCaretPosition(length(get(hObject,'string'))); 
