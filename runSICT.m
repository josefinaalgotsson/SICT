% Close UI figure
delete(findall(0,'Type','figure'))
close all 
clear all
set(0,'DefaultFigureWindowStyle','normal')
% 
%Old path
a=path;
%Old path as cell array
oldPath=strread(a,'%s','delimiter',':'); 
% Current directory
a=cd;
% Rows to remove since they could contain Transient_Transport.m files
removerows=strmatch(a(1:end-2),oldPath);
% Remove paths 

addpath(genpath(cd))

% Run SICT
Transient_Transport

global hdls
