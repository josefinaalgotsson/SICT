function showtrackers
global hdls
global tvector tend

% Update clock showing simulated time
set(hdls.editClock,'string',num2str(tvector(end)))
% Update clock showing simulation time (how long it has taken to run
% simulation)
set(hdls.editSimulationTime,'string',num2str(tend))

