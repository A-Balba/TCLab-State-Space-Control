% RUN_STATE_FEEDBACK
% Run the observer-based TCLab state-feedback experiment.
% The hardware section requires the external TCLab interface used in the lab.

run(fullfile(fileparts(mfilename('fullpath')), '..', 'control', 'observer_state_feedback.m'));
