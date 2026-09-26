% RUN_SYSTEM_IDENTIFICATION
% Run linear and nonlinear TCLab system identification using supplied data.
% Place data.mat in the MATLAB working directory before running.

run(fullfile(fileparts(mfilename('fullpath')), '..', 'identification', 'd2m1_system_identification.m'));
