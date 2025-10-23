% ------------------------------------------------------------
% Hypernatremia Example - Nephron Model Simulation
% ------------------------------------------------------------
% Description:
% This example simulates a nephron model response to hypernatremia,
% defined here as an elevated plasma sodium concentration (155 mM).
% All other solute concentrations are held at baseline values.
% ------------------------------------------------------------

clc; clear; close all;

% Define hypernatremic plasma sodium concentration
hyper_Na = 0.155;  % 155 mM = 0.155 mol/L

% Call nephron model function:
% nephronModel3(condition, conc_Na, conc_K, conc_HCO3, conc_Urea, conc_Cl, conc_Glucose)
[streams_hyper] = nephronModel32('Hypernatremia', hyper_Na, 0.00425, 0.024, 0.00475, 0.101, 0.07);

% Pause to prevent figure window from closing immediately (optional)
pause;

% Notes
% - You can adjust hyper_Na to explore different severities.
% - The baseline concentrations are approximate physiological values:
%   Na+: 0.140 M, K+: 0.004 M, HCO3-: 0.024 M, Urea: 0.005 M,
%   Cl-: 0.100 M, Glucose: 0.007 M
% - Output "streams_hyper" contains all nephron segment data.

% Example for healthy baseline run
% [streams_healthy] = nephronModel32('Baseline', 0.140, 0.00425, 0.024, 0.00475, 0.101, 0.07);
