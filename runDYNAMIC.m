clc;
clear;
close all;

% Normal healthy plasma concentrations (mM)
healthy_Na      = 140;
healthy_K       = 4.25;
healthy_HCO3    = 24;
healthy_Urea    = 4.75;
healthy_Cl      = 101;
healthy_Glucose = 5;
healthy_Creat   = 0.09; % typical plasma creatinine in mM (~1 mg/dL)

[healthy_streams] = nephronModel3DYNAMIC('Healthy State', ...
    healthy_Na, healthy_K, healthy_HCO3, healthy_Urea, ...
    healthy_Cl, healthy_Glucose, healthy_Creat);

pause;

