clc;
clear;
close all;

% Normal healthy concentrations
healthy_Na    = 140;
healthy_K     = 4.25;
healthy_HCO3  = 24;
healthy_Urea  = 4.75;
healthy_Cl    = 101;
healthy_Glucose = 5; % example glucose value because we don't have values!

% Call the nephron model (you need to do this first to run the whole code) 
[healthy_streams] = nephronModel3('Healthy State', ...
    healthy_Na, healthy_K, healthy_HCO3, healthy_Urea, healthy_Cl, healthy_Glucose);

pause;
