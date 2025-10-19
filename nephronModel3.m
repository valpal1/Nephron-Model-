function [streams] = nephronModel3(scenarioName, conc_Na, conc_K, conc_HCO3, conc_Urea, conc_Cl, conc_Glucose)

% This model simulates nephron solute transport and flow of a healthy
% kidney I hope 

%% Baseline GFR -> indicates good health 
GFR_L_per_min = 0.125;
GFR = GFR_L_per_min * 60; % L/hr

% Initialize a 7x8 matrix for 7 key streams and 8 species ig (that's wrong
% idk how to explain this matrix) 
streams = zeros(7, 8);

% Stream 1: entering PCT
streams(1, 2) = conc_Na * GFR / 1000;      
streams(1, 3) = conc_K * GFR / 1000;       
streams(1, 4) = conc_HCO3 * GFR / 1000;    
streams(1, 5) = conc_Urea * GFR / 1000;    
streams(1, 6) = conc_Cl * GFR / 1000;      
streams(1, 7) = conc_Glucose * GFR / 1000; 
streams(1, 8) = (1000 * GFR) / 18;         
streams(1, 1) = sum(streams(1, 2:8));

% You can run the next parts by section if needed in the Editor 

%% Segment reabsorption fractions
reab_pct      = [0.65; 0.65; 0.85; 0.50; 0.60; 0.00; 0.66];
reab_desc     = [0; 0; 0; 0.15; 0; 0; 0.15];
reab_asc      = [0.25; 0.20; 0; 0; 0.45; 0; 0];
reab_dct      = [0.075; 0.025; 0.085; 0; 0.075; 0; 0.075];
reab_cort_cd  = [0.035; 0.05; 0.045; 0.025; 0.035; 0; 0.075];
reab_med_cd   = [0.03; 0.045; 0; 0.225; 0.015; 0; 0.04];

remaining = {1-reab_pct, 1-reab_desc, 1-reab_asc, 1-reab_dct, 1-reab_cort_cd, 1-reab_med_cd};

% Compute all downstream flows
for i = 1:6
    streams(i+1, 2:8) = streams(i, 2:8) .* remaining{i}';
    streams(i+1, 1) = sum(streams(i+1, 2:8));
end

%% Compute concentrations
stream_labels = {'1 (PCT In)','2 (Desc In)','3 (Asc In)','4 (DCT In)', ...
                 '5 (Cort. CD In)','6 (Med. CD In)','7 (Final Urine)'};
volume_L = streams(:,8) * 18 / 1000;
concentrations = streams(:, 2:7) ./ volume_L;

%% Plot combined solute flow rates
stream_indices_for_plotting = 1:7;
species_names = {'Na+','K+','HCO3-','Urea','Cl-','Glucose'};

figure('Name', [scenarioName, ': Solute Flow Rates']);
semilogy(stream_indices_for_plotting, streams(:,2:7), 'LineWidth', 2);
title([scenarioName, ': Solute Flow Rates']);
xlabel('Stream Number'); ylabel('Molar Flow Rate (mol/hr)');
legend(species_names, 'Location', 'best');
grid on; xticks(stream_indices_for_plotting); xticklabels(stream_labels); xtickangle(45);

%% Plot combined solute concentrations
figure('Name', [scenarioName, ': Solute Concentrations']);
semilogy(stream_indices_for_plotting, concentrations, 'LineWidth', 2);
title([scenarioName, ': Solute Concentrations']);
xlabel('Stream Number'); ylabel('Concentration (mol/L)');
legend(species_names, 'Location', 'best');
grid on; xticks(stream_indices_for_plotting); xticklabels(stream_labels); xtickangle(45);

%% Plot individual solute graphs
for i = 1:length(species_names)
    % Flow rate
    figure('Name', [scenarioName, ': ', species_names{i}, ' Flow Rate']);
    semilogy(stream_indices_for_plotting, streams(:, i+1), '-o', 'LineWidth', 2);
    title([scenarioName, ': ', species_names{i}, ' Molar Flow Rate']);
    xlabel('Stream Number'); ylabel('Molar Flow Rate (mol/hr)');
    grid on; xticks(stream_indices_for_plotting); xticklabels(stream_labels); xtickangle(45);

    % Concentration
    figure('Name', [scenarioName, ': ', species_names{i}, ' Concentration']);
    semilogy(stream_indices_for_plotting, concentrations(:, i), '-o', 'LineWidth', 2);
    title([scenarioName, ': ', species_names{i}, ' Concentration']);
    xlabel('Stream Number'); ylabel('Concentration (mol/L)');
    grid on; xticks(stream_indices_for_plotting); xticklabels(stream_labels); xtickangle(45);
end
end
