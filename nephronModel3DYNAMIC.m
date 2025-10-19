function [streams] = nephronModel3DYNAMIC(scenarioName, conc_Na, conc_K, conc_HCO3, conc_Urea, conc_Cl, conc_Glucose, conc_Creat)
% nephronModel3DYNAMIC
% Simulates solute and water transport through the nephron of a healthy kidney.

% Tracks Na+, K+, HCO3-, Urea, Cl-, Glucose, Creatinine, and Water
% through nephron segments (PCT → Desc → Asc → DCT → Cortical CD → Medullary CD).

%% Baseline constants 
GFR_L_per_min = 0.125;        % Typical per-nephron flow, L/min (~125 mL/min total)
GFR = GFR_L_per_min * 60;     % L/hr

% Stream matrix:
% Columns: [TotalFlow Na K HCO3 Urea Cl Glucose Creat Water]
streams = zeros(7, 9);

% Stream 1: entering PCT (glomerular filtrate)
streams(1, 2) = conc_Na * GFR / 1000;      
streams(1, 3) = conc_K * GFR / 1000;       
streams(1, 4) = conc_HCO3 * GFR / 1000;    
streams(1, 5) = conc_Urea * GFR / 1000;    
streams(1, 6) = conc_Cl * GFR / 1000;      
streams(1, 7) = conc_Glucose * GFR / 1000; 
streams(1, 8) = conc_Creat * GFR / 1000;   % Creatinine filtered freely
streams(1, 9) = (1000 * GFR) / 18;         % water molar flow
streams(1, 1) = sum(streams(1, 2:9));

%% Segment fractional reabsorption/secretion 
% Order: [Na  K  HCO3  Urea  Cl  Glucose  Creat  Water]
% Note: Negative values indicate net secretion into tubule.

PCT       = [0.65  0.65  0.85  0.50  0.65  1.00  -0.10  0.65];
DescLoop  = [0.00  0.00  0.00  0.15  0.00  0.00   0.00  0.15];
AscLoop   = [0.25  0.20  0.00  0.00  0.45  0.00   0.00  0.00];
DCT       = [0.075 0.025 0.085 0.00  0.075 0.00   0.00  0.075];
CortCD    = [0.035 0.05  0.045 0.025 0.035 0.00   0.00  0.075];
MedCD     = [0.03  0.045 0.00  0.225 0.015 0.00   0.00  0.04];

reab_segments = [reab_PCT; reab_Desc; reab_Asc; reab_DCT; reab_CortCD; reab_MedCD];

%% Compute downstream flows
for i = 1:6
    reab = reab_segments(i,:);
    streams(i+1, 2:9) = streams(i, 2:9) .* (1 - reab);
    streams(i+1, 1) = sum(streams(i+1, 2:9));
end

%% Compute concentrations (mol/L) 
stream_labels = {'1 (PCT In)','2 (Desc In)','3 (Asc In)','4 (DCT In)', ...
                 '5 (Cort CD In)','6 (Med CD In)','7 (Final Urine)'};
volume_L = streams(:,9) * 18 / 1000;  % convert water molar flow → volume flow
concentrations = streams(:, 2:8) ./ volume_L;

%% Plot solute flow rates 
species_names = {'Na+','K+','HCO3-','Urea','Cl-','Glucose','Creatinine'};
stream_indices = 1:7;

figure('Name', [scenarioName, ': Solute Flow Rates']);
semilogy(stream_indices, streams(:,2:8), 'LineWidth', 2);
title([scenarioName, ': Solute Flow Rates']);
xlabel('Stream Number'); ylabel('Molar Flow Rate (mol/hr)');
legend(species_names, 'Location', 'best');
grid on; xticks(stream_indices); xticklabels(stream_labels); xtickangle(45);

%% Plot solute concentrations 
figure('Name', [scenarioName, ': Solute Concentrations']);
semilogy(stream_indices, concentrations, 'LineWidth', 2);
title([scenarioName, ': Solute Concentrations']);
xlabel('Stream Number'); ylabel('Concentration (mol/L)');
legend(species_names, 'Location', 'best');
grid on; xticks(stream_indices); xticklabels(stream_labels); xtickangle(45);

%% Plot individual solute profiles 
for i = 1:length(species_names)
    figure('Name', [scenarioName, ': ', species_names{i}, ' Flow Rate']);
    semilogy(stream_indices, streams(:, i+1), '-o', 'LineWidth', 2);
    title([scenarioName, ': ', species_names{i}, ' Flow Rate']);
    xlabel('Stream'); ylabel('Molar Flow Rate (mol/hr)');
    grid on; xticks(stream_indices); xticklabels(stream_labels); xtickangle(45);

    figure('Name', [scenarioName, ': ', species_names{i}, ' Concentration']);
    semilogy(stream_indices, concentrations(:, i), '-o', 'LineWidth', 2);
    title([scenarioName, ': ', species_names{i}, ' Concentration']);
    xlabel('Stream'); ylabel('Concentration (mol/L)');
    grid on; xticks(stream_indices); xticklabels(stream_labels); xtickangle(45);
end

%% Summary Table 
summaryTable = table(stream_labels', ...
    streams(:,2), streams(:,3), streams(:,4), streams(:,5), streams(:,6), streams(:,7), streams(:,8), ...
    'VariableNames', {'Segment','Na_mol_hr','K_mol_hr','HCO3_mol_hr','Urea_mol_hr','Cl_mol_hr','Glucose_mol_hr','Creat_mol_hr'});
disp(summaryTable);

%% Final homeostasis comparison: Plasma vs Urine composition 
plasma = [conc_Na, conc_K, conc_HCO3, conc_Urea, conc_Cl, conc_Glucose, conc_Creat];
urine = concentrations(end, :);

percent_diff = ((urine - plasma) ./ plasma) * 100;

figure('Name', [scenarioName, ': Urine vs Plasma Composition']);
bar(categorical(species_names), percent_diff);
ylabel('% Difference (Urine vs Plasma)');
title([scenarioName, ': Urine vs Plasma Concentration Comparison']);
grid on;
yline(0, '--k');
text(1:length(species_names), percent_diff, ...
    compose('%.1f%%', percent_diff), ...
    'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10);
ylim([-100 1000]); % to visualize large ratios like creatinine and urea
xlabel('Solute');
ylabel('% Difference from Plasma');

end
