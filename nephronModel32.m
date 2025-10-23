function [streams] = nephronModel32(condition, conc_Na, conc_K, conc_HCO3, conc_Urea, conc_Cl, conc_Glucose)
% nephronModel3: Dynamic, solute-driven nephron model (full).
%
% Inputs:
%   condition     - string label (e.g., 'Healthy', 'Hypernatremia')
%   conc_X        - input plasma concentrations (mol/L)
%       conc_Na, conc_K, conc_HCO3, conc_Urea, conc_Cl, conc_Glucose
% Output:
%   streams       - 7x8 numeric matrix of stream molar flows (mol/hr)
%
% Overview:
%   - Starling's forces (P_GC, P_BS, pi_GC, pi_BS), together with Kf,
%     determine a dynamic GFR (mL/min).
%   - GFR is clamped to a physiological range to avoid unrealistic values.
%   - Segmental reabsorption fractions (PCT, Desc, Asc, DCT, Cortical CD,
%     Medullary CD) are defined as baseline values 
%   - Those baseline reabsorption fractions are scaled by a function of the
%     new GFR so that: higher GFR -> reduced fractional reabsorption,
%     lower GFR -> increased fractional reabsorption (physiological:
%     less transit time reduces fractional reabsorption).
%   - Each solute has its own sensitivity to GFR-driven reabsorption changes
%     (e.g., glucose reabsorption saturable vs Na reabsorption more linear).
%
% Notes on units & conversions:
%   - Concentrations (inputs) are mol/L (e.g., 0.140 mol/L = 140 mM Na+)
%   - GFR is computed in mL/min, converted to mL/hr when computing mol/hr
%   - We report molar flow rates (mol/hr) per stream 
%
% Model limitations / assumptions:
%   - pi_BS (oncotic pressure in Bowman's space) is ~0 in healthy kidneys
%   - We use simple proportional relations linking solute deviation to P_GC
%     and pi_GC. These are phenomenological and tunable (alpha, beta).
%   - Reabsorption modulation by GFR is a simplified number
%

%% 0) Basic print
fprintf('\n=== nephronModel32 : %s ===\n', condition);


%% 1) STARLING FORCES & DYNAMIC GFR (solutes drive pressures)
% Baseline Starling parameters (typical textbook-ish values)
P_GC0  = 55;   % mmHg, glomerular capillary hydrostatic pressure baseline
P_BS0  = 15;   % mmHg, Bowman capsule hydrostatic pressure baseline
pi_GC0 = 30;   % mmHg, glomerular oncotic pressure baseline (plasma proteins)
pi_BS0 = 0;    % mmHg, Bowman space oncotic (negligible normally)
Kf = 12.5; % mL / min / mmHg, whole-kidney filtration coefficient (typical)

% Baseline concentrations for computing relative deviation (mol/L)
baseline = struct( ...
    'Na', 0.140, ...    % 140 mM
    'K', 0.00425, ...   % 4.25 mM
    'HCO3', 0.024, ...  % 24 mM
    'Urea', 0.00475, ...% 4.75 mM
    'Cl', 0.101, ...    % 101 mM
    'Glucose', 0.07);   % 70 mM? (you previously used 0.07; keep consistent)

% Compute fractional deviations from baseline (dimensionless)
% (positive -> above baseline; negative -> below baseline)
Na_frac  = (conc_Na - baseline.Na) / baseline.Na;
K_frac   = (conc_K  - baseline.K)  / baseline.K;
Cl_frac  = (conc_Cl - baseline.Cl) / baseline.Cl;
HCO3_frac= (conc_HCO3 - baseline.HCO3)/baseline.HCO3;
Urea_frac= (conc_Urea - baseline.Urea)/baseline.Urea;
Glu_frac = (conc_Glucose - baseline.Glucose)/baseline.Glucose;

% Sensitivities: how much a solute deviation influences pressures
% These are phenomenological tuning parameters. They should be adjusted
% if you calibrate the model to experimental data.
alpha_Na = 0.20;   % Na+ effect on hydrostatic pressure (P_GC)
alpha_K  = 0.05;   % K+ effect on hydrostatic pressure

% extra small sensitivities so computed fractions are actually used
alpha_HCO3 = 0.02; % small HCO3 effect on hydrostatic pressure
alpha_Urea = 0.01; % small urea effect on hydrostatic pressure
alpha_Glu  = 0.03; % small glucose effect on hydrostatic pressure

% Oncotic sensitivity: changes in plasma ionic strength/protein affec
beta_protein = 0.15; % Cl or plasma-proxy effect on oncotic pressure

% small additional oncotic contribution from urea (used as a minor proxy)
beta_urea = 0.05;

% Adjust pressures 
% Rationale:
%  - Increased Na and K (reflecting increased effective circulating volume/osmolality)
%    tend to increase intravascular hydrostatic pressures at the renal level.
%  - Increased 'protein proxy' (we use Cl as a rough proxy just to allow a
%    variable) increases oncotic pressure.
P_GC  = P_GC0 * (1 + alpha_Na*Na_frac + alpha_K*K_frac + alpha_HCO3*HCO3_frac + alpha_Urea*Urea_frac + alpha_Glu*Glu_frac);
pi_GC = pi_GC0 * (1 + beta_protein*Cl_frac + beta_urea*Urea_frac);

% Bowman space pressures left constant for now 
P_BS  = P_BS0;
pi_BS = pi_BS0;

% Compute GFR via Starling (mL/min)
GFR_dynamic = Kf * ((P_GC - P_BS) - (pi_GC - pi_BS));

% Physiological clamp: keep GFR inside reasonable bounds to avoid blow-ups
GFR_min = 50;   % mL/min, allow slightly low (hypofiltration)
GFR_max = 200;  % mL/min, allow hyperfiltration for simulation
GFR_mL_per_min = max(min(GFR_dynamic, GFR_max), GFR_min);

% Useful baseline GFR for scaling reabsorption (choose physiological)
GFR_baseline = 125; % mL/min: used as "normal" reference for scaling since this is a healthy number 

fprintf('\n--- Starling summary ---\n');
fprintf('P_GC = %.2f mmHg | P_BS = %.2f mmHg | pi_GC = %.2f mmHg | pi_BS=%.2f\n', P_GC, P_BS, pi_GC, pi_BS);
fprintf('GFR (calc) = %.2f mL/min  (clamped to [%.0f, %.0f])\n', GFR_mL_per_min, GFR_min, GFR_max);

%% 2) Convert GFR to volumetric flow used by your tubular calculations
% In your previous code you used GFR (mL/min) -> GFR (mL/hr) -> mol/hr flows.
GFR_mL_per_hr = GFR_mL_per_min * 60; % mL/hr
% You compute water molar flow as GFR/18 (because 1 mL water ~1 g, molar mass 18 g/mol)
% so GFR (mL/hr)/18 -> mol H2O / hr. We'll keep same convention.

%% 3) Segmental baseline reabsorption fractions (your original numbers)
% The order of solutes in streams/blood arrays: [total, Na, K, HCO3, Urea, Cl, Glucose, H2O]
% baseline fraction reabsorbed in each segment for each solute (7 rows -> segments)
% segments = [PCT; Desc; Asc; DCT; Cortical CD; Medullary CD]
reab_pct_baseline = [0.65; 0.5707; 0.85; 0.25; 0.60; 1.00; 0.66]; % for PCT
% 1=Na,2=K,3=HCO3,4=Urea,5=Cl,6=Glucose,7=H2O  (we'll store without total column)

baseline_reab = zeros(6,7); % 6 tubular segments (PCT, Desc, Asc, DCT, CortCD, MedCD)
% PCT (segment 1 -> streams(2) output in your original)
baseline_reab(1,:) = reab_pct_baseline'; % use the reab_pct_baseline vector for the PCT row
% Descending limb (segment 2)
baseline_reab(2,:) = [0, 0, 0, 0.1, 0, 0, (15/34)];
% Ascending limb (segment 3)
baseline_reab(3,:) = [(5/7), 0.409, 0, (1/9), (5/8), 0, 0];
% DCT (segment 4)
baseline_reab(4,:) = [0.5, 0.0865, (17/30), 0, 0.5, 0, 15/38];
% Cortical collecting duct (segment 5)
baseline_reab(5,:) = [0.4, 0.1894, (9/13), 0.03125, (7/15), 0, 0.075];
% Medullary collecting duct (segment 6)
baseline_reab(6,:) = [0.5, 0.2103, 0.25, (9/47), (3/8), 0, (8/13)];

% NOTE: baseline_reab rows correspond to outputs of each segment (i.e., fraction
% reabsorbed from the incoming stream by that segment).

%% 4) GFR-dependent modulation of reabsorption
% Rationale:
%  - Reabsorption often depends on transit time and available driving forces.
%  - Higher GFR -> shorter transit time -> less can be reabsorbed (fraction falls).
%  - Lower GFR -> longer transit -> more fractional reabsorption (up to a max).
% Implementation:
%  - For each solute we use a different sensitivity parameter (gamma). The
%    reabsorption at the new GFR is: reab_new = baseline_reab * modulation_factor(solute)
%  - modulation_factor = 1 / (1 + gamma * (GFR_new / GFR_baseline - 1))
%    -> If GFR_new > baseline, factor < 1 (reabsorption reduced).
%    -> If GFR_new < baseline, factor > 1 (reab increases), but we clamp to < 0.99.

gamma = struct('Na', 1.0, 'K', 0.6, 'HCO3', 0.8, 'Urea', 0.3, 'Cl', 0.7, 'Glucose', 2.0, 'H2O', 1.0);
GFR_ratio = GFR_mL_per_min / GFR_baseline;  % >1 => hyperfiltration; <1 => hypofiltration

% compute modulation factors per solute
mod_Na   = 1./(1 + gamma.Na   * (GFR_ratio - 1));
mod_K    = 1./(1 + gamma.K    * (GFR_ratio - 1));
mod_HCO3 = 1./(1 + gamma.HCO3 * (GFR_ratio - 1));
mod_Urea = 1./(1 + gamma.Urea * (GFR_ratio - 1));
mod_Cl   = 1./(1 + gamma.Cl   * (GFR_ratio - 1));
mod_Glu  = 1./(1 + gamma.Glucose*(GFR_ratio - 1));
mod_H2O  = 1./(1 + gamma.H2O  * (GFR_ratio - 1));

% pack mods into vector aligned with solute columns
mod_vec = [mod_Na, mod_K, mod_HCO3, mod_Urea, mod_Cl, mod_Glu, mod_H2O];

% Apply modulation to baseline_reab to get actual reabsorption fractions for this run
reab_actual = baseline_reab .* (ones(size(baseline_reab,1),1) * mod_vec);

% Clamp reabsorption fractions between 0 and 0.99 (avoid >1 or negative)
reab_actual = max(min(reab_actual, 0.99), 0);

fprintf('\n--- Reabsorption modulation (per-solute factors) ---\n');
fprintf('mod_Na=%.3f, mod_K=%.3f, mod_HCO3=%.3f, mod_Urea=%.3f, mod_Cl=%.3f, mod_Glu=%.3f\n', ...
    mod_Na, mod_K, mod_HCO3, mod_Urea, mod_Cl, mod_Glu);

%% 5) Build stream/blood arrays and compute flows (mol/hr)
% streams: 7 x 8 matrix (rows = streams 1..7, cols = [total, Na, K, HCO3, Urea, Cl, Glucose, H2O])
streams = zeros(7,8);
blood   = zeros(7,8);

% Stream 1 = fluid entering PCT (from glomerular filtrate)
% convert GFR (mL/hr) into mol/hr using conc (mol/L) and volume conversions:
% GFR_mL_per_hr / 1000 = L/hr ; conc (mol/L) * (L/hr) = mol/hr
streams(1,2) = conc_Na   * (GFR_mL_per_hr/1000);
streams(1,3) = conc_K    * (GFR_mL_per_hr/1000);
streams(1,4) = conc_HCO3 * (GFR_mL_per_hr/1000);
streams(1,5) = conc_Urea * (GFR_mL_per_hr/1000);
streams(1,6) = conc_Cl   * (GFR_mL_per_hr/1000);
streams(1,7) = conc_Glucose * (GFR_mL_per_hr/1000);
% water molar flow: approximate 1 mL water = 1 g = 1/18 mol
streams(1,8) = (GFR_mL_per_hr) / 18;
streams(1,1) = sum(streams(1,2:8));

% Blood (plasma leaving filtration) approximate: total plasma flow - filtered amount
% We approximate plasma flow as 1100 mL/min * 60 = 66000 mL/hr (this is rough)
% Keep your prior blood calculation pattern but safer: RPF_approx in mL/hr
RPF_approx_mL_per_min = 600; % typical renal plasma flow (mL/min)
RPF_approx_mL_per_hr  = RPF_approx_mL_per_min * 60;
remaining_plasma_L_per_hr = (RPF_approx_mL_per_hr - GFR_mL_per_hr)/1000; % L/hr
% For mol/hr in blood: conc * L/hr
blood(1,2) = conc_Na * remaining_plasma_L_per_hr;
blood(1,3) = conc_K  * remaining_plasma_L_per_hr;
blood(1,4) = conc_HCO3 * remaining_plasma_L_per_hr;
blood(1,5) = conc_Urea * remaining_plasma_L_per_hr;
blood(1,6) = conc_Cl * remaining_plasma_L_per_hr;
blood(1,7) = conc_Glucose * remaining_plasma_L_per_hr;
blood(1,8) = remaining_plasma_L_per_hr * 1000 / 18; % convert L/hr to mL/hr then to mol/hr for H2O approx
blood(1,1) = sum(blood(1,2:8));

%% 6) Process each tubular segment sequentially using reab_actual
% Segment mapping:
% input streams -> segment -> output streams
% Stream 1 -> PCT -> Stream 2
% Stream 2 -> Desc -> Stream 3
% Stream 3 -> Asc -> Stream 4
% Stream 4 -> DCT -> Stream 5
% Stream 5 -> Cort CD -> Stream 6
% Stream 6 -> Med CD -> Stream 7

% For indexing convenience, solute columns in streams: [total, Na, K, HCO3, Urea, Cl, Glucose, H2O]
% reab_actual rows correspond to segments 1..6 and columns to solutes 1..7 (Na..H2O)
for seg = 1:6
    % incoming stream index = seg (1..6)
    in_idx = seg;
    out_idx = seg + 1;
    % get reabsorption fractions for this segment (1x7)
    reab_frac = reab_actual(seg, :);
    % compute reabsorbed mol/hr for each solute column
    % map reab_frac columns to streams columns (Na->2, K->3, ... H2O->8)
    reab_vector = zeros(1,8);
    reab_vector(2:8) = streams(in_idx, 2:8) .* reab_frac; % solute-specific reabsorption (mol/hr)
    % output stream mol/hr = input - reabsorbed
    streams(out_idx, 2:8) = streams(in_idx, 2:8) - reab_vector(2:8);
    streams(out_idx,1) = sum(streams(out_idx, 2:8));
    % update blood (add reabsorbed amounts)
    blood(out_idx, 2:8) = blood(in_idx, 2:8) + reab_vector(2:8);
    blood(out_idx,1) = sum(blood(out_idx, 2:8));
end

%% 7) Concentrations and tables 
% convert water molar flow to L for concentration denominator:
volume_L = streams(:,8) * 18 / 1000; % streams(:,8) in mol/hr -> convert back to L/hr for concentration
% Note: earlier you used streams(:,8) * 18 /1000 which is reverse of mol->L,
% Keeping consistent with your prior code to avoid breaking table layout.
concentrations = zeros(7,6);
% avoid division by zero
for i = 1:7
    if volume_L(i) <= 0
        concentrations(i,:) = NaN(1,6);
    else
        concentrations(i,1) = streams(i,2) / volume_L(i); % C_Na
        concentrations(i,2) = streams(i,3) / volume_L(i); % C_K
        concentrations(i,3) = streams(i,4) / volume_L(i); % C_HCO3
        concentrations(i,4) = streams(i,5) / volume_L(i); % C_Urea
        concentrations(i,5) = streams(i,6) / volume_L(i); % C_Cl
        concentrations(i,6) = streams(i,7) / volume_L(i); % C_Glucose
    end
end

%% 8) Print tables (molar flows and concentrations)
stream_labels = {'1 (PCT In)', '2 (Desc In)', '3 (Asc In)', '4 (DCT In)', '5 (Cort. CD In)', '6 (Med. CD In)', '7 (Final Urine)'};
species_labels_n = {'n_Na+', 'n_K+', 'n_HCO3-', 'n_Urea', 'n_Cl-', 'n_Glucose', 'n_H2O'};
species_labels_c = {'C_Na+', 'C_K+', 'C_HCO3-', 'C_Urea', 'C_Cl-', 'C_Glucose'};

fprintf('\n\n TABLE 1: Molar Flow Rates (mol/hr) \n');
fprintf('%-16s', 'Stream');
for j = 1:length(species_labels_n), fprintf('\t%s', species_labels_n{j}); end, fprintf('\n');
for i = 1:7
    fprintf('%-16s', stream_labels{i});
    fprintf('\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.4f\n', streams(i, 2:8), streams(i,1));
end

fprintf('\n\n TABLE 2: Solute Concentrations (mol/L) \n');
fprintf('%-16s', 'Stream');
for j = 1:length(species_labels_c), fprintf('\t%s', species_labels_c{j}); end, fprintf('\n');
for i = 1:7
    fprintf('%-16s', stream_labels{i});
    if any(isnan(concentrations(i,:)))
        fprintf('\t   NaN\t   NaN\t   NaN\t   NaN\t   NaN\t   NaN\n');
    else
        fprintf('\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\n', concentrations(i,:));
    end
end


%% 9) more plotting 
stream_indices_for_plotting = 1:7;
figure('Name', [condition, ': Solute Flow Rates (Log)']);
semilogy(stream_indices_for_plotting, streams(:,2), '-s', 'LineWidth', 2, 'DisplayName', 'Na+'); hold on;
semilogy(stream_indices_for_plotting, streams(:,6), '-^', 'LineWidth', 2, 'DisplayName', 'Cl-');
semilogy(stream_indices_for_plotting, streams(:,5), '-p', 'LineWidth', 2, 'DisplayName', 'Urea');
semilogy(stream_indices_for_plotting, streams(:,3), '-d', 'LineWidth', 2, 'DisplayName', 'K+');
semilogy(stream_indices_for_plotting, streams(:,4), '-h', 'LineWidth', 2, 'DisplayName', 'HCO3-');
semilogy(stream_indices_for_plotting, streams(:,7), '-x', 'LineWidth', 2, 'DisplayName', 'Glucose');
hold off;
title([condition, ': Solute Flow Rates (Log Scale)']);
xlabel('Stream Number (Input to Segment)'); ylabel('Molar Flow Rate (mol/hr)');
legend('show', 'Location', 'southwest'); grid on; xticks(stream_indices_for_plotting); xticklabels(stream_labels); xtickangle(45);

figure('Name', [condition, ': Solute Concentrations']);
plot(stream_indices_for_plotting, concentrations(:,1), '-s', 'LineWidth', 2, 'DisplayName', 'Na+'); hold on;
plot(stream_indices_for_plotting, concentrations(:,5), '-^', 'LineWidth', 2, 'DisplayName', 'Cl-');
plot(stream_indices_for_plotting, concentrations(:,4), '-p', 'LineWidth', 2, 'DisplayName', 'Urea');
plot(stream_indices_for_plotting, concentrations(:,2), '-d', 'LineWidth', 2, 'DisplayName', 'K+');
plot(stream_indices_for_plotting, concentrations(:,3), '-h', 'LineWidth', 2, 'DisplayName', 'HCO3-');
plot(stream_indices_for_plotting, concentrations(:,6), '-x', 'LineWidth', 2, 'DisplayName', 'Glucose');
hold off;
title([condition, ': Solute Concentrations Along the Nephron']);
xlabel('Stream Number (Input to Segment)'); ylabel('Concentration (mol/L)');
legend('show', 'Location', 'best'); grid on; xticks(stream_indices_for_plotting); xticklabels(stream_labels); xtickangle(45);

%% 10) Return streams matrix for downstream use
% streams kept in same shape as prior code for compatibility
% and also returned so you can further process in calling scripts
end

% -------------------------
% Example runs below (uncomment to run directly from this file)
% -------------------------
% Healthy baselfunction [streams] = nephronModel32(condition, conc_Na, conc_K, conc_HCO3, conc_Urea, conc_Cl, conc_Glucose)

% Hypernatremia example (ready to run):
% Here we use an elevated plasma sodium concentration (e.g., 155 mM -> 0.155 mol/L).
% Other solutes kept at baseline; you can vary them to explore model responses.
% To run the hypernatremia example, uncomment the lines below and execute the file.
% clc; clear; close all;
% hyper_Na = 0.155; % 155 mM -> 0.155 mol/L
% [streams_hyper] = nephronModel3('Hypernatremia', hyper_Na, 0.00425, 0.024, 0.00475, 0.101, 0.07);
% pause;
