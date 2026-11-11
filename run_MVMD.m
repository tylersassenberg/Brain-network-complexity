%% run_MVMD.m
% Subsets participants, extracts functional timeseries acquired using
% individualized parcellation (GPIP), decomposes them with multivariate
% variational mode decomposition (MVMD) across parcels, averages mode
% functions within the 17 Schaefer/Yeo networks, and saves relative energy
% and central frequency information.
%
% Requires MVMD.m (Rehman & Aftab, 2019) on the MATLAB path.
%
% Outputs (in OUT_DIR):
%   subject_ids.csv                         - subject IDs used by all later scripts
%   MVMD/<id>/<id>_mode<k>.mat, _omega.mat   - parcel-level mode functions
%   network_timeseries/<id>_mode<k>.csv     - network-averaged mode timeseries
%   MVMD_network_relative_energy.csv        - relative energy per network and mode
%   MVMD_central_frequencies.csv            - central frequency of each mode

%%user settings
PARCEL_FILE = 'HCPRest_run1_indiv_RSFC_matrices_v6.mat'; % contains correlation_matrices.id
ANNOT_DIR   = fullfile('data', 'gpip_annot');           % per-subject GPIP timeseries (.dat)
OUT_DIR     = 'output';
MVMD_DIR    = fullfile(OUT_DIR, 'MVMD');
NET_TS_DIR  = fullfile(OUT_DIR, 'network_timeseries');

K = 10;                % number of MVMD modes
ANALYSIS_MODES = 2:5;  % modes carried forward to RQA/MRN analyses

%MVMD parameters: MVMD(signal, alpha, tau, K, DC, init, tol)
ALPHA = 1000; TAU = 0; DC = 0; INIT = 0; TOL = 0;

%%Network definitions (Schaefer 400 parcels, 17 networks; LH and RH parcels)
net_names = {'VisCent','VisPeri','SomMotA','SomMotB','DorsAttnA','DorsAttnB', ...
             'SalVentAttnA','SalVentAttnB','LimbicB','LimbicA','ContA','ContB', ...
             'ContC','DefaultA','DefaultB','DefaultC','TempPar'};
nets = {[1:12,201:212], [13:24,213:223], [25:43,224:243], [44:59,244:258], ...
        [60:72,259:272], [73:85,273:284], [86:100,285:303], [101:108,304:312], ...
        [109:113,313:318], [114:120,319:324], [121:133,325:335], [134:143,336:350], ...
        [144:148,351:357], [149:166,358:373], [167:187,374:384], [188:194,385:390], ...
        [195:200,391:400]};
n_nets = numel(nets);

%%Subject
load(PARCEL_FILE, 'correlation_matrices');

%subjects without missing/poor quality GPIP
subs = [1:7,9:35,37:47,49:60,62:64,66,68,70:80,83:87,89:98,101,103:106,108,109,111,112,114:139,141:145,147,149:161,163:189,191:203,206:210,212:223,225:230,232:248,251:253,255:280,282:302,305:325,327:337,339:342,344:346,348:349,352,354,355,357:361,363,365:388,390:399,401:415,418:428,430:456,458:471,474:475,477:490,492:498,500:502,504:527,529:532,534:566,568:575,577:603,605:634,636:657,659:663,665:670,672:688,690:728,730:736,738:753,755:770,773:796,798:802,804:830,832,834:840,842:850,853:883,886:892,894:915,917:926,928:936,938:947,949:960,962:967,969:991,993:1046,1048:1056,1058:1068,1071:1084];

subjlist = correlation_matrices.id(subs);
subjlist = subjlist(:);
n_subs = numel(subjlist);

if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end
writetable(table(subjlist, 'VariableNames', {'subject_id'}), fullfile(OUT_DIR, 'subject_ids.csv'));

%%MVMD decomposition
for n = 1:n_subs
    sid = subjlist{n};
    fprintf('--- subject %d/%d: %s ---\n', n, n_subs, sid);

    subj_dir = fullfile(MVMD_DIR, sid);
    if ~exist(subj_dir, 'dir'), mkdir(subj_dir); end

    lh_ts = load(fullfile(ANNOT_DIR, sid, [sid, '_run-01_gpip_avgts_lh.dat']));
    rh_ts = load(fullfile(ANNOT_DIR, sid, [sid, '_run-01_gpip_avgts_rh.dat']));
    ts_avg = [lh_ts, rh_ts]';  % parcels x timepoints

    [modes, ~, omega] = MVMD(ts_avg, ALPHA, TAU, K, DC, INIT, TOL);  % modes: K x timepoints x parcels

    %save each mode as a parcels x timepoints matrix in a variable named mode<k>
    for k = 1:K
        mode_name = sprintf('mode%d', k);
        s = struct(mode_name, squeeze(modes(k, :, :))');
        save(fullfile(subj_dir, [sid, '_', mode_name, '.mat']), '-struct', 's');
    end
    save(fullfile(subj_dir, [sid, '_omega.mat']), 'omega');
end

%%Network-averaged timeseries for analysis modes
if ~exist(NET_TS_DIR, 'dir'), mkdir(NET_TS_DIR); end

for k = ANALYSIS_MODES
    for n = 1:n_subs
        sid = subjlist{n};
        mode_ts = load_mode(MVMD_DIR, sid, k);  % timepoints x parcels

        ts_data = zeros(size(mode_ts, 1), n_nets);
        for c = 1:n_nets
            ts_data(:, c) = mean(mode_ts(:, nets{c}), 2);
        end

        writematrix(ts_data, fullfile(NET_TS_DIR, sprintf('%s_mode%d.csv', sid, k)));
    end
end

%%Mode energies
energy_data = zeros(n_subs, K, n_nets);

for n = 1:n_subs
    sid = subjlist{n};
    mode_ts = cell(1, K);
    for k = 1:K
        mode_ts{k} = load_mode(MVMD_DIR, sid, k);
    end

    for c = 1:n_nets
        %network-averaged timeseries for every mode (timepoints x K)
        imfs = cell2mat(cellfun(@(x) mean(x(:, nets{c}), 2), mode_ts, 'UniformOutput', false));

        %energy values
        abs_energy = var(imfs);
        total_net_energy = sum(abs_energy);  % only use with all K modes
        energy_data(n, :, c) = abs_energy / total_net_energy;
    end
end

%save as subjects x (network, mode) table with columns <Network>_RE_mode<k>
energy_table = reshape(energy_data, n_subs, K * n_nets);  %network-major column order
[mode_idx, net_idx] = ndgrid(1:K, 1:n_nets);
energy_names = arrayfun(@(c, k) sprintf('%s_RE_mode%d', net_names{c}, k), net_idx(:), mode_idx(:), ...
                        'UniformOutput', false);
energy_table = [table(subjlist, 'VariableNames', {'subject_id'}), ...
                array2table(energy_table, 'VariableNames', energy_names')];
writetable(energy_table, fullfile(OUT_DIR, 'MVMD_network_relative_energy.csv'));

%%central frequencies
freq_data = zeros(n_subs, K);
for n = 1:n_subs
    sid = subjlist{n};
    omegas = load(fullfile(MVMD_DIR, sid, [sid, '_omega.mat']));
    freq_data(n, :) = omegas.omega(end, :);  %final MVMD iteration (iteration 500 with TOL = 0)
end

freq_table = [table(subjlist, 'VariableNames', {'subject_id'}), ...
              array2table(freq_data, 'VariableNames', arrayfun(@(k) sprintf('Mode%d', k), 1:K, 'UniformOutput', false))];
writetable(freq_table, fullfile(OUT_DIR, 'MVMD_central_frequencies.csv'));

%%local functions
function mode_ts = load_mode(mvmd_dir, sid, k)
%load mode k for a subject as a timepoints x parcels matrix, with NaN/Inf set to 0.
    mode_name = sprintf('mode%d', k);
    s = load(fullfile(mvmd_dir, sid, [sid, '_', mode_name, '.mat']), mode_name);
    mode_ts = s.(mode_name)';
    mode_ts(~isfinite(mode_ts)) = 0;
end
