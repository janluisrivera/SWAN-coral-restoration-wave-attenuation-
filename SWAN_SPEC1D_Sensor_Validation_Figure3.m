%% ========================================================================
% SWAN SPEC1D vs Sensor Validation
% TIMESERIES + SCATTER FIGURES
% ------------------------------------------------------------------------
% What this matlab version does:
% - Reads BOTH T1P3 and T1P5 sensors used for this research
% - Skips NODATA blocks in SWAN SPEC1D
% - Computes band-limited bulk parameters in 0.05-0.35 Hz for model + sensor
% - Hs = 4*sqrt(m0)
% - Tp from 3-point parabolic peak fit
% - Tm01 = m0/m1
% - Uses PER = Tm01 and labels it as energy period, T_e
% - Keeps Measured = blue and SWAN = magenta
% - Forces consistent limits between timeseries y-axes and scatter x/y axes
% - Adds r, RMSE, and Bias to each scatter panel 
% - Keeps full Bias/RMSE/Corr values in a separate summary table (.csv)
% - BONUS: It can also export a separate Peak Period "Tp" validation figure if desired
% - Allows manual validation START time and END time
%% ========================================================================

clear; close all; clc

%% -----------------------
% USER INPUTS
% -----------------------
rho = 1025;
g   = 9.81;

fa = 0.05;
fb = 0.35;

retime_method = 'linear';      % 'linear' or 'mean'
ASSUME_SWAN_IS_Jm2Hz = true;   % SWAN header commonly says J/m2/Hz
MIN_M0 = 1e-6;                 % minimum band energy threshold

% ----- Trim controls -----
trimStartHours = 0;
trimEndHours   = 0;

% Optional: hard-start the plotted/validated window at a specific start time
useManualValidationStartTime = true;
manualValidationStartTime    = datetime(2025,5,14,21,0,0,'TimeZone','UTC');

% Optional: hard-stop the plotted/validated window at a specific end time
useManualValidationEndTime = true;
manualValidationEndTime    = datetime(2025,7,10,12,0,0,'TimeZone','UTC');

% Optional: reject very small Hs values in Hs stats
minHsForStats = 0.05;

% ----- Choose what the second paper figure shows -----
% 'PER' = Average Period = Tm01
% 'Tp'  = Peak Period
paperPeriodMetric = 'PER';

% ----- BONUS -----
% If true, also export a separate Tp validation figure even when the
% main paper-period figure is set to PER.
saveBonusTpFigure = false;

% ----- Figure layout -----
% 'combined4rows' = one paper-style figure with 4 rows:
%                   H_s(T1P3), H_s(T1P5), T_e(T1P3), T_e(T1P5)
% 'twoFigures'    = two separate figures:
%                   one for H_s and one for T_e / T_p
figureLayoutMode = 'twoFigures';

% ----- Save outputs -----
savePaperFigurePNG = true;
saveStatsCSV       = true;
saveStatsXLSX      = false;

% ----- Line/figure style -----
colSensor  = [0 0 1];          % blue
colSWAN    = [1 0 1];          % magenta
colScatter = colSWAN;          % scatter-point color
lwSensor   = 1.4;
lwSWAN     = 1.4;
fsTitle    = 16;               % row titles
fsAxis     = 14;               % axis labels + ticks
fsPanel    = 18;               % panel letters
fsMain     = 18;               % main figure title (if enabled)

% ----- Paired timeseries + scatter style -----
showMainFigureTitle    = true;
scatterMarkerSize      = 18;
scatterMarkerFaceAlpha = 0.35;
scatterMarkerEdgeAlpha = 0.35;

% Extra spacing requested between the upper Hs panels (a-d) and the lower
% Te panels (e-h). 'loose' adds more breathing room than 'compact'.
tileSpacingMode   = 'loose';
tilePaddingMode   = 'loose';
figureHeightPerRow = 320;  % larger value = more vertical space

% ----- Stats text box placement in scatter panels -----
statsX = 0.97;   % move left/right (0 to 1)
statsY = 0.25;   % move up/down    (0 to 1)

% ----- Cases -----
cases(1).name       = 'T1P3';
cases(1).spc1d_file = "D:\Paper_Thesis\Validation\New_1-30MINS\spec_T1P3_20250514_20250521.spc1d";
cases(1).sensor_mat = "D:\LPT1P320250513Pb_waves_035_heff_13m.mat";

cases(2).name       = 'T1P5';
cases(2).spc1d_file = "D:\Paper_Thesis\Validation\New_1-30MINS\spec_T1P5_20250514_20250521.spc1d";
cases(2).sensor_mat = "D:\T1P520250513Pb_waves_035_heff_13m.mat";

nCases = numel(cases);

%% -----------------------
% STORAGE
%% -----------------------
R = struct([]);

%% -----------------------
% PROCESS EACH CASE
%% -----------------------
for ic = 1:nCases

    fprintf('\n============================================================\n');
    fprintf('Processing case: %s\n', cases(ic).name);
    fprintf('============================================================\n');

    %% -----------------------
    % Load Sensor
    %% -----------------------
    S = load(cases(ic).sensor_mat, 'Pwav');
    Pwav = S.Pwav;

    f_obs = Pwav.f(:);
    E_obs = Pwav.E;      % [Nt x Nf] in m^2/Hz
    t_obs = Pwav.tt(:);

    try
        t_obs.TimeZone = 'UTC';
    catch
    end

    if ~isempty(t_obs)
        fprintf('First sensor time:      %s\n', char(string(t_obs(1))));
        fprintf('Last sensor time:       %s\n', char(string(t_obs(end))));
    end

    idxBand = find(f_obs >= fa & f_obs <= fb);
    if isempty(idxBand)
        error('No sensor frequencies in band %.2f-%.2f Hz for %s', fa, fb, cases(ic).name);
    end
    fB = f_obs(idxBand);

    Nt_obs   = numel(t_obs);
    Hs_obs   = nan(Nt_obs,1);
    Tp_obs   = nan(Nt_obs,1);
    Tm01_obs = nan(Nt_obs,1);

    for it = 1:Nt_obs
        SB = double(E_obs(it,idxBand)).';
        SB(~isfinite(SB) | SB < 0) = 0;

        m0 = trapz(fB, SB);
        if m0 < MIN_M0
            continue
        end

        Hs_obs(it) = 4 * sqrt(max(m0,0));

        fp = peak_freq_parabolic(fB, SB);
        if isfinite(fp) && fp > 0
            Tp_obs(it) = 1 / fp;
        end

        m1 = trapz(fB, fB .* SB);
        if m1 > 0
            Tm01_obs(it) = m0 / m1;
        end
    end

    %% -----------------------
    % Load SWAN SPEC1D
    %% -----------------------
    [f_swan, t_swan, EnDens_swan, nSkipped] = read_swan_spc1d_skip_nodata(char(cases(ic).spc1d_file));

    fprintf('Valid SWAN spectra read: %d\n', numel(t_swan));
    fprintf('NODATA blocks skipped:   %d\n', nSkipped);
    fprintf('First valid SWAN time:   %s\n', char(string(t_swan(1))));
    fprintf('Last valid SWAN time:    %s\n', char(string(t_swan(end))));

    if ~isempty(t_obs) && ~isempty(t_swan)
        if isempty(t_obs(1).TimeZone) && ~isempty(t_swan(1).TimeZone)
            t_swan.TimeZone = '';
        end
    end

    if ASSUME_SWAN_IS_Jm2Hz
        Seta_swan = EnDens_swan ./ (rho*g);
    else
        Seta_swan = EnDens_swan;
    end

    %% -----------------------
    % Interpolate SWAN spectra to sensor frequency grid
    %% -----------------------
    Nt_mod = numel(t_swan);
    Nf_obs = numel(f_obs);

    Seta_swan_interp = nan(Nt_mod, Nf_obs);

    for it = 1:Nt_mod
        Si = double(Seta_swan(it,:));
        goodf = isfinite(Si) & isfinite(f_swan(:)');

        if nnz(goodf) < 2
            Seta_swan_interp(it,:) = 0;
        else
            Seta_swan_interp(it,:) = interp1(f_swan(goodf), Si(goodf), f_obs, 'linear', 0);
        end
    end

    %% -----------------------
    % Model band-limited bulk params
    %% -----------------------
    Hs_model   = nan(Nt_mod,1);
    Tp_model   = nan(Nt_mod,1);
    Tm01_model = nan(Nt_mod,1);

    for it = 1:Nt_mod
        SB = double(Seta_swan_interp(it,idxBand)).';
        SB(~isfinite(SB) | SB < 0) = 0;

        m0 = trapz(fB, SB);
        if m0 < MIN_M0
            continue
        end

        Hs_model(it) = 4 * sqrt(max(m0,0));

        fp = peak_freq_parabolic(fB, SB);
        if isfinite(fp) && fp > 0
            Tp_model(it) = 1 / fp;
        end

        m1 = trapz(fB, fB .* SB);
        if m1 > 0
            Tm01_model(it) = m0 / m1;
        end
    end

    %% -----------------------
    % Align sensor to model times
    %% -----------------------
    sensTT = timetable(t_obs, Hs_obs, Tp_obs, Tm01_obs, ...
        'VariableNames', {'Hs_obs','Tp_obs','Tm01_obs'});

    TT_obs_at_mod = retime(sensTT, t_swan, retime_method);

    t_common = t_swan(:);
    Hs_sens_common   = TT_obs_at_mod.Hs_obs(:);
    Tp_sens_common   = TT_obs_at_mod.Tp_obs(:);
    Tm01_sens_common = TT_obs_at_mod.Tm01_obs(:);

    %% -----------------------
    % Trim start/end + optional hard start/hard stop
    %% -----------------------
    trimMask = true(size(t_common));

    if ~isempty(t_common)
        t0 = t_common(1) + hours(trimStartHours);
        t1 = t_common(end) - hours(trimEndHours);

        if useManualValidationStartTime
            tStartManual = manualValidationStartTime;

            if isempty(t_common(1).TimeZone)
                try
                    tStartManual.TimeZone = '';
                catch
                end
            elseif isempty(tStartManual.TimeZone)
                try
                    tStartManual.TimeZone = t_common(1).TimeZone;
                catch
                end
            end

            if t_common(1) > tStartManual
                warning('%s starts at %s, which is later than the requested start time %s. Check that you are using the full input files.', ...
                    cases(ic).name, char(string(t_common(1))), char(string(tStartManual)));
            end

            if tStartManual > t0
                t0 = tStartManual;
            end
        end

        if useManualValidationEndTime
            tEndManual = manualValidationEndTime;

            if isempty(t_common(1).TimeZone)
                try
                    tEndManual.TimeZone = '';
                catch
                end
            elseif isempty(tEndManual.TimeZone)
                try
                    tEndManual.TimeZone = t_common(1).TimeZone;
                catch
                end
            end

            if t_common(end) < tEndManual
                warning('%s ends at %s, which is earlier than the requested end time %s. Check that you are using the full input files.', ...
                    cases(ic).name, char(string(t_common(end))), char(string(tEndManual)));
            end

            if tEndManual < t1
                t1 = tEndManual;
            end
        end

        if t1 < t0
            error('After applying manual start/end times, start is later than end for %s.', cases(ic).name);
        end

        trimMask = (t_common >= t0) & (t_common <= t1);
    end

    t_plot = t_common(trimMask);

    Hs_model_plot   = Hs_model(trimMask);
    Tp_model_plot   = Tp_model(trimMask);
    Tm01_model_plot = Tm01_model(trimMask);

    Hs_sens_plot   = Hs_sens_common(trimMask);
    Tp_sens_plot   = Tp_sens_common(trimMask);
    Tm01_sens_plot = Tm01_sens_common(trimMask);

    if isempty(t_plot)
        error('After trimming, no data remain for %s. Reduce trimStartHours and/or check manual validation start/end times.', cases(ic).name);
    end

    fprintf('First plotted time:      %s\n', char(string(t_plot(1))));
    fprintf('Last plotted time:       %s\n', char(string(t_plot(end))));

    PER_model_plot = Tm01_model_plot;
    PER_sens_plot  = Tm01_sens_plot;

    %% -----------------------
    % Metrics
    %% -----------------------
    HsStats  = calc_metrics(Hs_model_plot,  Hs_sens_plot,  minHsForStats, 0);
    TpStats  = calc_metrics(Tp_model_plot,  Tp_sens_plot,  0, 0);
    PERStats = calc_metrics(PER_model_plot, PER_sens_plot, 0, 0);

    fprintf('Hs  Bias = %.6f m\n', HsStats.bias);
    fprintf('Hs  RMSE = %.6f m\n', HsStats.rmse);
    fprintf('Hs  Corr = %.6f\n',  HsStats.corr);

    fprintf('Tp  Bias = %.6f s\n', TpStats.bias);
    fprintf('Tp  RMSE = %.6f s\n', TpStats.rmse);
    fprintf('Tp  Corr = %.6f\n',  TpStats.corr);

    fprintf('PER Bias = %.6f s\n', PERStats.bias);
    fprintf('PER RMSE = %.6f s\n', PERStats.rmse);
    fprintf('PER Corr = %.6f\n',  PERStats.corr);

    R(ic).name        = cases(ic).name;
    R(ic).t           = t_plot;
    R(ic).tStart      = t_plot(1);
    R(ic).tEnd        = t_plot(end);

    R(ic).Hs_model    = Hs_model_plot;
    R(ic).Tp_model    = Tp_model_plot;
    R(ic).Tm01_model  = Tm01_model_plot;
    R(ic).PER_model   = PER_model_plot;

    R(ic).Hs_sensor   = Hs_sens_plot;
    R(ic).Tp_sensor   = Tp_sens_plot;
    R(ic).Tm01_sensor = Tm01_sens_plot;
    R(ic).PER_sensor  = PER_sens_plot;

    R(ic).HsStats     = HsStats;
    R(ic).TpStats     = TpStats;
    R(ic).PERStats    = PERStats;

    R(ic).nSkipped    = nSkipped;
end

%% -----------------------
% Choose which period-like variable goes into the second paper figure/table
%% -----------------------
switch upper(string(paperPeriodMetric))
    case "TP"
        paperPeriodMetric = 'Tp';
        paperPeriodTitleText = 'Peak Period';
        paperPeriodYLabel    = 'Peak Period (s)';
        paperPeriodShortName = 'Tp';
        statsFileCSV         = 'Validation_Stats_Hs_Tp.csv';
        statsFileXLSX        = 'Validation_Stats_Hs_Tp.xlsx';

        for ic = 1:nCases
            R(ic).PaperPeriod_model  = R(ic).Tp_model;
            R(ic).PaperPeriod_sensor = R(ic).Tp_sensor;
            R(ic).PaperPeriodStats   = R(ic).TpStats;
        end

    otherwise
        paperPeriodMetric = 'PER';
        paperPeriodTitleText = 'Energy Period';
        paperPeriodYLabel    = 'T_e (s)';
        paperPeriodShortName = 'Te';
        statsFileCSV         = 'Validation_Stats_Hs_PER.csv';
        statsFileXLSX        = 'Validation_Stats_Hs_PER.xlsx';

        for ic = 1:nCases
            R(ic).PaperPeriod_model  = R(ic).PER_model;
            R(ic).PaperPeriod_sensor = R(ic).PER_sensor;
            R(ic).PaperPeriodStats   = R(ic).PERStats;
        end
end

%% -----------------------
% Global x-limits for consistent panels
%% -----------------------
globalStart = R(1).tStart;
globalEnd   = R(1).tEnd;
for ic = 2:nCases
    if R(ic).tStart < globalStart
        globalStart = R(ic).tStart;
    end
    if R(ic).tEnd > globalEnd
        globalEnd = R(ic).tEnd;
    end
end

HsSharedLims = compute_common_limits({R(1).Hs_sensor, R(1).Hs_model, ...
                                      R(2).Hs_sensor, R(2).Hs_model}, true);

PeriodSharedLims = compute_common_limits({R(1).PaperPeriod_sensor, R(1).PaperPeriod_model, ...
                                          R(2).PaperPeriod_sensor, R(2).PaperPeriod_model}, false);

TpSharedLims = compute_common_limits({R(1).Tp_sensor, R(1).Tp_model, ...
                                      R(2).Tp_sensor, R(2).Tp_model}, false);

%% ========================================================================
% BUILD ROW DEFINITIONS FOR MAIN PAPER FIGURES
%% ========================================================================
if strcmpi(paperPeriodMetric, 'Tp')
    periodScatterXLabel = 'Measured T_p (s)';
    periodScatterYLabel = 'Model T_p (s)';
else
    periodScatterXLabel = 'Measured T_e (s)';
    periodScatterYLabel = 'Model T_e (s)';
end

rowDef = struct([]);

rowDef(1).t             = R(1).t;
rowDef(1).obs           = R(1).Hs_sensor;
rowDef(1).mod           = R(1).Hs_model;
rowDef(1).station       = R(1).name;
rowDef(1).yLabel        = 'H_s (m)';
rowDef(1).scatterXLabel = 'Measured H_s (m)';
rowDef(1).scatterYLabel = 'Model H_s (m)';
rowDef(1).lims          = HsSharedLims;
rowDef(1).corrText      = R(1).HsStats.corr;
rowDef(1).rmseText      = R(1).HsStats.rmse;
rowDef(1).biasText      = R(1).HsStats.bias;
rowDef(1).statUnits     = 'm';

rowDef(2).t             = R(2).t;
rowDef(2).obs           = R(2).Hs_sensor;
rowDef(2).mod           = R(2).Hs_model;
rowDef(2).station       = R(2).name;
rowDef(2).yLabel        = 'H_s (m)';
rowDef(2).scatterXLabel = 'Measured H_s (m)';
rowDef(2).scatterYLabel = 'Model H_s (m)';
rowDef(2).lims          = HsSharedLims;
rowDef(2).corrText      = R(2).HsStats.corr;
rowDef(2).rmseText      = R(2).HsStats.rmse;
rowDef(2).biasText      = R(2).HsStats.bias;
rowDef(2).statUnits     = 'm';

rowDef(3).t             = R(1).t;
rowDef(3).obs           = R(1).PaperPeriod_sensor;
rowDef(3).mod           = R(1).PaperPeriod_model;
rowDef(3).station       = R(1).name;
rowDef(3).yLabel        = paperPeriodYLabel;
rowDef(3).scatterXLabel = periodScatterXLabel;
rowDef(3).scatterYLabel = periodScatterYLabel;
rowDef(3).lims          = PeriodSharedLims;
rowDef(3).corrText      = R(1).PaperPeriodStats.corr;
rowDef(3).rmseText      = R(1).PaperPeriodStats.rmse;
rowDef(3).biasText      = R(1).PaperPeriodStats.bias;
rowDef(3).statUnits     = 's';

rowDef(4).t             = R(2).t;
rowDef(4).obs           = R(2).PaperPeriod_sensor;
rowDef(4).mod           = R(2).PaperPeriod_model;
rowDef(4).station       = R(2).name;
rowDef(4).yLabel        = paperPeriodYLabel;
rowDef(4).scatterXLabel = periodScatterXLabel;
rowDef(4).scatterYLabel = periodScatterYLabel;
rowDef(4).lims          = PeriodSharedLims;
rowDef(4).corrText      = R(2).PaperPeriodStats.corr;
rowDef(4).rmseText      = R(2).PaperPeriodStats.rmse;
rowDef(4).biasText      = R(2).PaperPeriodStats.bias;
rowDef(4).statUnits     = 's';

combinedFigureFile     = ['Fig_Validation_Hs_', paperPeriodShortName, '_Combined_TimeseriesScatter.png'];
paperHsScatterFile     = 'Figure3_Top_Hs_T1P3_T1P5.png';
paperPeriodScatterFile = 'Figure3_Bottom_Te_T1P3_T1P5.png';

switch lower(char(string(figureLayoutMode)))
    case 'combined4rows'
        if showMainFigureTitle
            combinedTitle = ['Model Validation: H_s and ', paperPeriodTitleText];
        else
            combinedTitle = '';
        end

        make_timeseries_scatter_paper_figure( ...
            rowDef, [globalStart globalEnd], combinedTitle, combinedFigureFile, ...
            savePaperFigurePNG, colSensor, colSWAN, colScatter, lwSensor, lwSWAN, ...
            fsAxis, fsTitle, fsPanel, fsMain, 'A', ...
            scatterMarkerSize, scatterMarkerFaceAlpha, scatterMarkerEdgeAlpha, ...
            statsX, statsY, tileSpacingMode, tilePaddingMode, figureHeightPerRow);

    case 'twofigures'
        if showMainFigureTitle
            titleHs = 'Model Validation: Significant Wave Height';
        else
            titleHs = '';
        end

        make_timeseries_scatter_paper_figure( ...
            rowDef(1:2), [globalStart globalEnd], titleHs, paperHsScatterFile, ...
            savePaperFigurePNG, colSensor, colSWAN, colScatter, lwSensor, lwSWAN, ...
            fsAxis, fsTitle, fsPanel, fsMain, 'A', ...
            scatterMarkerSize, scatterMarkerFaceAlpha, scatterMarkerEdgeAlpha, ...
            statsX, statsY, tileSpacingMode, tilePaddingMode, figureHeightPerRow);

        if showMainFigureTitle
            titlePer = ['Model Validation: ', paperPeriodTitleText];
        else
            titlePer = '';
        end

        make_timeseries_scatter_paper_figure( ...
            rowDef(3:4), [globalStart globalEnd], titlePer, paperPeriodScatterFile, ...
            savePaperFigurePNG, colSensor, colSWAN, colScatter, lwSensor, lwSWAN, ...
            fsAxis, fsTitle, fsPanel, fsMain, 'E', ...
            scatterMarkerSize, scatterMarkerFaceAlpha, scatterMarkerEdgeAlpha, ...
            statsX, statsY, tileSpacingMode, tilePaddingMode, figureHeightPerRow);

    otherwise
        error('Unknown figureLayoutMode = %s. Use ''combined4rows'' or ''twoFigures''.', figureLayoutMode);
end

%% ========================================================================
% BONUS FIGURE: SEPARATE TP VALIDATION FIGURE
%% ========================================================================
if saveBonusTpFigure
    rowDefTp = struct([]);

    rowDefTp(1).t             = R(1).t;
    rowDefTp(1).obs           = R(1).Tp_sensor;
    rowDefTp(1).mod           = R(1).Tp_model;
    rowDefTp(1).station       = R(1).name;
    rowDefTp(1).yLabel        = 'T_p (s)';
    rowDefTp(1).scatterXLabel = 'Measured T_p (s)';
    rowDefTp(1).scatterYLabel = 'Model T_p (s)';
    rowDefTp(1).lims          = TpSharedLims;
    rowDefTp(1).corrText      = R(1).TpStats.corr;
    rowDefTp(1).rmseText      = R(1).TpStats.rmse;
    rowDefTp(1).biasText      = R(1).TpStats.bias;
    rowDefTp(1).statUnits     = 's';

    rowDefTp(2).t             = R(2).t;
    rowDefTp(2).obs           = R(2).Tp_sensor;
    rowDefTp(2).mod           = R(2).Tp_model;
    rowDefTp(2).station       = R(2).name;
    rowDefTp(2).yLabel        = 'T_p (s)';
    rowDefTp(2).scatterXLabel = 'Measured T_p (s)';
    rowDefTp(2).scatterYLabel = 'Model T_p (s)';
    rowDefTp(2).lims          = TpSharedLims;
    rowDefTp(2).corrText      = R(2).TpStats.corr;
    rowDefTp(2).rmseText      = R(2).TpStats.rmse;
    rowDefTp(2).biasText      = R(2).TpStats.bias;
    rowDefTp(2).statUnits     = 's';

    bonusTpFile = 'Fig_Validation_Tp_T1P3_T1P5_TimeseriesScatter.png';

    if showMainFigureTitle
        titleTp = 'Model Validation: Peak Period';
    else
        titleTp = '';
    end

    make_timeseries_scatter_paper_figure( ...
        rowDefTp, [globalStart globalEnd], titleTp, bonusTpFile, ...
        savePaperFigurePNG, colSensor, colSWAN, colScatter, lwSensor, lwSWAN, ...
        fsAxis, fsTitle, fsPanel, fsMain, 'I', ...
        scatterMarkerSize, scatterMarkerFaceAlpha, scatterMarkerEdgeAlpha, ...
        statsX, statsY);
end

%% ========================================================================
% SUMMARY TABLE FOR PAPER
%% ========================================================================
Station   = strings(3*nCases,1);
Parameter = strings(3*nCases,1);
Units     = strings(3*nCases,1);
Bias      = nan(3*nCases,1);
RMSE      = nan(3*nCases,1);
Corr      = nan(3*nCases,1);
N         = nan(3*nCases,1);

row = 0;
for ic = 1:nCases
    row = row + 1;
    Station(row)   = string(R(ic).name);
    Parameter(row) = "Hs";
    Units(row)     = "m";
    Bias(row)      = R(ic).HsStats.bias;
    RMSE(row)      = R(ic).HsStats.rmse;
    Corr(row)      = R(ic).HsStats.corr;
    N(row)         = R(ic).HsStats.n;

    row = row + 1;
    Station(row)   = string(R(ic).name);
    Parameter(row) = "Tp";
    Units(row)     = "s";
    Bias(row)      = R(ic).TpStats.bias;
    RMSE(row)      = R(ic).TpStats.rmse;
    Corr(row)      = R(ic).TpStats.corr;
    N(row)         = R(ic).TpStats.n;

    row = row + 1;
    Station(row)   = string(R(ic).name);
    Parameter(row) = "PER";
    Units(row)     = "s";
    Bias(row)      = R(ic).PERStats.bias;
    RMSE(row)      = R(ic).PERStats.rmse;
    Corr(row)      = R(ic).PERStats.corr;
    N(row)         = R(ic).PERStats.n;
end

StatsTable = table(Station, Parameter, Units, Bias, RMSE, Corr, N);

fprintf('\n==================== PAPER SUMMARY TABLE ====================\n');
disp(StatsTable)

if saveStatsCSV
    writetable(StatsTable, statsFileCSV);
end

if saveStatsXLSX
    try
        writetable(StatsTable, statsFileXLSX);
    catch ME
        warning('Could not write %s: %s', statsFileXLSX, ME.message);
    end
end

fprintf('\n==================== FINAL SUMMARY ====================\n');
for ic = 1:nCases
    fprintf('%s:\n', R(ic).name);
    fprintf('   Plotted start = %s\n', char(string(R(ic).tStart)));
    fprintf('   Plotted end   = %s\n', char(string(R(ic).tEnd)));
    fprintf('   Hs  Bias = %.6f m | RMSE = %.6f m | Corr = %.6f\n', ...
        R(ic).HsStats.bias, R(ic).HsStats.rmse, R(ic).HsStats.corr);
    fprintf('   Tp  Bias = %.6f s | RMSE = %.6f s | Corr = %.6f\n', ...
        R(ic).TpStats.bias, R(ic).TpStats.rmse, R(ic).TpStats.corr);
    fprintf('   PER Bias = %.6f s | RMSE = %.6f s | Corr = %.6f\n', ...
        R(ic).PERStats.bias, R(ic).PERStats.rmse, R(ic).PERStats.corr);
    fprintf('   NODATA skipped = %d\n', R(ic).nSkipped);
end

%% ========================================================================
% Helper: parabolic peak frequency estimate
%% ========================================================================
function fp = peak_freq_parabolic(f, S)
    fp = NaN;

    if isempty(f) || isempty(S) || all(~isfinite(S))
        return
    end

    [~, k] = max(S);

    if isempty(k) || ~isfinite(S(k)) || S(k) <= 0
        return
    end

    if k <= 1 || k >= numel(S) || any(S([k-1 k k+1]) <= 0)
        fp = f(k);
        return
    end

    y1 = log(S(k-1)); y2 = log(S(k)); y3 = log(S(k+1));
    x1 = f(k-1);      x2 = f(k);      x3 = f(k+1);

    denom = (x1-x2)*(x1-x3)*(x2-x3);
    if abs(denom) < eps
        fp = f(k);
        return
    end

    a = (x3*(y2-y1) + x2*(y1-y3) + x1*(y3-y2)) / denom;
    b = (x3^2*(y1-y2) + x2^2*(y3-y1) + x1^2*(y2-y3)) / denom;

    if ~isfinite(a) || ~isfinite(b) || a == 0
        fp = f(k);
        return
    end

    x_peak = -b/(2*a);
    x_peak = max(min(x_peak, x3), x1);

    fp = x_peak;
end

%% ========================================================================
% Helper: validation metrics
%% ========================================================================
function S = calc_metrics(modelVec, obsVec, minModel, minObs)

    if nargin < 3 || isempty(minModel)
        minModel = -inf;
    end
    if nargin < 4 || isempty(minObs)
        minObs = -inf;
    end

    mask = isfinite(modelVec) & isfinite(obsVec) & ...
           modelVec > minModel & obsVec > minObs;

    x = modelVec(mask);
    y = obsVec(mask);

    S = struct('bias', NaN, 'rmse', NaN, 'corr', NaN, 'n', numel(x));

    if isempty(x)
        return
    end

    S.bias = mean(x - y);
    S.rmse = sqrt(mean((x - y).^2));

    xm = mean(x);
    ym = mean(y);
    denom = sqrt(sum((x - xm).^2) * sum((y - ym).^2));

    if denom > 0
        S.corr = sum((x - xm) .* (y - ym)) / denom;
    else
        S.corr = NaN;
    end

    S.n = numel(x);
end

%% ========================================================================
% Helper: paper-style figure with wide timeseries and narrow scatter panel
%% ========================================================================
function fig = make_timeseries_scatter_paper_figure(rowDef, xRange, mainTitleStr, outFile, savePNG, ...
    colObs, colMod, colScatter, lwObs, lwMod, fsAxis, fsTitle, fsPanel, fsMain, startPanelChar, ...
    scatterMarkerSize, scatterMarkerFaceAlpha, scatterMarkerEdgeAlpha, statsX, statsY, ...
    tileSpacingMode, tilePaddingMode, figureHeightPerRow)

    nRows = numel(rowDef);
    panelChars = char(startPanelChar + (0:(2*nRows - 1)));

    % Increased height and loose tile spacing provide additional separation
    % between panels (c,d) and (e,f) in the combined four-row figure.
    figHeight = max(800, figureHeightPerRow*nRows + 150);
    fig = figure('Color','w', 'Position', [60 30 1500 figHeight]);

    tl = tiledlayout(fig, nRows, 4, ...
        'TileSpacing', tileSpacingMode, ...
        'Padding', tilePaddingMode);

    tsAxes = gobjects(nRows,1);
    ip = 0;

    for ir = 1:nRows
        tsAxes(ir) = nexttile(tl, [1 3]);
        [hObs, hMod] = plot_timeseries_panel_with_limits( ...
            tsAxes(ir), rowDef(ir).t, rowDef(ir).obs, rowDef(ir).mod, ...
            rowDef(ir).yLabel, rowDef(ir).station, ...
            colObs, colMod, lwObs, lwMod, fsAxis, fsTitle, ...
            (ir == nRows), xRange, rowDef(ir).lims);

        if ir == 1
            try
                legend(tsAxes(ir), [hObs hMod], {'Measured','Model'}, ...
                    'Location','northeast', 'Box','on', 'FontSize', max(fsAxis-1,10), ...
                    'AutoUpdate','off');
            catch
                legend(tsAxes(ir), [hObs hMod], {'Measured','Model'}, ...
                    'Location','northeast', 'Box','on', 'FontSize', max(fsAxis-1,10));
            end
        end

        ip = ip + 1;
        add_panel_label_paper(tsAxes(ir), panelChars(ip), fsPanel);

        axSc = nexttile(tl);
        plot_scatter_panel_1to1( ...
            axSc, rowDef(ir).obs, rowDef(ir).mod, ...
            rowDef(ir).scatterXLabel, rowDef(ir).scatterYLabel, ...
            rowDef(ir).lims, rowDef(ir).corrText, rowDef(ir).rmseText, ...
            rowDef(ir).biasText, rowDef(ir).statUnits, ...
            colScatter, scatterMarkerSize, scatterMarkerFaceAlpha, ...
            scatterMarkerEdgeAlpha, fsAxis, statsX, statsY);

        ip = ip + 1;
        add_panel_label_paper(axSc, panelChars(ip), fsPanel);
    end

    linkaxes(tsAxes, 'x');

    for k = 1:numel(tsAxes)
        try
            tsAxes(k).XAxis.TickLabelFormat = 'MMM dd';
        catch
        end
    end

    if ~isempty(mainTitleStr)
        title(tl, mainTitleStr, 'FontWeight','bold', 'FontSize', fsMain);
    end

    if savePNG
        exportgraphics(fig, outFile, 'Resolution', 600);
    end
end

%% ========================================================================
% Helper: left-side timeseries panel
%% ========================================================================
function [hObs, hMod] = plot_timeseries_panel_with_limits(ax, t, yObs, yMod, yLabelStr, titleStr, ...
    colObs, colMod, lwObs, lwMod, fsAxis, fsTitle, showXLabel, xRange, yLims)

    axes(ax); %#ok<LAXES>
    hold(ax, 'on')

    hObs = plot(ax, t, yObs, '-', 'Color', colObs, 'LineWidth', lwObs);
    hMod = plot(ax, t, yMod, '-', 'Color', colMod, 'LineWidth', lwMod);

    grid(ax, 'off');
    box(ax, 'on');
    set(ax, 'FontSize', fsAxis, ...
            'LineWidth', 1.0, ...
            'TickDir', 'out', ...
            'Layer', 'top');

    ylabel(ax, yLabelStr, 'FontSize', fsAxis);
    title(ax, titleStr, 'FontWeight','bold', 'FontSize', fsTitle);

    if nargin >= 14 && ~isempty(xRange)
        xlim(ax, xRange);
    elseif ~isempty(t)
        xlim(ax, [t(1) t(end)]);
    end

    if nargin >= 15 && ~isempty(yLims)
        ylim(ax, yLims);
    end

    if showXLabel
        xlabel(ax, 'Date', 'FontSize', fsAxis);
    else
        ax.XTickLabel = [];
    end
end

%% ========================================================================
% Helper: right-side scatter panel with 1:1 line
%% ========================================================================
function plot_scatter_panel_1to1(ax, xObs, yMod, xlab, ylab, lims, rVal, rmseVal, biasVal, statUnits, ...
    markerColor, markerSize, markerFaceAlpha, markerEdgeAlpha, fsAxis, statsX, statsY)

    xObs = xObs(:);
    yMod = yMod(:);

    mask = isfinite(xObs) & isfinite(yMod);
    x = xObs(mask);
    y = yMod(mask);

    axes(ax); %#ok<LAXES>
    hold(ax, 'on')
    box(ax, 'on')
    grid(ax, 'on')
    set(ax, 'FontSize', fsAxis, ...
            'LineWidth', 1.0, ...
            'TickDir', 'out', ...
            'Layer', 'top');

    plot(ax, lims, lims, 'k-', 'LineWidth', 1.2);

    if ~isempty(x)
        try
            scatter(ax, x, y, markerSize, ...
                'MarkerFaceColor', markerColor, ...
                'MarkerEdgeColor', markerColor, ...
                'MarkerFaceAlpha', markerFaceAlpha, ...
                'MarkerEdgeAlpha', markerEdgeAlpha, ...
                'Marker', 'o', 'LineWidth', 0.5);
        catch
            scatter(ax, x, y, markerSize, ...
                'MarkerFaceColor', markerColor, ...
                'MarkerEdgeColor', markerColor, ...
                'Marker', 'o', 'LineWidth', 0.5);
        end
    else
        text(ax, 0.5, 0.5, 'No valid paired data', ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle');
    end

    xlim(ax, lims);
    ylim(ax, lims);
    try
        pbaspect(ax, [1 1 1]);
    catch
        axis(ax, 'square');
    end

    if isfinite(rVal)
        rStr = sprintf('r = %.2f', rVal);
    else
        rStr = 'r = NaN';
    end

    if isfinite(rmseVal)
        rmseStr = sprintf('RMSE = %.2f %s', rmseVal, statUnits);
    else
        rmseStr = sprintf('RMSE = NaN %s', statUnits);
    end

    if isfinite(biasVal)
        biasStr = sprintf('Bias = %.2f %s', biasVal, statUnits);
    else
        biasStr = sprintf('Bias = NaN %s', statUnits);
    end

    txt = sprintf('%s\n%s\n%s', rStr, rmseStr, biasStr);

    text(ax, statsX, statsY, txt, ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', ...
        'FontWeight', 'normal', ...
        'FontSize', max(fsAxis-4,9), ...
        'BackgroundColor', 'w', ...
        'Margin', 2);

    xlabel(ax, xlab, 'FontSize', fsAxis);
    ylabel(ax, ylab, 'FontSize', fsAxis);
end

%% ========================================================================
% Helper: compute shared y-limits (also reused for scatter x/y limits)
%% ========================================================================
function lims = compute_common_limits(vecCell, zeroFloor)

    yAll = [];

    for i = 1:numel(vecCell)
        yi = vecCell{i};

        if isempty(yi)
            continue
        end

        yi = yi(:);
        yi = yi(isfinite(yi));

        if ~isempty(yi)
            yAll = [yAll; yi]; %#ok<AGROW>
        end
    end

    if isempty(yAll)
        lims = [0 1];
        return
    end

    if zeroFloor
        yAll = yAll(yAll >= 0);

        if isempty(yAll)
            lims = [0 1];
            return
        end

        ymax = max(yAll);
        lims = [0, max(0.1, 1.12*ymax)];

    else
        ymin = min(yAll);
        ymax = max(yAll);

        if ymax <= ymin
            pad = max(0.10, 0.05*max(abs(ymax),1));
        else
            pad = 0.12 * (ymax - ymin);
        end

        lims = [max(0, ymin - pad), ymax + pad];

        if lims(2) <= lims(1)
            lims = [max(0, ymin - 0.5), ymax + 0.5];
        end
    end
end

%% ========================================================================
% Helper: panel label A/B/C/...
%% ========================================================================
function add_panel_label_paper(ax, labelChar, fsPanel)
    labelStr = sprintf('(%s)', lower(char(labelChar)));

    text(ax, 0.015, 0.98, labelStr, ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', ...
        'FontWeight', 'bold', ...
        'FontSize', fsPanel);
end

%% ========================================================================
% Local function: read SWAN SPEC1D and skip NODATA timestamps
%% ========================================================================
function [f, t, EnDens, nSkipped] = read_swan_spc1d_skip_nodata(filename)

    if ~isfile(filename)
        error('File not found: %s', filename);
    end

    fid = fopen(filename,'r');
    if fid < 0
        error('Could not open: %s', filename);
    end
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

    C = {};
    while true
        ln = fgetl(fid);
        if ~ischar(ln)
            break
        end
        C{end+1,1} = ln; %#ok<AGROW>
    end

    trimU = @(s) upper(strtrim(s));

    iAF = [];
    for i = 1:numel(C)
        if startsWith(trimU(C{i}), 'AFREQ')
            iAF = i;
            break
        end
    end
    if isempty(iAF)
        error('AFREQ block not found.');
    end

    Nf = sscanf(C{iAF+1}, '%d');
    if isempty(Nf) || Nf < 2
        error('Could not read Nf after AFREQ.');
    end

    f = nan(1,Nf);
    for k = 1:Nf
        tmp = sscanf(C{iAF+1+k}, '%f');
        if isempty(tmp)
            error('Could not parse AFREQ row %d', k);
        end
        f(k) = tmp(1);
    end

    iQ = [];
    for i = 1:numel(C)
        if startsWith(trimU(C{i}), 'QUANT')
            iQ = i;
            break
        end
    end
    if isempty(iQ)
        error('QUANT block not found.');
    end

    Nq = sscanf(C{iQ+1}, '%d');
    if isempty(Nq) || Nq < 1
        error('Could not read Nq after QUANT.');
    end

    exc = nan(Nq,1);
    for q = 1:Nq
        excLine = C{iQ + 1 + (q-1)*3 + 3};
        tmp = sscanf(excLine, '%f');
        if ~isempty(tmp)
            exc(q) = tmp(1);
        end
    end

    t = NaT(0,1);
    try
        t.TimeZone = 'UTC';
    catch
    end

    EnDens   = nan(0,Nf);
    nSkipped = 0;

    tsRE = '^\s*(\d{8}\.\d{6})';
    i = iQ + 1 + Nq*3 + 1;

    while i <= numel(C)

        ln = strtrim(C{i});
        m = regexp(ln, tsRE, 'tokens', 'once');

        if isempty(m)
            i = i + 1;
            continue
        end

        ts = m{1};

        try
            tt = datetime(ts, 'InputFormat','yyyyMMdd.HHmmss', 'TimeZone','UTC');
        catch
            tt = datetime(ts, 'InputFormat','yyyyMMdd.HHmmss');
        end

        iLoc = i + 1;
        while iLoc <= numel(C) && isempty(strtrim(C{iLoc}))
            iLoc = iLoc + 1;
        end

        if iLoc > numel(C)
            break
        end

        nextLine = trimU(C{iLoc});

        if startsWith(nextLine, 'NODATA')
            nSkipped = nSkipped + 1;
            i = iLoc + 1;
            continue
        end

        if ~startsWith(nextLine, 'LOCATION')
            error('Expected LOCATION or NODATA after time %s, got: %s', ts, C{iLoc});
        end

        rowStart = iLoc + 1;
        e = nan(1,Nf);

        for k = 1:Nf
            thisRow = rowStart + (k-1);

            if thisRow > numel(C)
                error('Unexpected end of file while reading spectra at time %s', ts);
            end

            nums = sscanf(C{thisRow}, '%f');
            if isempty(nums)
                error('Could not read spectral row at time %s (row %d). Line: %s', ts, k, C{thisRow});
            end

            e(k) = nums(1);
        end

        if isfinite(exc(1))
            e(abs(e - exc(1)) < 1e-6) = NaN;
        end

        t(end+1,1) = tt; %#ok<AGROW>
        EnDens(end+1,:) = e; %#ok<AGROW>

        i = rowStart + Nf;
    end

    if isempty(t)
        error('No valid spectral time blocks parsed.');
    end
end
