%% Acropora palmata sensitivity analysis
% Generates a three-panel figure for ±25% perturbations in bottom friction, coral drag, and incident wave direction.

clear; clc; close all;

%% Settings

rootDir = 'D:\Paper_Thesis\Sensitivity_Acropora_TW';
outputDir = fullfile(rootDir, 'Sensitivity_Figure_Output');

yearsPlot = [0 1 5 10 15 20];
timePickMode = 'peakMeanHs';

HsMin = 0.05;
clipNegativeAttenuation = true;

saveFigurePNG = true;
saveFigurePDF = false;

figureNamePNG = fullfile(outputDir, ...
    'Acropora_Sensitivity_Cf_CD_Direction_GRIDMAX.png');
figureNamePDF = fullfile(outputDir, ...
    'Acropora_Sensitivity_Cf_CD_Direction_GRIDMAX.pdf');

if ~isfolder(outputDir)
    mkdir(outputDir);
end

panels(1).title = 'Sensitivity to bottom friction C_f';
panels(1).cases(1) = makeLineCase('C_f -25%', 'Cf_minus25', '--', 'o');
panels(1).cases(2) = makeLineCase('Baseline', 'baseline', '-', 'o');
panels(1).cases(3) = makeLineCase('C_f +25%', 'Cf_plus25', '--', 'o');

panels(2).title = 'Sensitivity to coral drag C_D';
panels(2).cases(1) = makeLineCase('C_D -25%', 'CD_minus25', '--', 's');
panels(2).cases(2) = makeLineCase('Baseline', 'baseline', '-', 's');
panels(2).cases(3) = makeLineCase('C_D +25%', 'CD_plus25', '--', 's');

panels(3).title = 'Sensitivity to wave direction';
panels(3).cases(1) = makeLineCase( ...
    'Direction -25%', 'DIR_minus25percent', '--', '^');
panels(3).cases(2) = makeLineCase( ...
    'Baseline', 'baseline', '-', '^');
panels(3).cases(3) = makeLineCase( ...
    'Direction +25%', 'DIR_plus25percent', '--', '^');

colors = [
    hex2rgb_local('#91bfdb')
    hex2rgb_local('#ffffbf')
    hex2rgb_local('#fc8d59')
];

for p = 1:numel(panels)
    for c = 1:numel(panels(p).cases)
        panels(p).cases(c).color = colors(c,:);
    end
end

lineWidth = 2.5;
markerSize = 6;
fontAxes = 11;
fontTitle = 12;
fontPanel = 13;
panelLabels = {'(a)', '(b)', '(c)'};

commonYLim = [0 40];

%% Calculate attenuation

nPanels = numel(panels);
nYears = numel(yearsPlot);
allY = [];

for p = 1:nPanels
    for c = 1:numel(panels(p).cases)
        caseDir = fullfile(rootDir, panels(p).cases(c).folder);

        if ~isfolder(caseDir)
            error('Sensitivity case folder not found: %s', caseDir);
        end

        [Hs0Cell, t0] = loadHsOnly_local(fullfile(caseDir, '0'));
        targetIndex = chooseTimeIndex_local(Hs0Cell, timePickMode);
        targetTime = t0(targetIndex);
        Hs0 = getAtTime_local(Hs0Cell, t0, targetTime);

        attenuation = nan(1, nYears);

        for k = 1:nYears
            yearDir = fullfile(caseDir, num2str(yearsPlot(k)));
            matFile = fullfile(yearDir, 'hsig.mat');

            if ~isfile(matFile)
                warning('Missing file: %s', matFile);
                continue;
            end

            [HsYearCell, tYear] = loadHsOnly_local(yearDir);
            HsYear = getAtTime_local(HsYearCell, tYear, targetTime);

            attenuationField = calcAttField_local( ...
                Hs0, HsYear, HsMin, clipNegativeAttenuation);

            attenuation(k) = max(attenuationField(:), [], 'omitnan');
        end

        panels(p).cases(c).attenuation = attenuation;
        allY = [allY attenuation]; %#ok<AGROW>
    end
end

%% Create figure

fig = figure( ...
    'Color', 'w', ...
    'Name', 'Acropora Sensitivity Analysis', ...
    'Units', 'inches', ...
    'Position', [0.7 1.0 13.0 4.7]);

layout = tiledlayout(fig, 1, 3, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

axesHandles = gobjects(1, nPanels);

for p = 1:nPanels
    ax = nexttile(layout, p);
    axesHandles(p) = ax;

    hold(ax, 'on');
    styleAxes_local(ax, fontAxes);

    title(ax, panels(p).title, ...
        'FontSize', fontTitle, ...
        'FontWeight', 'bold', ...
        'Interpreter', 'tex');

    xlabel(ax, 'Restoration Year', 'FontSize', fontAxes);

    if p == 1
        ylabel(ax, 'Peak wave-height attenuation (%)', ...
            'FontSize', fontAxes);
    end

    xlim(ax, [min(yearsPlot) max(yearsPlot)]);
    xticks(ax, yearsPlot);

    legendHandles = gobjects(1, numel(panels(p).cases));

    for c = 1:numel(panels(p).cases)
        caseData = panels(p).cases(c);

        legendHandles(c) = plot(ax, yearsPlot, caseData.attenuation, ...
            caseData.lineStyle, ...
            'Color', caseData.color, ...
            'LineWidth', lineWidth, ...
            'Marker', caseData.marker, ...
            'MarkerSize', markerSize, ...
            'MarkerFaceColor', caseData.color, ...
            'MarkerEdgeColor', caseData.color, ...
            'DisplayName', caseData.label);
    end

    legend(ax, legendHandles, ...
        'Location', 'northwest', ...
        'Box', 'off', ...
        'Interpreter', 'tex', ...
        'FontSize', fontAxes - 1);

    addPanelLabel_local(ax, panelLabels{p}, fontPanel);
end

if isempty(commonYLim)
    yTop = niceYTop_local(max(allY, [], 'omitnan'));
    commonYLim = [0 yTop];
end

for p = 1:nPanels
    ylim(axesHandles(p), commonYLim);
end

%% Export figure

if saveFigurePNG
    exportgraphics(fig, figureNamePNG, 'Resolution', 600);
end

if saveFigurePDF
    exportgraphics(fig, figureNamePDF, 'ContentType', 'vector');
end

%% Local functions

function caseInfo = makeLineCase(label, folder, lineStyle, marker)
    caseInfo.label = label;
    caseInfo.folder = folder;
    caseInfo.lineStyle = lineStyle;
    caseInfo.marker = marker;
    caseInfo.color = [0 0 0];
    caseInfo.attenuation = [];
end

function attenuation = calcAttField_local(Hs0, HsYear, HsMin, clipNegative)
    attenuation = 100 .* (Hs0 - HsYear) ./ Hs0;

    if clipNegative
        attenuation(attenuation < 0) = 0;
    end

    attenuation(Hs0 < HsMin) = NaN;
    attenuation(~isfinite(attenuation)) = NaN;
end

function targetIndex = chooseTimeIndex_local(HsCell, timePickMode)
    switch lower(timePickMode)
        case 'first'
            targetIndex = 1;

        case 'peakmeanhs'
            meanHs = nan(numel(HsCell), 1);

            for i = 1:numel(HsCell)
                values = HsCell{i};
                meanHs(i) = mean(values(:), 'omitnan');
            end

            [~, targetIndex] = max(meanHs);

        otherwise
            error('Unknown timePickMode: %s', timePickMode);
    end
end

function [HsCell, times] = loadHsOnly_local(runPath)
    matFile = fullfile(runPath, 'hsig.mat');

    if ~isfile(matFile)
        error('Missing hsig.mat: %s', matFile);
    end

    data = load(matFile);
    names = fieldnames(data);
    values = struct2cell(data);

    keep = cellfun(@(x) isnumeric(x) && ~isscalar(x), values);
    names = names(keep);
    values = values(keep);

    if isempty(values)
        error('No numeric wave-height arrays found in %s', matFile);
    end

    times = parseTimesFromNames_local(names);
    [times, order] = sort(times);
    HsCell = values(order);
end

function Hs = getAtTime_local(HsCell, times, targetTime)
    [found, index] = ismember(targetTime, times);

    if ~found
        [~, index] = min(abs(times - targetTime));
    end

    Hs = squeeze(HsCell{index});

    if ndims(Hs) > 2
        Hs = Hs(:,:,1);
    end
end

function times = parseTimesFromNames_local(names)
    nNames = numel(names);
    timeStrings = cell(nNames, 1);
    hasTimestamps = false;

    for i = 1:nNames
        match = regexp(names{i}, '\d{8}_\d{6}', 'match', 'once');

        if ~isempty(match)
            timeStrings{i} = match;
            hasTimestamps = true;
        end
    end

    if hasTimestamps
        for i = 1:nNames
            if isempty(timeStrings{i})
                timeStrings{i} = sprintf('19000101_%06d', i);
            end
        end

        times = datetime(timeStrings, ...
            'InputFormat', 'yyyyMMdd_HHmmss');
    else
        times = (1:nNames)';
    end
end

function styleAxes_local(ax, fontSize)
    set(ax, ...
        'FontSize', fontSize, ...
        'LineWidth', 1.0, ...
        'Layer', 'top', ...
        'Box', 'off', ...
        'Color', 'w', ...
        'XColor', 'k', ...
        'YColor', 'k', ...
        'TickDir', 'out');

    grid(ax, 'off');
end

function addPanelLabel_local(ax, label, fontSize)
    text(ax, 0.00, 1.04, label, ...
        'Units', 'normalized', ...
        'FontSize', fontSize, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'bottom', ...
        'Clipping', 'off');
end

function yTop = niceYTop_local(yMax)
    if ~isfinite(yMax) || yMax <= 0
        yTop = 1;
        return;
    end

    yMax = 1.10 * yMax;

    if yMax <= 1
        step = 0.2;
    elseif yMax <= 3
        step = 0.5;
    elseif yMax <= 10
        step = 1;
    elseif yMax <= 30
        step = 5;
    else
        step = 10;
    end

    yTop = ceil(yMax / step) * step;
end

function rgb = hex2rgb_local(hexColor)
    hexColor = char(hexColor);
    hexColor = strrep(hexColor, '#', '');

    if numel(hexColor) ~= 6
        error('Hex color must contain six characters.');
    end

    rgb = [
        hex2dec(hexColor(1:2))
        hex2dec(hexColor(3:4))
        hex2dec(hexColor(5:6))
    ]' ./ 255;
end
