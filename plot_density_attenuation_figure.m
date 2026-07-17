%% SWAN peak wave-height attenuation by density scenario
% Generates Comparison figure for Acropora palmata and
% Orbicella faveolata under trade-wind, long-period, and hurricane forcing.

clear; clc; close all;

%% Paths and scenarios

dataRoot = 'D:\Paper_Thesis\Revisions';
outputDir = fullfile(dataRoot, 'DensityComparison_Outputs_GRIDMAX');

if ~isfolder(outputDir)
    mkdir(outputDir);
end

densityFolders = {'1', '0.666667', '0.333333', '0.111111'};
densityLabels = {
    '9 stems per 9 m^{2}'
    '6 stems per 9 m^{2}'
    '3 stems per 9 m^{2}'
    '1 stem per 9 m^{2}'
};
survivalPercent = [100 67 33 11];
markers = {'o', 'd', 's', '^'};

for d = 1:numel(densityFolders)
    dens(d).root = fullfile(dataRoot, densityFolders{d});
    dens(d).label = densityLabels{d};
    dens(d).survival = survivalPercent(d);
    dens(d).marker = markers{d};
end

acroporaHex = {'#fed98e', '#fe9929', '#d95f0e', '#993404'};
orbicellaHex = {'#cbc9e2', '#9e9ac8', '#756bb1', '#54278f'};

for d = 1:numel(dens)
    dens(d).colorA = hex2rgb_local(acroporaHex{d});
    dens(d).colorO = hex2rgb_local(orbicellaHex{d});
end

scenDir = {'TW', 'LP', 'H'};
scenTitle = {'Trade-Wind', 'Long-Period', 'Hurricane'};

runs(1).mean = 'A.palmata_3.0';
runs(1).lower = 'A.palmata_3.0_Lower';
runs(1).upper = 'A.palmata_3.0_Upper';
runs(1).plotLabel = '\itAcropora palmata';

runs(2).mean = 'O.faveolata';
runs(2).lower = 'O.faveolata_lower';
runs(2).upper = 'O.faveolata_Upper';
runs(2).plotLabel = '\itOrbicella faveolata';

yearsPlot = 0:20;
timePickMode = 'peakMeanHs';
HsMin = 0.05;
clipNegativeAttenuation = true;

applyCoverCutoff = true;
firstYearCoverGE100 = [
     8, 10, 14, NaN
   NaN, NaN, NaN, NaN
];

plotUncertaintyBands = true;
bandAlpha = 0.14;

autoSpeciesYLim = false;
manualSpeciesYLim = {
    [0 14]
    [0 1.0]
};

saveFigurePNG = true;
saveFigurePDF = false;

figureNamePNG = fullfile(outputDir, ...
    'Density_Attenuation_GRIDMAX_2x3_9_6_3_1_100cover_cutoff.png');
figureNamePDF = fullfile(outputDir, ...
    'Density_Attenuation_GRIDMAX_2x3_9_6_3_1_100cover_cutoff.pdf');

lineWidth = 3.0;
markerSize = 6;
fontAxes = 13;
fontTitle = 15;
fontPanel = 15;
fontLegend = 12;
fontLegendTitle = 13;
panelLabels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};

nDens = numel(dens);
nScen = numel(scenDir);
nRuns = numel(runs);
nYears = numel(yearsPlot);

for d = 1:nDens
    if ~isfolder(dens(d).root)
        warning('Density root not found: %s', dens(d).root);
    end
end

%% Calculate attenuation

results = struct();

for r = 1:nRuns
    for s = 1:nScen
        for d = 1:nDens

            baselinePath = fullfile( ...
                dens(d).root, scenDir{s}, runs(r).mean, '0');

            [Hs0Cell, tt0] = loadHsOnly_local(baselinePath);
            targetIndex = chooseTimeIndex_local(Hs0Cell, timePickMode);
            targetTime = tt0(targetIndex);
            Hs0 = getAtTime_local(Hs0Cell, tt0, targetTime);

            attMean = nan(1, nYears);

            for k = 1:nYears
                yearPath = fullfile( ...
                    dens(d).root, scenDir{s}, runs(r).mean, ...
                    num2str(yearsPlot(k)));

                if ~isfile(fullfile(yearPath, 'hsig.mat'))
                    warning('Missing file: %s', fullfile(yearPath, 'hsig.mat'));
                    continue;
                end

                [HsYearCell, ttYear] = loadHsOnly_local(yearPath);
                HsYear = getAtTime_local(HsYearCell, ttYear, targetTime);
                attenuation = calcAttField_local( ...
                    Hs0, HsYear, HsMin, clipNegativeAttenuation);

                attMean(k) = max(attenuation(:), [], 'omitnan');
            end

            attLow = nan(1, nYears);
            attHigh = nan(1, nYears);

            lowerPath = fullfile( ...
                dens(d).root, scenDir{s}, runs(r).lower);
            upperPath = fullfile( ...
                dens(d).root, scenDir{s}, runs(r).upper);

            if plotUncertaintyBands && isfolder(lowerPath) && isfolder(upperPath)
                lowerSeries = calcGridMaxSeries_local( ...
                    dens(d).root, scenDir{s}, runs(r).lower, ...
                    yearsPlot, targetTime, Hs0, HsMin, ...
                    clipNegativeAttenuation);

                upperSeries = calcGridMaxSeries_local( ...
                    dens(d).root, scenDir{s}, runs(r).upper, ...
                    yearsPlot, targetTime, Hs0, HsMin, ...
                    clipNegativeAttenuation);

                attLow = min(lowerSeries, upperSeries);
                attHigh = max(lowerSeries, upperSeries);
            end

            results(r, s, d).attMean = attMean;
            results(r, s, d).attLow = attLow;
            results(r, s, d).attHigh = attHigh;
        end
    end
end

%% Create figure

fig = figure( ...
    'Color', 'w', ...
    'Name', 'Density Scenario Attenuation Comparison', ...
    'Units', 'inches', ...
    'Position', [0.7 0.7 15.0 8.5]);

layout = tiledlayout(fig, nRuns, nScen, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax = gobjects(nRuns, nScen);

for r = 1:nRuns
    for s = 1:nScen

        tileIndex = (r - 1) * nScen + s;
        ax(r, s) = nexttile(layout, tileIndex);
        hold(ax(r, s), 'on');
        styleAxes_local(ax(r, s), fontAxes);

        title(ax(r, s), ...
            sprintf('%s: %s', runs(r).plotLabel, scenTitle{s}), ...
            'FontSize', fontTitle, ...
            'FontWeight', 'bold', ...
            'Interpreter', 'tex');

        xlabel(ax(r, s), 'Restoration Year', 'FontSize', fontAxes);

        if s == 1
            ylabel(ax(r, s), ...
                'Peak wave-height attenuation (%)', ...
                'FontSize', fontAxes);
        end

        xlim(ax(r, s), [min(yearsPlot) max(yearsPlot)]);
        xticks(ax(r, s), 0:5:20);

        legendHandles = gobjects(nDens, 1);

        if plotUncertaintyBands
            for d = 1:nDens
                curveColor = pickDensityColor_local(r, dens(d));

                [xBand, yLow] = applyCoverCutoff_local( ...
                    yearsPlot, results(r, s, d).attLow, ...
                    applyCoverCutoff, firstYearCoverGE100(r, d));

                [~, yHigh] = applyCoverCutoff_local( ...
                    yearsPlot, results(r, s, d).attHigh, ...
                    applyCoverCutoff, firstYearCoverGE100(r, d));

                if any(isfinite(yLow)) && any(isfinite(yHigh))
                    fillBand_local( ...
                        ax(r, s), xBand, yLow, yHigh, ...
                        curveColor, bandAlpha);
                end
            end
        end

        for d = 1:nDens
            curveColor = pickDensityColor_local(r, dens(d));

            [xPlot, yPlot] = applyCoverCutoff_local( ...
                yearsPlot, results(r, s, d).attMean, ...
                applyCoverCutoff, firstYearCoverGE100(r, d));

            legendHandles(d) = plot(ax(r, s), xPlot, yPlot, '-', ...
                'Color', curveColor, ...
                'LineWidth', lineWidth, ...
                'Marker', dens(d).marker, ...
                'MarkerSize', markerSize, ...
                'MarkerFaceColor', 'w', ...
                'DisplayName', sprintf( ...
                    '%d%% survival (%s)', ...
                    dens(d).survival, dens(d).label));
        end

        if s == 1
            lgd = legend(ax(r, s), legendHandles, ...
                'Location', 'northwest', ...
                'Interpreter', 'tex');

            set(lgd, ...
                'Box', 'off', ...
                'FontSize', fontLegend, ...
                'ItemTokenSize', [18 10]);

            try
                lgdTitle = title( ...
                    lgd, runs(r).plotLabel, 'Interpreter', 'tex');
                set(lgdTitle, ...
                    'FontSize', fontLegendTitle, ...
                    'FontWeight', 'bold');
            catch
            end
        end

        addPanelLabel_local( ...
            ax(r, s), panelLabels{tileIndex}, fontPanel);
    end
end

if autoSpeciesYLim
    for r = 1:nRuns
        rowMaximum = 0;

        for s = 1:nScen
            for d = 1:nDens
                rowMaximum = max( ...
                    rowMaximum, ...
                    max(results(r, s, d).attMean, [], 'omitnan'));

                if plotUncertaintyBands
                    rowMaximum = max( ...
                        rowMaximum, ...
                        max(results(r, s, d).attHigh, [], 'omitnan'));
                end
            end
        end

        rowLimit = niceYTop_local(rowMaximum);

        for s = 1:nScen
            ylim(ax(r, s), [0 rowLimit]);
        end
    end
else
    for r = 1:nRuns
        for s = 1:nScen
            ylim(ax(r, s), manualSpeciesYLim{r});
        end
    end
end

if saveFigurePNG
    exportgraphics(fig, figureNamePNG, 'Resolution', 600);
end

if saveFigurePDF
    exportgraphics(fig, figureNamePDF, 'ContentType', 'vector');
end

%% Local functions

function attSeries = calcGridMaxSeries_local( ...
    baseRoot, scenarioName, runName, yearsPlot, targetTime, ...
    Hs0, HsMin, clipNegative)

    attSeries = nan(size(yearsPlot));

    for k = 1:numel(yearsPlot)
        yearPath = fullfile( ...
            baseRoot, scenarioName, runName, num2str(yearsPlot(k)));

        if ~isfile(fullfile(yearPath, 'hsig.mat'))
            warning('Missing file: %s', fullfile(yearPath, 'hsig.mat'));
            continue;
        end

        [HsYearCell, ttYear] = loadHsOnly_local(yearPath);
        HsYear = getAtTime_local(HsYearCell, ttYear, targetTime);
        attenuation = calcAttField_local( ...
            Hs0, HsYear, HsMin, clipNegative);

        attSeries(k) = max(attenuation(:), [], 'omitnan');
    end
end

function attenuation = calcAttField_local( ...
    Hs0, HsYear, HsMin, clipNegative)

    attenuation = 100 .* (Hs0 - HsYear) ./ Hs0;

    if clipNegative
        attenuation(attenuation < 0) = 0;
    end

    attenuation(Hs0 < HsMin) = NaN;
end

function targetIndex = chooseTimeIndex_local(HsCell, timePickMode)

    switch lower(timePickMode)
        case 'first'
            targetIndex = 1;

        case 'peakmeanhs'
            meanHs = nan(numel(HsCell), 1);

            for i = 1:numel(HsCell)
                field = HsCell{i};
                meanHs(i) = mean(field(:), 'omitnan');
            end

            [~, targetIndex] = max(meanHs);

        otherwise
            error('Unknown timePickMode: %s', timePickMode);
    end
end

function [HsCell, timeValues] = loadHsOnly_local(runPath)

    matFile = fullfile(runPath, 'hsig.mat');

    if ~isfile(matFile)
        error('Missing hsig.mat at: %s', matFile);
    end

    data = load(matFile);
    fieldNames = fieldnames(data);
    fieldValues = struct2cell(data);

    keep = cellfun( ...
        @(value) isnumeric(value) && ~isscalar(value), ...
        fieldValues);

    fieldNames = fieldNames(keep);
    fieldValues = fieldValues(keep);

    if isempty(fieldValues)
        error('No numeric wave-height arrays found in %s', matFile);
    end

    timeValues = parseTimesFromNames_local(fieldNames);
    [timeValues, order] = sort(timeValues);
    HsCell = fieldValues(order);
end

function Hs = getAtTime_local(HsCell, timeValues, targetTime)

    [found, index] = ismember(targetTime, timeValues);

    if ~found
        [~, index] = min(abs(timeValues - targetTime));
    end

    Hs = HsCell{index};

    if ndims(Hs) > 2
        Hs = Hs(:, :, 1);
    end
end

function timeValues = parseTimesFromNames_local(fieldNames)

    nFields = numel(fieldNames);
    timeStrings = cell(nFields, 1);
    foundTimestamp = false;

    for i = 1:nFields
        match = regexp( ...
            fieldNames{i}, '\d{8}_\d{6}', 'match', 'once');

        if ~isempty(match)
            timeStrings{i} = match;
            foundTimestamp = true;
        end
    end

    if foundTimestamp
        for i = 1:nFields
            if isempty(timeStrings{i})
                timeStrings{i} = sprintf('19000101_%06d', i);
            end
        end

        timeValues = datetime( ...
            timeStrings, 'InputFormat', 'yyyyMMdd_HHmmss');
    else
        timeValues = (1:nFields)';
    end
end

function fillBand_local(ax, x, yLow, yHigh, colorValue, alphaValue)

    x = x(:)';
    yLow = yLow(:)';
    yHigh = yHigh(:)';

    valid = ~(isnan(x) | isnan(yLow) | isnan(yHigh));

    if ~any(valid)
        return;
    end

    validIndex = find(valid);
    segmentBreaks = [ ...
        1, find(diff(validIndex) > 1) + 1, numel(validIndex) + 1];

    for b = 1:numel(segmentBreaks) - 1
        segment = validIndex( ...
            segmentBreaks(b):segmentBreaks(b + 1) - 1);

        xSegment = x(segment);
        lowSegment = yLow(segment);
        highSegment = yHigh(segment);

        patchHandle = fill(ax, ...
            [xSegment fliplr(xSegment)], ...
            [lowSegment fliplr(highSegment)], ...
            colorValue, ...
            'FaceAlpha', alphaValue, ...
            'EdgeColor', 'none', ...
            'HandleVisibility', 'off');

        try
            patchHandle.Annotation.LegendInformation.IconDisplayStyle = 'off';
        catch
        end
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

function addPanelLabel_local(ax, labelText, fontSize)

    text(ax, 0.00, 1.04, labelText, ...
        'Units', 'normalized', ...
        'FontSize', fontSize, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'bottom', ...
        'Clipping', 'off');
end

function yTop = niceYTop_local(yMaximum)

    if ~isfinite(yMaximum) || yMaximum <= 0
        yTop = 1;
        return;
    end

    yMaximum = 1.10 * yMaximum;

    if yMaximum <= 1
        step = 0.2;
    elseif yMaximum <= 3
        step = 0.5;
    elseif yMaximum <= 10
        step = 1;
    elseif yMaximum <= 30
        step = 5;
    else
        step = 10;
    end

    yTop = ceil(yMaximum / step) * step;
end

function rgb = hex2rgb_local(hexValue)

    hexValue = strrep(char(hexValue), '#', '');

    if numel(hexValue) ~= 6
        error('Hex color must contain six characters.');
    end

    rgb = [
        hex2dec(hexValue(1:2))
        hex2dec(hexValue(3:4))
        hex2dec(hexValue(5:6))
    ]' ./ 255;
end

function colorValue = pickDensityColor_local(speciesIndex, density)

    if speciesIndex == 1
        colorValue = density.colorA;
    else
        colorValue = density.colorO;
    end
end

function [xOut, yOut] = applyCoverCutoff_local( ...
    xIn, yIn, applyCutoff, firstYearCoverGE100)

    xOut = xIn(:)';
    yOut = yIn(:)';

    if applyCutoff && isfinite(firstYearCoverGE100)
        keep = xOut < firstYearCoverGE100;
        xOut = xOut(keep);
        yOut = yOut(keep);
    end
end
