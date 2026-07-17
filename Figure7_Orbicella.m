%% Figure 7: Orbicella faveolata spatial wave attenuation
% Generates the Year 0, Year 20, and DeltaHs figure.
%
% Author: Janluis Rivera Ramos

clear; clc; close all;

%% Paths

dataRoot = 'D:\Paper_Thesis\Revisions';
maskRoot = 'E:\VegetationPatch_SWAN\Mound_6m_Vegetation_Polygon';

cmapFolder = fullfile(getenv('USERPROFILE'), ...
    'OneDrive', 'Documents', 'MATLAB', 'Colormaps');

if isfolder(cmapFolder)
    addpath(genpath(cmapFolder));
end

outputDir = fullfile(dataRoot, 'SpatialComparison_Outputs_GRIDMAX');

if ~isfolder(outputDir)
    mkdir(outputDir);
end

%% Settings

densityRoots = { ...
    fullfile(dataRoot, '1'), ...
    fullfile(dataRoot, '0.666667'), ...
    fullfile(dataRoot, '0.333333'), ...
    fullfile(dataRoot, '0.111111')};

densityLabels = { ...
    '9 stems per 9 m^{2}', ...
    '6 stems per 9 m^{2}', ...
    '3 stems per 9 m^{2}', ...
    '1 stem per 9 m^{2}'};

densityShortNames = { ...
    '9_stems_per_9m2', ...
    '6_stems_per_9m2', ...
    '3_stems_per_9m2', ...
    '1_stem_per_9m2'};

survivorship = [100 67 33 11];
selectedDensityIndex = 4;

year0 = 0;
yearRestored = 20;

scenarioFolders = {'TW', 'LP', 'H'};
scenarioTitles = {'Trade wind', 'Long-period', 'Hurricane'};

speciesFolder = 'O.faveolata';
speciesLabel = '\itOrbicella faveolata';
figureNumber = 7;

timePickMode = 'peakMeanHs';
minimumHs = 0.05;
onlyShowReductions = true;

restorationMaskFile = fullfile(maskRoot, ...
    'Orbicella', 'orbicella_veg_year51.dat');

maskOutsideValue = NaN;
maskDeltaHsOutsidePolygon = true;

applyDisplayThreshold = false;
displayThreshold = 0.005;

hsLimitsTWLP = [0 3];
useHurricanePercentile = true;
hurricanePercentile = 98;
deltaHsLimits = [0 0.003];

arrowSpacing = 15;
arrowScale = 0.5;
arrowLineWidth = 1.0;

bathyLevels = [0.5 0.5];
bathyColor = [0.7 0.7 0.7];
bathyLineWidth = 0.7;

outlineColor = [1 0 1];
colormapName = 'fake_parula';

fontAxis = 10;
fontLabel = 11;
fontTitle = 12;
fontPanel = 12;
fontRow = 12;

savePNG = true;
savePDF = false;

%% Load model output

densityIndex = selectedDensityIndex;

if densityIndex < 1 || densityIndex > numel(densityRoots)
    error('selectedDensityIndex must be 1, 2, 3, or 4.');
end

selectedRoot = densityRoots{densityIndex};

if ~isfolder(selectedRoot)
    error('Density folder not found:\n%s', selectedRoot);
end

scenario = struct();

for s = 1:numel(scenarioFolders)
    scenario(s).year0 = fullfile(selectedRoot, scenarioFolders{s}, ...
        speciesFolder, num2str(year0));

    scenario(s).yearRestored = fullfile(selectedRoot, scenarioFolders{s}, ...
        speciesFolder, num2str(yearRestored));

    requireFolder(scenario(s).year0);
    requireFolder(scenario(s).yearRestored);
end

[X, Y] = loadGridXY(scenario(1).year0);
[nRows, nCols] = size(X);

utmProjection = projcrs(32619);
[latitude, longitude] = projinv(utmProjection, X, Y);

longitudeLimits = [min(longitude(:)) max(longitude(:))];
latitudeLimits = [min(latitude(:)) max(latitude(:))];

bathymetry = [];
bathymetryFile = fullfile(scenario(1).year0, 'depth.mat');

if isfile(bathymetryFile)
    bathymetry = loadGridMatrix( ...
        bathymetryFile, nRows, nCols, X);
end

restorationMask = false(nRows, nCols);

if isfile(restorationMaskFile)
    restorationMask = loadMask( ...
        restorationMaskFile, nRows, nCols, X);
end

selected = struct();

for s = 1:numel(scenarioFolders)
    [hsYear0, directionYear0, timeYear0] = ...
        loadHsDirection(scenario(s).year0);

    [hsRestored, directionRestored, timeRestored] = ...
        loadHsDirection(scenario(s).yearRestored);

    targetIndex = chooseTimeIndex(hsYear0, timePickMode);
    targetTime = timeYear0(targetIndex);

    selected(s).hsYear0 = ...
        getAtTime(hsYear0, timeYear0, targetTime);

    selected(s).directionYear0 = ...
        getAtTime(directionYear0, timeYear0, targetTime);

    selected(s).hsRestored = ...
        getAtTime(hsRestored, timeRestored, targetTime);

    selected(s).directionRestored = ...
        getAtTime(directionRestored, timeRestored, targetTime);

    deltaHs = selected(s).hsYear0 - selected(s).hsRestored;

    if onlyShowReductions
        deltaHs(deltaHs < 0) = 0;
    end

    deltaHs(~isfinite(deltaHs)) = NaN;
    deltaHs(selected(s).hsYear0 < minimumHs) = NaN;

    if applyDisplayThreshold
        deltaHs(deltaHs < displayThreshold) = 0;
    end

    if maskDeltaHsOutsidePolygon && any(restorationMask(:))
        deltaHs(~restorationMask) = maskOutsideValue;
    end

    selected(s).deltaHs = deltaHs;
end

%% Create figure

fig = figure( ...
    'Color', 'w', ...
    'Units', 'pixels', ...
    'Position', [50 50 1650 900]);

leftMargin = 0.105;
rightMargin = 0.105;
bottomMargin = 0.090;
topMargin = 0.085;
columnGap = 0.025;
rowGap = 0.040;

axisWidth = ...
    (1 - leftMargin - rightMargin - 2 * columnGap) / 3;

axisHeight = ...
    (1 - bottomMargin - topMargin - 2 * rowGap) / 3;

axesGrid = gobjects(3, 3);
panelLetters = {'a','d','g'; 'b','e','h'; 'c','f','i'};

for s = 1:numel(scenarioFolders)
    if s < 3
        hsLimits = hsLimitsTWLP;
    elseif useHurricanePercentile
        hsLimits = hurricaneColorLimits( ...
            selected(s).hsYear0, ...
            selected(s).hsRestored, ...
            hurricanePercentile);
    else
        hsLimits = [0 4];
    end

    for row = 1:3
        xPosition = leftMargin + ...
            (s - 1) * (axisWidth + columnGap);

        yPosition = bottomMargin + ...
            (3 - row) * (axisHeight + rowGap);

        ax = axes( ...
            'Parent', fig, ...
            'Position', ...
            [xPosition yPosition axisWidth axisHeight]);

        axesGrid(row, s) = ax;

        showLongitude = row == 3;
        showLatitude = s == 1;
        showColorbar = s == 3;

        switch row
            case 1
                plotWaveField( ...
                    ax, longitude, latitude, ...
                    selected(s).hsYear0, ...
                    selected(s).directionYear0, ...
                    restorationMask, bathymetry, ...
                    longitudeLimits, latitudeLimits, ...
                    hsLimits, arrowSpacing, arrowScale, ...
                    arrowLineWidth, bathyLevels, ...
                    bathyColor, bathyLineWidth, ...
                    false, showLatitude, ...
                    fontAxis, fontLabel, ...
                    outlineColor, colormapName);

                title(ax, scenarioTitles{s}, ...
                    'FontSize', fontTitle, ...
                    'FontWeight', 'bold');

                if showColorbar
                    addColorbar(ax, 'H_s (m)', fontAxis);
                end

            case 2
                plotWaveField( ...
                    ax, longitude, latitude, ...
                    selected(s).hsRestored, ...
                    selected(s).directionRestored, ...
                    restorationMask, bathymetry, ...
                    longitudeLimits, latitudeLimits, ...
                    hsLimits, arrowSpacing, arrowScale, ...
                    arrowLineWidth, bathyLevels, ...
                    bathyColor, bathyLineWidth, ...
                    false, showLatitude, ...
                    fontAxis, fontLabel, ...
                    outlineColor, colormapName);

                if showColorbar
                    addColorbar(ax, 'H_s (m)', fontAxis);
                end

            case 3
                plotDifferenceField( ...
                    ax, longitude, latitude, ...
                    selected(s).deltaHs, ...
                    restorationMask, bathymetry, ...
                    longitudeLimits, latitudeLimits, ...
                    deltaHsLimits, bathyLevels, ...
                    bathyColor, bathyLineWidth, ...
                    showLongitude, showLatitude, ...
                    fontAxis, fontLabel, outlineColor);

                if showColorbar
                    addColorbar(ax, '\DeltaH_s (m)', fontAxis);
                end
        end

        addPanelLabel(ax, panelLetters{row, s}, fontPanel);
    end
end

addRowLabel(fig, axesGrid(1,1), ...
    'Restoration Year 0', fontRow);

addRowLabel(fig, axesGrid(2,1), ...
    sprintf('Restoration Year %d', yearRestored), fontRow);

addRowLabel(fig, axesGrid(3,1), ...
    sprintf('\\DeltaH_s (Year %d - Year %d)', ...
    year0, yearRestored), fontRow);

sgtitle(fig, sprintf('%s spatial response: %s (%d%% survival)', ...
    speciesLabel, ...
    densityLabels{densityIndex}, ...
    survivorship(densityIndex)), ...
    'FontSize', 14, ...
    'FontWeight', 'bold', ...
    'Interpreter', 'tex');

drawnow;

outputBase = sprintf( ...
    'Figure%d_Orbicella_spatial_masked_%s_Year%d_GRIDMAX', ...
    figureNumber, ...
    densityShortNames{densityIndex}, ...
    yearRestored);

if savePNG
    exportgraphics( ...
        fig, ...
        fullfile(outputDir, [outputBase '.png']), ...
        'Resolution', 600);
end

if savePDF
    exportgraphics( ...
        fig, ...
        fullfile(outputDir, [outputBase '.pdf']), ...
        'ContentType', 'vector');
end

%% Local functions

function requireFolder(folderPath)
    if ~isfolder(folderPath)
        error('Folder not found:\n%s', folderPath);
    end
end

function [X, Y] = loadGridXY(runPath)
    xData = load(fullfile(runPath, 'xp.mat'));
    yData = load(fullfile(runPath, 'yp.mat'));

    xFields = fieldnames(xData);
    yFields = fieldnames(yData);

    X = xData.(xFields{1});
    Y = yData.(yFields{1});
end

function [hsCell, directionCell, commonTimes] = ...
        loadHsDirection(runPath)

    hsFile = fullfile(runPath, 'hsig.mat');
    directionFile = fullfile(runPath, 'wdir.mat');

    if ~isfile(hsFile)
        error('Missing hsig.mat:\n%s', hsFile);
    end

    if ~isfile(directionFile)
        error('Missing wdir.mat:\n%s', directionFile);
    end

    hsData = load(hsFile);
    directionData = load(directionFile);

    hsNames = fieldnames(hsData);
    directionNames = fieldnames(directionData);

    hsValues = struct2cell(hsData);
    directionValues = struct2cell(directionData);

    validHs = cellfun( ...
        @(x) isnumeric(x) && ~isscalar(x), hsValues);

    validDirection = cellfun( ...
        @(x) isnumeric(x) && ~isscalar(x), directionValues);

    hsNames = hsNames(validHs);
    hsValues = hsValues(validHs);

    directionNames = directionNames(validDirection);
    directionValues = directionValues(validDirection);

    hsTimes = parseTimesFromNames(hsNames);
    directionTimes = parseTimesFromNames(directionNames);

    [hsTimes, hsOrder] = sort(hsTimes);
    [directionTimes, directionOrder] = sort(directionTimes);

    hsValues = hsValues(hsOrder);
    directionValues = directionValues(directionOrder);

    [commonTimes, hsIndex, directionIndex] = ...
        intersect(hsTimes, directionTimes);

    if isempty(commonTimes)
        error(['No common timestamps between hsig.mat ' ...
            'and wdir.mat in:\n%s'], runPath);
    end

    hsCell = hsValues(hsIndex);
    directionCell = directionValues(directionIndex);
end

function times = parseTimesFromNames(names)
    count = numel(names);
    timeStrings = cell(count, 1);
    foundTimestamp = false;

    for i = 1:count
        match = regexp( ...
            names{i}, '\d{8}_\d{6}', ...
            'match', 'once');

        if ~isempty(match)
            timeStrings{i} = match;
            foundTimestamp = true;
        end
    end

    if foundTimestamp
        for i = 1:count
            if isempty(timeStrings{i})
                timeStrings{i} = ...
                    sprintf('19000101_%06d', i);
            end
        end

        times = datetime( ...
            timeStrings, ...
            'InputFormat', 'yyyyMMdd_HHmmss');
    else
        times = (1:count)';
    end
end

function targetIndex = chooseTimeIndex(hsCell, mode)
    switch lower(mode)
        case 'first'
            targetIndex = 1;

        case 'peakmeanhs'
            meanHs = nan(numel(hsCell), 1);

            for i = 1:numel(hsCell)
                values = hsCell{i};

                if ndims(values) > 2
                    values = values(:,:,1);
                end

                meanHs(i) = mean(values(:), 'omitnan');
            end

            [~, targetIndex] = max(meanHs);

        otherwise
            error('Unknown timePickMode: %s', mode);
    end
end

function values = getAtTime(valueCell, times, targetTime)
    [found, index] = ismember(targetTime, times);

    if ~found
        if isdatetime(times)
            [~, index] = min(abs(times - targetTime));
        else
            [~, index] = min( ...
                abs(double(times) - double(targetTime)));
        end
    end

    values = squeeze(valueCell{index});

    if ndims(values) > 2
        values = values(:,:,1);
    end
end

function mask = loadMask(filePath, nRows, nCols, referenceGrid)
    fileId = fopen(filePath, 'r');

    if fileId == -1
        error('Cannot open vegetation file:\n%s', filePath);
    end

    cleanup = onCleanup(@() fclose(fileId)); %#ok<NASGU>
    rawValues = fscanf(fileId, '%f');

    expectedCount = nRows * nCols;

    if numel(rawValues) < expectedCount
        error('Mask file has fewer values than expected:\n%s', ...
            filePath);
    end

    rawValues = rawValues(1:expectedCount);

    try
        mask = reshape(rawValues, [nCols, nRows])';
    catch
        mask = reshape(rawValues, [nRows, nCols]);
    end

    mask = logical(mask);

    if ~isequal(size(mask), size(referenceGrid))
        mask = rot90(mask, -1);
    end
end

function values = loadGridMatrix( ...
        filePath, nRows, nCols, referenceGrid)

    data = load(filePath);
    fields = fieldnames(data);
    values = squeeze(data.(fields{1}));

    if isvector(values)
        try
            values = reshape(values, [nCols, nRows])';
        catch
            values = reshape(values, [nRows, nCols]);
        end
    end

    if ~isequal(size(values), size(referenceGrid))
        values = rot90(values, -1);
    end

    values(~isfinite(values)) = NaN;
end

function plotWaveField( ...
        ax, longitude, latitude, hs, direction, ...
        mask, bathymetry, longitudeLimits, latitudeLimits, ...
        colorLimits, arrowSpacing, arrowScale, ...
        arrowLineWidth, bathyLevels, bathyColor, ...
        bathyLineWidth, showLongitude, showLatitude, ...
        fontAxis, fontLabel, outlineColor, colormapName)

    pcolor(ax, longitude, latitude, hs);
    shading(ax, 'flat');
    axis(ax, 'equal');
    axis(ax, 'tight');

    applyColormap(ax, colormapName);
    clim(ax, colorLimits);
    hold(ax, 'on');

    if ~isempty(bathymetry)
        contour(ax, longitude, latitude, bathymetry, ...
            bathyLevels, ...
            'Color', bathyColor, ...
            'LineWidth', bathyLineWidth);
    end

    if any(mask(:))
        contour(ax, longitude, latitude, mask, [1 1], ...
            'Color', outlineColor, ...
            'LineWidth', 2.0);
    end

    directionRadians = deg2rad(270 - direction);
    u = cos(directionRadians);
    v = sin(directionRadians);

    quiver(ax, ...
        longitude(1:arrowSpacing:end, 1:arrowSpacing:end), ...
        latitude(1:arrowSpacing:end, 1:arrowSpacing:end), ...
        u(1:arrowSpacing:end, 1:arrowSpacing:end), ...
        v(1:arrowSpacing:end, 1:arrowSpacing:end), ...
        arrowScale, 'k', ...
        'LineWidth', arrowLineWidth);

    formatMapAxis( ...
        ax, longitudeLimits, latitudeLimits, ...
        showLongitude, showLatitude, ...
        fontAxis, fontLabel);
end

function plotDifferenceField( ...
        ax, longitude, latitude, differenceField, ...
        mask, bathymetry, longitudeLimits, latitudeLimits, ...
        colorLimits, bathyLevels, bathyColor, ...
        bathyLineWidth, showLongitude, showLatitude, ...
        fontAxis, fontLabel, outlineColor)

    pcolor(ax, longitude, latitude, differenceField);
    shading(ax, 'flat');
    axis(ax, 'equal');
    axis(ax, 'tight');

    colormap(ax, whiteRed(256));
    clim(ax, colorLimits);
    hold(ax, 'on');

    if ~isempty(bathymetry)
        contour(ax, longitude, latitude, bathymetry, ...
            bathyLevels, ...
            'Color', bathyColor, ...
            'LineWidth', bathyLineWidth);
    end

    if any(mask(:))
        contour(ax, longitude, latitude, mask, [1 1], ...
            'Color', outlineColor, ...
            'LineWidth', 2.0);
    end

    formatMapAxis( ...
        ax, longitudeLimits, latitudeLimits, ...
        showLongitude, showLatitude, ...
        fontAxis, fontLabel);
end

function formatMapAxis( ...
        ax, longitudeLimits, latitudeLimits, ...
        showLongitude, showLatitude, fontAxis, fontLabel)

    xlim(ax, longitudeLimits);
    ylim(ax, latitudeLimits);

    ax.FontSize = fontAxis;
    ax.Box = 'on';
    ax.Layer = 'top';

    if showLongitude
        xlabel(ax, 'Longitude (°)', ...
            'FontSize', fontLabel);
    else
        ax.XTickLabel = {};
    end

    if showLatitude
        ylabel(ax, 'Latitude (°)', ...
            'FontSize', fontLabel);
    else
        ax.YTickLabel = {};
    end
end

function applyColormap(ax, name)
    try
        colormapFunction = str2func(name);
        colormap(ax, colormapFunction(256));
    catch
        colormap(ax, parula(256));
    end
end

function colorLimits = hurricaneColorLimits( ...
        hsYear0, hsRestored, percentileValue)

    valuesYear0 = hsYear0(isfinite(hsYear0));
    valuesRestored = hsRestored(isfinite(hsRestored));

    if isempty(valuesYear0) || isempty(valuesRestored)
        colorLimits = [0 2.5];
        return;
    end

    upperLimit = max( ...
        prctile(valuesYear0, percentileValue), ...
        prctile(valuesRestored, percentileValue));

    if ~isfinite(upperLimit) || upperLimit <= 0
        upperLimit = 2.5;
    end

    colorLimits = [0 upperLimit];
end

function addPanelLabel(ax, letter, fontSize)
    text(ax, 0.02, 0.98, sprintf('(%s)', letter), ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', ...
        'FontSize', fontSize, ...
        'FontWeight', 'bold', ...
        'BackgroundColor', 'w', ...
        'Margin', 1);
end

function addColorbar(ax, labelText, fontSize)
    axisPosition = ax.Position;

    colorbarHandle = colorbar(ax);
    colorbarHandle.Label.String = labelText;
    colorbarHandle.FontSize = fontSize;

    ax.Position = axisPosition;

    colorbarHandle.Position = [ ...
        axisPosition(1) + axisPosition(3) + 0.006, ...
        axisPosition(2), ...
        0.012, ...
        axisPosition(4)];
end

function addRowLabel(figHandle, referenceAxis, labelText, fontSize)
    axisPosition = referenceAxis.Position;

    labelAxis = axes( ...
        'Parent', figHandle, ...
        'Position', ...
        [0.020 axisPosition(2) 0.045 axisPosition(4)], ...
        'Visible', 'off');

    text(labelAxis, 0.5, 0.5, labelText, ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'Rotation', 90, ...
        'FontSize', fontSize, ...
        'FontWeight', 'bold', ...
        'Interpreter', 'tex');
end

function colormapValues = whiteRed(count)
    if nargin < 1
        count = 256;
    end

    colormapValues = [ ...
        ones(count, 1), ...
        linspace(1, 0, count)', ...
        linspace(1, 0, count)'];
end
