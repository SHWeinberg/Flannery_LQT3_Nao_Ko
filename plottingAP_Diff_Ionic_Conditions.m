
%% WT-normal, MT- normal, MT-HypoNa, MT-HyperNa, MT-HypoK, MT-HyperK

cellIndices = [1, 55, 2656, 888, 4088, 14]; % examples of cells in the population that are susceptible under each example condition
% 

titleLabels = {... 
    'Susceptible to None', ...
    'Only Susceptible to HypoNa', ...
    'Only Susceptible to HyperNa', ...
    'Only Susceptible to HypoK', ...
    'Only Susceptible to HyperK', ...
    'Susceptible to All'};

mutant_type = 'dKPQ';
cellType = 'Torord_endo';
beats = 101; ignoreFirst = 99;

ion_conditions = [120, 5; 160, 5; 140, 3; 140, 7; 140, 5];
labels_mutant = {'LQT3 - HypoNa (120 mM)', 'LQT3 - HyperNa (160 mM)', ...
                 'LQT3 - HypoK (3 mM)', 'LQT3 - HyperK (7 mM)', 'LQT3 - Normal (140/5 mM)'};
colors = lines(4);      
colors(5,:) = [0 0 0];  % Black for normal

multiplierFieldNames = {'INa', 'Ito', 'ICaL', 'IKr', 'IKs', 'IK1', ...
                    'INaCa_i', 'INaCa_ss', 'INaK', 'IKb', ...
                    'INab', 'ICab', 'IpCa'};

multiplierFieldLabels = {'Na', 'to', 'CaL', 'Kr', 'Ks', 'K1', ...
                    'NaCa_i', 'NaCa_{ss}', 'NaK', 'Kb', ...
                    'Nab', 'Cab', 'pCa'};


for subplotIdx = 1:length(cellIndices)
    cellIdx = cellIndices(subplotIdx);
    scalingVec = passScalingFactors(cellIdx,:);

    Vm_all_mutant = cell(1, size(ion_conditions,1));
    time_all_mutant = cell(1, size(ion_conditions,1));

    for i = 1:size(ion_conditions, 1)
        nao = ion_conditions(i,1);
        ko  = ion_conditions(i,2);

        param = struct();
        param.bcl = 1000;
        param.model = @model_Torord_Clancy;
        param.verbose = false;
        param.nao = nao;
        param.ko = ko;
        param.cao = 1.8;
        param.mt1 = mutant_type;
        param.mt2 = mutant_type;
        param.mf1 = 1;
        param.mf2 = 0;

        for k = 1:length(multiplierFieldNames)
            param.([multiplierFieldNames{k} '_Multiplier']) = scalingVec(k);
        end

        X0 = getStartingState(cellType);
        INa_SS1 = determine_INa_MC_steadystate(X0(1), mutant_type, 1);
        INa_SS2 = determine_INa_MC_steadystate(X0(1), mutant_type, 0);
        X0 = [X0, INa_SS1, INa_SS1, INa_SS2, INa_SS2];

        [time_raw, X] = modelRunner_Clancy(X0, [], param, beats, ignoreFirst);
        currents = getCurrentsStructure_Clancy(time_raw, X, param, ignoreFirst);
        Vm_all_mutant{i} = currents.V;
        time_all_mutant{i} = currents.time - currents.time(1);
    end

    % WT normal
    param = struct();
    param.bcl = 1000;
    param.model = @model_Torord_Clancy;
    param.verbose = false;
    param.nao = 140; param.ko = 5; param.cao = 1.8;
    param.mt1 = mutant_type;
    param.mt2 = mutant_type;
    param.mf1 = 0;
    param.mf2 = 0;

    for k = 1:length(multiplierFieldNames)
        param.([multiplierFieldNames{k} '_Multiplier']) = scalingVec(k);
    end

    X0 = getStartingState(cellType);
    INa_SS1 = determine_INa_MC_steadystate(X0(1), mutant_type, 0);
    INa_SS2 = determine_INa_MC_steadystate(X0(1), mutant_type, 0);
    X0 = [X0, INa_SS1, INa_SS1, INa_SS2, INa_SS2];

    [timeWT_raw, XWT] = modelRunner_Clancy(X0, [], param, beats, ignoreFirst);
    currentsWT = getCurrentsStructure_Clancy(timeWT_raw, XWT, param, ignoreFirst);
    VmWT = currentsWT.V;
    timeWT = currentsWT.time - currentsWT.time(1);

    % === PLOTTING ===
    figWidth = 500;  % narrower width
    figHeight = 600;
    fig = figure('Name', sprintf('Cell %d - %s', cellIdx, titleLabels{subplotIdx}), ...
            'Position', [400 200 figWidth figHeight]);

    % Top (Vm) plot – square
    ax1 = axes(fig);
    hold(ax1, 'on');
        legendEntries = cell(1, length(Vm_all_mutant) + 1);
       legendEntries = cell(1, length(Vm_all_mutant) + 1);
    for i = 1:length(Vm_all_mutant)
        plot(ax1, time_all_mutant{i}, Vm_all_mutant{i}, 'LineWidth', 1.3, 'Color', colors(i,:));
        legendEntries{i} = labels_mutant{i};
    end
    plot(ax1, timeWT, VmWT, 'k:', 'LineWidth', 1.8);
    legendEntries{end} = 'WT - Normal (140/5 mM)';

    if subplotIdx == 1
        lgd = legend(ax1, legendEntries, 'Location', 'northeast', 'FontSize', 13);
        lgd.Box = 'off';
    end

    % === Remove x-axis and add thick black bar aligned to bottom of y-axis ===
    ax1.XTick = [];
    ax1.XColor = 'none';
    barY = -110;  % bottom of y-axis
    blackLine = line(ax1, [1000 1500], [barY barY], 'Color', 'k', 'LineWidth', 4);
    blackLine.Annotation.LegendInformation.IconDisplayStyle = 'off';  % exclude from legend
    
    % Add label above bar for the last plot only
    if subplotIdx == 6
        text(ax1, 250, barY + 5, '500 ms', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', 18, 'FontWeight', 'normal');
    end


    % === Remove y-axis for all but subplotIdx 1 and 4 ===
    if subplotIdx ~= 1 && subplotIdx ~= 4
        ax1.YTick = [];
        ax1.YColor = 'none';
        ylabel(ax1, '');
    end


    title(ax1, titleLabels{subplotIdx}, 'FontSize', 16);
    ylabel(ax1, 'V_m (mV)', 'FontSize', 16);
    xlim(ax1, [950,2000]);
    ylim(ax1, [-110, 100]);
    % grid(ax1, 'on');
    axis(ax1, 'square');

    
    ax1.Units = 'normalized';
    ax1.Position = [0.15, 0.42, 0.7, 0.5];  % lowered slightly from 0.48 → 0.42

    
    % Bottom (bar) plot – shorter and same width
    ax2 = axes(fig);
    logVals = log2(scalingVec);
    bar(ax2, logVals, 'FaceColor', [0.2 0.4 1]);
    yline(ax2, 0, '--k');
    latexLabels = arrayfun(@(s) ['$' s{1} '$'], multiplierFieldLabels, 'UniformOutput', false);
    set(ax2, 'XTick', 1:13, ...
         'XTickLabel', multiplierFieldLabels, ...
         'XTickLabelRotation', 45, ...
         'TickLabelInterpreter', 'tex', 'FontSize', 16);  % default font, subscripts with {}

    %ylabel(ax2, 'Conductance');
    ylim(ax2, [-2.2 2.2]);
    ax2.YTick = [-2 -1 0 1 2];
    ax2.YTickLabel = {'0.25', '0.5', '1', '2', '4'};
    % grid(ax2, 'on');
    box(ax2, 'on');

    % Position bar plot closer to top plot
    ax2.Units = 'normalized';
    ax2.Position = [ax1.Position(1), 0.24, ax1.Position(3), 0.14];  % raised from 0.08 → 0.24

    ax1.YLabel.FontSize = 16;
    ax1.YAxis.FontSize = 16;
    ax2.YLabel.FontSize = 16;
    ax2.YAxis.FontSize = 16;


end

%% WT- Normal, WT - HypoNa, WT - HyperNa, WT - HypoK, WT - HyperK

cellIndices = [1, 55, 2656, 888, 4088, 14]; % examples of cells in the population that are susceptible under each example condition
titleLabels = {
    'Susceptible to None', ...
    'Only Susceptible to HypoNa', ...
    'Only Susceptible to HyperNa', ...
    'Only Susceptible to HypoK', ...
    'Only Susceptible to HyperK', ...
    'Susceptible to All'};

cellType = 'Torord_endo';
beats = 101; ignoreFirst = 99;

ion_conditions = [120, 5; 160, 5; 140, 3; 140, 7; 140, 5];
labels_WT = {'WT - HypoNa (120 mM)', 'WT - HyperNa (160 mM)', ...
             'WT - HypoK (3 mM)', 'WT - HyperK (7 mM)', 'WT - Normal (140/5 mM)'};
colors = lines(4);
colors(5,:) = [0 0 0];  % Black for normal

multiplierFieldNames = {'INa', 'Ito', 'ICaL', 'IKr', 'IKs', 'IK1', ...
                        'INaCa_i', 'INaCa_ss', 'INaK', 'IKb', ...
                        'INab', 'ICab', 'IpCa'};
multiplierFieldLabels = {'Na', 'to', 'CaL', 'Kr', 'Ks', 'K1', ...
                         'NaCa_i', 'NaCa_{ss}', 'NaK', 'Kb', ...
                         'Nab', 'Cab', 'pCa'};

for subplotIdx = 1:length(cellIndices)
    cellIdx = cellIndices(subplotIdx);
    scalingVec = passScalingFactors(cellIdx,:);

    Vm_all_WT = cell(1, size(ion_conditions,1));
    time_all_WT = cell(1, size(ion_conditions,1));

    for i = 1:size(ion_conditions, 1)
        nao = ion_conditions(i,1);
        ko  = ion_conditions(i,2);

        param = struct();
        param.bcl = 1000;
        param.model = @model_Torord_Clancy;
        param.verbose = false;
        param.nao = nao;
        param.ko = ko;
        param.cao = 1.8;
        param.mt1 = 'dKPQ';
        param.mt2 = 'dKPQ';
        param.mf1 = 0;
        param.mf2 = 0;


        for k = 1:length(multiplierFieldNames)
            param.([multiplierFieldNames{k} '_Multiplier']) = scalingVec(k);
        end

        X0 = getStartingState(cellType);
        INa_SS1 = determine_INa_MC_steadystate(X0(1), 'dKPQ', 0);
        INa_SS2 = determine_INa_MC_steadystate(X0(1), 'dKPQ', 0);
        X0 = [X0, INa_SS1, INa_SS1, INa_SS2, INa_SS2];

        [time_raw, X] = modelRunner_Clancy(X0, [], param, beats, ignoreFirst);
        currents = getCurrentsStructure_Clancy(time_raw, X, param, ignoreFirst);
        Vm_all_WT{i} = currents.V;
        time_all_WT{i} = currents.time - currents.time(1);
    end

    % === PLOTTING ===
    fig = figure('Name', sprintf('Cell %d - %s', cellIdx, titleLabels{subplotIdx}), ...
            'Position', [400 200 500 600]);

    ax1 = axes(fig); hold(ax1, 'on');
    legendEntries = cell(1, length(Vm_all_WT));

    for i = 1:length(Vm_all_WT)
        plot(ax1, time_all_WT{i}, Vm_all_WT{i}, 'LineWidth', 1.3, 'Color', colors(i,:));
        legendEntries{i} = labels_WT{i};
    end

    if subplotIdx == 1
        lgd = legend(ax1, legendEntries, 'Location', 'northeast', 'FontSize', 16);
        lgd.Box = 'off';
    end

    ax1.XTick = [];
    ax1.XColor = 'none';
    barY = -110;
    blackLine = line(ax1, [1000 1500], [barY barY], 'Color', 'k', 'LineWidth', 4);
    blackLine.Annotation.LegendInformation.IconDisplayStyle = 'off';

    if subplotIdx == 6
        text(ax1, 250, barY + 5, '500 ms', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', 18);
    end

    if subplotIdx ~= 1 && subplotIdx ~= 4
        ax1.YTick = [];
        ax1.YColor = 'none';
        ylabel(ax1, '');
    end

    title(ax1, titleLabels{subplotIdx}, 'FontSize', 18);
    ylabel(ax1, 'V_m (mV)');
    xlim(ax1, [950,2000]);
    ylim(ax1, [-110, 100]);
    set(ax1, 'FontSize', 16)
    axis(ax1, 'square');
    ax1.Units = 'normalized';
    ax1.Position = [0.15, 0.42, 0.7, 0.5];

    % Bar plot
    ax2 = axes(fig);
    logVals = log2(scalingVec);
    bar(ax2, logVals, 'FaceColor', [0.2 0.4 1]);
    yline(ax2, 0, '--k');
    set(ax2, 'XTick', 1:13, ...
         'XTickLabel', multiplierFieldLabels, ...
         'XTickLabelRotation', 45, ...
         'TickLabelInterpreter', 'tex', 'FontSize', 16);
    ylim(ax2, [-2.2 2.2]);
    ax2.YTick = [-2 -1 0 1 2];
    ax2.YTickLabel = {'0.25', '0.5', '1', '2', '4'};
    set(ax2, 'FontSize', 16);
    box(ax2, 'on');
    ax2.Units = 'normalized';
    ax2.Position = [ax1.Position(1), 0.24, ax1.Position(3), 0.14];

    ax1.YLabel.FontSize = 16;
    ax1.YAxis.FontSize = 16;
    ax2.YLabel.FontSize = 16;
    ax2.YAxis.FontSize = 16;
end

% %% 
% rowSums = sum(sus_matrix, 2);              % Sum across each row
% targetRows = find(rowSums == 1);           % Rows with exactly one '1'
% 
% for i = 1:length(targetRows)
%     rowIdx = targetRows(i);
%     colIdx = find(sus_matrix(rowIdx, :) == 1);
%     fprintf('Row %d has exactly one 1 at column %d.\n', rowIdx, colIdx);
% end


%% WT- Normal, WT - HypoNa, WT - HyperNa, WT - HypoK, WT - HyperK

cellIndices = [1 , 55, 2656, 888, 4088, 14]; %same indices as above, wild-type traces
titleLabels = {
    'Susceptible to None', ...
    'Only Susceptible to HypoNa', ...
    'Only Susceptible to HyperNa', ...
    'Only Susceptible to HypoK', ...
    'Only Susceptible to HyperK', ...
    'Susceptible to All'};

cellType = 'Torord_endo';
beats = 101; ignoreFirst = 100;

ion_conditions = [120, 5; 160, 5; 140, 3; 140, 7; 140, 5];
labels_WT = {'WT - HypoNa (120 mM)', 'WT - HyperNa (160 mM)', ...
             'WT - HypoK (3 mM)', 'WT - HyperK (7 mM)', 'WT - Normal (140/5 mM)'};
colors = lines(4);
colors(5,:) = [0 0 0];  % Black for normal

multiplierFieldNames = {'INa', 'Ito', 'ICaL', 'IKr', 'IKs', 'IK1', ...
                        'INaCa_i', 'INaCa_ss', 'INaK', 'IKb', ...
                        'INab', 'ICab', 'IpCa'};
multiplierFieldLabels = {'Na', 'to', 'CaL', 'Kr', 'Ks', 'K1', ...
                         'NaCa_i', 'NaCa_{ss}', 'NaK', 'Kb', ...
                         'Nab', 'Cab', 'pCa'};

for subplotIdx = 1:length(cellIndices)
    cellIdx = cellIndices(subplotIdx);
    scalingVec = passScalingFactors(cellIdx,:);

    Vm_all_WT = cell(1, size(ion_conditions,1));
    time_all_WT = cell(1, size(ion_conditions,1));

    for i = 1:size(ion_conditions, 1)
        nao = ion_conditions(i,1);
        ko  = ion_conditions(i,2);

        param = struct();
        param.bcl = 1000;
        param.model = @model_Torord_Clancy;
        param.verbose = false;
        param.nao = nao;
        param.ko = ko;
        param.cao = 1.8;
        param.MT1 = 'dKPQ';
        param.MT2 = 'dKPQ';
        param.mf1 = 0;
        param.mf2 = 0;


        for k = 1:length(multiplierFieldNames)
            param.([multiplierFieldNames{k} '_Multiplier']) = scalingVec(k);
        end

        X0 = getStartingState(cellType);
        INa_SS1 = determine_INa_MC_steadystate(X0(1), 'dKPQ', 0);
        INa_SS2 = determine_INa_MC_steadystate(X0(1), 'dKPQ', 0);
        X0 = [X0, INa_SS1, INa_SS1, INa_SS2, INa_SS2];

        [time_raw, X] = modelRunner_Clancy(X0, [], param, beats, ignoreFirst);
        currents = getCurrentsStructure_Clancy(time_raw, X, param, ignoreFirst);
        Vm_all_WT{i} = currents.V;
        time_all_WT{i} = currents.time - currents.time(1);
    end

    % === PLOTTING ===
    fig = figure('Name', sprintf('Cell %d - %s', cellIdx, titleLabels{subplotIdx}), ...
            'Position', [400 200 500 600]);

    ax1 = axes(fig); hold(ax1, 'on');
    legendEntries = cell(1, length(Vm_all_WT));

    for i = 1:length(Vm_all_WT)
        plot(ax1, time_all_WT{i}, Vm_all_WT{i}, 'LineWidth', 1.3, 'Color', colors(i,:));
        legendEntries{i} = labels_WT{i};
    end

    if subplotIdx == 1
        lgd = legend(ax1, legendEntries, 'Location', 'northeast', 'FontSize', 13);
        lgd.Box = 'off';
    end

    ax1.XTick = [];
    ax1.XColor = 'none';
    barY = -110;
    blackLine = line(ax1, [0 500], [barY barY], 'Color', 'k', 'LineWidth', 4);
    blackLine.Annotation.LegendInformation.IconDisplayStyle = 'off';

    if subplotIdx == 6
        text(ax1, 250, barY + 5, '500 ms', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', 18);
    end

    if subplotIdx ~= 1 && subplotIdx ~= 4
        ax1.YTick = [];
        ax1.YColor = 'none';
        ylabel(ax1, '');
    end

    title(ax1, titleLabels{subplotIdx}, 'FontSize', 18);
    ylabel(ax1, 'V_m (mV)');
    xlim(ax1, [-20, max(cellfun(@(t) max(t), time_all_WT))]);
    ylim(ax1, [-110, 100]);
    set(ax1, 'FontSize', 16)
    axis(ax1, 'square');
    ax1.Units = 'normalized';
    ax1.Position = [0.15, 0.42, 0.7, 0.5];

    % Bar plot
    ax2 = axes(fig);
    logVals = log2(scalingVec);
    bar(ax2, logVals, 'FaceColor', [0.2 0.4 1]);
    yline(ax2, 0, '--k');
    set(ax2, 'XTick', 1:13, ...
         'XTickLabel', multiplierFieldLabels, ...
         'XTickLabelRotation', 45, ...
         'TickLabelInterpreter', 'tex', 'FontSize', 16);
    ylim(ax2, [-2.2 2.2]);
    ax2.YTick = [-2 -1 0 1 2];
    ax2.YTickLabel = {'0.25', '0.5', '1', '2', '4'};
    set(ax2, 'FontSize', 16);
    box(ax2, 'on');
    ax2.Units = 'normalized';
    ax2.Position = [ax1.Position(1), 0.24, ax1.Position(3), 0.14];

    ax1.YLabel.FontSize = 16;
    ax1.YAxis.FontSize = 16;
    ax2.YLabel.FontSize = 16;
    ax2.YAxis.FontSize = 16;
end

%% 
rowSums = sum(sus_matrix, 2);              % Sum across each row
targetRows = find(rowSums == 1);           % Rows with exactly one '1'

for i = 1:length(targetRows)
    rowIdx = targetRows(i);
    colIdx = find(sus_matrix(rowIdx, :) == 1);
    fprintf('Row %d has exactly one 1 at column %d.\n', rowIdx, colIdx);
end
