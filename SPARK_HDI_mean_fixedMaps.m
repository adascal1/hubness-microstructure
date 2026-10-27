function [HDI,HEI, error_HDI, error_HEI, x_nn, y_nn, R2]=SPARK_HDI_mean_fixedMaps(patient_data, healthy_data, group, computationType, idx, varargin);
%% This script will calculate the HDI and HEI values in surface or volume, and write the plotted points to a table
% SPARK_HDI_mean_fixedMaps(patient_data, healthy_data, group)
% Mandatory Parameters: 
% 'patient_data' : The group of patient data along the cortex, for the whole cortex.
% 'healthy_data' : The group of healthy data along the cortex, for the whole cortex.
% 'Group'        : Name of Subject Group
% 'computationType': Whether you want to use cohen's d or full data values for the calculation
% 'idx;
% Optional Parameters: 
% 'label'     : the name of the region being studied
% 'writeTable': whether we want the plotted points to be saved as a .csv file. Default set to 'true'
% 'toPlot'    : whether we want to plot the the graph. Default set to 'true'
% 'analysis'  : options are 'surf' for surface based computation or 'vol' for volume based computation. Default set to surface.

    % Preparing Function
    opt = inputParser;

    % Adding parameters with default values and validation
    addParameter(opt, 'label', '', @ischar); % no default, just validation
    addParameter(opt, 'writeTable', true, @(x) islogical(x) || isnumeric(x)); % default is true
    addParameter(opt, 'toPlot', true, @(x) islogical(x) || isnumeric(x)); % default is true
    addParameter(opt, 'xLim', false, @(x) islogical(x) || isnumeric(x));
    addParameter(opt, 'yLim', false, @(x) islogical(x) || isnumeric(x));


    % parse input
    parse(opt, varargin{:});

    % Extract Parsed results
    label      = opt.Results.label;
    writeTable = logical(opt.Results.writeTable);
    toPlot     = logical(opt.Results.toPlot);
    xLim =  opt.Results.xLim;
    yLim =  opt.Results.yLim;


    % Checking we have proper input types
    if ~(strcmp(computationType, 'cohend') || strcmp(computationType, 'originalMethod'))
        error('Input must be cohend or originalMethod')
    end

    % Loading Data and only keeping the ROIs
    patient_data = patient_data(:, idx);
    healthy_data = healthy_data(:, idx);

    if strcmp(computationType, 'originalMethod')
    % assigning data to correct axis
        x_axis_data = mean(healthy_data)./std(healthy_data);
        y_axis_data = (mean(patient_data) - mean(healthy_data))./std(healthy_data);
        nameXaxis = 'Normalized Mean of Hemispheric Average of Controls'; nameYaxis = {'Normalized Difference Between Average Epilepsy',...
                                                                                        'and Hemispheric Average of Controls'};
    elseif strcmp(computationType, 'cohend')
        x_axis_data = mean(healthy_data) ./ std(healthy_data);
        y_axis_data = patient_data;
        nameXaxis = 'Normalized Mean of Hemispheric Average of Controls'; nameYaxis = "Cohen's D";

    end
   
    % removing any x-y pairs which have either a nan or an inf 
    x_nn=x_axis_data((~isnan(x_axis_data)) & (~isinf(x_axis_data)) & (~isnan(y_axis_data)) & (~isinf(y_axis_data)));
    y_nn=y_axis_data((~isnan(x_axis_data)) & (~isinf(x_axis_data)) & (~isnan(y_axis_data)) & (~isinf(y_axis_data)));
    
    % generating the linear regression model
    mdl = fitlm(x_nn, y_nn);
    
    % Saving HEI and HDI Values
    HEI = mdl.Coefficients.Estimate(1);
    HDI = mdl.Coefficients.Estimate(2);
        
    error_HEI = mdl.Coefficients.SE(1);
    error_HDI = mdl.Coefficients.SE(2);
        
    % saving p-value of the slope
    p_value_HDI = mdl.Coefficients.pValue(2);
    p_value_HEI = mdl.Coefficients.pValue(1);
    
    R2 = mdl.Rsquared.Ordinary;
    
    if toPlot
    
        % plotting
        fig = figure('Position', [680 302 806 576]);
        g  = gramm('x', x_nn, 'y', y_nn);
        g.geom_point();
        g.stat_glm('disp_fit', true);
        g.set_names('x', nameXaxis, 'y', nameYaxis);
        g.set_title(sprintf('Hub Disruption Index \n %s', label));
        g.set_color_options('map',[0 0 0], 'n_color',1,'n_lightness',1);
        g.set_text_options('font', 'Gillius ADF', 'base_size', 16)
        if isnumeric(xLim)
            g.axe_property('XLim',xLim, 'YLim', yLim)
        end
        g.draw();
        set([g.results.stat_glm.area_handle],...
            'FaceColor',[0.4 0.4 0.4]);
        
        set(g.results.stat_glm.text_handle, ...
                'String', sprintf('HDI = %0.3f \nHEI = %0.3f \nR2  = %0.3f', HDI, HEI,  mdl.Rsquared.Ordinary), ...
                'Position', [ 0.7 0.975 1], ...
                'FontSize', 15, ...
                'FontName', 'Gillius ADF');     
   
        set(g.results.geom_point_handle, ... 
                                'MarkerEdgeColor', [ 0 0 0 ], ...
                                'SizeData', 40, ...
                                'LineWidth', 0.3, ...
                                'MarkerFaceColor', [0.4 0.4 0.4]);
         % set(g.facet_axes_handles, ...
         %                        'XTick',g.facet_axes_handles.XTick([1,length(g.facet_axes_handles.XTick)]), ...
         %                        'YTick',g.facet_axes_handles.YTick([1,length(g.facet_axes_handles.YTick)]), ...
         %                        'TickLength', [0 0]);

         set(g.facet_axes_handles, ...
                                'TickLength', [0 0]);
        % saving
        saveas(fig, sprintf('/mfip/mfip1/arielle/MICs_dataset/05_SPARK_ANALYSIS/figures/o2e_hubDisruptionPlot_method-%s_ROI-%s_group-%s.svg',computationType, strrep(label, ' ', ''), group))
                
        clear g;
    end
  
    % Writing to table of points
    if writeTable
        T = table(x_nn', y_nn','VariableNames', {'x', 'y'});
        writetable(T, sprintf('/mfip/mfip1/arielle/MICs_dataset/05_SPARK_ANALYSIS/tables/o2e_HubDisruptionPlot_method-%s_ROI-%s_group-%s.csv', computationType, strrep(label, ' ', ''), group))
    end
end

