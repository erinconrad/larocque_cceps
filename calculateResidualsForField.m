%{
Written by Joshua LaRocque 2025
%}



function all_vars_reg = calculateResidualsForField(all_vars_sort, fieldName, threshold, outlier_cfg, do_plot)
%% Function to add field to existing table with the residuals of subject by subject regression   
% For the 
% Filter data based on a threshold for the specified field (if threshold
% exists - otherwise, do not threshold
do_outlier = outlier_cfg.do_outliers;
if isfield(outlier_cfg, 'fancy')
    fancy = outlier_cfg.fancy;
else
    fancy = 0;
end
%thr_fac = outlier_cfg.thr_fac;

if do_outlier 
    if issorted(all_vars_sort.Distance)
        all_vars_reg = all_vars_sort(all_vars_sort.(fieldName) >= threshold, :);
    else
        all_vars_sort = sortrows(all_vars_sort, 'Distance');
        if ~isempty(threshold)
            all_vars_reg = all_vars_sort(all_vars_sort.(fieldName) >= threshold, :);
        else
            all_vars_reg = all_vars_sort;
        end
    end
else
    if ~isempty(threshold)
        all_vars_reg = all_vars_sort(all_vars_sort.(fieldName) >= threshold, :);
    else
        all_vars_reg = all_vars_sort;
    end
end
    
    % Assuming you have your data in a table named 'all_vars'
   
    % Initialize the residuals column with NaN values
    residualsFieldName = [fieldName '_residuals'];
    all_vars_reg.(residualsFieldName) = NaN(height(all_vars_reg), 1);

    a = unique(all_vars_reg.SOZ_stim); %using only non-SOZ electrodes
    b = unique(all_vars_reg.SOZ_resp); 
    
    regress_idx = all_vars_reg.SOZ_resp == b(1) & all_vars_reg.SOZ_stim == a(1);

        % Extract the variables for the current subject
        Distance = all_vars_reg.Distance;
        var2use = all_vars_reg.(fieldName);
        
        % Perform the regression analysis for the current subject
        [xData, yData] = prepareCurveData(Distance(regress_idx), var2use(regress_idx));

        % Set up fittype and options.
        ft = fittype('rat01');
        opts = fitoptions('Method', 'NonlinearLeastSquares');
        opts.Display = 'Off';
        opts.Lower = [-Inf -1];
        opts.Robust = 'Bisquare';
        opts.StartPoint = [0.1 0.1];

        % Fit model to data.
        [mdl, gof, ~] = fit(xData, yData, ft, opts);

        % Calculate the residuals for the current subject
        residuals = var2use - feval(mdl, Distance);
        
        % Store the residuals in the corresponding rows of the residuals field
        all_vars_reg.(residualsFieldName) = residuals;

        % Optionally, you can plot the data and the regression curve for each subject
        if do_plot
            figure;
            % Plot the red points
            plot(Distance(~regress_idx), var2use(~regress_idx), '.', 'Color', [1 0 0]); hold on
            % Plot the blue points used for fitting
            blue_points = plot(Distance(regress_idx), var2use(regress_idx), '.', 'Color', [0 0 1]);
            % Plot the line of best fit (in green)
            h = plot(mdl, Distance(regress_idx), var2use(regress_idx));
            h(2).LineWidth = 2; h(2).Color = [0 1 0];
            
            % Add legend entries
            legend({'SOZ', 'non-SOZ', 'Best Fit'}, 'Location', 'Best');
            hold off;
            
            xlabel('Distance (mm)')
            title(['Regression for ' fieldName ', all subjects: ' sprintf(', P1 = %0.2f, Q1 = %0.2f, R2 = %0.2f; threshold %0.1f', mdl.p1, mdl.q1, gof.rsquare, threshold)], 'Interpreter', 'none', 'FontWeight', 'bold');
            
            if fancy==1
                set(gca, 'Box', 'off')
                set(gca, 'LineWidth', 2)
                set(gca, 'FontSize', 12)
            end

        end
    
end
