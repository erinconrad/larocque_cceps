%{
Written by Joshua LaRocque 2025
%}

function tbl_out = bootstrap_RMS_residuals(tbl, subject_col, rms_fitted_col)
    % Check if required columns exist
    requiredVars = {subject_col, rms_fitted_col, 'residuals'};
    if ~all(ismember(requiredVars, tbl.Properties.VariableNames))
        error("Input table must contain specified Subject and RMS_fitted columns, and residuals");
    end
    
    % Ensure Subject column is treated as categorical for grouping
    if isnumeric(tbl.(subject_col)) || isstring(tbl.(subject_col))
        tbl.(subject_col) = string(tbl.(subject_col)); % Convert to string for consistency
    end
    
    % Preallocate output column
    bootstrapped = nan(height(tbl), 1);
    
    % Get unique subjects
    uniqueSubjects = unique(tbl.(subject_col));
    
    for i = 1:numel(uniqueSubjects)
        subj = uniqueSubjects(i);
        idx = tbl.(subject_col) == subj;
        
        % Extract rows for the subject
        try
            fitted_vals = tbl.(rms_fitted_col)(idx);
        catch
            
        end

        residual_vals = tbl.residuals(idx);
        
        % Bootstrap: Add each RMS_fitted to a randomly sampled residual
        bootstrapped(idx) = fitted_vals + residual_vals(randi(length(residual_vals), size(fitted_vals)));
    end
    
    % Append the new column to the table
    tbl_out = tbl;
    tbl_out.bootstrapped = bootstrapped;
end

