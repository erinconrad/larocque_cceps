%{
Written by Joshua LaRocque 2025
%}


rng("default")
fprintf('\nInitialized RNG\n');
nboot = 1e3;

file2analyze = '../data/all_vars_peds_BL_06022025.mat'; % now we can use HUP267 - previously case mismatch error
load(file2analyze)
remove_WM_resp = 1;
init_threshold = 1;

% If old field names, update to new ones
if ismember('Before_stim_BL',all_vars_store.Properties.VariableNames)
    all_vars_store = renamevars(all_vars_store,["Before_stim_BL","rms_scaled","rms_unscaled","rms_BL"],["trueBL","RMS","RMS_unscaled","RMS_BL"]);
end
all_vars_store.RMS_trueBL = all_vars_store.RMS_unscaled./all_vars_store.trueBL;

% adds field for string-formatted subject code - dont' want regression to
% treat it as a number
all_vars_store.substr = num2str(all_vars_store.subject);
all_vars_store.sub_resp_ch = strcat(all_vars_store.substr,'_',all_vars_store.resp_ch);
all_vars_store.sub_stim_ch = strcat(all_vars_store.substr,'_',all_vars_store.stim_ch);

if remove_WM_resp == 1
    all_vars_reg2 = all_vars_store(ismember(all_vars_store.resp_matter,[1 3]),:);%(ismember(all_vars_reg2.stim_matter,[2 4]) & ismember(all_vars_reg2.resp_matter,[1 3]),:);
else
    all_vars_reg2 = all_vars_store;
end
all_vars_reg2 = all_vars_reg2(~all_vars_reg2.stim_sz,:); %takes out CCEPS which caused a seizure

outlier_cfg.do_outliers = 0;
outlier_cfg.thr_fac = 10;
outlier_cfg.fancy=1;
do_plot = 1; threshold = 1; % at first, threshold above 1

all_vars_reg2.BL_mean_corr_abs = abs(all_vars_reg2.BL_mean_corr);
all_vars_reg2 = all_vars_reg2(~isnan(all_vars_reg2.Distance),:);

all_vars_reg2 = calculateResidualsForField(all_vars_reg2, 'RMS', threshold, outlier_cfg, do_plot);


c=1;
col_names = {'RMS_residuals'};%,'RMS_unscaled_residuals','RMS_BL_residuals','RMS_trueBL_residuals'}; % + BL_mean_corr_abs_residuals
col_name = col_names{c};

small_tab = all_vars_reg2(:,{'RMS_residuals','Distance','SOZ_resp','SOZ_stim','substr','sub_stim_ch','sub_resp_ch','subject'});

clear all_vars_store all_vars_reg2
glme = fitglme(small_tab,[col_name ' ~ 1 + Distance + SOZ_resp + SOZ_stim + SOZ_resp:SOZ_stim + (1|sub_resp_ch) + (1|sub_stim_ch) + (1|substr)']);


%{
% Erin check to see if the 0/1 numerical matters - it does not
% Ensure grouping vars for random effects are categorical
small_tab2 = small_tab;
small_tab2.sub_resp_ch = categorical(string(small_tab2.sub_resp_ch));
small_tab2.sub_stim_ch = categorical(string(small_tab2.sub_stim_ch));
small_tab2.substr      = categorical(string(small_tab2.substr));

% Convert SOZ indicators to categorical with explicit reference level
% (reference = 'nonSOZ', so the intercept is nonSOZ/nonSOZ)
small_tab2.SOZ_resp = categorical(small_tab2.SOZ_resp, [0 1], {'nonSOZ','SOZ'});
small_tab2.SOZ_stim = categorical(small_tab2.SOZ_stim, [0 1], {'nonSOZ','SOZ'});

% (If your SOZ columns might be logical or strings, you can coerce first:)
% small_tab.SOZ_resp = categorical(double(logical(small_tab.SOZ_resp)), [0 1], {'nonSOZ','SOZ'});
% small_tab.SOZ_stim = categorical(double(logical(small_tab.SOZ_stim)), [0 1], {'nonSOZ','SOZ'});

% Fit (no 'CategoricalVars' arg needed)
glme2 = fitglme(small_tab2, ...
  'RMS_residuals ~ 1 + Distance + SOZ_resp + SOZ_stim + SOZ_resp:SOZ_stim + (1|sub_resp_ch) + (1|sub_stim_ch) + (1|substr)', ...
  'DummyVarCoding','reference');


%}

small_tab.RMS_fitted = fitted(glme);
small_tab.residuals = residuals(glme);
store_coefs_randfx = table([],[],[],[],[],'VariableNames',glme.CoefficientNames);
for i = 1:nboot
    i
    newtab = bootstrap_RMS_residuals(small_tab,'subject','RMS_fitted');
    glme = fitglme(newtab,['bootstrapped ~ 1 + Distance + SOZ_resp + SOZ_stim + SOZ_resp:SOZ_stim + (1|sub_resp_ch) + (1|sub_stim_ch) + (1|substr)']);
    store_coefs_randfx(i,:) = array2table(fixedEffects(glme)');
    if mod(i,100)==0
        save(['../output/RMS_reg_residuals_CHRAND_' 'main_peds' '_bootstrapping_prog_0825.mat'],'store_coefs_randfx','newtab','glme','-mat')
    end
end

