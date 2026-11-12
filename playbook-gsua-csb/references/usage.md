# GSUA-CSB: real usage outside the toolbox

Mechanical scan (read-only) of .m files under: `C:\Users\Lanre\MATLAB Drive`, `D:\MATLAB Drive.source-of-truth-20260718T095207`, `C:\Users\Lanre\Documents`, `D:\Other-projects\Semillero Matemáticas`, `D:\Other-projects\generalizedTransmission`, `D:\Other-projects\VIH`. Comments skipped; files mirrored under C: and D: counted once. Data for the pack author, not conclusions.

| function | calls | files |
|---|---|---|
| `gsua_dataprep` | 137 | 128 |
| `gsua_eval` | 93 | 52 |
| `gsua_pe` | 81 | 62 |
| `gsua_dmatrix` | 71 | 52 |
| `gsua_deval` | 34 | 28 |
| `gsua_ia` | 33 | 30 |
| `gsua_costf` | 20 | 19 |
| `gsua_sa` | 18 | 13 |
| `gsua_ua` | 12 | 12 |
| `gsua_plot` | 10 | 7 |
| `gsua_covmetric` | 10 | 10 |
| `gsua_csb` | 6 | 6 |
| `gsua_depth` | 6 | 3 |
| `gsua_dia` | 5 | 3 |
| `gsua_oatr` | 4 | 4 |
| `gsua_likelihood` | 3 | 3 |
| `gsua_corr2omitnan` | 3 | 3 |
| `gsua_load` | 3 | 3 |
| `gsua_costfMulti` | 3 | 2 |
| `gsua_intrp` | 2 | 2 |
| `gsua_save` | 2 | 1 |
| `gsua_pardeval` | 2 | 2 |
| `gsua_dpmat` | 1 | 1 |
| `gsua_bfm` | 1 | 1 |
| `gsua_noisefloor` | 1 | 1 |
| `gsua_rcostf` | 1 | 1 |
| `gsua_MCF` | 1 | 1 |
| `gsua_bootstrap` | 1 | 1 |

## Examples (up to 3 per function)

### `gsua_dataprep`
- `C:\Users\Lanre\MATLAB Drive\chimerita\create_model.m`:32: `[T, ~] = gsua_dataprep(ModelName, RangeT, 'domain', Domain, 'names', FullNames, 'out_names', OutNames);`
- `C:\Users\Lanre\MATLAB Drive\CHIMERITApaper\Daniel\exprun.m`:13: `[T,~] = gsua_dataprep("JCHIMERA_SIR", RangeT, 'domain', Domain, 'names',...`
- `C:\Users\Lanre\MATLAB Drive\enviar\Examples\default_example.m`:10: `[T,~] = gsua_dataprep(odes,vars,[0 264],'SIR','range',[0 0; 100 3000; 1 500; 0 1; 0 1],'output',3);`

### `gsua_eval`
- `C:\Users\Lanre\MATLAB Drive\ChimeraAedes\Chimera SEIR Vectorial\5 eqs\Otras localidades\EstimationLocalities.m`:21: `gsua_eval(EstimatedParams(1,:)', T,0:D-1,ydata);`
- `C:\Users\Lanre\MATLAB Drive\CHIMERITApaper\Daniel\exprun.m`:45: `gsua_eval(Parameters(:,1), T, xdata,ydata(1:2,:));`
- `C:\Users\Lanre\MATLAB Drive\Gaming place\izi.m`:15: `gsua_eval(prop.Solver,T.gaming,T,xdata,y_nom);`

### `gsua_pe`
- `C:\Users\Lanre\MATLAB Drive\ChimeraAedes\Chimera SEIR Vectorial\5 eqs\Otras localidades\EstimationLocalities.m`:13: `[T_est, ResTemp] = gsua_pe(T,xdata, ydata,'solver',Solver,'N',1,'opt',Opt,'save',false);`
- `C:\Users\Lanre\MATLAB Drive\CHIMERITApaper\Daniel\exprun.m`:31: `[T_aux, Res_aux] = gsua_pe(T, xdata, ydata(1:2,:),...`
- `C:\Users\Lanre\MATLAB Drive\estimationTrouble\classic mortality\main_china.m`:106: `[T_est,res,IP]=gsua_pe(T,xdata,ydata,'solver',solver,'N',1,'opt',opt,'save',false);`

### `gsua_dmatrix`
- `C:\Users\Lanre\MATLAB Drive\Henry-working\testR.m`:24: `M = gsua_dmatrix(T_oat,5000);`
- `C:\Users\Lanre\MATLAB Drive\Henry-working\testR.m`:46: `M = gsua_dmatrix(T_oat,5000);`
- `C:\Users\Lanre\MATLAB Drive\Henry-working\testR.m`:68: `M = gsua_dmatrix(T_oat,5000);`

### `gsua_deval`
- `C:\Users\Lanre\MATLAB Drive\Gaming place\izi.m`:9: `y_nom=gsua_deval(T.Nominal',prop,xdata);`
- `C:\Users\Lanre\MATLAB Drive\Henry-working\lossFunction\script_covid.m`:50: `yrawdata = gsua_deval(paramCurves(:,w)',T2,xdata);`
- `C:\Users\Lanre\MATLAB Drive\Henry-working\lossFunction\script_exps.m`:50: `yrawdata = gsua_deval(paramCurves(:,w)',T2,xdata);`

### `gsua_ia`
- `C:\Users\Lanre\MATLAB Drive\Trabajo con Daniel\ModeloBasadoEnAgentesVectoresHumanos\TOOLBOX_SVH\f25_gsua_pilot.m`:72: `Tia = gsua_ia(Te, est(:, acc), false, false, false);`
- `C:\Users\Lanre\MATLAB Drive\Trabajo con Daniel\Wolbachia-ThermalMismatch\model\calibration\endemicPool.m`:48: `Tia=gsua_ia(T,pool,true,false,false);`
- `C:\Users\Lanre\MATLAB Drive\Trabajo con Daniel\Wolbachia-ThermalMismatch\model\calibration\redoCaliPhase2.m`:119: `Tia = gsua_ia(Te, EstFiltered, true, false, false);`

### `gsua_costf`
- `C:\Users\Lanre\MATLAB Drive\Trabajo con Daniel\ModeloBasadoEnAgentesVectoresHumanos\TOOLBOX_SVH\f25_estimate.m`:67: `xi = @(yy, f) gsua_costf(1, mean((yy - omega * yy).^2), numel(yy), yy, f, costAlpha);`
- `C:\Users\Lanre\MATLAB Drive\Trabajo con Daniel\ModeloBasadoEnAgentesVectoresHumanos\TOOLBOX_SVH\f25_gsua_pilot.m`:44: `xi = @(y, f) gsua_costf(1, mean((y - omega * y).^2), numel(y), y, f, costAlpha);`
- `C:\Users\Lanre\MATLAB Drive\Trabajo con Daniel\ModeloBasadoEnAgentesVectoresHumanos\TOOLBOX_SVH\f25_n3_srange.m`:46: `xi = gsua_costf(1, mean((yref - 1.3 * yref).^2), numel(yref), yref, yy, 2);`

### `gsua_sa`
- `C:\Users\Lanre\MATLAB Drive\Henry-working\testR.m`:25: `T_sa = gsua_sa(M,T_oat,'parallel',false);`
- `C:\Users\Lanre\MATLAB Drive\Henry-working\testR.m`:47: `T_sa = gsua_sa(M,T_oat,'parallel',false);`
- `C:\Users\Lanre\MATLAB Drive\Henry-working\testR.m`:69: `T_sa = gsua_sa(M,T_oat,'parallel',false);`

### `gsua_ua`
- `C:\Users\Lanre\MATLAB Drive\SC\IETS\Daniel\SemanalReborn\UA.m`:268: `y = gsua_ua(Experiments',T,'parallel',true);`
- `C:\Users\Lanre\MATLAB Drive\Trabajo con Daniel\Wolbachia-ThermalMismatch\model\calibration\redoCaliPhase2.m`:122: `Y_diag = gsua_ua(M_diag, T, 'xdata', xdata, 'ynom', ywork2, 'parallel', true);`
- `C:\Users\Lanre\MATLAB Drive\Trabajo con Daniel\Wolbachia-ThermalMismatch\model\calibration\runCaliCounterfactual.m`:174: `Y_diag = gsua_ua(M_diag, T, 'xdata', xdata, 'ynom', ywork2, 'parallel', true);`

### `gsua_plot`
- `C:\Users\Lanre\MATLAB Drive\Henry-working\testR.m`:28: `gsua_plot('Pie',T_sa,T_sa.STi)`
- `C:\Users\Lanre\MATLAB Drive\Henry-working\testR.m`:50: `gsua_plot('Pie',T_sa,T_sa.STi)`
- `C:\Users\Lanre\MATLAB Drive\Henry-working\testR.m`:72: `gsua_plot('Pie',T_sa,T_sa.STi)`

### `gsua_covmetric`
- `C:\Users\Lanre\MATLAB Drive\Trabajo con Daniel\Wolbachia-ThermalMismatch\model\calibration\redoCaliPhase2.m`:123: `[cost_data, cost_band] = gsua_covmetric(Y_diag, ywork2, 'margin', 0.1);`
- `C:\Users\Lanre\MATLAB Drive\Trabajo con Daniel\Wolbachia-ThermalMismatch\model\calibration\runCaliCounterfactual.m`:175: `[cost_data, cost_band] = gsua_covmetric(Y_diag, ywork2, 'margin', 0.1);`
- `C:\Users\Lanre\MATLAB Drive\Trabajo con Daniel\Wolbachia-ThermalMismatch\model\calibration\runStageA.m`:130: `[cost_data, cost_band] = gsua_covmetric(Y_diag, ywork2, 'margin', 0.1);`

### `gsua_csb`
- `C:\Users\Lanre\MATLAB Drive\Henry-working\uniquenessExp.m`:11: `[tcsb,a,b,c,d] = gsua_csb(T_oat,N,'method','regression');`
- `C:\Users\Lanre\MATLAB Drive\Henry-working\ongoingfuncs\uniquenessExp.m`:11: `[tcsb,a,b,c,d] = gsua_csb(T_oat,N);`
- `C:\Users\Lanre\MATLAB Drive\Henry-working\ongoingfuncs\definitive tests\CSBPrueba.m`:78: `[ParU, ~, ~, ~,sup, cont] = gsua_csb(T_oat, N, 'reps', reps, 'protect', true)`

### `gsua_depth`
- `C:\Users\Lanre\MATLAB Drive\Trabajo con Daniel\Wolbachia IETS\Modelo\Refinamiento\plot_dengue_forecast.m`:144: `[depthY2_sub, idxY2_sub] = gsua_depth(xDepthY2, parallel);   %#ok<ASGLU>`
- `C:\Users\Lanre\MATLAB Drive\Trabajo con Daniel\Wolbachia IETS\Modelo\Refinamiento\plot_dengue_forecast.m`:153: `[depthY_sub, idxY_sub] = gsua_depth(xDepthY, parallel);      %#ok<ASGLU>`
- `C:\Users\Lanre\MATLAB Drive\Trabajo con Daniel\Wolbachia IETS\Modelo\Ryoiki Tenkai\plot_dengue_forecast.m`:144: `[depthY2_sub, idxY2_sub] = gsua_depth(xDepthY2, parallel);   %#ok<ASGLU>`

### `gsua_dia`
- `C:\Users\Lanre\MATLAB Drive\ValidacionChimera\ValidationFunctions\vc_report.m`:125: `T2 = gsua_dia(T, estimations2, false);`
- `C:\Users\Lanre\MATLAB Drive\ValidacionChimera\ValidationFunctions\vc_report.m`:183: `T = gsua_dia(T, T.Est);%Identifiability table`
- `C:\Users\Lanre\MATLAB Drive\ValidacionChimera\ValidationFunctions\vc_validation.m`:35: `T = gsua_dia(T, estimations);                   % Identifiability table`

### `gsua_oatr`
- `C:\Users\Lanre\MATLAB Drive\Henry-working\ongoingfuncs\definitive tests\CSBPrueba.m`:64: `T_oat = gsua_oatr(T);`
- `C:\Users\Lanre\MATLAB Drive\Intervals article\multiple_uci\procesamiento.m`:4: `T=gsua_oatr(T);`
- `C:\Users\Lanre\MATLAB Drive\ValidacionChimera\ValidationFunctions\vc_validation.m`:61: `Toat = gsua_oatr(Toat, 'parallel',parallel, 'correct', true);`

### `gsua_likelihood`
- `C:\Users\Lanre\MATLAB Drive\Trabajo con Daniel\Wolbachia-ThermalMismatch\model\calibration\profileLikelihood_nu8nu9.m`:165: `Tfinal = gsua_likelihood(T, xdata2, ywork2, 0.95, 0.25, 1.1, 0.05, 0.01, 12, 3, ...`
- `C:\Users\Lanre\MATLAB Drive\Trabajo con Daniel\Wolbachia-ThermalMismatch-v2\model\calibration\profileLikelihood_nu8nu9.m`:165: `Tfinal = gsua_likelihood(T, xdata2, ywork2, 0.95, 0.25, 1.1, 0.05, 0.01, 12, 3, ...`
- `D:\Other-projects\VIH\Estimacion_Medellin\estimacion_VIH_Medellin.m`:169: `Tperfil = gsua_likelihood(Tpe, xdata, ydata, 0.95, 0.05, 1.1, 0.01, 0.01, 15, 1, false, usarParalelo, false, [], false);`

### `gsua_corr2omitnan`
- `C:\Users\Lanre\MATLAB Drive\Trabajo con Daniel\Wolbachia-ThermalMismatch-v2\model\calibration\priorUA.m`:62: `ci=(2-gsua_corr2omitnan(pin,introData))*mean((pin(vI)-introData(vI)).^2)/regI;`
- `C:\Users\Lanre\MATLAB Drive\Trabajo con Daniel\Wolbachia-ThermalMismatch-v2\model\calibration\priorUA_golden2024_lowEfum005.m`:161: `ci = (2-gsua_corr2omitnan(pin,introData)) * mean((pin(vI)-introData(vI)).^2) / regI;`
- `C:\Users\Lanre\MATLAB Drive\Trabajo con Daniel\Wolbachia-ThermalMismatch-v2\model\calibration\recreate_priorUA_earlyfum_lowEfum.m`:120: `ci = (2-gsua_corr2omitnan(pin,introData)) * mean((pin(vI)-introData(vI)).^2) / regI;`

### `gsua_load`
- `C:\Users\Lanre\MATLAB Drive\ValidacionChimera\main.m`:130: `T = gsua_load(T);`
- `C:\Users\Lanre\MATLAB Drive\ValidacionChimera\ValidationFunctions\vc_report.m`:180: `T = gsua_load(T);`
- `C:\Users\Lanre\MATLAB Drive\ValidacionChimera\ValidationFunctions\WriteDataJSON.m`:11: `T = gsua_load(T);`

### `gsua_costfMulti`
- `C:\Users\Lanre\MATLAB Drive\ValidacionChimera\ValidationFunctions\vc_report.m`:200: `cos=gsua_costfMulti(Y_test(:,Cut(1):Cut(2),:),ycut(:,Cut(1):Cut(2)),margin,false);`
- `C:\Users\Lanre\MATLAB Drive\ValidacionChimera\ValidationFunctions\vc_validation.m`:49: `activator = gsua_costfMulti(ymedian, ybest, 0.3, false);`
- `C:\Users\Lanre\MATLAB Drive\ValidacionChimera\ValidationFunctions\vc_validation.m`:81: `costf = gsua_costfMulti(ydata, yfun, margin, parallel);% UA cost function`

### `gsua_intrp`
- `C:\Users\Lanre\MATLAB Drive\SC\IETS\Daniel\SemanalReborn\UA.m`:280: `y=squeeze(gsua_intrp(fun(M(i,:),domain,opt),xdata,out));`
- `C:\Users\Lanre\MATLAB Drive\SC\IETS\Daniel\SemanalReborn\CodigosGeneral\UA.m`:280: `y=squeeze(gsua_intrp(fun(M(i,:),domain,opt),xdata,out));`

### `gsua_save`
- `C:\Users\Lanre\MATLAB Drive\ValidacionChimera\main.m`:108: `gsua_save(T)`
- `C:\Users\Lanre\MATLAB Drive\ValidacionChimera\main.m`:177: `gsua_save(T)`

### `gsua_pardeval`
- `C:\Users\Lanre\MATLAB Drive\ValidacionChimera\ValidationFunctions\vc_report.m`:197: `Y_test=gsua_pardeval(M,T,xdata,false);`
- `C:\Users\Lanre\MATLAB Drive\ValidacionChimera\ValidationFunctions\vc_validation.m`:79: `ydata = gsua_pardeval(M, T, xdata, parallel);   % UA data`

### `gsua_dpmat`
- `C:\Users\Lanre\MATLAB Drive\Henry-working\ongoingfuncs\definitive tests\CSBPrueba.m`:58: `[T,~] = gsua_dpmat(odes,vars,[0 39],'7m','output',2,'opt',opts,'Range',Range);`

### `gsua_bfm`
- `C:\Users\Lanre\MATLAB Drive\Intervals article\multiple_uci\procesamiento.m`:7: `[New_range,T3,J_test,Y_test,sup]= gsua_bfm(T,1000);`

### `gsua_noisefloor`
- `C:\Users\Lanre\MATLAB Drive\Trabajo con Daniel\ModeloBasadoEnAgentesVectoresHumanos\TOOLBOX_SVH\f25_gsua_pilot.m`:64: `nf = gsua_noisefloor(Te, xdata, y, est, res, 'B', 500, 'seed', 1);`

### `gsua_rcostf`
- `C:\Users\Lanre\MATLAB Drive\ValidacionChimera\ValidationFunctions\CheckingUpgrade.m`:13: `PartialResidual = gsua_rcostf(ydataCut', ysimCut', bound);`

### `gsua_MCF`
- `C:\Users\Lanre\Documents\Profesor\2026-2\Modelacion3Semillero\03_Incertidumbre_MonteCarlo\matlab\demo_03_incertidumbre.m`:117: `[Jsup, Jinf] = gsua_MCF(T, M, Y, ydata);`

### `gsua_bootstrap`
- `D:\Other-projects\VIH\Estimacion_Medellin\estimacion_VIH_Medellin.m`:158: `[estBootstrap, costosBootstrap] = gsua_bootstrap(Tpe, xdata, ydata, Tpe.Estfmincon, costos, ...`

## Octopus lessons and failures mentioning GSUA

- (modelobasadoenagentesvectoreshumanos-be5d220643) E5 methodology redesign (2026-07-24) — supersedes fixed nEarly/nLate windows; STATUS marker for other agents
- (modelobasadoenagentesvectoreshumanos-be5d220643) CURRENT PLAN (2026-10-09): 4 determinant factors -> 4 intermediate quantities; tracked in TOOLBOX_SVH/PLAN_FACTORES_TRANSMISION.md
- (modelobasadoenagentesvectoreshumanos-be5d220643) Step 7 result (2026-10-13): reduced model + closure reproduces ABM trajectories
- (modelobasadoenagentesvectoreshumanos-be5d220643) Plan de fase 2 aprobado (2026-10-13): constantes del cierre no congeladas
- (modelobasadoenagentesvectoreshumanos-be5d220643) fmincon interior-point abandons starts placed on a bound
- (sc-82233348c4) Cómo ejecutar MATLAB y Python en SC sin fallos
- (semillero-matem-ticas-743f12f23b) Correr MATLAB: el MCP no conecta, usar matlab -batch
- (toolbox-684dd703ef) GSUA-CSB toolbox architecture, workflow, and current doc-pass state
- (toolbox-684dd703ef) Real Python interpreter is Anaconda, not the PATH stub
- (toolbox-684dd703ef) MATLAB ode45 vs scipy solve_ivp RK45 diverge at default tolerance on sensitive ODEs
- (toolbox-684dd703ef) GSUA-CSB toolbox-presentation paper: two-track state + playbook mapping
- (toolbox-684dd703ef) New Octopus playbook redaccion-software (software/toolbox papers)
- (vih-4e66e1ec43) Precondiciones de calibración implementadas
- (vih-4e66e1ec43) Estado (2026-11-05): piloto Medellín listo para correr el ciclo completo en GSUA-CSB
- (vih-4e66e1ec43) Resultados preliminares Medellín (corrida corta, 2026-11-05) y documentos para supervisores
- (wolbachia-iets-aa3b83cbe1) Calibration acceptance standard settled: containment-and-shape, covmetric as diagnostic
- (wolbachia-iets-aa3b83cbe1) 2026-07-24 session: model verification complete, ready to start Bello/Itagui runs
- (wolbachia-iets-aa3b83cbe1) Cali methodology confirmed: single-stage counterfactual, no release data
- (wolbachia-iets-aa3b83cbe1) Bello/Itagui Stage A + Cali runs complete, reviewed and accepted (2026-07-26)
- (wolbachia-iets-aa3b83cbe1) Cluster-of-stories result: Wolbachia loss robust across 3 explanatory families; mechanism NOT identifiable (2026-09-04)
- (wolbachia-iets-aa3b83cbe1) MaxEnt golden-expansion pipeline validated + launched (7000 starts, ~245 golden expected) (2026-09-04)
- (wolbachia-iets-aa3b83cbe1) MaxEnt golden harvest RESULT: 162 new golden (4.6% yield), golden set now ~194; harvest stalled at block 8 on a pathological start (2026-09-04)
- (wolbachia-thermalmismatch-v2-9e259e175f) runStageA smoke-size bug fixed; gsua_ia still unsuitable for tiny pools
- (wolbachia-thermalmismatch-v2-9e259e175f) Paper figure provenance regenerated 2026-10-10 (commits a1ddb04, 674cfab)
- (wolbachia-thermalmismatch-v2-9e259e175f) Local runs 2026-10-10: GSA rerun, unbiased band, pool/NB sensitivity; F12/F13 (commits to 0783974)
- (wolbachia-thermalmismatch-v2-9e259e175f) Canonical conditional ensemble = tmax; condEnsemble.m switch; outputs regenerated
- (wolbachia-thermalmismatch-v2-9e259e175f) F12 re-run launched: runStageB_onL (8 candidates x Bello, Medellin)
- (wolbachia-thermalmismatch-v2-9e259e175f) Round prep 2026-10-12: clusters, S1 Table, abstract, form pass
