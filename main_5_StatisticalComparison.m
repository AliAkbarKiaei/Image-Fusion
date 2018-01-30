clc
clear
close all
addpath(strcat(pwd,'/functions'));
addpath(strcat(pwd,'/images'));
addpath(strcat(pwd,'/fusion metrics'));

warning('OFF', 'MATLAB:xlswrite:AddSheet');
%% Defining Metrics, Fusion Methods, Source Images, excel file that you want to write the results in
%filter='_allMetrics'; metricNames = {'metric_MI', 'metric_MI2', 'metric_ssim', 'std2', 'entropy', 'metric_Edge_Intensity', 'metric_PSNR', 'metric_AverageGradient',...
%  'metric_Qabf', 'metricXydeas', 'metricCvejic', 'metricPeilla', 'metricChen', 'metricChenBlum', 'metricWang', 'metricYang', 'metricZhao'}' % , 'metricZheng', 'metricPWW', 'metricHossny'}'

filter='_1_12'; metricNames = {'entropy', 'metric_PSNR', 'metric_ssim', 'metric_MI2', 'metric_Qabf', 'metricChen','metricWang', 'metricPeilla', 'metricZhao'}'

choice = questdlg('Which experiments do you like to Evaluate?', ...
    'Menu', ...
    'expr 1-4','expr 5-8', 'expr 1-12', 'expr 1-4');
% Handle response
switch choice
    case 'expr 1-4'
        fusionMethodNames = {'PCNN_NSCT', 'm_PCNN', 'SCM_F', 'NSCT', 'NSCT_SR', 'SCM_M', 'Del_PCA', 'Del_max', 'Del_weighted'}'
        sourceImages = {'G01' 'CT' 'MR_T2' ; 'G02' 'CT' 'MR_T2' ; 'G03' 'MR_T1' 'MR_T2' ; 'G04' 'MR_T1' 'MR_T2'}
        xlsFile = strcat('Qc_1_4',filter,'.xlsx') %write the results to this file
    case 'expr 5-8'
        fusionMethodNames = {'CST','MFDF_NSST','NNSST','ST_NSST', 'FMSAP', 'Del_PCA', 'Del_max', 'Del_weighted'}'
        sourceImages = {'G05' 'CT' 'MR' ; 'G06' 'CT' 'MR' ; 'G07' 'CT' 'MR' ; 'G08' 'CT' 'MR' }
        xlsFile = strcat('Qc_5_8',filter,'.xlsx') %write the results to this file
    case 'expr 1-12'
        %fusionMethodNames = {'FSD', 'GP','DWT','RP', 'MDP', 'PCA', 'LP', 'SIDWT', 'Del_PCA', 'Del_max', 'Del_weighted'}'
        fusionMethodNames = {'FSD', 'GP','DWT','RP', 'MDP', 'LP', 'SIDWT', 'Del_PCA', 'Del_max', 'Del_weighted'}'
        sourceImages = {'G01' 'CT' 'MR_T2' ; 'G02' 'CT' 'MR_T2' ; 'G03' 'MR_T1' 'MR_T2' ; 'G04' 'MR_T1' 'MR_T2' ; ...
            'G05' 'CT' 'MR' ; 'G06' 'CT' 'MR' ; 'G07' 'CT' 'MR' ; 'G08' 'CT' 'MR' ; ...
            'G09' 'CT' 'MR_T2' ; 'G10' 'MR_PD' 'MR_T2' ; 'G11' 'CT' 'MR_GAD' ; 'G12' 'CT' 'MR_T1'}
        xlsFile = strcat('Qc_1_12',filter,'.xlsx') %write the results to this file
        
end
%% initialization
fileName = split(xlsFile,'.xls'); fileName = char(fileName(1));
%%
fn_QC(metricNames,fusionMethodNames,sourceImages, xlsFile);
fn_ANOVA_Friedman(fileName, metricNames,fusionMethodNames,sourceImages)
fn_Friedman_PostHocs(fileName)
fn_metricsCorrelations(fileName, metricNames,fusionMethodNames,sourceImages)
%%
disp('-------------------');
disp(['the Quantitative Comparisons are written in this excel file: ', xlsFile]);
disp(['the ANOVA results are written in this excel file: ', strcat(fileName,'_ANOVA.xlsx')]);
disp(['the Friedman results are written in this excel file: ', strcat(fileName,'_Friedman.xlsx')]);
