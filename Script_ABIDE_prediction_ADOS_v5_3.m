% Compare static-FC and/or FC-variability to predict ADOS score

% Dependencies:
% MVPA-Light: https://github.com/treder/MVPA-Light
% LIBSVM: https://github.com/cjlin1/libsvm
% Permutation ANOVA: https://www.mathworks.com/matlabcentral/fileexchange/44307-randanova1?s_tid=srchtitle
% ShadedErrorBar: https://www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar

clear all
close all
clc

%% Set Environment
% Set the number of cores to use
Core_max = feature('numcores')
% maxNumCompThreads(Core_max-3);
% parpool('local',Core_max-4);

% add path
addpath(pwd);

% set path
cd('../../..');
Dirworking = pwd;
Dirlog = [Dirworking filesep 'log'];
Dirdocu1 = [Dirworking filesep 'documents' filesep 'abide1']; % folder whose ABIDE1 document
Dirdocu2 = [Dirworking filesep 'documents' filesep 'abide2']; % folder whose ABIDE2 document
Dirmask = [Dirworking filesep 'mask'];
Dirdata = [Dirworking filesep 'analysis' filesep 'abide' filesep 'data']; % where extracted time-series (txt files) located
Dirresult = [Dirworking filesep 'analysis' filesep 'MSSD_classification'];

%% Parameters Settings
cd(Dirworking);

% Set parameters
option.method = 'svm'; %'kernel_fda', 'svm'
option.kernel = 'linear';
option.preprocess = {'oversample'};
option.statmetric = 'accuracy';

option.windowsize = 20;         % current TR parameter (20/30/40 in ABIDE dataset)
option.slidingwindow = 'robust';  % current sliding window type
option.signaltype = 'Sig';  % Raw Dfc signals
option.metricindex = 4;     % 1 is Raw Dfc, 2 is mean, 3 is SD, 4 is MSSD,...


option.atlas = 4;   % 1: aal2bilateral / 2: aal2unilateral / 3: LAIRD / 4: Schaefer200 / 5: Gordon

switch option.atlas
    case 1
        option.allROI = 1:64; % bilateral AAL2 atals number
        option.outROI = [48:64]; % cerebellum areas
        Idx_data = 8;
        Label_atlas = 'aal2bilateral';
    case 2
        option.allROI = 1:120;
        option.outROI = 95:120;
        Idx_data = 7;
        Label_atlas = 'aal2unilateral';
    case 3
        option.allROI = 1:20;
        option.outROI = 19:20; % noise networks
        Idx_data = 9;
        Label_atlas = 'LAIRD';
    case 4
        option.allROI = 1:200;
        option.outROI = [];
        Label_atlas = 'Schaefer200'; % In Schaefer200, no cerebellum
end
option.incROI = setdiff(option.allROI, option.outROI);
option.nROI = length(option.incROI);

option.lassoiter = 5;

% Downsampling rate
TR_down = 3;

%% Make an Index matrix (lower triangle of connect matrix)
% This will be used to find where are significant features
TMP_mat1 = 1:(option.nROI)^2;
TMP_mat1 = reshape(TMP_mat1, option.nROI,option.nROI);
TMP_mat1 = tril(TMP_mat1,-1);
TMP_vec1 = TMP_mat1(TMP_mat1~=0);
for nComp = 1:length(TMP_vec1), Idxref_tril(nComp) = find(TMP_mat1 == TMP_vec1(nComp)); end

TMP_mat2 = TMP_mat1 + (option.nROI)^2;
TMP_mat2 = reshape(TMP_mat2, option.nROI,option.nROI);
TMP_mat2 = tril(TMP_mat2,-1);
TMP_mat2 = [TMP_mat1 TMP_mat2];
TMP_vec2 = TMP_mat2(TMP_mat2~=0);
for nComp = 1:length(TMP_vec2), Idxref_tril_combine(nComp) = find(TMP_mat2 == TMP_vec2(nComp)); end

%% Get subjects Information
cd(Dirresult)
if exist('BasicInfo_ABIDE.mat', 'file') == 2
    load BasicInfo_ABIDE;
    
else

    Script_BasicInfo_ABIDE12; % Merge documents of ABIDE1 & ABIDE2

    cd(Dirresult);
    save BasicInfo_ABIDE12 BasicInfo

end

DeleteList = isnan(BasicInfo.ADOS_Total);

%% Select Subjects for the further analysis
Idx_pat = find(DeleteList == 0 & BasicInfo.ADOS_Total ~= 0); % Only patients have ADOS scores

%% Load data
mkdir(Dirresult); cd(Dirresult);
Label_rawdata = sprintf('Sig_%s_allsubj',Label_atlas);

if exist([Label_rawdata, '.mat'], 'file') == 2
    load(Label_rawdata);
else   
    
    
    List_subj = BasicInfo.subID;
    
    % Load ABIDE ===============
    cd(Dirdata);
    for nSubj = 1:size(BasicInfo,1)
        % example: load Schaefer200 time series
        cd([List_subj{nSubj}]);
        List_txt = dir('ts_Schaefer200*');
        TMP_signal = load(List_txt.name); % each column is each brain region

        Flip_TMP_signal = [flip(TMP_signal); TMP_signal; flip(TMP_signal)];

        % Downsampling to TR = 3sec
        TR_orig = BasicInfo.TRSec(nSubj);
        L = ceil(size(TMP_signal,1)*(floor(TR_orig*1000))/(floor(TR_down*1000)));
        TMP_signal = resample(Flip_TMP_signal,floor(TR_orig*1000),floor(TR_down*1000));
        signal = TMP_signal([L:2*L-1],:);

        Sig_allsubj{nSubj,1} = signal;

        clear LAIRD_signal TMP_LAIRD_signal
        cd(Dirdata);
    end
    
   
    % Merge & Save =========================
    cd(Dirresult);
    Filename = ['Sig_' Label_atlas '_allsubj.mat'];
    eval(['save ',Filename,' Sig_allsubj']);
       
       
end

%% Calculate Dynamic connectivty using sliding window (varying window size)
if strcmp(option.slidingwindow,'robust')
    Name_fol = sprintf('Connectivity_%s_%s_%dTR',Label_atlas,option.slidingwindow,option.windowsize);
end
mkdir(Name_fol); cd(Name_fol);

if exist('Dfc_metrics', 'var') == 1
    fprintf('Load Dfc_metrics success \n');
elseif exist('Dfc_metrics', 'var') == 0

    if exist('Dfc_metrics.mat', 'file') ~= 2
        
    windowsize = option.windowsize;

    parfor nSubj = 1:size(Sig_allsubj,1)  % When you use 'DCC' change to for instead of parfor

    Name_file = strcat('Connectivity_', BasicInfo.subID(nSubj), '.mat');
    Name_file = Name_file{1};

    if exist(Name_file, 'file') ~= 2
        if strcmp(option.slidingwindow,'robust')
            [TMP_Dfc] = Script_Dfc(Sig_allsubj{nSubj,1},windowsize,4);
        end
        parsave(Name_file, TMP_Dfc);
    end
        TMP_Dfc = [];
    end

    end
    
    %% Calculate Metrics
    parfor nSubj = 1:size(Sig_allsubj,1)

    Name_file = strcat('Connectivity_', BasicInfo.subID(nSubj), '.mat');
    Name_file = Name_file{1};

    Mat_connect = parload(Name_file);

    [Mat_MEAN Mat_SD Mat_MSSD Mat_VSD] = Script_calMSSD(Mat_connect);

    Mat_MEAN = tril(Mat_MEAN,-1);
    Mat_SD = tril(Mat_SD,-1);
    Mat_MSSD = tril(Mat_MSSD,-1);
    Mat_VSD = tril(Mat_VSD,-1);

    Dfc_metrics1{nSubj,1} = Mat_connect;
    Dfc_metrics2{nSubj,1} = Mat_MEAN(Mat_MEAN~=0);
    Dfc_metrics3{nSubj,1} = Mat_SD(Mat_SD~=0);
    Dfc_metrics4{nSubj,1} = Mat_MSSD(Mat_MSSD~=0);
    Dfc_metrics5{nSubj,1} = Mat_VSD(Mat_VSD~=0);

    Mat_MEAN=[]; Mat_SD=[]; Mat_MSSD=[]; Mat_VSD=[];


    end

    % Summarize to one variable
    for nIdx = 1:5
        eval(['Dfc_metrics(:,' num2str(nIdx) ') = Dfc_metrics' num2str(nIdx) '(:,1);']);
    end

    clear Dfc_metrics1 Dfc_metrics2 Dfc_metrics3 Dfc_metrics4 Dfc_metrics5

    save Dfc_metrics Dfc_metrics

end

%% Calculating static connectivity strength
if exist('fc_static.mat','file') == 2
    load('fc_static.mat');
else

    if strcmp(option.slidingwindow,'robust')
        parfor nSubj = 1:size(Sig_allsubj,1)
            TMP_sig = Sig_allsubj{nSubj,1};
            [T,p] = size(TMP_sig);
            TMP_corr = NaN(p,p);

            for j = 1:size(TMP_sig,2)
                for k = 1:size(TMP_sig,2)
                    data1 = TMP_sig(:,j); data2 = TMP_sig(:,k);
                    [~,~,rob_corrw] = andlab_robustfit( data1, data2);
                    TMP_corr(j,k) = atanh(rob_corrw);
                end
            end

            TMP_corr = tril(TMP_corr,-1);
            fc_static{nSubj,1} = TMP_corr(TMP_corr~=0);
            TMP_sig = [];        
        end        
        save fc_static fc_static;
    end

end

%% Reshape data matrix
parfor nSubj = 1:size(fc_static,1)
    Data_fc = fc_static{nSubj};
    Mat_fc_concat(nSubj,:) = zscore(Data_fc);

    Data_metric = Dfc_metrics{nSubj,option.metricindex};
    Mat_Dfc_concat(nSubj,:) = zscore(Data_metric);

    Mat_combine(nSubj,:) = cat(1,zscore(Data_fc),zscore(Data_metric));
    Data_fc = []; Data_metric = [];
end
fprintf('Reshape is finished \n');

Mat_fc_whole = Mat_fc_concat(Idx_pat,:);
Mat_Dfc_whole = Mat_Dfc_concat(Idx_pat,:);
Mat_combine_whole = Mat_combine(Idx_pat,:);

% Select a specific score (ADOS)
TMP_Score = BasicInfo.ADOS_Total;   % ADOS_Total, SRS_RawTotal
Array_Score = TMP_Score(Idx_pat,:);
clear TMP_Score*

%% Dimension Reduction
% LASSO regression
Fol_save = sprintf('Lasso_%diter_based_n%d_%s_predict_ADOS',option.lassoiter,numel(Idx_pat),option.slidingwindow);
mkdir(Fol_save); cd(Fol_save);

if exist('Cell_B_fc.mat') > 1 && exist('Cell_B_Dfc.mat') > 1 && exist('Cell_B_combine.mat') > 1
    load Cell_B_fc; load Cell_B_Dfc; load Cell_B_combine;
    load Cell_FitInfo_fc; load Cell_FitInfo_Dfc; load Cell_FitInfo_combine;
    load Devi_fc; load Devi_Dfc; load Devi_combine;
else
    for nIter = 1:option.lassoiter
        nIter

        % Lasso regression on whole data (all subject not limited to age group)
        [B_fc,FitInfo_fc] = lassoglm(Mat_fc_whole,Array_Score,'normal','CV',5,'NumLambda',200,'Options',statset('UseParallel',true));
        [B_Dfc,FitInfo_Dfc] = lassoglm(Mat_Dfc_whole,Array_Score,'normal','CV',5,'NumLambda',200,'Options',statset('UseParallel',true));
        [B_combine,FitInfo_combine] = lassoglm(Mat_combine_whole,Array_Score,'normal','CV',5,'NumLambda',200,'Options',statset('UseParallel',true)); % In case of Elasticnet, add: 'Alpha',0.5
        Cell_B_fc{nIter} = B_fc; Cell_FitInfo_fc{nIter} = FitInfo_fc;
        Cell_B_Dfc{nIter} = B_Dfc; Cell_FitInfo_Dfc{nIter} = FitInfo_Dfc;
        Cell_B_combine{nIter} = B_combine; Cell_FitInfo_combine{nIter} = FitInfo_combine;

        Devi_fc(nIter,:) = FitInfo_fc.Deviance;
        Devi_Dfc(nIter,:) = FitInfo_Dfc.Deviance;
        Devi_combine(nIter,:) = FitInfo_combine.Deviance;

        clear B_fc B_Dfc B_combine FitInfo_fc FitInfo_Dfc FitInfo_combine
    end
    
    save Cell_B_fc Cell_B_fc; save Cell_FitInfo_fc Cell_FitInfo_fc;
    save Cell_B_Dfc Cell_B_Dfc; save Cell_FitInfo_Dfc Cell_FitInfo_Dfc;
    save Cell_B_combine Cell_B_combine; save Cell_FitInfo_combine Cell_FitInfo_combine;
    save Devi_fc Devi_fc; save Devi_Dfc Devi_Dfc; save Devi_combine Devi_combine;
end
  
Ave_Devi_fc = mean(Devi_fc,1);
Ave_Devi_Dfc = mean(Devi_Dfc,1);
Ave_Devi_combine = mean(Devi_combine,1);

Idx_min_Devi_fc = find(Ave_Devi_fc == min(Ave_Devi_fc));
Idx_min_Devi_Dfc = find(Ave_Devi_Dfc == min(Ave_Devi_Dfc));
Idx_min_Devi_combine = find(Ave_Devi_combine == min(Ave_Devi_combine));

for nGrandIter = 1:2
    
option.featuretype = nGrandIter; % 1: for a specific lambda / 2: smallest lambda
clearvars Opt_B_* Predictor_* Idxref_Predictor_* Cell_MSE_* Cell_svr_* TXT_* Area_wave_* Vec_MSE Idx_type Summary_perm Mean_MSE_* Std_MSE_*


% Select the optimal Lambda
if option.featuretype == 1 % optimal lambda
    
    for nIter = 1:size(Cell_B_fc,2)
       Opt_B_fc(:,nIter) = Cell_B_fc{nIter}(:,Idx_min_Devi_fc);
       Opt_B_Dfc(:,nIter) = Cell_B_Dfc{nIter}(:,Idx_min_Devi_Dfc); 
       Opt_B_combine(:,nIter) = Cell_B_combine{nIter}(:,Idx_min_Devi_combine);          
    end
    
elseif option.featuretype == 2 % liberal lambda
    for nIter = 1:size(Cell_B_fc,2)
       Opt_B_fc(:,nIter) = Cell_B_fc{nIter}(:,1);
       Opt_B_Dfc(:,nIter) = Cell_B_Dfc{nIter}(:,1);
       Opt_B_combine(:,nIter) = Cell_B_combine{nIter}(:,1);
    end    
   
end

Opt_B_fc = mean(Opt_B_fc,2);
Opt_B_Dfc = mean(Opt_B_Dfc,2);
Opt_B_combine = mean(Opt_B_combine,2);

Predictor_fc = find(Opt_B_fc(:,1)~=0);
Predictor_Dfc = find(Opt_B_Dfc(:,1)~=0);
Predictor_combine = find(Opt_B_combine(:,1)~=0);
    
% Find the locations of the selected feature within the original
% connectivity matrix (n X n)
TMP_Idxref_Predictor_fc = Idxref_tril(Predictor_fc);
for nComp = 1:length(TMP_Idxref_Predictor_fc), Idxref_Predictor_fc(nComp) = find(TMP_mat1 == TMP_Idxref_Predictor_fc(nComp)); end

TMP_Idxref_Predictor_Dfc = Idxref_tril(Predictor_Dfc);
for nComp = 1:length(TMP_Idxref_Predictor_Dfc), Idxref_Predictor_Dfc(nComp) = find(TMP_mat1 == TMP_Idxref_Predictor_Dfc(nComp)); end

TMP_Idxref_Predictor_combine = Idxref_tril_combine(Predictor_combine);
for nComp = 1:length(TMP_Idxref_Predictor_combine), Idxref_Predictor_combine(nComp) = find(TMP_mat2 == TMP_Idxref_Predictor_combine(nComp)); end

clear TMP_Idxref*

%% Machine Learning
% SVR Parameter
cfg = [];
cfg = mv_get_hyperparameter('libsvm');
cfg.svm_type   = 3;
cfg.kernel     = 'linear';
num_cv         = 5;
num_iter       = 100; % for statistical testing  
    
if option.featuretype == 1
    n_feature = 1;    
elseif option.featuretype == 2
    n_feature = min([size(Predictor_fc,1) size(Predictor_Dfc,1) size(Predictor_combine,1)]);
end

for ifeature = 1:n_feature

    fprintf('%d/%d feature processing... \n',ifeature,n_feature)
    
    if option.featuretype == 2
        [~,Predictor_fc] = maxk(abs(Opt_B_fc(:,1)), ifeature);
        [~,Predictor_Dfc] = maxk(abs(Opt_B_Dfc(:,1)), ifeature);
        [~,Predictor_combine] = maxk(abs(Opt_B_combine(:,1)), ifeature);
    end
    % Make sparse matrices
    i_Mat_fc_whole = Mat_fc_whole(:,Predictor_fc);
    i_Mat_Dfc_whole = Mat_Dfc_whole(:,Predictor_Dfc);
    i_Mat_combine_whole = Mat_combine_whole(:,Predictor_combine);


    %% SVR run 
    parfor nIter = 1:num_iter
    c = cvpartition(Array_Score, 'KFold',num_cv);
    
    for nCV = 1:num_cv
        % Partitioning for Cross-validation
        Idx_train = training(c,nCV);
        Idx_test = test(c,nCV);
        
        Score_true = Array_Score(Idx_test,:);
        model_fc = train_libsvm(cfg, i_Mat_fc_whole(Idx_train,:), Array_Score(Idx_train));
        [Pred_fc] = test_libsvm(model_fc, i_Mat_fc_whole(Idx_test,:));
        MSE_fc(nCV,1,nIter) = mv_calculate_performance('mean_squared_error', [], Pred_fc, Score_true);
        w_fc(nCV,:,nIter) = (model_fc.model.SVs' * model_fc.model.sv_coef)'; 
        b_fc(nCV,1,nIter) = -model_fc.model.rho;

        model_Dfc = train_libsvm(cfg, i_Mat_Dfc_whole(Idx_train,:), Array_Score(Idx_train));
        [Pred_Dfc] = test_libsvm(model_Dfc, i_Mat_Dfc_whole(Idx_test,:));
        MSE_Dfc(nCV,1,nIter) = mv_calculate_performance('mean_squared_error', [], Pred_Dfc, Score_true);
        w_Dfc(nCV,:,nIter) = (model_Dfc.model.SVs' * model_Dfc.model.sv_coef)';
        b_Dfc(nCV,1,nIter) = -model_Dfc.model.rho; 


        model_combine = train_libsvm(cfg, i_Mat_combine_whole(Idx_train,:), Array_Score(Idx_train));
        [Pred_combine] = test_libsvm(model_combine, i_Mat_combine_whole(Idx_test,:));
        MSE_combine(nCV,1,nIter) = mv_calculate_performance('mean_squared_error', [], Pred_combine, Score_true);
        w_combine(nCV,:,nIter) = (model_combine.model.SVs' * model_combine.model.sv_coef)';
        b_combine(nCV,1,nIter) = -model_combine.model.rho;        
        
    end
    Idx_train = []; Idx_test = []; Score_true = []; c = [];
           
    end
     
    % resolve error when no feature selected. 10 is a just arbitrary number
    if numel(Predictor_fc) == 0, w_fc = zeros(num_cv,10,num_iter); end
    if numel(Predictor_Dfc) == 0, w_Dfc = zeros(num_cv,10,num_iter); end
    if numel(Predictor_combine) == 0, w_combine = zeros(num_cv,10,num_iter); end
 
    % Save
    Cell_MSE_fc{1,ifeature} = squeeze(mean(MSE_fc,1)); 
    Cell_MSE_Dfc{1,ifeature} = squeeze(mean(MSE_Dfc,1)); 
    Cell_MSE_combine{1,ifeature} = squeeze(mean(MSE_combine,1)); 
    Cell_svr_w_fc{1,ifeature} = squeeze(mean(w_fc,1))';
    Cell_svr_w_Dfc{1,ifeature} = squeeze(mean(w_Dfc,1))';
    Cell_svr_w_combine{1,ifeature} = squeeze(mean(w_combine,1))';
    Cell_svr_b_fc{1,ifeature} = squeeze(mean(b_fc,1));
    Cell_svr_b_Dfc{1,ifeature} = squeeze(mean(b_Dfc,1));
    Cell_svr_b_combine{1,ifeature} = squeeze(mean(b_combine,1));

    clear i_Mat_fc_whole i_Mat_Dfc_whole i_Mat_combine_whole MSE_* w_* b_*
end

%% Save SVR results
if ifeature == 1
    Dir_name = sprintf('SVR_setlambda_lasso_n%d_%s_%diter_predict_ADOS_v5',numel(Idx_pat),option.slidingwindow,num_iter);
else
    Dir_name = sprintf('SVR_feature%d_lasso_n%d_%s_%diter_predict_ADOS_v5',ifeature,numel(Idx_pat),option.slidingwindow,num_iter);
end

mkdir(Dir_name); cd(Dir_name);
save Cell_MSE_fc Cell_MSE_fc
save Cell_MSE_Dfc Cell_MSE_Dfc
save Cell_MSE_combine Cell_MSE_combine

save Cell_svr_b_fc Cell_svr_b_fc
save Cell_svr_b_Dfc Cell_svr_b_Dfc
save Cell_svr_b_combine Cell_svr_b_combine
save Cell_svr_w_fc Cell_svr_w_fc
save Cell_svr_w_Dfc Cell_svr_w_Dfc
save Cell_svr_w_combine Cell_svr_w_combine

% Save as text
for nIter = 1:size(Cell_MSE_fc,2)
    TXT_MSE_fc(:,nIter) = Cell_MSE_fc{nIter};
    TXT_MSE_Dfc(:,nIter) = Cell_MSE_Dfc{nIter};
    TXT_MSE_combine(:,nIter) = Cell_MSE_combine{nIter};
end

csvwrite('TXT_MSE_fc.txt', TXT_MSE_fc);
csvwrite('TXT_MSE_Dfc.txt', TXT_MSE_Dfc);
csvwrite('TXT_MSE_combine.txt', TXT_MSE_combine);
csvwrite('TXT_MSE_fc.csv', TXT_MSE_fc);
csvwrite('TXT_MSE_Dfc.csv', TXT_MSE_Dfc);
csvwrite('TXT_MSE_combine.csv', TXT_MSE_combine);

if option.featuretype == 1
   Num_selected_feature = [length(Predictor_fc) length(Predictor_Dfc) length(Predictor_combine)]; 
elseif option.featuretype == 2
    Num_selected_feature = n_feature;
end
csvwrite('Num_selected_feature.txt', Num_selected_feature);

if option.featuretype == 2
    Value_maximum = max([max(max(TXT_MSE_fc)) max(max(TXT_MSE_Dfc)) max(max(TXT_MSE_combine))]); % To calculate Area under curve, we need normalization
    TXT_MSE_fc = TXT_MSE_fc./Value_maximum;
    TXT_MSE_Dfc = TXT_MSE_Dfc./Value_maximum;
    TXT_MSE_combine = TXT_MSE_combine./Value_maximum;
    
    Area_wave_fc = sum(TXT_MSE_fc,2)./size(TXT_MSE_fc,2);
    Area_wave_Dfc = sum(TXT_MSE_Dfc,2)./size(TXT_MSE_Dfc,2);
    Area_wave_combine = sum(TXT_MSE_combine,2)./size(TXT_MSE_combine,2);
    
    csvwrite('TXT_Area_wave_fc.txt',Area_wave_fc);
    csvwrite('TXT_Area_wave_fc.csv',Area_wave_fc);
    csvwrite('TXT_Area_wave_Dfc.txt',Area_wave_Dfc);
    csvwrite('TXT_Area_wave_Dfc.csv',Area_wave_Dfc);
    csvwrite('TXT_Area_wave_combine.txt',Area_wave_combine);
    csvwrite('TXT_Area_wave_combine.csv',Area_wave_combine);
end

%% Permutation Testing on MSE
Test_Perm = {'pval'; 'Factual';'Fdist'};

if option.featuretype == 1
    % Overall Permutation Test
    Vec_MSE = [Cell_MSE_fc{1}' Cell_MSE_Dfc{1}' Cell_MSE_combine{1}'];
    Idx_type = [ones(1,num_iter) ones(1,num_iter)*2 ones(1,num_iter)*3];
    [Test_Perm{1,2},Test_Perm{2,2},Test_Perm{3,2}] = randanova1(Vec_MSE,Idx_type,5000);   % 5000 iterations; outputs are pval / Factual / Fdist
    
    % Testing static-FC & MSSD (1 & 2)
    Vec_MSE = [Cell_MSE_fc{1}' Cell_MSE_Dfc{1}'];
    Idx_type = [ones(1,num_iter) ones(1,num_iter)*2];
    [Test_Perm{1,3},Test_Perm{2,3},Test_Perm{3,3}] = randanova1(Vec_MSE,Idx_type,5000);   % 5000 iterations
    
     % Testing static-FC & Combined (1 & 3)
    Vec_MSE = [Cell_MSE_fc{1}' Cell_MSE_combine{1}'];
    Idx_type = [ones(1,num_iter) ones(1,num_iter)*2];
    [Test_Perm{1,4},Test_Perm{2,4},Test_Perm{3,4}] = randanova1(Vec_MSE,Idx_type,5000);   % 5000 iterations
    
     % Testing MSSD * Combined (2 & 3)
    Vec_MSE = [Cell_MSE_Dfc{1}' Cell_MSE_combine{1}'];
    Idx_type = [ones(1,num_iter) ones(1,num_iter)*2];
    [Test_Perm{1,5},Test_Perm{2,5},Test_Perm{3,5}] = randanova1(Vec_MSE,Idx_type,5000);   % 5000 iterations
   
    Summary_Perm = cell2table(Test_Perm);
    Summary_Perm.Properties.VariableNames = {'Stat','Overall','static_MSSD','static_combine','MSSD_combine'};
    save Summary_Perm Summary_Perm
    
    Summary_Perm(3,:) = []; % Remove Fdist for saving a subsequent txt file
    writetable(Summary_Perm,'Summary_Perm.csv');
    
elseif option.featuretype == 2
    % Overall Permutation Test
    Vec_MSE = [Area_wave_fc' Area_wave_Dfc' Area_wave_combine'];
    Idx_type = [ones(1,num_iter) ones(1,num_iter)*2 ones(1,num_iter)*3];
    [Test_Perm{1,2},Test_Perm{2,2},Test_Perm{3,2}] = randanova1(Vec_MSE,Idx_type,5000);   % 5000 iterations; outputs are pval / Factual / Fdist
    
    % Testing static-FC & MSSD (1 & 2)
    Vec_MSE = [Area_wave_fc' Area_wave_Dfc'];
    Idx_type = [ones(1,num_iter) ones(1,num_iter)*2];
    [Test_Perm{1,3},Test_Perm{2,3},Test_Perm{3,3}] = randanova1(Vec_MSE,Idx_type,5000);   % 5000 iterations
    
     % Testing static-FC & Combined (1 & 3)
    Vec_MSE = [Area_wave_fc' Area_wave_combine'];
    Idx_type = [ones(1,num_iter) ones(1,num_iter)*2];
    [Test_Perm{1,4},Test_Perm{2,4},Test_Perm{3,4}] = randanova1(Vec_MSE,Idx_type,5000);   % 5000 iterations
    
     % Testing MSSD * Combined (2 & 3)
    Vec_MSE = [Area_wave_Dfc' Area_wave_combine'];
    Idx_type = [ones(1,num_iter) ones(1,num_iter)*2];
    [Test_Perm{1,5},Test_Perm{2,5},Test_Perm{3,5}] = randanova1(Vec_MSE,Idx_type,5000);   % 5000 iterations
   
    Summary_Perm = cell2table(Test_Perm);
    Summary_Perm.Properties.VariableNames = {'Stat','Overall','static_MSSD','static_combine','MSSD_combine'};
    save Summary_Perm Summary_Perm
    
    Summary_Perm(3,:) = []; % Remove Fdist for saving a subsequent txt file
    writetable(Summary_Perm,'Summary_Perm.csv');
            
end


%% Plotting a matrix for which connectivity (or MSSD) features influence the classification
% Summarize mean and std
for nIter = 1:n_feature
    Mean_MSE_fc(1,nIter) = mean(Cell_MSE_fc{nIter});
    Std_MSE_fc(1,nIter) = std(Cell_MSE_fc{nIter});
    Mean_MSE_Dfc(1,nIter) = mean(Cell_MSE_Dfc{nIter});
    Std_MSE_Dfc(1,nIter) = std(Cell_MSE_Dfc{nIter});
    Mean_MSE_combine(1,nIter) = mean(Cell_MSE_combine{nIter});
    Std_MSE_combine(1,nIter) = std(Cell_MSE_combine{nIter});
end

if option.featuretype == 1
    %% Bar graph
    Values_mean = [Mean_MSE_fc; Mean_MSE_Dfc; Mean_MSE_combine];
    Values_std = [Std_MSE_fc; Std_MSE_Dfc; Std_MSE_combine];
    
    b = bar(1:3, Values_mean, 'FaceColor','Flat'); 
    b.CData(1,:) = [0.7 0.7 0.7];
    b.CData(2,:) = [0 0 0.8];
    b.CData(3,:) = [0 0.8 0];
    hold on;
    er = errorbar(1:3, Values_mean, Values_std);
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    set(gca,'XTick', 1:3, 'XTickLabels', {'Static','MSSD','Combine'})
    ylabel('MSE');
    saveas(gcf,[Dir_name '.png']);
    close all
    
    %% Weights matrix
    % Calculate the average of weights
    Ave_svr_w_fc = mean(Cell_svr_w_fc{1},1);
    Ave_svr_w_Dfc = mean(Cell_svr_w_Dfc{1},1);
    Ave_svr_w_combine = mean(Cell_svr_w_combine{1},1);

    % Preparing plotting
    clim = [-0.5, 0.5];  % colorbar range
    Range_x = [-1 1];
    Range_y = [-1 1];
    
    for nPlot = 1:2
        if nPlot == 1
            dx = diff(Range_x)/(option.nROI-1);
            xg = linspace(Range_x(1)-dx/2,Range_x(2)+dx/2, option.nROI+1);
        elseif nPlot == 2
            dx = diff(Range_x)/(2*option.nROI-1);
            xg = linspace(Range_x(1)-dx/2,Range_x(2)+dx/2, 2*option.nROI+1);
        end
        dy = diff(Range_y)/(option.nROI-1);
        yg = linspace(Range_y(1)-dy/2,Range_y(2)+dy/2, option.nROI+1);

        if nPlot == 1
            Plot_svr_w_fc = zeros(option.nROI);
            for nFeature = 1:length(Ave_svr_w_fc), Plot_svr_w_fc(Idxref_Predictor_fc(nFeature)) = Ave_svr_w_fc(nFeature); end
            figure(1); Plot_color = imagesc(Range_x,Range_y,Plot_svr_w_fc); axis off; hold on
            Plot_grid = mesh(xg,yg,zeros([option.nROI+1, option.nROI+1]))
            Plot_grid.FaceColor = 'none'; Plot_grid.EdgeColor = 'k';
            set(gca, 'Ydir','reverse','clim',clim); colormap jet; pbaspect([1 1 1]);
            saveas(gcf,'Plot_svr_w_static.png');
            close all
            save Plot_svr_w_fc Plot_svr_w_fc

            Plot_svr_w_Dfc = zeros(option.nROI);
            for nFeature = 1:length(Ave_svr_w_Dfc), Plot_svr_w_Dfc(Idxref_Predictor_Dfc(nFeature)) = Ave_svr_w_Dfc(nFeature); end
            figure(1); Plot_color = imagesc(Range_x,Range_y,Plot_svr_w_Dfc); axis off; hold on
            Plot_grid = mesh(xg,yg,zeros([option.nROI+1, option.nROI+1]))
            Plot_grid.FaceColor = 'none'; Plot_grid.EdgeColor = 'k';
            set(gca, 'Ydir','reverse','clim',clim); colormap jet; pbaspect([1 1 1]);
            saveas(gcf,'Plot_svr_w_mssd.png');
            close all
            save Plot_svr_w_Dfc Plot_svr_w_Dfc

        elseif nPlot == 2
            Plot_svr_w_combine = zeros(option.nROI, option.nROI*2);
            for nFeature = 1:length(Ave_svr_w_combine), Plot_svr_w_combine(Idxref_Predictor_combine(nFeature)) = Ave_svr_w_combine(nFeature); end
            figure(1); Plot_color = imagesc(Range_x,Range_y,Plot_svr_w_combine); axis off; hold on
            Plot_grid = mesh(xg,yg,zeros([option.nROI+1, 2*option.nROI+1]))
            Plot_grid.FaceColor = 'none'; Plot_grid.EdgeColor = 'k';
            set(gca, 'Ydir','reverse','clim',[-0.5,0.5]); colormap jet; pbaspect([2 1 1]);
            saveas(gcf,'Plot_svr_w_combine.png');
            close all
            save Plot_svr_w_combine Plot_svr_w_combine
        end
    end
    
elseif option.featuretype == 2  
    %% Bar graph
    Values_mean = [mean(Area_wave_fc); mean(Area_wave_Dfc); mean(Area_wave_combine)];
    Values_std = [std(Area_wave_fc); std(Area_wave_Dfc); std(Area_wave_combine)];
    
    b = bar(1:3, Values_mean, 'FaceColor','Flat'); 
    b.CData(1,:) = [0.7 0.7 0.7];
    b.CData(2,:) = [0 0 0.8];
    b.CData(3,:) = [0 0.8 0];
    hold on;
    er = errorbar(1:3, Values_mean, Values_std);
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    set(gca,'XTick', 1:3, 'XTickLabels', {'Static','MSSD','Combine'})
    ylabel('AUC');
    saveas(gcf,'Bargraph_Area_wave.png');
    close all
    
    %% Wave graph
    % Plotting wave graph
    shadedErrorBar(1:size(Mean_MSE_fc,2),Mean_MSE_fc,Std_MSE_fc,{'k','markerfacecolor','r'}); hold on
    shadedErrorBar(1:size(Mean_MSE_Dfc,2),Mean_MSE_Dfc,Std_MSE_Dfc,{'b','markerfacecolor','r'}); hold on
    shadedErrorBar(1:size(Mean_MSE_combine,2),Mean_MSE_combine,Std_MSE_combine,{'g','markerfacecolor','r'}); 
    legend({'Static','MSSD','Combined'},'Location','southeast');
    %set(gca,'Ylim',[0.5 1]);
    ylabel('MSE');
    saveas(gcf,[Dir_name '.png']);
    close all

end

    cd('..');
end

