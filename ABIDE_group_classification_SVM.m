% Compare whole connectivity and/or MSSD to classiciation ASD/CON

% Dependencies:
% MVPA-Light: https://github.com/treder/MVPA-Light
% LIBSVM: https://github.com/cjlin1/libsvm
% Permutation ANOVA: https://www.mathworks.com/matlabcentral/fileexchange/44307-randanova1?s_tid=srchtitle
% ShadedErrorBar: https://www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar

clc;
clear all;
close all;

%% Set Environments
% Set the number of cores to use
Core_max = feature('numcores')
LASTN = maxNumCompThreads(30); % number of cores to use

% add path
addpath(pwd);

% set path
Dirworking = ['project_folder'];
Dirlog = [Dirworking filesep 'log'];
Dirdocu1 = [Dirworking filesep 'documents' filesep 'abide1']; % folder whose ABIDE1 document
Dirdocu2 = [Dirworking filesep 'documents' filesep 'abide2']; % folder whose ABIDE2 document
Dirmask = [Dirworking filesep 'mask'];
Dirdata = [Dirworking filesep 'analysis' filesep 'ts']; % where extracted time-series (txt files) located
Dirresult = [Dirworking filesep 'analysis' filesep 'result'];

%% Parameters settings
cd(Dirworking);

% Set parameters
option.method = 'svm'; %'kernel_fda', 'svm'
option.kernel = 'linear';
option.metric = {'auc'};
option.preprocess = {'oversample'};
option.statmetric = 'accuracy';

option.windowsize = 20;         % current TR parameter (In ABIDE, 20/30/40)
option.slidingwindow = 'sw';  % current sliding window type 'sw'. 'tsw'
option.signaltype = 'Sig';  % Raw Dfc signals
option.metricindex = 4;     % 1 is Raw Dfc, 2 is mean, 3 is SD, 4 is MSSD,...


option.atlas = 4;   % 1: aal2bilateral / 2: aal2unilateral / 3: LAIRD / 4: Schaefer200

switch option.atlas
    case 1
        option.allROI = 1:64; % bilateral AAL2 atals number
        option.outROI = [48:64]; % cerebellum areas
        Label_atlas = 'aal2bilateral';
    case 2
        option.allROI = 1:120;
        option.outROI = 95:120;
        Label_atlas = 'aal2unilateral';
    case 3
        option.allROI = 1:20;
        option.outROI = [14,19:20]; % Exclude cerebellum network
        Label_atlas = 'LAIRD';
    case 4
        option.allROI = 1:200;
        option.outROI = [];
        Label_atlas = 'Schaefer200'; % In Schaefer200, no cerebellum
    case 5
        option.allROI = 1:333;
        option.outROI = [];
        Label_atlas = 'Gordon';
end
option.incROI = setdiff(option.allROI, option.outROI);
option.nROI = length(option.incROI);

option.whichlambda = 2; % 2 = using smallest lambda (i.e., least stringent) at the feature selection stage
option.lassoiter = 5; % number of Elastic-net (or Lasso) iteration for feature selection

% Downsampling rate
TR_down = 3; % Due to inhomogenity of TRs of dataset, downsampled to 0.333 Hz sampling rate

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
if exist('BasicInfo_ABIDE12.mat', 'file') == 2
    load BasicInfo_ABIDE12; load BasicInfo1; load BasicInfo2;
    
else

    Script_BasicInfo_ABIDE12; % Merge documents of ABIDE1 & ABIDE2

    cd(Dirresult);
    save BasicInfo_ABIDE12 BasicInfo

end

%% Select Subjects for the further analysis
Idx_con = find(BasicInfo.Sub_group==1);
Idx_pat = find(BasicInfo.Sub_group==2);

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

        clear signal TMP_signal TR_orig L Flip_TMP_signal
        cd(Dirdata);
    end
    
    % Merge & Save =========================
    cd(Dirresult);
    Filename = ['Sig_' Label_atlas '_allsubj.mat'];
    eval(['save ',Filename,' Sig_allsubj']);
       
end

%% Calculate Dynamic connectivty using sliding window (varying window size)
if strcmp(option.slidingwindow,'sw') || strcmp(option.slidingwindow,'tsw')
    Name_fol = sprintf('Connectivity_%s_%s_%dTR',Label_atlas,option.slidingwindow,option.windowsize);
end
mkdir(Name_fol); cd(Name_fol);

if exist('Dfc_metrics.mat', 'file') == 2
    load('Dfc_metrics.mat');
else
    windowsize = option.windowsize;

    parfor nSubj = 1:size(Sig_allsubj,1)  % When you use 'DCC' change to for instead of parfor

    Name_file = strcat('Connectivity_', BasicInfo.subID(nSubj), '.mat');
    Name_file = Name_file{1};

    if exist(Name_file, 'file') ~= 2
        if strcmp(option.slidingwindow,'sw')
            [TMP_Dfc] = Script_Dfc(Sig_allsubj{nSubj,1},windowsize,1);
        elseif strcmp(option.slidingwindow,'tsw')
            [TMP_Dfc] = Script_Dfc(Sig_allsubj{nSubj,1},windowsize,2);
        elseif strcmp(option.slidingwindow,'DCC')
            [TMP_Dfc] = Script_Dfc(Sig_allsubj{nSubj,1},windowsize,3);
        end
        parsave(Name_file, TMP_Dfc);
    end
        TMP_Dfc = [];
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

    if strcmp(option.slidingwindow,'sw') || strcmp(option.slidingwindow,'tsw') || strcmp(option.slidingwindow,'DCC')
        parfor nSubj = 1:size(Sig_allsubj,1)
            TMP_sig = Sig_allsubj{nSubj,1};
            TMP_corr = tril( atanh(corr(TMP_sig)),-1 );
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

Mat_con = Mat_fc_concat(Idx_con,:);
Mat_pat = Mat_fc_concat(Idx_pat,:);
Mat_fc_whole = cat(1,Mat_con,Mat_pat);

Mat_con = Mat_Dfc_concat(Idx_con,:);
Mat_pat = Mat_Dfc_concat(Idx_pat,:);
Mat_Dfc_whole = cat(1,Mat_con,Mat_pat);

Mat_con = Mat_combine(Idx_con,:);
Mat_pat = Mat_combine(Idx_pat,:);
Mat_combine_whole = cat(1,Mat_con,Mat_pat);

clear Mat_fc_concat Mat_Dfc_concat Mat_combine
Fol_save = sprintf('Elastic_%diter_based_n%d_%s_groupclassify',option.lassoiter,numel(Idx_con)+numel(Idx_pat),option.slidingwindow);
mkdir(Fol_save); cd(Fol_save);

%% ========== Machine learning Parts ==========
% Parameter setting (Using MVPA-Light)
cfg = [];
cfg = mv_get_hyperparameter('svm');
cfg.classifier = option.method;
cfg.preprocess = option.preprocess;
cfg.kernel = option.kernel;
cfg.metric     = option.metric;
num_cv         = 5;
num_iter       = 100;
cfg.design     = [ones(size(Mat_con,1),1); 2*ones(size(Mat_pat,1),1)]; % Group index (1 = Control, 2 = Autism)


%% Feature selection (nested)
% Note: although it looks like feature selection is not embedded in the
% main n-fold cross-validation, it is.
% By using the same 'cvpartition', the calculated Elasticnet results are
% used in the main CV. (i.e., nested n-fold CV)
% This is for utilizing all available parallel cores.

% if there is already calculated the feature selection results, load them.
if exist('Cell_cv.mat') && exist('Cell_B_fc.mat') > 1 && exist('Cell_B_Dfc.mat') > 1 && exist('Cell_B_combine.mat') > 1
    fprintf('\n\n Feature Selection is loaded.... \n');
    load Cell_cv; 
    load Cell_B_fc; load Cell_B_Dfc; load Cell_B_combine;
    load Cell_FitInfo_fc; load Cell_FitInfo_Dfc; load Cell_FitInfo_combine;
else

% ===== doing feature selection.. =====
fprintf('\n\n Feature Selection is on going.... \n');

for tmp_nIter = 1:num_iter
    c = cvpartition(cfg.design, 'KFold',num_cv); 
    Cell_cv{tmp_nIter} = c;
end
save Cell_cv Cell_cv; clear tmp_nIter c

% mkdir('Feature_selection'); cd('Feature_selection');
for nCV = 1:num_cv
    for nLasso = 1:option.lassoiter
        parfor nIter = 1:num_iter
            c = Cell_cv{nIter};
    
            % clear variances
            Mat_fc_train = []; Mat_Dfc_train = []; Mat_combine_train = []; Mat_fc_test = []; Mat_Dfc_test = []; Mat_combine_test = [];
        
            % Partitioning for Cross-validation
            Idx_train = training(c,nCV);
            Idx_test = test(c,nCV);
            Label_train = cfg.design(Idx_train,:); 
            Label_test = cfg.design(Idx_test,:);
        
            % Elastic-net
            Mat_fc_train = Mat_fc_whole(Idx_train,:);
            Mat_Dfc_train = Mat_Dfc_whole(Idx_train,:);
            Mat_combine_train = Mat_combine_whole(Idx_train,:);
            
            Mat_fc_test = Mat_fc_whole(Idx_test,:);
            Mat_Dfc_test = Mat_Dfc_whole(Idx_test,:);
            Mat_combine_test = Mat_combine_whole(Idx_test,:);
            
            fprintf('%dth iter Elastic-net \n',nIter);
            % Lasso regression on whole data (all subject not limited to age group)
            [B_fc,FitInfo_fc] = lassoglm(Mat_fc_train,Label_train-1,'binomial','CV',5,'NumLambda',200,'Options',statset('UseParallel',true),'Alpha',0.5); % '-1' to make binary array (0 or 1). In other words (1 or 2 -> 0 or 1)
            [B_Dfc,FitInfo_Dfc] = lassoglm(Mat_Dfc_train,Label_train-1,'binomial','CV',5,'NumLambda',200,'Options',statset('UseParallel',true),'Alpha',0.5);
            [B_combine,FitInfo_combine] = lassoglm(Mat_combine_train,Label_train-1,'binomial','CV',5,'NumLambda',200,'Options',statset('UseParallel',true),'Alpha',0.5);
            
            Cell_B_fc{nCV,nLasso,nIter} = B_fc; Cell_FitInfo_fc{nCV,nLasso,nIter} = FitInfo_fc;
            Cell_B_Dfc{nCV,nLasso,nIter} = B_Dfc; Cell_FitInfo_Dfc{nCV,nLasso,nIter} = FitInfo_Dfc;
            Cell_B_combine{nCV,nLasso,nIter} = B_combine; Cell_FitInfo_combine{nCV,nLasso,nIter} = FitInfo_combine;
        end
    end    
       
end

% Save Lasso (Elastic-net) results for the current n-fold cross-validation
save Cell_B_fc Cell_B_fc; save Cell_FitInfo_fc Cell_FitInfo_fc;
save Cell_B_Dfc Cell_B_Dfc; save Cell_FitInfo_Dfc Cell_FitInfo_Dfc;
save Cell_B_combine Cell_B_combine; save Cell_FitInfo_combine Cell_FitInfo_combine;

end

%% Start main n-fold cross-validation
for nIter = 1:num_iter
    fprintf('\n\n ========== %dth iteration ========== \n',nIter);
    c = Cell_cv{nIter};

for nCV = 1:num_cv
    clear Mat_fc_train Mat_Dfc_train Mat_combine_train Mat_fc_test Mat_Dfc_test Mat_combine_test
    clear Devi_fc Devi_Dfc Devi_combine Opt_B_* i_Mat_*
    
    fprintf('\n %dth cross-validation \n',nCV);
    % Partitioning for Cross-validation
    Idx_train = training(c,nCV);
    Idx_test = test(c,nCV);
    Label_train = cfg.design(Idx_train,:); 
    Label_test = cfg.design(Idx_test,:);

    %% Feature selection (LASSO regression or Elastic-net)
    Mat_fc_train = Mat_fc_whole(Idx_train,:);
    Mat_Dfc_train = Mat_Dfc_whole(Idx_train,:);
    Mat_combine_train = Mat_combine_whole(Idx_train,:);
    
    Mat_fc_test = Mat_fc_whole(Idx_test,:);
    Mat_Dfc_test = Mat_Dfc_whole(Idx_test,:);
    Mat_combine_test = Mat_combine_whole(Idx_test,:);
    
    for nLasso = 1:option.lassoiter % To run 'parfor' above, these lines were added
        Devi_fc(nLasso,:) = Cell_FitInfo_fc{nCV,nLasso,nIter}.Deviance;
        Devi_Dfc(nLasso,:) = Cell_FitInfo_Dfc{nCV,nLasso,nIter}.Deviance;
        Devi_combine(nLasso,:) = Cell_FitInfo_combine{nCV,nLasso,nIter}.Deviance;
    end
    Ave_Devi_fc = mean(Devi_fc,1);
    Ave_Devi_Dfc = mean(Devi_Dfc,1);
    Ave_Devi_combine = mean(Devi_combine,1);
    clear TMP_Devi_*
    
    Idx_min_Devi_fc = find(Ave_Devi_fc == min(Ave_Devi_fc));
    Idx_min_Devi_Dfc = find(Ave_Devi_Dfc == min(Ave_Devi_Dfc));
    Idx_min_Devi_combine = find(Ave_Devi_combine == min(Ave_Devi_combine));

    % When several lambdas were selected, choose one
    if numel(Idx_min_Devi_fc) > 1, Idx_min_Devi_fc = Idx_min_Devi_fc(1); end
    if numel(Idx_min_Devi_Dfc) > 1, Idx_min_Devi_Dfc = Idx_min_Devi_Dfc(1); end
    if numel(Idx_min_Devi_combine) > 1, Idx_min_Devi_combine = Idx_min_Devi_combine(1); end
    
    % Select the Lambda
    if option.whichlambda == 1 % optimal lambda
        for nLasso = 1:option.lassoiter
           Opt_B_fc(:,nLasso) = Cell_B_fc{nCV,nLasso,nIter}(:,Idx_min_Devi_fc);
           Opt_B_Dfc(:,nLasso) = Cell_B_Dfc{nCV,nLasso,nIter}(:,Idx_min_Devi_Dfc); 
           Opt_B_combine(:,nLasso) = Cell_B_combine{nCV,nLasso,nIter}(:,Idx_min_Devi_combine);          
        end   
    elseif option.whichlambda == 2 % smallest lambda
        for nLasso = 1:option.lassoiter
           Opt_B_fc(:,nLasso) = Cell_B_fc{nCV,nLasso,nIter}(:,1);
           Opt_B_Dfc(:,nLasso) = Cell_B_Dfc{nCV,nLasso,nIter}(:,1);
           Opt_B_combine(:,nLasso) = Cell_B_combine{nCV,nLasso,nIter}(:,1);
        end 
    end
  
    Opt_B_fc = mean(Opt_B_fc,2); 
    Opt_B_Dfc = mean(Opt_B_Dfc,2);
    Opt_B_combine = mean(Opt_B_combine,2);

    Predictor_fc = find(Opt_B_fc(:,1)~=0); % it will use for the optimal lambda condition. In the smallest lambda condition, the Predictor_* variables will calculate again later
    Predictor_Dfc = find(Opt_B_Dfc(:,1)~=0);
    Predictor_combine = find(Opt_B_combine(:,1)~=0);

    % Find the locations of the selected feature within the original
    % connectivity matrix (n X n). (For subsequent SVM/SVR weights matrix outputs)
    TMP_Idxref_Predictor_fc = Idxref_tril(Predictor_fc);
    for nComp = 1:length(TMP_Idxref_Predictor_fc), Idxref_Predictor_fc{nCV,nIter}(nComp) = find(TMP_mat1 == TMP_Idxref_Predictor_fc(nComp)); end
    TMP_Idxref_Predictor_Dfc = Idxref_tril(Predictor_Dfc);
    for nComp = 1:length(TMP_Idxref_Predictor_Dfc), Idxref_Predictor_Dfc{nCV,nIter}(nComp) = find(TMP_mat1 == TMP_Idxref_Predictor_Dfc(nComp)); end
    TMP_Idxref_Predictor_combine = Idxref_tril_combine(Predictor_combine);
    for nComp = 1:length(TMP_Idxref_Predictor_combine), Idxref_Predictor_combine{nCV,nIter}(nComp) = find(TMP_mat2 == TMP_Idxref_Predictor_combine(nComp)); end
    clear TMP_Idxref*
    
    if option.whichlambda == 1 
        n_feature = 1;   
    elseif option.whichlambda == 2 % smallest lambda (i.e., least stringent)
        n_feature = min([size(Predictor_fc,1) size(Predictor_Dfc,1) size(Predictor_combine,1)]);
    end
 
    %% SVM / SVR parts
    for ifeature = 1:n_feature   
  
        if option.whichlambda == 2 % smallest lambda; Put feature with the order of feature selection coefficient
            Predictor_fc = []; Predictor_Dfc = []; Predictor_combine = [];
            [~,Predictor_fc] = maxk(abs(Opt_B_fc(:,1)), ifeature);
            [~,Predictor_Dfc] = maxk(abs(Opt_B_Dfc(:,1)), ifeature);
            [~,Predictor_combine] = maxk(abs(Opt_B_combine(:,1)), ifeature);
        end
     
        % Make sparse matrices
        i_Mat_fc_train{ifeature} = Mat_fc_train(:,Predictor_fc);
        i_Mat_Dfc_train{ifeature} = Mat_Dfc_train(:,Predictor_Dfc);
        i_Mat_combine_train{ifeature} = Mat_combine_train(:,Predictor_combine);
        
        i_Mat_fc_test{ifeature} = Mat_fc_test(:,Predictor_fc);
        i_Mat_Dfc_test{ifeature} = Mat_Dfc_test(:,Predictor_Dfc);
        i_Mat_combine_test{ifeature} = Mat_combine_test(:,Predictor_combine);        
        

        % Stop when all features were null
        if numel(i_Mat_fc_train) == 0 && numel(i_Mat_Dfc_train) == 0 && numel(i_Mat_combine_train) == 0
            break;
        end
    end
    
    % ===== Train & Test =====
    fprintf('Following feature number is processing... up to %d \n',n_feature);
    parfor ifeature = 1:n_feature  
        fprintf('%d ',ifeature)
        
        if numel(Predictor_fc) > 0
            model_fc = train_svm(cfg, i_Mat_fc_train{ifeature}, Label_train);
            [Pred_fc, dval] = test_svm(model_fc, i_Mat_fc_test{ifeature});
            AUC_fc = mv_calculate_performance('auc', 'dval', dval, Label_test); 
            TMP_Confusion_fc = mv_calculate_performance('confusion', 'clabel', Pred_fc, Label_test);
            Sensitivity_fc = TMP_Confusion_fc(1,1) / (TMP_Confusion_fc(1,1) + TMP_Confusion_fc(1,2));
            Specificity_fc = TMP_Confusion_fc(2,2) / (TMP_Confusion_fc(2,1) + TMP_Confusion_fc(2,2));
            w_fc = model_fc.w'; 
            b_fc = model_fc.b;
        elseif numel(Predictor_fc) == 0 % When there is no selected features (almost impossible to occur; Just in case)
            AUC_fc = 0;
            Sensitivity_fc = 0;
            Specificity_fc = 0;
        end

        if numel(Predictor_Dfc) > 0
            model_Dfc = train_svm(cfg, i_Mat_Dfc_train{ifeature}, Label_train);
            [Pred_Dfc, dval] = test_svm(model_Dfc, i_Mat_Dfc_test{ifeature});
            AUC_Dfc = mv_calculate_performance('auc', 'dval', dval, Label_test);
            TMP_Confusion_Dfc = mv_calculate_performance('confusion', 'clabel', Pred_Dfc, Label_test);
            Sensitivity_Dfc = TMP_Confusion_Dfc(1,1) / (TMP_Confusion_Dfc(1,1) + TMP_Confusion_Dfc(1,2));
            Specificity_Dfc = TMP_Confusion_Dfc(2,2) / (TMP_Confusion_Dfc(2,1) + TMP_Confusion_Dfc(2,2));
            w_Dfc = model_Dfc.w'; 
            b_Dfc = model_Dfc.b; 
        elseif numel(Predictor_Dfc) == 0
            AUC_Dfc = 0;
            Sensitivity_Dfc = 0;
            Specificity_Dfc = 0;    
        end

        if numel(Predictor_combine) > 0
            model_combine = train_svm(cfg, i_Mat_combine_train{ifeature}, Label_train);
            [Pred_combine, dval] = test_svm(model_combine, i_Mat_combine_test{ifeature});
            AUC_combine = mv_calculate_performance('auc', 'dval', dval, Label_test);
            TMP_Confusion_combine = mv_calculate_performance('confusion', 'clabel', Pred_combine, Label_test);
            Sensitivity_combine = TMP_Confusion_combine(1,1) / (TMP_Confusion_combine(1,1) + TMP_Confusion_combine(1,2));
            Specificity_combine = TMP_Confusion_combine(2,2) / (TMP_Confusion_combine(2,1) + TMP_Confusion_combine(2,2));
            w_combine = model_combine.w'; 
            b_combine = model_combine.b;   
        elseif numel(Predictor_combine) == 0
            AUC_combine = 0;
            Sensitivity_combine = 0;
            Specificity_combine = 0;
        end

        
        % resolve error when no feature selected. 10 is a just arbitrary number
        if numel(Predictor_fc) == 0, w_fc = []; b_fc = []; end
        if numel(Predictor_Dfc) == 0, w_Dfc = []; b_Dfc = []; end
        if numel(Predictor_combine) == 0, w_combine = []; b_combine = []; end

        % Save
        TMP_AUC_fc{1,ifeature}(nCV,nIter) = AUC_fc; 
        TMP_AUC_Dfc{1,ifeature}(nCV,nIter) = AUC_Dfc; 
        TMP_AUC_combine{1,ifeature}(nCV,nIter) = AUC_combine; 
        TMP_Sensitivity_fc{1,ifeature}(nCV,nIter) = Sensitivity_fc; 
        TMP_Sensitivity_Dfc{1,ifeature}(nCV,nIter) = Sensitivity_Dfc; 
        TMP_Sensitivity_combine{1,ifeature}(nCV,nIter) = Sensitivity_combine; 
        TMP_Specificity_fc{1,ifeature}(nCV,nIter) = Specificity_fc; 
        TMP_Specificity_Dfc{1,ifeature}(nCV,nIter) = Specificity_Dfc; 
        TMP_Specificity_combine{1,ifeature}(nCV,nIter) = Specificity_combine; 

        Cell_svm_w_fc{nCV,ifeature,nIter} = w_fc;
        Cell_svm_w_Dfc{nCV,ifeature,nIter} = w_Dfc;
        Cell_svm_w_combine{nCV,ifeature,nIter} = w_combine;
        Cell_svm_b_fc{nCV,ifeature,nIter} = b_fc;
        Cell_svm_b_Dfc{nCV,ifeature,nIter} = b_Dfc;
        Cell_svm_b_combine{nCV,ifeature,nIter} = b_combine;

   
    end
    
end

end

%% Save results
% Save SVM/SVR results
if option.whichlambda == 1
    Dir_name = sprintf('SVM_setlambda_lasso_n%d_%s_%diter',numel(Idx_con)+numel(Idx_pat),option.slidingwindow,num_iter);
elseif option.whichlambda == 2
    Dir_name = sprintf('SVM_feature%d_lasso_n%d_%s_%diter',n_feature,numel(Idx_con)+numel(Idx_pat),option.slidingwindow,num_iter);
end
mkdir(Dir_name); cd(Dir_name);

for ifeature = 1:n_feature
    Cell_AUC_fc{1,ifeature} = mean(TMP_AUC_fc{1,ifeature},1)'; 
    Cell_AUC_Dfc{1,ifeature} = mean(TMP_AUC_Dfc{1,ifeature},1)'; 
    Cell_AUC_combine{1,ifeature} = mean(TMP_AUC_combine{1,ifeature},1)'; 
    Cell_Sensitivity_fc{1,ifeature} = mean(TMP_Sensitivity_fc{1,ifeature},1)'; 
    Cell_Sensitivity_Dfc{1,ifeature} = mean(TMP_Sensitivity_Dfc{1,ifeature},1)'; 
    Cell_Sensitivity_combine{1,ifeature} = mean(TMP_Sensitivity_combine{1,ifeature},1)'; 
    Cell_Specificity_fc{1,ifeature} = mean(TMP_Specificity_fc{1,ifeature},1)'; 
    Cell_Specificity_Dfc{1,ifeature} = mean(TMP_Specificity_Dfc{1,ifeature},1)'; 
    Cell_Specificity_combine{1,ifeature} = mean(TMP_Specificity_combine{1,ifeature},1)'; 
end
save Cell_AUC_fc Cell_AUC_fc
save Cell_AUC_Dfc Cell_AUC_Dfc
save Cell_AUC_combine Cell_AUC_combine
save Cell_Sensitivity_fc Cell_Sensitivity_fc
save Cell_Sensitivity_Dfc Cell_Sensitivity_Dfc
save Cell_Sensitivity_combine Cell_Sensitivity_combine
save Cell_Specificity_fc Cell_Specificity_fc
save Cell_Specificity_Dfc Cell_Specificity_Dfc
save Cell_Specificity_combine Cell_Specificity_combine

save Cell_svm_b_fc Cell_svm_b_fc
save Cell_svm_b_Dfc Cell_svm_b_Dfc
save Cell_svm_b_combine Cell_svm_b_combine
save Cell_svm_w_fc Cell_svm_w_fc
save Cell_svm_w_Dfc Cell_svm_w_Dfc
save Cell_svm_w_combine Cell_svm_w_combine

% Save as text
for nIter = 1:size(Cell_AUC_fc,2)   % Each Column means 'num of features'
    TXT_AUC_fc(:,nIter) = Cell_AUC_fc{nIter};
    TXT_AUC_Dfc(:,nIter) = Cell_AUC_Dfc{nIter};
    TXT_AUC_combine(:,nIter) = Cell_AUC_combine{nIter};
end

csvwrite('TXT_AUC_fc.txt', TXT_AUC_fc);
csvwrite('TXT_AUC_Dfc.txt', TXT_AUC_Dfc);
csvwrite('TXT_AUC_combine.txt', TXT_AUC_combine);
csvwrite('TXT_AUC_fc.csv', TXT_AUC_fc);
csvwrite('TXT_AUC_Dfc.csv', TXT_AUC_Dfc);
csvwrite('TXT_AUC_combine.csv', TXT_AUC_combine);

if option.whichlambda == 1
   Num_selected_feature = [length(Predictor_fc) length(Predictor_Dfc) length(Predictor_combine)]; 
elseif option.whichlambda == 2
    Num_selected_feature = n_feature;
end
csvwrite('Num_selected_feature.txt', Num_selected_feature);

if option.whichlambda == 2
    Area_wave_fc = sum(TXT_AUC_fc,2)./size(TXT_AUC_fc,2);
    Area_wave_Dfc = sum(TXT_AUC_Dfc,2)./size(TXT_AUC_Dfc,2);
    Area_wave_combine = sum(TXT_AUC_combine,2)./size(TXT_AUC_combine,2);
    
    csvwrite('TXT_Area_wave_fc.txt',Area_wave_fc);
    csvwrite('TXT_Area_wave_fc.csv',Area_wave_fc);
    csvwrite('TXT_Area_wave_Dfc.txt',Area_wave_Dfc);
    csvwrite('TXT_Area_wave_Dfc.csv',Area_wave_Dfc);
    csvwrite('TXT_Area_wave_combine.txt',Area_wave_combine);
    csvwrite('TXT_Area_wave_combine.csv',Area_wave_combine);
end

%% Permutation Testing on AUC
Test_Perm = {'pval'; 'Factual';'Fdist'};

if option.whichlambda == 1
    % Overall Permutation Test
    Vec_AUC = [Cell_AUC_fc{1}' Cell_AUC_Dfc{1}' Cell_AUC_combine{1}'];
    Idx_type = [ones(1,num_iter) ones(1,num_iter)*2 ones(1,num_iter)*3];
    [Test_Perm{1,2},Test_Perm{2,2},Test_Perm{3,2}] = randanova1(Vec_AUC,Idx_type,5000);   % 5000 iterations; outputs are pval / Factual / Fdist
    
    % Testing static-FC & MSSD (1 & 2)
    Vec_AUC = [Cell_AUC_fc{1}' Cell_AUC_Dfc{1}'];
    Idx_type = [ones(1,num_iter) ones(1,num_iter)*2];
    [Test_Perm{1,3},Test_Perm{2,3},Test_Perm{3,3}] = randanova1(Vec_AUC,Idx_type,5000);   % 5000 iterations
    
     % Testing static-FC & Combined (1 & 3)
    Vec_AUC = [Cell_AUC_fc{1}' Cell_AUC_combine{1}'];
    Idx_type = [ones(1,num_iter) ones(1,num_iter)*2];
    [Test_Perm{1,4},Test_Perm{2,4},Test_Perm{3,4}] = randanova1(Vec_AUC,Idx_type,5000);   % 5000 iterations
    
     % Testing MSSD * Combined (2 & 3)
    Vec_AUC = [Cell_AUC_Dfc{1}' Cell_AUC_combine{1}'];
    Idx_type = [ones(1,num_iter) ones(1,num_iter)*2];
    [Test_Perm{1,5},Test_Perm{2,5},Test_Perm{3,5}] = randanova1(Vec_AUC,Idx_type,5000);   % 5000 iterations
   
    Summary_Perm = cell2table(Test_Perm);
    Summary_Perm.Properties.VariableNames = {'Stat','Overall','static_MSSD','static_combine','MSSD_combine'};
    save Summary_Perm Summary_Perm
    
    Summary_Perm(3,:) = []; % Remove Fdist for saving a subsequent txt file
    writetable(Summary_Perm,'Summary_Perm.csv');
    
elseif option.whichlambda == 2
    % Overall Permutation Test
    Vec_AUC = [Area_wave_fc' Area_wave_Dfc' Area_wave_combine'];
    Idx_type = [ones(1,num_iter) ones(1,num_iter)*2 ones(1,num_iter)*3];
    [Test_Perm{1,2},Test_Perm{2,2},Test_Perm{3,2}] = randanova1(Vec_AUC,Idx_type,5000);   % 5000 iterations; outputs are pval / Factual / Fdist
    
    % Testing static-FC & MSSD (1 & 2)
    Vec_AUC = [Area_wave_fc' Area_wave_Dfc'];
    Idx_type = [ones(1,num_iter) ones(1,num_iter)*2];
    [Test_Perm{1,3},Test_Perm{2,3},Test_Perm{3,3}] = randanova1(Vec_AUC,Idx_type,5000);   % 5000 iterations
    
     % Testing static-FC & Combined (1 & 3)
    Vec_AUC = [Area_wave_fc' Area_wave_combine'];
    Idx_type = [ones(1,num_iter) ones(1,num_iter)*2];
    [Test_Perm{1,4},Test_Perm{2,4},Test_Perm{3,4}] = randanova1(Vec_AUC,Idx_type,5000);   % 5000 iterations
    
     % Testing MSSD * Combined (2 & 3)
    Vec_AUC = [Area_wave_Dfc' Area_wave_combine'];
    Idx_type = [ones(1,num_iter) ones(1,num_iter)*2];
    [Test_Perm{1,5},Test_Perm{2,5},Test_Perm{3,5}] = randanova1(Vec_AUC,Idx_type,5000);   % 5000 iterations
   
    Summary_Perm = cell2table(Test_Perm);
    Summary_Perm.Properties.VariableNames = {'Stat','Overall','static_MSSD','static_combine','MSSD_combine'};
    save Summary_Perm Summary_Perm
    
    Summary_Perm(3,:) = []; % Remove Fdist for saving a subsequent txt file
    writetable(Summary_Perm,'Summary_Perm.csv');
    
    % FDR correction was adopted on the result txt files subsequently
    
end


%% Plotting a matrix for which connectivity (or MSSD) features influence the classification
% Summarize mean and std
for nIter = 1:n_feature
    Mean_AUC_fc(1,nIter) = mean(Cell_AUC_fc{nIter});
    Std_AUC_fc(1,nIter) = std(Cell_AUC_fc{nIter});
    Mean_AUC_Dfc(1,nIter) = mean(Cell_AUC_Dfc{nIter});
    Std_AUC_Dfc(1,nIter) = std(Cell_AUC_Dfc{nIter});
    Mean_AUC_combine(1,nIter) = mean(Cell_AUC_combine{nIter});
    Std_AUC_combine(1,nIter) = std(Cell_AUC_combine{nIter});
end

if option.whichlambda == 1
    %% Bar graph
    Values_mean = [Mean_AUC_fc; Mean_AUC_Dfc; Mean_AUC_combine];
    Values_std = [Std_AUC_fc; Std_AUC_Dfc; Std_AUC_combine];
    
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
    saveas(gcf,[Dir_name '.png']);
    close all
    
    %% Weights matrix
    % Calculate the average of weights
    Template_weights1 = zeros(option.nROI);
    Template_weights2 = zeros(option.nROI, option.nROI*2);
    
    for nIter = 1:num_iter  % Cell_svm_w_fc{nCV,ifeature,nIter}
        for nCV = 1:num_cv
            Cur_weight_fc = Cell_svm_w_fc{nCV,1,nIter};
            Cur_weight_Dfc = Cell_svm_w_Dfc{nCV,1,nIter};
            Cur_weight_combine = Cell_svm_w_combine{nCV,1,nIter};
            
            % Case: FC (static-fc)
            if numel(Cur_weight_fc) == 0 % when no feature was survived
                Mat_svm_w_fc(:,:,nCV,nIter) = Template_weights1;
            elseif numel(Cur_weight_fc) > 0
                tmp_svm_w_fc = Template_weights1;
                for nFeature = 1:length(Cur_weight_fc), tmp_svm_w_fc(Idxref_Predictor_fc{nCV,nIter}(nFeature)) = Cur_weight_fc(nFeature); end
                Mat_svm_w_fc(:,:,nCV,nIter) = tmp_svm_w_fc; clear tmp_svm_w_fc
            end

            % Case: FC-variability (MSSD)
            if numel(Cur_weight_Dfc) == 0 % when no feature was survived
                Mat_svm_w_Dfc(:,:,nCV,nIter) = Template_weights1;
            elseif numel(Cur_weight_Dfc) > 0
                tmp_svm_w_Dfc = Template_weights1;
                for nFeature = 1:length(Cur_weight_Dfc), tmp_svm_w_Dfc(Idxref_Predictor_Dfc{nCV,nIter}(nFeature)) = Cur_weight_Dfc(nFeature); end
                Mat_svm_w_Dfc(:,:,nCV,nIter) = tmp_svm_w_Dfc; clear tmp_svm_w_Dfc
            end
            
            % Case: Combine
            if numel(Cur_weight_combine) == 0 % when no feature was survived
                Mat_svm_w_combine(:,:,nCV,nIter) = Template_weights2;
            elseif numel(Cur_weight_combine) > 0
                tmp_svm_w_combine = Template_weights2;
                for nFeature = 1:length(Cur_weight_combine), tmp_svm_w_combine(Idxref_Predictor_combine{nCV,nIter}(nFeature)) = Cur_weight_combine(nFeature); end
                Mat_svm_w_combine(:,:,nCV,nIter) = tmp_svm_w_combine;
            end
            
        end
    end

    
    % Average over 3rd and 4th dimension (nCV, nIter) to calculate the
    % averaged feature matrix
    Plot_svm_w_fc = mean(Mat_svm_w_fc,4); Plot_svm_w_fc = mean(Plot_svm_w_fc,3);
    Plot_svm_w_Dfc = mean(Mat_svm_w_Dfc,4); Plot_svm_w_Dfc = mean(Plot_svm_w_Dfc,3);
    Plot_svm_w_combine = mean(Mat_svm_w_combine,4); Plot_svm_w_combine = mean(Plot_svm_w_combine,3);

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
            if numel(Predictor_fc) > 0
                figure(1); Plot_color = imagesc(Range_x,Range_y,Plot_svm_w_fc); axis off; hold on
                Plot_grid = mesh(xg,yg,zeros([option.nROI+1, option.nROI+1]));
                Plot_grid.FaceColor = 'none'; Plot_grid.EdgeColor = 'k';
                set(gca, 'Ydir','reverse','clim',clim); colormap jet; pbaspect([1 1 1]);
                saveas(gcf,'Plot_svm_w_static.png');
                close all
                save Plot_svm_w_fc Plot_svm_w_fc
            end

            if numel(Predictor_Dfc) > 0
                figure(1); Plot_color = imagesc(Range_x,Range_y,Plot_svm_w_Dfc); axis off; hold on
                Plot_grid = mesh(xg,yg,zeros([option.nROI+1, option.nROI+1]));
                Plot_grid.FaceColor = 'none'; Plot_grid.EdgeColor = 'k';
                set(gca, 'Ydir','reverse','clim',clim); colormap jet; pbaspect([1 1 1]);
                saveas(gcf,'Plot_svm_w_mssd.png');
                close all
                save Plot_svm_w_Dfc Plot_svm_w_Dfc
            end
            
        elseif nPlot == 2
            if numel(Predictor_combine) > 0
                figure(1); Plot_color = imagesc(Range_x,Range_y,Plot_svm_w_combine); axis off; hold on
                Plot_grid = mesh(xg,yg,zeros([option.nROI+1, 2*option.nROI+1]));
                Plot_grid.FaceColor = 'none'; Plot_grid.EdgeColor = 'k';
                set(gca, 'Ydir','reverse','clim',[-0.5,0.5]); colormap jet; pbaspect([2 1 1]);
                saveas(gcf,'Plot_svm_w_combine.png');
                close all
                save Plot_svm_w_combine Plot_svm_w_combine
            end
        end
    end
    
elseif option.whichlambda == 2    
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
    shadedErrorBar(1:n_feature,Mean_AUC_fc,Std_AUC_fc,{'k','markerfacecolor','r'}); hold on
    shadedErrorBar(1:n_feature,Mean_AUC_Dfc,Std_AUC_Dfc,{'b','markerfacecolor','r'}); hold on
    shadedErrorBar(1:n_feature,Mean_AUC_combine,Std_AUC_combine,{'g','markerfacecolor','r'}); 
    legend({'Static','MSSD','Combined'},'Location','southeast');
    ylabel('AUC');
    saveas(gcf,[Dir_name '.png']);
    close all

    
end





