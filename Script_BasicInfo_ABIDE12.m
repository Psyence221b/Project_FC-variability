%% ABIDE 1 part

cd(Dirlog);
MRIInfo1 = readtable('BasicInfo_onlyRunData.csv');
% Change wrong variable name
MRIInfo1.Properties.VariableNames{'SamlpingRate'} = 'SamplingRate';

FDInfo1 = readtable('final_FDmotion_before.csv');

cd(Dirdocu1);
score1 = readtable('Andlab_ABIDE_summary.xlsx');
score1(1,:) = [];

% Trim SamplingRate
List06 = MRIInfo1.SamplingRate >= 0.59 & MRIInfo1.SamplingRate <= 0.61;
MRIInfo1.SamplingRate(List06) = 0.6;
List05 = MRIInfo1.SamplingRate >= 0.49 & MRIInfo1.SamplingRate <= 0.55;
MRIInfo1.SamplingRate(List05) = 0.5;

% Remove Last row
MRIInfo1(end,:) = [];

%% Load demographic information
cd(Dirdocu1);
DemoInfo1 = readtable('Andlab_ABIDE_basic_subject_info.csv');

if strcmp(DemoInfo1.Properties.VariableNames{1}, 'x___Lab_ID')
    DemoInfo1.Properties.VariableNames{'x___Lab_ID'} = 'Lab_ID';
end


TMP_table = table();
for nTable1 = 1:size(MRIInfo1,1)
    for nTable2 = 1:size(DemoInfo1,1)
        
        if strcmp(MRIInfo1.subID(nTable1), DemoInfo1.Lab_ID(nTable2))
            TMP_table(nTable1,:) = DemoInfo1(nTable2,:);
        end
    end
end
TMP_table = removevars(TMP_table, 'Lab_ID');

BasicInfo1 = [MRIInfo1 TMP_table FDInfo1(:,2)];
BasicInfo1.SubjectType{131} = 'CONTROL';
BasicInfo1.Properties.VariableNames{'Var2'} = 'FD';

%% Remove subjects
% Remove subjects if SamplingRate is too low
ListDelete = BasicInfo1.VolumeN <= 100;
BasicInfo1(ListDelete,:) = [];

% Remove subjects if other information is missing
ListDelete = BasicInfo1.Sex ~= 1 & BasicInfo1.Sex ~= 2;
BasicInfo1(ListDelete,:) = [];

% Remove subjects if he/she has lots of movement
ListDelete = BasicInfo1.FD >= 0.25;
BasicInfo1(ListDelete,:) = [];

%% Add ADOS, SRS, SQT scores
% find subject indices corresponding to score variable
Idx_subj = [];
parfor nSubj = 1:size(BasicInfo1,1)
    Idx_subj(nSubj) = find(strcmp(BasicInfo1.subID(nSubj),score1.Var3));    
end

L_ADOS = score1(Idx_subj,[34:36,44:46,66:69]);
for i = 1:size(L_ADOS,2), L_ADOS.(i) = str2double(L_ADOS{:,i}); end % change to number instead of character

L_ADI = score1(Idx_subj,[37:41]);
for i = 1:size(L_ADI,2), L_ADI.(i) = str2double(L_ADI{:,i}); end

L_SRS = score1(Idx_subj,[48:49,73:77]);
for i = 1:size(L_SRS,2), L_SRS.(i) = str2double(L_SRS{:,i}); end

L_VABS = score1(Idx_subj,[51:65]);
for i = 1:size(L_VABS,2), L_VABS.(i) = str2double(L_VABS{:,i}); end

% remove -9999 values
Table_score = [L_ADOS L_ADI L_SRS L_VABS];
for i = 1:size(Table_score,2)
    Idx_NaN = find(Table_score{:,i} == -9999);
    Table_score(Idx_NaN,i) = {NaN};
end

% Add to the BasicInfo table

BasicInfo1 = [BasicInfo1 Table_score];

%%
clear TMP_table nTable1 nTable2 List05 List06 i Idx_subj L_ADOS L_ADI L_SRS L_VABS Table_score


%% ABIDE 2 Part ============================================

cd(Dirlog);
MRIInfo2 = readtable('BasicInfo_onlyRunData_abide2.csv');
FDInfo2 = readtable('final_FDmotion_before_abide2.csv');

cd(Dirdocu2);
score2 = readtable('Andlab_ABIDE_summary_abide2.xlsx');

%% Load demographic information
cd(Dirdocu2);
DemoInfo2 = readtable('Andlab_ABIDE_basic_subject_info_abide2.csv');

if strcmp(DemoInfo2.Properties.VariableNames{1}, 'x___Lab_ID')
    DemoInfo2.Properties.VariableNames{'x___Lab_ID'} = 'Lab_ID';
end


TMP_table = table();
for nTable1 = 1:size(MRIInfo2,1)
    for nTable2 = 1:size(DemoInfo2,1)
        
        if strcmp(MRIInfo2.subID(nTable1), DemoInfo2.Lab_ID(nTable2))
            TMP_table(nTable1,:) = DemoInfo2(nTable2,:);
        end
    end
end
TMP_table = removevars(TMP_table, 'Lab_ID');

BasicInfo2 = [MRIInfo2 TMP_table FDInfo2(:,2)];
BasicInfo2.Properties.VariableNames{'Var2'} = 'FD';



%% Remove subjects
% Remove subjects if SamplingRate is too low
ListDelete = BasicInfo2.VolumeN <= 100;
BasicInfo2(ListDelete,:) = [];

% Remove subjects if other information is missing
ListDelete = BasicInfo2.Sex ~= 1 & BasicInfo2.Sex ~= 2;
BasicInfo2(ListDelete,:) = [];

% Remove subjects if he/she has lots of movement
ListDelete = BasicInfo2.FD >= 0.25;
BasicInfo2(ListDelete,:) = [];

%% Add ADOS, SRS, SQT scores
% find subject indices corresponding to score variable
Idx_subj = [];
parfor nSubj = 1:size(BasicInfo2,1)
    Idx_subj(nSubj) = find(strcmp(BasicInfo2.subID(nSubj),score2.Lab_ID));    
end

L_ADOS = score2(Idx_subj,[19,27,21:26,20,28]);
% for i = 1:size(L_ADOS,2), L_ADOS.(i) = str2double(L_ADOS{:,i}); end % change to number instead of character

L_ADI = score2(Idx_subj,[31:35]);
% for i = 1:size(L_ADI,2), L_ADI.(i) = str2double(L_ADI{:,i}); end

L_SRS = score2(Idx_subj,[36:37,40:44]);
% for i = 1:size(L_SRS,2), L_SRS.(i) = str2double(L_SRS{:,i}); end

L_VABS = score2(Idx_subj,[62:76]);
% for i = 1:size(L_VABS,2), L_VABS.(i) = str2double(L_VABS{:,i}); end

% remove -9999 values
Table_score = [L_ADOS L_ADI L_SRS L_VABS];
for i = 1:size(Table_score,2)
    Idx_NaN = find(Table_score{:,i} == -9999);
    Table_score(Idx_NaN,i) = {NaN};
end

% Add to the BasicInfo table

BasicInfo2 = [BasicInfo2 Table_score];

%% Calculate quartile of age to choose a criterion of age group
Tot_AgeAtScan = [DemoInfo1.AgeAtScan; DemoInfo2.AgeAtScan];

Age_quartile = quantile(Tot_AgeAtScan,3);

%% Make a merged file
% Match varaible names
BasicInfo2.Properties.VariableNames = BasicInfo1.Properties.VariableNames;

BasicInfo = vertcat(BasicInfo1, BasicInfo2);

%%
clear DemoInfo* MRIInfo* TMP_table nTable1 nTable2 List05 List06 FDInfo* i Idx_subj L_ADOS L_ADI L_SRS L_VABS Table_score

