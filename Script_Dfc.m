function [TMP_Dfc] = Script_Dfc(Sig_input,windowsize,option)
% Calculate Dfc
% if option = 1 sliding window / option = 2 tapered sliding window
if option == 1
    [TMP_Dfc] = sliding_window(Sig_input, windowsize); % sliding window
elseif option == 2
    [TMP_Dfc] = tapered_sliding_window(Sig_input,windowsize,3); % tapered sliding window
elseif option == 3
    [TMP_Dfc] = DCC(Sig_input, 'whiten'); % DCC
elseif option == 4                      % Robust fit
    [T,p] = size(Sig_input);

    Ct = NaN(p,p,T);

    
        for j = 1:size(Sig_input,2)
            for k = 1:size(Sig_input,2)
                parfor i=(windowsize):T
                    data1 = Sig_input(:,j); data2 = Sig_input(:,k);
                    [~,~,rob_corrw] = andlab_robustfit( data1((i-windowsize+1):i),data2((i-windowsize+1):i) );
                    Ct(j,k,i) = rob_corrw;
                end
            end
        end
    [TMP_Dfc] = Ct;
    
elseif option == 5 % edge-centric time series
    TMP_Sig_input = normalize(Sig_input,1); % normalize (z-transform) each column; each column is each brain region
    for nRow = 1:size(Sig_input,2)
        for nCol = 1:size(Sig_input,2)
           TMP_Dfc(nRow,nCol,:) = TMP_Sig_input(:,nRow) .* TMP_Sig_input(:,nCol); 
        end
    end
    
end

    Idx_nan = isnan(TMP_Dfc);
    TMP_Dfc(Idx_nan) = [];
    TMP_Dfc = reshape(TMP_Dfc,size(Sig_input,2),size(Sig_input,2),[]);
    if option == 1 || option == 2 || option == 3 || option == 4
        TMP_Dfc = atanh(TMP_Dfc); % Fisher normalization r-to-z
    end
    Idx_nan = [];
end

