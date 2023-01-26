function [Outmatrix_MEAN Outmatrix_SD Outmatrix_MSSD1 Outmatrix_VSD] = Script_calMSSD(Data_Dfc)

for i = 1:size(Data_Dfc,1)
    parfor j = 1:size(Data_Dfc,2)
        
        DataInput = Data_Dfc(i,j,:);
        % Calculate mean, MSSD1, MSSD2, SD, VSD
        array1 = DataInput(1:end-1);
        array2 = DataInput(1+1:end);
        MEAN = mean(DataInput);
        MSSD1 = sqrt((sum((array2-array1).^2)/(length(array1)))).*1000;
        SD = std(DataInput);
    
        Data2 = DataInput + 1;
        array1 = Data2(1:end-1);
        array2 = Data2(1+1:end);    
        VSD = ( std(abs(array2-array1)) / (sum(Data2)/length(Data2)) ).*1000;
        
        Outmatrix_MEAN(i,j) = MEAN;
        Outmatrix_SD(i,j) = SD;
        Outmatrix_MSSD1(i,j) = MSSD1;

        Outmatrix_VSD(i,j) = VSD;
    end    
    clear array1 array2
end
    
end

