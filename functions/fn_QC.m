function fn_QC(metricNames,fusionMethodNames,sourceImages, xlsFileWrite)

figure('units','normalized','outerposition',[0 0 1 1])
disp('the results of Metrics are written in these sheets:')
for metricNumber = 1:length(metricNames)
    QC = zeros(length(fusionMethodNames),length(sourceImages));
    metricName = metricNames{metricNumber};
    %fusionMethodName names
    
    disp([num2str(metricNumber),': ',metricName]);
    
    for fusionMethodNumber = 1:length(fusionMethodNames)
        for sourceImageNumber=1:length(sourceImages)
            % fusion fusionMethodName
            fusionMethodName = fusionMethodNames{fusionMethodNumber};
            sourceImage = sourceImages{sourceImageNumber,1};
            modal_1 = sourceImages{sourceImageNumber,2};
            modal_2 = sourceImages{sourceImageNumber,3};
            %input images
            
            
            U1 = imread(strcat(sourceImage,'_',modal_1,'.png'));
            U2 = imread(strcat(sourceImage,'_',modal_2,'.png'));
            %strcat(sourceImage,'_fused_',fusionMethodName,'.png')
            U_fused = imread(strcat(sourceImage,'_fused_',fusionMethodName,'.png'));
            if (size(U1,3))>1, U1=rgb2gray(U1); end
            if (size(U2,3))>1, U2=rgb2gray(U2); end
            if (size(U_fused,3))>1, U_fused=rgb2gray(U_fused); end
            U1 = double(U1); U2 = double(U2); U_fused = double(U_fused);
            
            % U1 = histeq(uint8(U1)); U2 = histeq(uint8(U2)); U_fused = histeq(uint8(U_fused));
            % U1=double(U1); U2=double(U2); U_fused=double(U_fused);
            
            %U1 = im2std(U1); U2 = im2std(U2); U_fused = im2std(U_fused);
            
            U1 = AMMSE(U1, U_fused); U2 = AMMSE(U2, U_fused);
            
            if ( strcmp(metricName,'std2') || strcmp(metricName,'entropy') || strcmp(metricName,'metric_Edge_Intensity') || strcmp(metricName,'metric_AverageGradient') || strcmp(metricName,'metric_Shannon') )
                func = str2func(metricName);
                QC(fusionMethodNumber,sourceImageNumber) = func(U_fused);
            elseif (strcmp(metricName,'metric_MI') || strcmp(metricName,'metric_ssim') || strcmp(metricName,'metric_PSNR') )
                func = str2func(metricName);
                metric1 = func(U1, U_fused);
                metric2 = func(U2, U_fused);
                metricMean = (metric1 + metric2)/2;
                QC(fusionMethodNumber,sourceImageNumber) = metricMean;
            elseif ( strcmp(metricName,'metricCvejic') || strcmp(metricName,'metricPeilla'))
                func = str2func(metricName);
                QC(fusionMethodNumber,sourceImageNumber) = func(U1, U2,U_fused,2);
            else
                func = str2func(metricName);
                QC(fusionMethodNumber,sourceImageNumber) = func(U1, U2, U_fused);
            end
            
        end
    end
    
    firstRow = {metricName, fusionMethodNames{1:end}};
    experiments={};
    for exptNumbers=1:length(sourceImages)
        experiments={experiments{1:end}, strcat('expt',num2str(exptNumbers))};
    end
    firstCol =[experiments';'mean';'std'];
    QC = [QC';mean(QC');std(QC')];
    
    subplot(ceil(sqrt(length(metricNames))), ceil(sqrt(length(metricNames))), metricNumber);
    errorbar(1:size(QC,2),QC(end-1,:),QC(end,:)); 
    myTitle = regexprep(metricName,'metric_','','ignorecase');myTitle = regexprep(myTitle,'metric','','ignorecase'); myTitle = regexprep(myTitle,'2','','ignorecase'); title(myTitle);
    myXlabels = strrep(fusionMethodNames,'_','-');
    xticklabels(myXlabels); xtickangle(30); grid on;
    set(gca, 'FontSize', 12); set(gca, 'FontWeight', 'bold')
     
    xlswrite(strcat(pwd,'/',xlsFileWrite),firstRow,metricName);
    xlswrite(strcat(pwd,'/',xlsFileWrite),firstCol,metricName,'A2');
    xlswrite(strcat(pwd,'/',xlsFileWrite),QC,metricName,'B2');
    
end
%disp('the results are written in this file: QC_group1.xlsx');


    function outImg = AMMSE(inpImg, targetNorm)
        % Affine Min Mean Square Error
        [rI cI hI] = size(inpImg);
        outImg = inpImg;
        for ch = 1:hI
            stdNormalize = (inpImg(:,:,ch) - mean2(inpImg(:,:,ch)))/std2(inpImg(:,:,ch));
            ouImg(:,:,ch) = stdNormalize*std2(targetNorm(:,:,ch)) + mean2(targetNorm(:,:,ch));
        end
    end
end