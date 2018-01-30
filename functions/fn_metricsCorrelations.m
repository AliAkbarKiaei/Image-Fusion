function fn_metricCorrelations(fileName, metricNames,fusionMethodNames,sourceImages)

xlsFileRead = strcat(fileName,'.xlsx');

%% correlation between metrics
correlations={};
for metricNumber1 = 1:length(metricNames)-1
    for metricNumber2 = metricNumber1+1:length(metricNames)
        metricName1 = metricNames{metricNumber1};
        metricName2 = metricNames{metricNumber2};
        [QC1, methods1]=xlsread(xlsFileRead, metricName1);    QC1 = QC1(1:length(sourceImages),1:length(fusionMethodNames));
        [QC2, methods2]=xlsread(xlsFileRead, metricName2);    QC2 = QC2(1:length(sourceImages),1:length(fusionMethodNames));
        corrMetrics = corr2(QC1,QC2);
        %disp([metricName1,', ',metricName2,': ', num2str(corrMetrics) ]);
        correlations = [correlations; {metricName1,metricName2,corrMetrics}];
    end
end
[s, idx]=sort(cell2mat(correlations(:,end)));
sortedCorrelations = correlations(idx,:);
xlswrite(strcat(fileName,'_Friedman.xlsx'),correlations, 'corrMethods');
%% correlation with Friedman ranking
[FriedmanRanks, methods] = xlsread(strcat(fileName,'_Friedman.xlsx'),'FriedmanRanks');
FriedmanRanks = FriedmanRanks(end:-1:1); methods = methods(end:-1:1);
FriedRes = max(FriedmanRanks) - FriedmanRanks +1;
corrMtricFried={};
for metricNumber = 1:length(metricNames)
    metricName = metricNames{metricNumber};
    [QC, methods1]=xlsread(xlsFileRead, metricName);    QC = QC(1:length(sourceImages),1:length(fusionMethodNames));
    temp = corr2(mean(QC),FriedRes');
    % figure, plot(1: length(FriedRes), mean(QC)), title(metricName)
    % m = mean(QC); disp (m(end-2)-m(end-3))
    corrMtricFried = [corrMtricFried;{metricName,temp }];
end
[s, idx]=sort(cell2mat(corrMtricFried(:,end)));
sortedMtricFried= corrMtricFried(idx,:);
xlswrite(strcat(fileName,'_Friedman.xlsx'),sortedMtricFried(end:-1:1,:), 'bestMatch2Fried');
end