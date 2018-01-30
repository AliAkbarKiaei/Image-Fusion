function fn_ANOVA_Friedman(fileName, metricNames,fusionMethodNames,sourceImages)
%%
xlsFileRead = strcat(fileName,'.xlsx');
whole_ANOVA=[]; whole_Friedman=[];whole_Friedman_sort=[];
disp('the results of these Metrics are read to compute ANOVA and Friedman:')

for metricNumber = 1:length(metricNames)
    metricName = metricNames{metricNumber};

    disp([num2str(metricNumber),': ',metricName]);
    [QC, methods]=xlsread(xlsFileRead, metricName); methods = methods(1,2:end);%,'B1:L9');
    QC = QC(1:length(sourceImages),1:length(fusionMethodNames));
    [rQC, cQC] = size(QC);

    %% making Excel for ANOVA
    QC_ANOVA=[]; Exp={}; method={}; measure={};
    for c=1:cQC
       for r=1:rQC
           QC_ANOVA =  [QC_ANOVA; QC(r,c)];
           measure = [measure; char(metricName)];
           method = [method; char(methods(c))];
           Exp = [Exp; strcat('Exp', num2str(r))];
       end
    end
    metricNameANOVA=[measure, method, Exp, num2cell(QC_ANOVA)];
    whole_ANOVA = [whole_ANOVA; metricNameANOVA]; 
    xlswrite(strcat(pwd,'\',fileName,'_ANOVA.xlsx'),[method, Exp, num2cell(QC_ANOVA)],metricName);

    %% making exel for Friedman
    whole_Friedman=[whole_Friedman;QC];
        
    %%{
    [QCsort, idx] = sort(QC,2,'descend');
    for r=1:rQC
        for c=1:cQC
            x = idx(r,c);
         QC_Friedman(r,x) = c;
        end
    end
    %}
    whole_Friedman_sort = [whole_Friedman_sort; QC_Friedman];
    
    %xlswrite(strcat(pwd,'/QC_ANOVA.xlsx'),QC_ANOVA,metricName);
    %xlswrite(strcat(pwd,'/QC_Friedman.xlsx'),QC_Friedman,metricName);

    
end
xlswrite(strcat(pwd,'\',fileName,'_ANOVA.xlsx'),whole_ANOVA,'whole');
xlswrite(strcat(pwd,'\',fileName,'_Friedman.xlsx'),[methods; num2cell(whole_Friedman)],'whole');
xlswrite(strcat(pwd,'\',fileName,'_Friedman.xlsx'),[methods; num2cell(whole_Friedman_sort)],'whole_Sort');

end