function fn_Friedman_PostHocs(fileName)
%%
xlsFileRead = strcat(fileName,'_friedman.xlsx');
[I, methods]=xlsread(xlsFileRead,'whole_sort');
[N,k] = size(I); %N:#data     k:#methods
mI = mean(I);
ChiSquare = ((12*N)/(k*(k+1))) * ( sum(mI.^2) - (k*(k+1)^2)/4 )

Ff = ((N-1)*ChiSquare) / ( N*(k-1) - ChiSquare )

[mIs, idx] = sort(mI);
mI = mI(idx);
methods = methods(idx); methods=strrep(methods,'_',' ');
postHocs=[];
%% creating axis showing Friedman:
FriedmanRanks = [methods', num2cell(mI')]
xlswrite(xlsFileRead,FriedmanRanks,'FriedmanRanks')
figure,
rmI = mI(end)-mI; %reverse of mI to show better in axis (from smallest performance to largest one)
axis([rmI(end) rmI(1), 1 k]);
xticks(rmI(end:-1:1));
xticklabels(methods(end:-1:1));
xtickangle(45)
yticks(2);
yticklabels('Friedman Ranks');
ytickangle(90)
for idx=1:size(mI,2)
    %dif_mI = round(mI(idx) - mI(idx-1),2);
    x = rmI(idx);%.5*(rmI(idx-1)+rmI(idx));
    text(x,1.5,num2str(mI(idx)),'rotation',90);
end
%% creating axis showing distances:
figure, 
rmI = mI(end)-mI; %reverse of mI to show better in axis (from smallest performance to largest one)
axis([rmI(end) rmI(1), 1 k]);
xticks(rmI(end:-1:1));
xticklabels(methods(end:-1:1));
xtickangle(45)
yticks(2);
yticklabels('distance');
ytickangle(90)
for idx=2:size(mI,2)
    dif_mI = round(mI(idx) - mI(idx-1),2);
    x = .5*(rmI(idx-1)+rmI(idx));
    text(x,1.5,num2str(dif_mI),'rotation',90);
end
%% creating axis for Nemenyi tests:
figure,
rmI = mI(end)-mI; %reverse of mI to show better in axis (from smallest performance to largest one)
axis([rmI(end) rmI(1), 1 k]);
xticks(rmI(end:-1:1));
xticklabels(methods(end:-1:1));
xtickangle(45)
%yticks(mI(end:-1:1));
yticks(1:k);
yticklabels(methods(end:-1:1));
disp('-----------------------------------------------');
%title({'Nemenyi and Bonferroni-Dunn tests:';''});
title({'{\color[rgb]{.6 .6 .6}Nemenyi(\alpha:0.05)  \color[rgb]{0 0 0}Nemenyi(\alpha:0.1) }';''})
%t1 = text(mI(2) ,k-1,'Nemenyi: \alpha=0.05'); t1.BackgroundColor= [0.8 0.8 0.8];
%t2 = text(mI(2) ,k-2,'Nemenyi: \alpha=0.1'); t2.BackgroundColor= [0.6 0.6 0.6];
%t3 = text(mI(2) ,k-3,'Bonferroni: \alpha=0.05'); t3.BackgroundColor= [0.2 0.2 0.2];

%% Nemenyi test with alph = 0.05
q_Nem_5_11=3.22
CD = q_Nem_5_11 * sqrt( k*(k+1)/(6*N) )
disp ('NOT Signigficant differences by Nemenyi (alph=0.05):');
result=[];
for i=1:k
    for j=i+1:k
        if (mI(j)-mI(i) < CD)
            disp (['(',char(methods(i)) ,', ', char(methods(j)),')']);
            line([rmI(i) rmI(j)],[k-i+1 k-i+1],'LineWidth',12,'color',[.6 .6 .6]);
            res=j;
        end
    end
    result=[result;k-res+1];
end
postHocs=[postHocs,{'Nemenyi (\alpha=.05)'};result];

disp('-----------------------------------------------');

%% Nemenyi test with alph = 0.1
q_Nem_10_11=2.976
CD = q_Nem_10_11 * sqrt( k*(k+1)/(6*N) )
disp ('NOT Signigficant differences by Nemenyi (alph=0.05):');
result=[];
for i=1:k
    for j=i+1:k
        if (mI(j)-mI(i) < CD)
            disp (['(',char(methods(i)) ,', ', char(methods(j)),')']);
            line([rmI(i) rmI(j)],[k-i+1 k-i+1],'LineWidth',4,'color',[0 0 0])
            res=j;
        end
    end
    result=[result;k-res+1];
end
postHocs=[postHocs, [{'Nemenyi (\alpha=.1)'}; result]];
disp('-----------------------------------------------');

%% creating axis for Kruskal–Wallis and Bonferroni-Dunn tests:
figure,
rmI = mI(end)-mI; %reverse of mI to show better in axis (from smallest performance to largest one)
axis([rmI(end) rmI(1), 1 k]);
xticks(rmI(end:-1:1));
xticklabels(methods(end:-1:1));
xtickangle(45)
%yticks(mI(end:-1:1));
yticks(1:k);
yticklabels(methods(end:-1:1));
disp('-----------------------------------------------');
%title({'Kruskal–Wallis (\alpha=0.05) ';''});
title({'{\color[rgb]{.6 .6 .6}Kruskal–Wallis(\alpha:0.05)  \color[rgb]{0 0 0}Bonferroni-Dunn(\alpha:0.05) }';''})
%t1 = text(mI(2) ,k-1,'Nemenyi: \alpha=0.05'); t1.BackgroundColor= [0.8 0.8 0.8];
%t2 = text(mI(2) ,k-2,'Nemenyi: \alpha=0.1'); t2.BackgroundColor= [0.6 0.6 0.6];
%t3 = text(mI(2) ,k-3,'Bonferroni: \alpha=0.05'); t3.BackgroundColor= [0.2 0.2 0.2];

%% Kruskal–Wallis test with alph = 0.05
%z_KrsWal_5_11=3.25
alpha=0.05;
p_KrsWal = 1- alpha/((k-1)*(k));
z_KrsWal = norminv(p_KrsWal, 0 ,1);
CD = z_KrsWal * sqrt( k*(k+1)/(6*N) )
disp ('NOT Signigficant differences by Kruskal–Wallis (alph=0.05):');
result=[];
for i=1:k
    for j=i+1:k
        if (mI(j)-mI(i) < CD)
            disp (['(',char(methods(i)) ,', ', char(methods(j)),')']);
            line([rmI(i) rmI(j)],[k-i+1 k-i+1],'LineWidth',12,'color',[.6 .6 .6])
            res=j;
        end
    end
    result=[result;k-res+1];
end
postHocs=[postHocs, [{'Kruskal Wallis (\alpha=.05)'}; result]];
disp('-----------------------------------------------');

%% Bonferroni-Dunn test
%q_Bonf_5_11 = 2.772
alpha=0.05;
q_Bonf = finv(1-alpha/5,k,N)
CD = q_Bonf * sqrt( k*(k+1)/(6*N) )
disp ('NOT Signigficant differences by Nemenyi (alph=0.05):');
result=[];
for i=1:k
    for j=i+1:k
        if (mI(j)-mI(i) < CD)
            disp (['(',char(methods(i)) ,', ', char(methods(j)),')']);
            line([rmI(i) rmI(j)],[k-i+1 k-i+1],'LineWidth',4,'color',[0 0 0])
            res=j;
        end
    end
    result=[result;k-res+1];
end
postHocs=[postHocs, [{'Bonferroni (\alpha=.05)'}; result]];

disp('-----------------------------------------------');

%% creating axis for Holm and Hochberg
figure,
rmI = mI(end)-mI; %reverse of mI to show better in axis (from smallest performance to largest one)
axis([rmI(end) rmI(1), 1 k]);
xticks(rmI(end:-1:1));
xticklabels(methods(end:-1:1));
xtickangle(45)
%yticks(mI(end:-1:1));
yticks(2:3);
yticklabels({'Holm';'Hochberg'});
title('Holm and Hochberg')

%% Holm
alpha=0.05;
SE = sqrt( k*(k+1)/(6*N) )
Z = (mI(end:-1:2)-mI(1))/SE;
P = 1- normcdf(Z) + normcdf(-Z);
alpha_i =  alpha./(k-1:-1:1);
disp ('Signigficant differences by Holm (alph=0.05):');
i=1;
while (i<k && P(i)<alpha_i(i))
    disp(char(methods(k-i+1)));
    i=i+1;
end
if (i==1) i=2;end
line([rmI(k) rmI(k-i+2)],[2 2],'LineWidth',4,'color',[0 0 0])
%% Hochberg
i=k-1;
disp ('NOT signigficant differences by Hochberg (alph=0.05):');
while (P(i)>alpha_i(i) && i>0)
    disp(char(methods(k-i+1)));
    i=i-1;
end
if (i==k-1) i=k-2;end
line([rmI(1) rmI(k-i)],[3 3],'LineWidth',4,'color',[0 0 0])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% creating axis for Holm:
figure,
rmI = mI(end)-mI; %reverse of mI to show better in axis (from smallest performance to largest one)
axis([rmI(end) rmI(1), 1 k]);
xticks(rmI(end:-1:1));
xticklabels(methods(end:-1:1));
xtickangle(45)
%yticks(mI(end:-1:1));
yticks(1:k);
yticklabels(methods(end:-1:1));
disp('-----------------------------------------------');
title('Holm')
%title({'{\color[rgb]{.6 .6 .6}Kruskal–Wallis(\alpha:0.05)  \color[rgb]{0 0 0}Bonferroni-Dunn(\alpha:0.05) }';''})
%t1 = text(mI(2) ,k-1,'Nemenyi: \alpha=0.05'); t1.BackgroundColor= [0.8 0.8 0.8];
%t2 = text(mI(2) ,k-2,'Nemenyi: \alpha=0.1'); t2.BackgroundColor= [0.6 0.6 0.6];
%t3 = text(mI(2) ,k-3,'Bonferroni: \alpha=0.05'); t3.BackgroundColor= [0.2 0.2 0.2];

%% Holm
alpha=0.05;
SE = sqrt( k*(k+1)/(6*N) )
disp ('Signigficant differences by Holm (alph=0.05):');
result=[];
for idx=k:-1:2
    Z = (mI(end:-1:k-idx+2)-mI(k-idx+1))/SE;
    P = 1- normcdf(Z) + normcdf(-Z);
    
    alpha_i =  alpha./(idx-1:-1:1);
    i=1;
    res=1;
    while (i<idx && P(i)<alpha_i(i) )
        disp(['(',char(methods(k-idx+1)),', ',char(methods(k-i+1)),')']);
        i=i+1;
        res=i-1;
    end
    result=[result;res];
    if (i==1) i=2;end
    line([rmI(k) rmI(k-i+2)],[idx idx],'LineWidth',4,'color',[0 0 0])
    text(rmI(k-idx+1),idx,'\leftarrow','FontWeight','Bold');
    
end
postHocs=[postHocs, [{'Holm (\alpha=.05)'}; result]];

%% creating axis for Hochberg:
figure,
rmI = mI(end)-mI; %reverse of mI to show better in axis (from smallest performance to largest one)
axis([rmI(end) rmI(1), 1 k]);
xticks(rmI(end:-1:1));
xticklabels(methods(end:-1:1));
xtickangle(45)
%yticks(mI(end:-1:1));
yticks(1:k);
yticklabels(methods(end:-1:1));
disp('-----------------------------------------------');
title('Hochberg')
%title({'{\color[rgb]{.6 .6 .6}Kruskal–Wallis(\alpha:0.05)  \color[rgb]{0 0 0}Bonferroni-Dunn(\alpha:0.05) }';''})
%t1 = text(mI(2) ,k-1,'Nemenyi: \alpha=0.05'); t1.BackgroundColor= [0.8 0.8 0.8];
%t2 = text(mI(2) ,k-2,'Nemenyi: \alpha=0.1'); t2.BackgroundColor= [0.6 0.6 0.6];
%t3 = text(mI(2) ,k-3,'Bonferroni: \alpha=0.05'); t3.BackgroundColor= [0.2 0.2 0.2];

%% Hochberg
alpha=0.05;
SE = sqrt( k*(k+1)/(6*N) )
disp ('NOT signigficant differences by Hochberg (alph=0.05):');
result=[];
for idx=k:-1:2
    Z = (mI(end:-1:k-idx+2)-mI(k-idx+1))/SE;
    P = 1- normcdf(Z) + normcdf(-Z);
    alpha_i =  alpha./(idx-1:-1:1);
    i=idx-1;
    res=idx-1;
    while (i>0 && P(i)>alpha_i(i) )
        %disp(char(methods(k-i+1)));
        disp(['(',char(methods(k-idx+1)),', ',char(methods(k-i+1)),')']);
        i=i-1;
        res=i;
    end
    result=[result;res];
    if (i==k-1) i=k-2;end
    line([rmI(k-idx+1) rmI(k-i)],[idx idx],'LineWidth',4,'color',[0 0 0])
    text(rmI(k-idx+1),idx,'\leftarrow','FontWeight','Bold');
end
postHocs=[postHocs, [{'Hotchberg (\alpha=.05)'}; result]];

%% Results
figure,
rmI = mI(end)-mI; %reverse of mI to show better in axis (from smallest performance to largest one)
[rP, cP]=size(postHocs);
axis([rmI(end) rmI(1), 0 cP]);
xticks(rmI(end:-1:1));
xticklabels(methods(end:-1:1));
xtickangle(45)
yticks(1:cP);
yticklabels(postHocs(1,1:end)');
disp('-----------------------------------------------');
title({'{\color[rgb]{.6 .6 .6}Del_P_C_A  \color[rgb]{.4 .4 .4}Del_M_a_x  \color[rgb]{0 0 0}Del_w_e_i_g_h_t_e_d }';''})

for idx=1:cP
    result = cell2mat(postHocs(2,idx));
    line([rmI(k) rmI(k-result(3)+1)],[idx idx],'LineWidth',20,'color',[.8 .8 .8])
    line([rmI(k) rmI(k-result(2)+1)],[idx idx],'LineWidth',12,'color',[.6 .6 .6])
    line([rmI(k) rmI(k-result(1)+1)],[idx idx],'LineWidth',4,'color',[0 0 0])
    
end
