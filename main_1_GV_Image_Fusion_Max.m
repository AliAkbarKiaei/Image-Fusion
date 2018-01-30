clc
clear
close all
addpath(strcat(pwd,'/functions'));
addpath(strcat(pwd,'/images'));

%%
sourceImages = {'G01' 'CT' 'MR_T2' ; 'G02' 'CT' 'MR_T2' ; 'G03' 'MR_T1' 'MR_T2' ; 'G04' 'MR_T1' 'MR_T2' ; ...
    'G05' 'CT' 'MR' ; 'G06' 'CT' 'MR' ; 'G07' 'CT' 'MR' ; 'G08' 'CT' 'MR' ; ...
    'G09' 'CT' 'MR_T2' ; 'G10' 'MR_PD' 'MR_T2' ; 'G11' 'CT' 'MR_GAD' ; 'G12' 'CT' 'MR_T1'}
%%
for sourceNumber=1:length(sourceImages)
    
        sourceImage = char(sourceImages(sourceNumber,1));
        modal1 = char(sourceImages(sourceNumber,2));
        modal2 = char(sourceImages(sourceNumber,3));

        U1 = imread(strcat(sourceImage,'_',modal1,'.png'));
        U2 = imread(strcat(sourceImage,'_',modal2,'.png'));
        if (size(U1,3)>1), U1 = rgb2gray(U1);   end
        if (size(U2,3)>1), U2 = rgb2gray(U2);   end
        U1=double(U1); U2=double(U2);
    
    clc
    disp(['Del_Image_Fusion_Max: ', sourceImage]);
    %%
    close all;
    imshow(uint8(U1)); figure, imshow(uint8(U2));
    pos=0;
    if (pos==[0])
        pos=[1,1];
        [rU, cU, hU]=size(U1);
        ROI_BW=ones(rU,cU);% ROI_BW(2:end-1, 2:end-1)=1;
        ROI_U2= U2; %(2:end-1, 2:end-1,:);
    else
        [ROI_BW, ROI_U2] = fn_ROI(U2);   % find the Region Of Interest from U2
    end
    clone_Standard = fn_cloneStandard(U1, ROI_BW, ROI_U2, pos);   % cloning image by standard method
    %clone_Del = fn_clonePerez(U1, ROI_BW, ROI_U2, pos);   % cloning image by Perez method
    U_fused = uint8(fn_fusion_GV_max(U1, ROI_BW, ROI_U2, pos));   % cloning image by Perez method
    imwrite( U_fused , strcat(cd,'/images/',sourceImage,'_fused_Del_max.png'));
end
