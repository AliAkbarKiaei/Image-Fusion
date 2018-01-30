clc
clear
close all
addpath(strcat(cd,'/fusion methods'));
addpath(strcat(cd,'/functions'));
addpath(strcat(cd,'/images'));
%%
sourceImages = {'G01' 'CT' 'MR_T2' ; 'G02' 'CT' 'MR_T2' ; 'G03' 'MR_T1' 'MR_T2' ; 'G04' 'MR_T1' 'MR_T2' ; ...
    'G05' 'CT' 'MR' ; 'G06' 'CT' 'MR' ; 'G07' 'CT' 'MR' ; 'G08' 'CT' 'MR' ; ...
    'G09' 'CT' 'MR_T2' ; 'G10' 'MR_PD' 'MR_T2' ; 'G11' 'CT' 'MR_GAD' ; 'G12' 'CT' 'MR_T1'}

fusionMethods = {'fsd' 'FSD' ; 'gra' 'GP' ; 'dwb' 'DWT'; 'rat' 'RP' ; ...
    'mod' 'MDP' ; 'pca' 'PCA' ; 'lap' 'LP' ; 'sih' 'SIDWT'}

for sourceNumber=1:length(sourceImages)
    for methodNumber=1:length(fusionMethods)
        %% input images
        sourceImage = char(sourceImages(sourceNumber,1));
        modal1 = char(sourceImages(sourceNumber,2));
        modal2 = char(sourceImages(sourceNumber,3));
        
        U1 = imread(strcat(sourceImage,'_',modal1,'.png'));
        U2 = imread(strcat(sourceImage,'_',modal2,'.png'));
        if (size(U1,3)>1)
            U1 = rgb2gray(U1);
        end
        if (size(U2,3)>1)
            U2 = rgb2gray(U2);
        end
        U1=double(U1); U2=double(U2);
        %% fusion methods
        fusionMethodFile = char(fusionMethods(methodNumber,1));
        fusionMethodName = char(fusionMethods(methodNumber,2));
        args={U1,U2,4,1,3};
        if (fusionMethodFile == 'pca')
            args={U1,U2};
        end
        %% fuse based on images and methods
        fusedFunc = str2func(strcat('fuse_',fusionMethodFile));
        disp(['input source: ', sourceImage, '       method: fuse_', fusionMethodName]);
        switch size(args,2)
            case 2
                fusedImage = uint8(fusedFunc( cell2mat(args(1)), cell2mat(args(2)) ));
            case 5
                fusedImage = uint8(fusedFunc( cell2mat(args(1)), cell2mat(args(2)) , cell2mat(args(3)) , cell2mat(args(4))  , cell2mat(args(5)) ));
        end
        
        %imshow((fusedImage));
        imwrite( fusedImage, strcat(cd,'/images/',sourceImage,'_fused_',fusionMethodName,'.png'));
        %disp(['fused image is written in image:  ', fileName, '_',method,'.png']);
        
    end
    disp('-------------------------------------');
end
