function weightedFusionDel = fn_fusion_GV_max(U1, ROI_BW, ROI_U2, pos)

[rU, cU, hU] = size(U1);
[rROI, cROI] = size(ROI_BW); sizeROI = rROI*cROI*hU;
ROI_U1 = double(U1(pos(1):pos(1)+rROI-1, pos(2):pos(2)+cROI-1, :));
ROI_U2 = double(ROI_U2);

G_ROI_U1 = zeros(size(ROI_U1)); G_ROI_U2 = zeros(size(ROI_U2));
filterForward=[0 -1 1]; filterBackward = [-1 1 0];
for channel = 1:hU
    Gx_ROI_U1(:,:,channel) = conv2(ROI_U1(:,:,channel),filterForward(end:-1:1),'same');
    Gy_ROI_U1(:,:,channel) = conv2(ROI_U1(:,:,channel),filterForward(end:-1:1)','same');
    Gx_ROI_U2(:,:,channel) = conv2(ROI_U2(:,:,channel),filterForward(end:-1:1),'same');
    Gy_ROI_U2(:,:,channel) = conv2(ROI_U2(:,:,channel),filterForward(end:-1:1)','same');
end

%% fusion method

maxAbsX =  max(abs(Gx_ROI_U1) , abs(Gx_ROI_U2));
Gx_ROI_U = maxAbsX.*sign(Gx_ROI_U1).*(abs(Gx_ROI_U1)>abs(Gx_ROI_U2)) + maxAbsX.*sign(Gx_ROI_U2).*(abs(Gx_ROI_U1)<abs(Gx_ROI_U2));

maxAbsY =  max(abs(Gy_ROI_U1 ), abs(Gy_ROI_U2));
Gy_ROI_U = maxAbsY.*sign(Gy_ROI_U1).*(abs(Gy_ROI_U1)>abs(Gy_ROI_U2)) + maxAbsY.*sign(Gy_ROI_U2).*(abs(Gy_ROI_U1)<abs(Gy_ROI_U2));

%% making Laplace from gradients
for channel = 1:hU
    Gxx_ROI_U(:,:,channel) = conv2(Gx_ROI_U(:,:,channel),filterBackward(end:-1:1),'same');
    Gyy_ROI_U(:,:,channel) = conv2(Gy_ROI_U(:,:,channel),filterBackward(end:-1:1)','same');
end
Laplace_ROI_U = Gxx_ROI_U + Gyy_ROI_U;

%% Del iterative Solver
%errMax = fix(rROI*cROI/1000);
mseMax=1e-2;
MSE=mseMax+1;

%ROI_temp=ROI_U1;
%ROI_temp = max(ROI_U1,ROI_U2);
ROI_temp = (ROI_U1.^4 + ROI_U2.^4) ./ (ROI_U1.^3 + ROI_U2.^3+eps);

iteration=1;
disp('iteration:     MSE:     min:     max:');
while ((MSE>mseMax) && (iteration<33))
    for r=3:rROI-2
        for c=3:cROI-2
            if ROI_BW(r,c)==1
                
                ROI_temp(r,c,:) = 0.25*(  ROI_temp(r-1,c,:)+ROI_temp(r+1,c,:)+ROI_temp(r,c-1,:)+ROI_temp(r,c+1,:)  - Laplace_ROI_U(r,c,:) );
                
            end
        end
    end
    ROI_temp(ROI_temp>255)=255; ROI_temp(ROI_temp<0)=0;

    MSE = sum(sum(sum(ROI_temp - ROI_U1).^2))/sizeROI;
    disp([num2str(iteration), ':  ', num2str(MSE), '  ', num2str(min(min(min(ROI_temp)))), '  ', num2str(max(max(max(ROI_temp))))   ]);
    ROI_U1 = ROI_temp;
    iteration = iteration+1;
    if (log2(iteration)==fix(log2(iteration)))
        figure, imshow(uint8(ROI_U1),[]);title(['iteration: ',num2str(iteration)]);
        %imwrite( uint8(ROI_U1), strcat(cd,'/functions/it',num2str(iteration),'.png'));
    end
end
%ROI_U1 = uint8(ROI_U1);
figure, imshow(uint8(ROI_U1));
%% importing to the image
weightedFusionDel = U1;
for r=1:rROI
    for c=1:cROI
        for channel=1:hU %RGB
            if ROI_BW(r,c)==1
                weightedFusionDel(pos(1)+r-1, pos(2)+c-1, channel) = ROI_U1(r,c,channel);
            end
        end
    end
end
figure, imshow(uint8(weightedFusionDel)); title('weighted Fusion: Del edition');




end

