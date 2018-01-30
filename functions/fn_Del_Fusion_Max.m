function waterize_Del = fn_Del_Fusion_Max(U1, ROI_BW, ROI_U2, pos)

[rROI, cROI] = size(ROI_BW); sizeROI = rROI*cROI*3;
ROI_U1 = double(U1(pos(1):pos(1)+rROI-1, pos(2):pos(2)+cROI-1, :));
ROI_U2 = double(ROI_U2);

G_ROI_U1 = zeros(size(ROI_U1)); G_ROI_U2 = zeros(size(ROI_U2));
filter=[-1 0 1]
for channel = 1:3
    Gx_ROI_U1(:,:,channel) = conv2(ROI_U1(:,:,channel),filter(end:-1:1),'same')/sum(abs(filter));
    Gy_ROI_U1(:,:,channel) = conv2(ROI_U1(:,:,channel),filter(end:-1:1)','same')/sum(abs(filter));
    Gx_ROI_U2(:,:,channel) = conv2(ROI_U2(:,:,channel),filter(end:-1:1),'same')/sum(abs(filter));
    Gy_ROI_U2(:,:,channel) = conv2(ROI_U2(:,:,channel),filter(end:-1:1)','same')/sum(abs(filter));
end


    mag_G_U1 = Gx_ROI_U1.^2 + Gy_ROI_U1.^2;
    mag_G_U2 = Gx_ROI_U2.^2 + Gy_ROI_U2.^2;
    max_mag_G1 = max(mag_G_U1,[],3);
    max_mag_G2 = max(mag_G_U2,[],3);
    isG1 = (max_mag_G1>max_mag_G2);
    
    for channel =1:3
        Gx_ROI_U2(:,:,channel) = isG1.*Gx_ROI_U1(:,:,channel) + (1-isG1).*Gx_ROI_U2(:,:,channel);
        Gy_ROI_U2(:,:,channel) = isG1.*Gy_ROI_U1(:,:,channel) + (1-isG1).*Gy_ROI_U2(:,:,channel);
    end
    



%Gx_ROI_U2 = max(Gx_ROI_U2, Gx_ROI_U1);
%Gy_ROI_U2 = max(Gy_ROI_U2, Gy_ROI_U1);

%% Del iterative Solver
%errMax = fix(rROI*cROI/1000);  
mseMax=1e-4;
MSE=mseMax+1;    ROI_temp=ROI_U1; iteration=1;
disp('iteration:     MSE:     min:     max:');
while ((MSE>mseMax) && (iteration<1025))
    for r=3:rROI-2
        for c=3:cROI-2
            if ROI_BW(r,c)==1
                %ROI_temp(r, c, :) = (ROI_U1(r-1,c,:) + ROI_U1(r+1,c,:)+...
                %    ROI_U1(r,c-1,:) + ROI_U1(r,c+1,:))/4 - w*G_ROI_U2(r,c,:);
                
                %ROI_temp(r,c,:) = .25*( ROI_temp(r-2,c-2,:)+ROI_temp(r-2,c+2,:)+ROI_temp(r+2,c-2,:)+ROI_temp(r+2,c+2,:) )...
                %    + 0.25*( Gx_ROI_U2(r-2,c-1,:)+Gx_ROI_U2(r+2,c-1,:)-Gx_ROI_U2(r-2,c+1,:)-Gx_ROI_U2(r+2,c+1,:)+2*Gx_ROI_U2(r,c-1,:)-2*Gx_ROI_U2(r,c+1,:) )...
                %    + 0.25*( Gy_ROI_U2(r-1,c-2,:)+Gy_ROI_U2(r-1,c+2,:)-Gy_ROI_U2(r+1,c-2,:)-Gy_ROI_U2(r+1,c+2,:)+2*Gy_ROI_U2(r-1,c,:)-2*Gy_ROI_U2(r+1,c,:) );
                
                ROI_temp(r,c,:) = .25*( ROI_temp(r-2,c-2,:)+ROI_temp(r-2,c+2,:)+ROI_temp(r+2,c-2,:)+ROI_temp(r+2,c+2,:) )...
                
            end
        end
    end
    MSE = sum(sum(sum(ROI_temp - ROI_U1).^2))/sizeROI;
    disp([num2str(iteration), ':  ', num2str(MSE), '  ', num2str(min(min(min(ROI_temp)))), '  ', num2str(max(max(max(ROI_temp))))   ]);
    ROI_U1 = ROI_temp;
    %ROI_U1(ROI_U1>255)=255; ROI_U1(ROI_U1<0)=0;
    iteration = iteration+1;
    if (log2(iteration)==fix(log2(iteration)))
        figure, imshow(uint8(ROI_U1));title(['iteration: ',num2str(iteration)]);
    end
end
ROI_U1 = uint8(ROI_U1);
figure, imshow(ROI_U1);
%% importing to the image
waterize_Del = U1;
for r=1:rROI
    for c=1:cROI
        for channel=1:3 %RGB
            if ROI_BW(r,c)==1
                waterize_Del(pos(1)+r-1, pos(2)+c-1, channel) = ROI_U1(r,c,channel);
            end
        end
    end
end
figure, imshow(waterize_Del); title('waterize: Del edition');




end

