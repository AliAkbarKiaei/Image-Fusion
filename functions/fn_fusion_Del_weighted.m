function weightedFusionDel = fn_fusion_Del_weighted(U1, ROI_BW, ROI_U2, pos)

[rU, cU, hU] = size(U1);
[rROI, cROI] = size(ROI_BW); sizeROI = rROI*cROI*hU;
ROI_U1 = double(U1(pos(1):pos(1)+rROI-1, pos(2):pos(2)+cROI-1, :));
ROI_U2 = double(ROI_U2);

G_ROI_U1 = zeros(size(ROI_U1)); G_ROI_U2 = zeros(size(ROI_U2));
filter=[-1 0 1];
for channel = 1:hU
    Gx_ROI_U1(:,:,channel) = conv2(ROI_U1(:,:,channel),filter(end:-1:1),'same')/sum(abs(filter));
    Gy_ROI_U1(:,:,channel) = conv2(ROI_U1(:,:,channel),filter(end:-1:1)','same')/sum(abs(filter));
    Gx_ROI_U2(:,:,channel) = conv2(ROI_U2(:,:,channel),filter(end:-1:1),'same')/sum(abs(filter));
    Gy_ROI_U2(:,:,channel) = conv2(ROI_U2(:,:,channel),filter(end:-1:1)','same')/sum(abs(filter));
end

%% fusion method
[fx1, maxIterate] = myFunc(Gx_ROI_U1); fy1 = myFunc(Gy_ROI_U1);
fx2 = myFunc(Gx_ROI_U2); fy2 = myFunc(Gy_ROI_U2);
Gx_ROI_U = (fx1.*Gx_ROI_U1 + fx2.*Gx_ROI_U2)./(fx1+fx2);
Gy_ROI_U = (fy1.*Gy_ROI_U1 + fy2.*Gy_ROI_U2)./(fy1+fy2);

Gx_ROI_U(isnan(Gx_ROI_U))=0; Gx_ROI_U = min(Gx_ROI_U , max(Gx_ROI_U1,Gx_ROI_U2));  Gx_ROI_U = max(Gx_ROI_U , min(Gx_ROI_U1,Gx_ROI_U2));
Gy_ROI_U(isnan(Gy_ROI_U))=0; Gy_ROI_U = min(Gy_ROI_U , max(Gy_ROI_U1,Gy_ROI_U2));  Gy_ROI_U = max(Gy_ROI_U , min(Gy_ROI_U1,Gy_ROI_U2));

%% Del iterative Solver
%errMax = fix(rROI*cROI/1000);
mseMax=1e-2;
MSE=mseMax+1;    ROI_temp=ROI_U1; iteration=1;
disp('iteration:     MSE:     min:     max:');
while ((MSE>mseMax) && (iteration<maxIterate))
    for r=3:rROI-2
        for c=3:cROI-2
            if ROI_BW(r,c)==1
                
                if ROI_BW(r-1,c)==0,  ROI_temp(r-2,c)=ROI_temp(r-1,c); end
                if ROI_BW(r+1,c)==0,  ROI_temp(r+2,c)=ROI_temp(r+1,c); end
                if ROI_BW(r,c-1)==0,  ROI_temp(r,c-2)=ROI_temp(r,c-1); end
                if ROI_BW(r,c+1)==0,  ROI_temp(r,c+2)=ROI_temp(r,c+1); end
                
                
                ROI_temp(r,c,:) = .25*( ROI_temp(r-2,c-2,:)+ROI_temp(r-2,c+2,:)+ROI_temp(r+2,c-2,:)+ROI_temp(r+2,c+2,:) )...
                    + 0.25*( Gx_ROI_U(r-2,c-1,:)+Gx_ROI_U(r+2,c-1,:)-Gx_ROI_U(r-2,c+1,:)-Gx_ROI_U(r+2,c+1,:)+2*Gx_ROI_U(r,c-1,:)-2*Gx_ROI_U(r,c+1,:) )...
                    + 0.25*( Gy_ROI_U(r-1,c-2,:)+Gy_ROI_U(r-1,c+2,:)-Gy_ROI_U(r+1,c-2,:)-Gy_ROI_U(r+1,c+2,:)+2*Gy_ROI_U(r-1,c,:)-2*Gy_ROI_U(r+1,c,:) );
                
            end
        end
    end
    MSE = sum(sum(sum(ROI_temp - ROI_U1).^2))/sizeROI;
    disp([num2str(iteration), ':  ', num2str(MSE), '  ', num2str(min(min(min(ROI_temp)))), '  ', num2str(max(max(max(ROI_temp))))   ]);
    ROI_U1 = ROI_temp;
    %ROI_U1(ROI_U1>255)=255; ROI_U1(ROI_U1<0)=0;
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

function [f, maxIterate] = myFunc(G)
%f = nthroot(G,3); maxIterate=1025;
%f = log(abs(G)); if (f<0) f=(-1)./f; end; f=sign(G).*f; maxIterate=1025;
f=exp(abs(G)).*sign(G); maxIterate=33;
%f=1; maxIterate=1025;
%f=G.^(5); maxIterate=129;
%f=G./100; maxIterate=1025;
end

