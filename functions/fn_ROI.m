function [ROI_BW, ROI_U2]=fn_ROI(U2)
figure;
[rU, cU, hU] = size(U2);
[BW, y, x] = roipoly(U2);%   figure, imshow(BW);
xmin = single(min(x)); xmax = single(max(x));
ymin = single(min(y)); ymax = single(max(y));
ROI_BW = BW(xmin:xmax, ymin:ymax); ROI_BW(1:2,:)=0; ROI_BW(end-1:end,:)=0; ROI_BW(:,1:2)=0; ROI_BW(:,end-1:end)=0;
figure, imshow(ROI_BW); title('ROI BW');

ROI_U2 = U2(xmin:xmax, ymin:ymax, :);
ROI_temp = ROI_U2;
for channel = 1:hU   %RGB
    temp = ROI_U2(:,:,channel);
    temp(ROI_BW==0)=0;
    ROI_temp(:,:,channel) = temp;
end
figure, imshow(ROI_temp); title('ROI');
end
