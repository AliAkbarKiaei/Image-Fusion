function [clone_standard] = fn_cloneStandard(U1, ROI_BW, ROI_U2, pos)
[rU, cU, hU] = size(U1);
[row, col] = size(ROI_BW);
clone_standard = U1;
for r=1:row
    for c=1:col
        for channel=1:hU %RGB
            if ROI_BW(r,c)==1
                clone_standard(pos(1)+r-1, pos(2)+c-1, channel) = ROI_U2(r,c,channel);
            end
        end
    end
end
figure, imshow(uint8(clone_standard)); title('cloning: standard edition');