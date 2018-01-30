function theEI = metric_Edge_Intensity(I)
    I = double(I);     [rI cI hI] = size(I);    filter = fspecial('sobel');
    Ix = imfilter(I,filter,'replicate');    Iy = imfilter(I,filter','replicate');
    for r = 1 : rI
        for c = 1 : cI
            for h = 1 : hI
                EIs(r,c,h) = sqrt(Ix(r,c,h)*Ix(r,c,h) + Iy(r,c,h)*Iy(r,c,h));
            end
        end
    end
    theEI = mean(mean(mean(EIs)));
end