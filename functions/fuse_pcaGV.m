function [a1,a2] = fuse_pcaGV(M1, M2)
%Y = fuse_pca(M1, M2) image fusion with PCA method
%
%    M1 - input image #1
%    M2 - input image #2
%
%    Y  - fused image

%    (Oliver Rockinger 16.08.99)

% check inputs
[z1, s1] = size(M1);
[z2, s2] = size(M2);
if (z1 ~= z2) || (s1 ~= s2)
    error('Input images are not of same size');
end;

% compute, select & normalize eigenvalues
[V, D] = eig(cov([M1(:) M2(:)]));
if (D(1,1) > D(2,2))
    a = V(:,1)./sum(V(:,1));
else
    a = V(:,2)./sum(V(:,2));
end;
a1=a(1); a2=a(2);

% and fuse
Y = a(1)*M1+a(2)*M2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%{
if (sum(isnan(a))>0) || (sum(abs(a)==inf)>0) || (sum(abs(a))>1)
    mM1 = sum(M1(:).^2); mM2 = sum(M2(:).^2); 
    a1 = mM1/(mM1+mM2); a2 = mM2/(mM1+mM2); 
end
%}
%{
a(a>1)=1; a(a<0)=0;
a1=a(1); a2=a(2);
Y = a(1)*M1+a(2)*M2;

if (sum(isnan(Y))>0)
   % disp('NaN');
    a1 = 0.5;
    a2 = 0.5;
end
if (sum(abs(Y)==inf)>0)
    %disp('inf');
    a1 = 0.5;
    a2 = 0.5;
end

if (a1<0)
    a1=-a1;
    a2=a2-2*a1;
elseif (a2<0)
    a2=-a2;
    a1=a1-2*a2;
end
if (a1>1000) || (a1<-1000) || (a2>1000) || (a2<-1000)
   disp('exceed 1000') 
end
%}

end
