function thePSNR=metric_PSNR(I1,I2)
	M1=max(I1(:));
	M2=max(I2(:));
	M=max(M1,M2);
	s=mean2((I1-I2).^2);
	thePSNR=10*log((M.^2)./s);
end
