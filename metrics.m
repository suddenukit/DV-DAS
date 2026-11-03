fig1 = openfig('./images/DV15143/das256_[DV15143]_frame320_TTdTT_neg_Env_est.fig','invisible');
h1 = findobj(fig1,'Type','image');
I1 = h1.CData;      % first image array
close(fig1);

fig2 = openfig('./images/DV15143/das256_[DV15143]_frame320_TTdTT_neg_Env.fig','invisible');
h2 = findobj(fig2,'Type','image');
I2 = h2.CData;      % second image array
close(fig2);

%% 2. Build overlap mask from non-zero pixels in each
mask1 = (I1 ~= 0);
mask2 = (I2 ~= 0);
overlapMask = mask1 & mask2;  

%% 3. Vectorize the overlapping pixels
v1 = double(I1(overlapMask));
v2 = double(I2(overlapMask));

%% 4. Compute MSE
MSE = immse(v1,v2);

%% 5. Compute PSNR
MAX_I = max([v1; v2]);
PSNR = psnr(v1,v2,MAX_I);

%% 6. Compute SSIM over the overlapping region
SSIM = ssim(v1, v2);

%% 7. Display
fprintf('Overlap MSE   = %.4f\n', MSE);
fprintf('Overlap PSNR  = %.2f dB\n', PSNR);
fprintf('Overlap SSIM  = %.4f\n', SSIM);