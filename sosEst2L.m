function c2Est = sosEst2L(params,frontDist,xErr,zErr)
standoff = params.standoff;
layer2 = params.layer2;
RFdata = params.RFdata;
TXdelays = params.TXdelays;
Nt = params.Nt;
Ns = params.Ns;
Na = params.Na;
t_axis = params.t_axis;
dt = params.dt;

c0 = params.c0;
c1 = params.c1;
c2L = params.c2L;

filterSobel = fspecial('sobel');
minDepth = 0.25*standoff;
maxDepth = 1.25*standoff;
nPeaks = 1;
peakHeightFrac = 0.25;
filterSize = [3 3];

%% Data Warping -- prepare for further analysis
ia = (1+Na)/2;
Dts = RFdata(:,:,ia);
chEnvTX = abs(hilbert(Dts)); 
tExt = round(minDepth/(c0/2)/dt) + 1;
tStart = round(TXdelays(:,ia)/dt) + tExt;
tEnd = tStart + round(maxDepth/(c0/2)/dt);
maxLoc = zeros(1,Ns);
for is = 1:Ns
    tTrace = chEnvTX(tStart(is):tEnd(is),is);
    [~,loc] = findpeaks(tTrace,'NPeaks',1,'MinPeakHeight',peakHeightFrac*max(tTrace));
    maxLoc(is) = loc + tStart(is) - 1;
end
tMax = t_axis(maxLoc);
tMax = movmean(tMax,filterSize);

NtFFT = pow2(nextpow2(Nt+1)+2);     % 4X extra padding
F0 = (0:NtFFT/2-1)/dt/NtFFT;        % positive temporal frequencies
Dfs = fft(Dts,NtFFT,1);
Dfs = Dfs(1:NtFFT/2,:);
for is = 1:Ns
    Dfs(:,is) = Dfs(:,is).*exp(1i*2*pi*tMax(is)*F0');
end
Dts_shift = ifft(Dfs,NtFFT,1,'symmetric');
Dts_shift = hilbert(Dts_shift(1:Nt,:));
angle_t = mean(calculate_angles_from_delays(params, frontDist, zErr, xErr, ia));

%% Estimate SoS based on warped data
c2Est = c2L;
angle_r = asin(c2Est/c1*sin(angle_t)); 
% extract time section containing multiple backwall reflections
tIndRange = round(0.5*layer2*cos(angle_r)/(c2L/2)/dt):round(1.5*layer2*cos(angle_r)/(c2L/2)/dt);
Dts_cut = Dts_shift(tIndRange,:);
Dts_filtered = imfilter(Dts_cut,filterSobel,"replicate");
tTrace = abs(sum(Dts_filtered,2));
[~,pLoc] = findpeaks(tTrace,'NPeaks',nPeaks,'MinPeakHeight',peakHeightFrac*max(tTrace));
delay = (tIndRange(1)+pLoc-2)*dt;         % average reverberation delay
% Solve the equation for c2: c2 = 2b*cos(arcsin(a*c2))
a = sin(angle_t) / c1;
b = 2*layer2/delay; 
c2Est = round(b / sqrt(a^2*b^2 + 1));       % estimated longitudinal SoS
end