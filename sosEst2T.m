function c2Est = sosEst2T(params,frontDist,xErr,zErr)
standoff = params.standoff;
layer2 = params.layer2;
RFdata = params.RFdata;
TXdelays = params.TXdelays;
fc = params.fc;
Nt = params.Nt;
Ns = params.Ns;
Na = params.Na;
t_axis = params.t_axis;
dt = params.dt;

c0 = params.c0;
c1 = params.c1;
c2T = params.c2T;

filterSobel = fspecial('sobel');
minDepth = 0.25*standoff;
maxDepth = 1.25*standoff;
nPeaks = 1;
peakHeightFrac = 0.25;
filterSize = [3 3];
minPeakWidth = 1/fc/dt;

%% Data Warping -- prepare for further analysis
Dts_shift = cell(1,3);
angle_t = cell(1,3);
for ia = [1 Na]
    Dts = RFdata(:,:,ia);
    chEnvTX = abs(hilbert(Dts));     
    tExt = round(minDepth/(c0/2)/dt) + 1;
    tStart = round(TXdelays(:,ia)/dt) + tExt;
    tEnd = tStart + round(maxDepth/(c0/2)/dt);    
    maxLoc = ones(1,Ns);
    for is = 1:Ns
        tTrace = chEnvTX(tStart(is):tEnd(is),is);
        [vals, locs, ~, proms] = findpeaks(tTrace,'MinPeakHeight',peakHeightFrac*max(tTrace),'MinPeakWidth',minPeakWidth);
        [~, idxMostProminent] = max(proms);
        loc = locs(idxMostProminent);
        if (loc)
            maxLoc(is) = loc + tStart(is) - 1;
        end
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
    temp = ifft(Dfs,NtFFT,1,'symmetric');
    Dts_shift{ia} = hilbert(temp(1:Nt,:));
    angle_t{ia} = mean(calculate_angles_from_delays(params, frontDist, zErr, xErr, ia));
end

%% Estimate SoS based on warped data
c2Est = c2T;
sosSum = 0;
for ia = [1 Na]
    angle_r = asin(c2Est/c1*sin(angle_t{ia}));   
    if (~isreal(angle_r))
        fprintf('An error occurred: refraction angle out of the bound.\n');
        return;
    end
    % extract time section containing 1 backwall reflection
    tIndRange = round(0.5*layer2*cos(angle_r)/(c2Est/2)/dt):round(1.5*layer2*cos(angle_r)/(c2Est/2)/dt);        
    Dts_cut = Dts_shift{ia}(tIndRange,:);
    Dts_filtered = imfilter(Dts_cut,filterSobel,"replicate");
    tTrace = abs(sum(Dts_filtered,2));
    [~,pLoc] = findpeaks(tTrace,'NPeaks',nPeaks,'MinPeakHeight',peakHeightFrac*max(tTrace));
    delay = (tIndRange(1)+pLoc-2)*dt;         % reverberation delay
    % Solve the equation for c2: c2 = 2b*cos(arcsin(a*c2))
    a = sin(angle_t{ia}) / c1;
    b = layer2/delay; 
    c2 = 2 * b / sqrt(4*a^2*b^2 + 1);       % estimated transverse SoS
    sosSum = sosSum + c2;
end
c2Est = round(sosSum / 2);


end