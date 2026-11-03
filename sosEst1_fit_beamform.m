function [zErr, xErr, c1_est, flag] = sosEst1_fit_beamform(params)
standoff = params.standoff;
R = params.pipe_IR;
zOffset = params.zOffset;
RFdata = params.RFdata;
TXdelays = params.TXdelays;
Nt = params.Nt;
Ns = params.Ns;
Na = params.Na;
ds = params.ds;
xv = params.xv;
zv = params.zv;
t_axis = params.t_axis;
s_axis = params.s_axis;
dt = params.dt;
dx = params.dx;
c0 = params.c0;

minDepth = 0.25*standoff;
maxDepth = 1.25*standoff;
nPeaks = 1;
peakHeightFrac = 0.25;
filterSize = [3 3];
searchLen = 60;
thres = 1e6;                    % beamformed data energy threshold in misalignment est.


skip = 0;%round(0.25*Ns);
sIndRange = ((1+skip):(Ns-skip))';
Nss = numel(sIndRange);
tMax = cell(1,3);
delay = cell(1,3);
for ia = 1:Na
    delay{ia} = TXdelays(:,ia);
    Dts = RFdata(:,:,ia);
    chEnvTX = abs(hilbert(Dts)); 
    tStart = round(TXdelays(:,ia)/dt) + round(minDepth/(c0/2)/dt);
    tEnd = tStart + round(maxDepth/(c0/2)/dt);
    % Find out the time of the first receive 
    maxLoc = zeros(Ns,1);
    for is = 1:Ns
        tTrace = chEnvTX(tStart(is):tEnd(is),is);
        [~,loc] = findpeaks(tTrace,'NPeaks',nPeaks,'MinPeakHeight',peakHeightFrac*max(tTrace));
        maxLoc(is) = loc + tStart(is) - 1;
    end
    tMax{ia} = t_axis(maxLoc)';
    %tMax{ia} = movmean(tM,filterSize);
end

% Initial guess: [zerr xerr c1]
p0 = [0; 0; 1480];
% Set lower and upper bounds
lb = [ -0.01, -0.01, 1400 ];
ub = [  0.01,  0.01, 1560 ];
% Create the optimization problem for lsqcurvefit.
options = optimoptions('lsqcurvefit', ...
    'OptimalityTolerance', 1e-18, ...   
    'StepTolerance', 1e-12, ...         
    'FunctionTolerance', 1e-18, ...
    'Display','final'); 


tolRatioSoS = 0.005;  
tolRatiaMis = 0.05;
maxIter = 100; % Maximum number of iterations
p = p0;

for k = 1:maxIter
    p_prev = p;
    params.c1 = p(3);    
    params.xErr0 = p(2);               % initial lateral misalignment [m]
    params.zErr0 = p(1);               % initial axial misalignment [m]
    [p(2),p(1),flag] = f_misalignEst_IQ(params,thres);

    fun3 = @(p3,x) modelFun([p(1); p(2); p3], x); % fix z_error, x_error, fit c1
    p(3) = lsqcurvefit(fun3, p(3), s_axis, zeros(Nss,1), lb(3), ub(3), options);
    
    if abs(p(3)-p_prev(3)) < tolRatioSoS * c0 
        % abs(p(2)-p_prev(2)) < abs(tolRatiaMis * p_prev(2)) && ...
        % abs(p(1)-p_prev(1)) < abs(tolRatiaMis * p_prev(1))
        fprintf('Converged after %d iterations.\n', k);
        break
    end 
end

[p(2),p(1),flag] = f_misalignEst_IQ(params,thres);
assert(flag==true); 
zErr = p(1);
xErr = p(2);
c1_est = round(p(3));
fprintf('Final fit: z_error = %g, x_error = %g, c1 = %g\n', p);

function res = modelFun(p, x)
    z_error = p(1);
    x_error = p(2);
    c1      = p(3); 
    res = [];
    for ia=2
        curRes = zeros(Nss,1);        
        for isR = 1:Nss
            searchRange = (max(1,sIndRange(isR)-searchLen) : min(Ns,sIndRange(isR)+searchLen))';
            travelDist = zeros(numel(searchRange),1);
            for isT = 1:numel(searchRange)
                wallRange = min(x(searchRange(isT)),x(isR)) : dx : max(x(searchRange(isT)),x(isR));
                wallDist = sqrt(R^2-(wallRange+x_error).^2)-zOffset+z_error;
                TRDist = sqrt((wallRange-x(searchRange(isT))).^2 + wallDist.^2) + sqrt((wallRange-x(isR)).^2 + wallDist.^2);
                travelDist(isT) = min(TRDist);
            end
            minTravelTime = min(travelDist/c1 + delay{ia}(searchRange));

            curRes(isR) = minTravelTime - tMax{ia}(sIndRange(isR));
        end
        res = [res;curRes];
    end
end

end