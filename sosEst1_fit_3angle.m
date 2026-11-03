function [zErr,xErr,c1_est,flag,diagnostics] = sosEst1_fit_3angle(params)
% Jointly estimate layer-1 sound speed and probe misalignment from 3 angles.

standoff = params.standoff;
R = params.pipe_IR;
RFdata = double(params.RFdata);
TXdelays = params.TXdelays;
t_axis = params.t_axis(:);
s_axis = params.s_axis(:);
dt = params.dt;
dx = params.dx;
c0 = params.c0;
t0 = params.t0;

minDepth = 0.25*standoff;
maxDepth = 1.25*standoff;
peakHeightFrac = 0.25;

[tObserved,peakConfidence] = pick_arrival_times(RFdata,TXdelays, ...
    t_axis,dt,c0,minDepth,maxDepth,peakHeightFrac);

maxAngle = max(abs(params.TXangles(:)));
wallMargin = standoff*tan(maxAngle) + 2*mean(diff(s_axis));
wallX = (s_axis(1)-wallMargin:dx:s_axis(end)+wallMargin).';

% q is dimensionless so the optimizer sees similarly scaled variables.
pOffset = [0,0,c0];
pScale = [1e-3,1e-3,100];
p0 = [0,0,c0];
lb = [-10e-3,-10e-3,1400];
ub = [ 10e-3, 10e-3,1560];
q0 = (p0-pOffset)./pScale;
qlb = (lb-pOffset)./pScale;
qub = (ub-pOffset)./pScale;

valid = isfinite(tObserved) & peakConfidence > 0;
weights = sqrt(peakConfidence(valid));
residualFun = @(q) joint_residual(q,pOffset,pScale,wallX,R, ...
    standoff,s_axis,TXdelays,t0,tObserved,valid,weights,dt);

options = optimoptions('lsqnonlin', ...
    'Display','final', ...
    'MaxIterations',100, ...
    'MaxFunctionEvaluations',500, ...
    'FunctionTolerance',1e-10, ...
    'StepTolerance',1e-10, ...
    'OptimalityTolerance',1e-8);

[q,resnorm,~,exitflag,output,~,jacobian] = ...
    lsqnonlin(residualFun,q0,qlb,qub,options);
p = pOffset + q.*pScale;

xErr = p(1);
zErr = p(2);
c1_est = round(p(3));

tPredicted = predict_arrival_times(p,wallX,R,standoff, ...
    s_axis,TXdelays,t0);
timeResidual = tPredicted-tObserved;
rmsResidualSamples = sqrt(mean((timeResidual(valid)/dt).^2));

boundTolerance = [0.05e-3,0.05e-3,0.5];
atLowerBound = p-lb <= boundTolerance;
atUpperBound = ub-p <= boundTolerance;
boundaryHit = any(atLowerBound | atUpperBound);

singularValues = svd(full(jacobian),'econ');
if isempty(singularValues) || singularValues(end) == 0
    jacobianCondition = Inf;
else
    jacobianCondition = singularValues(1)/singularValues(end);
end

flag = exitflag > 0 && ~boundaryHit && isfinite(rmsResidualSamples);

diagnostics = struct();
diagnostics.parameters = p;
diagnostics.tObserved = tObserved;
diagnostics.tPredicted = tPredicted;
diagnostics.timeResidual = timeResidual;
diagnostics.peakConfidence = peakConfidence;
diagnostics.rmsResidualSamples = rmsResidualSamples;
diagnostics.resnorm = resnorm;
diagnostics.exitflag = exitflag;
diagnostics.output = output;
diagnostics.boundaryHit = boundaryHit;
diagnostics.jacobianSingularValues = singularValues;
diagnostics.jacobianCondition = jacobianCondition;
diagnostics.wallX = wallX;

fprintf(['Joint 3-angle fit: z_error = %.4f mm, x_error = %.4f mm, ' ...
    'c1 = %.3f m/s\n'],zErr*1e3,xErr*1e3,p(3));
fprintf('RMS arrival residual = %.3f samples, scaled Jacobian condition = %.3g\n', ...
    rmsResidualSamples,jacobianCondition);

end

function [tObserved,confidence] = pick_arrival_times(RFdata,TXdelays, ...
    t_axis,dt,c0,minDepth,maxDepth,peakHeightFrac)

[Nt,Ns,Na] = size(RFdata);
tObserved = nan(Ns,Na);
confidence = zeros(Ns,Na);

startOffset = round(2*minDepth/c0/dt);
windowLength = round(2*maxDepth/c0/dt);

for ia = 1:Na
    channelEnvelope = abs(hilbert(RFdata(:,:,ia)));
    for is = 1:Ns
        firstIndex = round(TXdelays(is,ia)/dt) + startOffset;
        firstIndex = max(1,min(Nt,firstIndex));
        lastIndex = min(Nt,firstIndex+windowLength);
        trace = channelEnvelope(firstIndex:lastIndex,is);
        traceMax = max(trace);
        if traceMax <= 0
            continue;
        end

        [peakValue,location] = findpeaks(trace,'NPeaks',1, ...
            'MinPeakHeight',peakHeightFrac*traceMax);
        if isempty(location)
            [peakValue,location] = max(trace);
        end

        sampleIndex = firstIndex+location(1)-1;
        tObserved(is,ia) = t_axis(sampleIndex);
        confidence(is,ia) = peakValue(1)/traceMax;
    end
end

end

function residual = joint_residual(q,pOffset,pScale,wallX,R,standoff, ...
    s_axis,TXdelays,t0,tObserved,valid,weights,dt)

p = pOffset + q.*pScale;
tPredicted = predict_arrival_times(p,wallX,R,standoff, ...
    s_axis,TXdelays,t0);
sampleResidual = (tPredicted(valid)-tObserved(valid))/dt;

% Squaring this transformed residual is equivalent to a Huber loss.
% huberDelta = 2;
% absResidual = abs(sampleResidual);
% residual = sampleResidual;
% large = absResidual > huberDelta;
% residual(large) = sign(sampleResidual(large)).* ...
%     sqrt(2*huberDelta*absResidual(large)-huberDelta^2);
residual = weights.*sampleResidual;

end

function tPredicted = predict_arrival_times(p,wallX,R,standoff, ...
    s_axis,TXdelays,t0)

xErr = p(1);
zErr = p(2);
c1 = p(3);

wallArgument = R^2-(wallX+xErr).^2;
if any(wallArgument <= 0)
    tPredicted = inf(numel(s_axis),size(TXdelays,2));
    return;
end

wallZ = sqrt(wallArgument)-R+standoff+zErr;
distance = hypot(wallX-s_axis.',wallZ);

Ns = numel(s_axis);
Na = size(TXdelays,2);
tPredicted = zeros(Ns,Na);

for ia = 1:Na
    transmitArrival = min(distance/c1+TXdelays(:,ia).',[],2);
    tPredicted(:,ia) = min(transmitArrival+distance/c1,[],1).' + t0;
end

end
