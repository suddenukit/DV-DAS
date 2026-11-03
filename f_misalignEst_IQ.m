%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transducer misalignment estimation using DAS-beamformed envelope data
% Ideal pipe center location in transducer coordinates:
%   (xc, zc) = (0, standoff-pipe_IR)
% Estimated pipe center location in transducer coordinates: 
%   (x_est, z_est) = (xc+xErr, zc+zErr)
% Inputs:
%   thres - energy threshold to flag missing frontwall
%   params - initialization and beamforming parameters
% Author(s): Daler Rakhmatov, Haoyun Chen
% Contact: daler@ece.uvic.ca, Date: 8 January 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[xErr,zErr,flag] = f_misalignEst_IQ(params,thres)

pipe_IR = params.pipe_IR;               % pipe's inner radius [m]
standoff = params.standoff;             % nominal standoff [m]
c1 = params.c1;                         % SoS in layer1 [m/s]

IQdata = params.IQdata;
fc = params.fc;
Na = params.Na;
TXangles = params.TXangles;
TXdelays = params.TXdelays;

Ns = params.Ns;                         % number of sensor channels
ds = params.ds;
s_axis = params.s_axis;

Nx = params.Nx;                         % number of lateral grid points
dx = params.dx;
x_axis = params.x_axis; 
xIndRange = params.xIndRange;


Nt = params.Nt;                         % number of temporal samples
dt = params.dt;
t0 = params.t0;                 

Nz = params.Nz;                         % number of axial grid points
dz = params.dz;
z_axis = params.z_axis;

zOffset = params.zOffset;

xErr0 = params.xErr0;
zErr0 = params.zErr0;
span = params.span;


%% Perform DAS beamforming

fprintf('\nEstimating transducer misalignment...');

% axial index range for beamforming
front_dist = sqrt(pipe_IR^2 - x_axis(xIndRange).^2) - pipe_IR + standoff;

z_min = min(front_dist) - span;
z_max = max(front_dist) + span;

zIndRange = (round(z_min/dz):round(z_max/dz));


NzBF = numel(zIndRange);
NxBF = numel(xIndRange);

BFdata = zeros(NzBF,NxBF);


for ia = 1:Na

    Dts = IQdata(1:Nt,1:Ns,ia);     % IQ channel data frame to process

    % Calculate TX/RX delays
    delTX = zeros(Nz,Nx);
    delRX = zeros(Nz,Nx,Ns);
    for iz = zIndRange
        for ix = xIndRange
            xDiff = s_axis - x_axis(ix);
            sDelay = sqrt(z_axis(iz)^2 + xDiff.^2)/c1;
            delTX(iz,ix) = min(sDelay + TXdelays(:,ia));
            delRX(iz,ix,:) = sDelay;
        end
    end

    % Apply delays and sum
    Dzx = zeros(Nz,Nx);
    for iz = zIndRange
        for ix = xIndRange
            snapshot = zeros(1,Ns);
            for is = 1:Ns
                tInd = (delTX(iz,ix) + delRX(iz,ix,is) + t0)/dt + 1;
                tau = delTX(iz,ix) + delRX(iz,ix,is) + t0;
                tInd_int = floor(tInd);
                tInd_frac = tInd - tInd_int;
                snapshot(is) = (Dts(tInd_int,is) + ...
                    (Dts(tInd_int+1,is) - Dts(tInd_int,is)) * tInd_frac) ...
                    * exp(1i*2*pi*fc*tau);
            end
            Dzx(iz,ix) = sum(snapshot)/Ns;
        end
    end

    BFdata = BFdata + Dzx(zIndRange,xIndRange);     % coherent compounding

end


%% Estimate errors of pipe's center coordinates

energy = sum(BFdata(:).*conj(BFdata(:)));
if energy > thres, flag = true;
else, flag = false; end

env = abs(BFdata);
[~,ind] = max(env);

z_data = z_axis(ind).' + z_axis(zIndRange(1));

% -- Apply nonlinear regression
beta0 = [zErr0 xErr0];
options = statset('RobustWgtFun','huber');
beta = nlinfit(x_axis(xIndRange), z_data, @frontwallCurve, beta0, options);
zErr = beta(1);
xErr = beta(2);


% fw = sqrt(pipe_IR^2 - (x_axis(xIndRange)+xErr).^2) - pipe_IR + (standoff+zErr);
% figure; imagesc(x_axis(xIndRange),z_axis(zIndRange),abs(BFdata)); grid on;
% hold on; plot(x_axis(xIndRange),fw,'k:','LineWidth',2); hold off;
% figure; imagesc(x_axis(xIndRange),z_axis(zIndRange),imag(BFdata)); grid on;
% hold on; plot(x_axis(xIndRange),fw,'k:','LineWidth',2); hold off;


fprintf('\nZ error est. = %.2f mm,  X error est. = %.2f mm\n',zErr*1e3,xErr*1e3);


%% Curvature model function
function z_model = frontwallCurve(coeff,x)
z_error = coeff(1);
x_error = coeff(2);
z_model = sqrt(pipe_IR^2 - (x+x_error).^2) - pipe_IR + (standoff+z_error);
end


end