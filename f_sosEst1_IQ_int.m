%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Speed of sound (SoS) estimation in liquid using DAS-beamformed data
% Authors: Daler Rakhmatov (daler@uvic.ca), Haoyun Chen
% Date: 9 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[c1_est] = f_sosEst1_IQ_int(params,c1_sweep)

pipe_IR = params.pipe_IR;               % pipe's inner radius [m]
standoff = params.standoff;             % nominal standoff [m]

IQdata = params.IQdata;
fc = params.fc;
Na = params.Na;
TXdelays = params.TXdelays;

Ns = params.Ns;                         % number of sensor channels
s_axis = params.s_axis;

Nx = params.Nx;                         % number of lateral grid points
x_axis = params.x_axis; 
xIndRange = params.xIndRange;

Nt = params.Nt;                         % number of temporal samples
dt = params.dt;
t0 = params.t0;                 

Nz = params.Nz;                         % number of axial grid points
dz = params.dz;
z_axis = params.z_axis;

span = params.span;


%% Perform DAS beamforming

fprintf('\nEstimating speed of sound in liquid...\n');

% axial index range for beamforming
front_dist = sqrt(pipe_IR^2 - x_axis(xIndRange).^2) - pipe_IR + standoff;

z_min = min(front_dist) - span;
z_max = max(front_dist) + span;
zIndRange = (round(z_min/dz):round(z_max/dz));


% Sweep across selected c1 values

Nc1 = numel(c1_sweep);                  % number of candidate SoS values
css1 = zeros(1,Nc1);                    % cosine similarity sum

for ic = 1:Nc1

    c1 = c1_sweep(ic);

    %-- Beamform zero-angle data

    ia = round((1+Na)/2);
    Dts = IQdata(1:Nt,1:Ns,ia);         % IQ channel data frame to process

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
    
    BFdata0 = Dzx(zIndRange,xIndRange);
    BFvec0 = imag(BFdata0(:));


    %-- Beamform negative/positive-angle data

    aIndSel = [1 Na];
    for ia = aIndSel

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

        BFdata = Dzx(zIndRange,xIndRange);
        BFvec = imag(BFdata(:));

        % update cosine similarity sum for vectorized (stacked) data
        css1(ic) = css1(ic) + ...
            dot(BFvec,BFvec0)/(norm(BFvec,2)*norm(BFvec0,2));

    end

    % normalize cosine similarity sum
    css1(ic) = css1(ic)/numel(aIndSel);     

end


% interpolate cosine similarity sum curve
c1_int = c1_sweep(1):1:c1_sweep(end);
css1_int = interp1(c1_sweep,css1,c1_int,'pchip');

% find SoS estimate
[~,loc] = max(css1_int);
c1_est = c1_int(loc);

figure; plot(c1_int,css1_int); grid on;
set(gca,'FontName','Courier','FontSize',14,'FontWeight','demi');
xlabel('c1 (m/s)','FontName','Courier','FontSize',14,'FontWeight','bold');
title('Cosine Similarity Sum','FontName','Courier','FontSize',14,'FontWeight','bold');


end