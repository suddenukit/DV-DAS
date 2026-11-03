%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Longitudinal SoS estimation in steel using DAS-beamformed data
% Authors: Daler Rakhmatov (daler@uvic.ca), Haoyun Chen
% Date: 10 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[c2L_est] = f_sosEst2L_IQ_int(params,c2L_sweep)

pipe_IR = params.pipe_IR;               % pipe's inner radius [m]
pipe_OR = params.pipe_OR;               % pipe's outer radius [m]
standoff = params.standoff;             % nominal standoff [m]

zErr = params.zErr;
xErr = params.xErr;

IQdata = params.IQdata;
c1 = params.c1;

Na = params.Na;
xv = params.xv;
zv = params.zv;
fc = params.fc;

Ns = params.Ns;                         % number of sensor channels
ds = params.ds;
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

zOffset = params.zOffset;

lambda2 = c2L_sweep(end)/fc; 

layer1 = standoff + zErr;

% frontwall and backwall boundaries
frontDist = sqrt(pipe_IR^2 - (x_axis+xErr).^2) - pipe_IR + layer1;
backDist = sqrt(pipe_OR^2 - (x_axis+xErr).^2) - pipe_IR + layer1;

alpha1 = asin((x_axis+xErr)./pipe_IR);  % frontwall normal angle [rad]


%% Perform DAS beamforming

fprintf('\nEstimating longitudinal speed of sound in steel...\n');

bwIndMin = round((backDist - lambda2/2)/dz);
bwIndMax = round((backDist + lambda2/2)/dz);

BWdata = zeros(Nz,Nx);
for ix = xIndRange
    BWdata(bwIndMin(ix):bwIndMax(ix),ix) = 1;
end
zIndMin = round((backDist - 2*lambda2)/dz);
zIndMax = round((backDist + 2*lambda2)/dz);

% Sweep across selected c2L values
Nc2L = numel(c2L_sweep);                % number of candidate SoS values
BWerr = zeros(1,Nc2L);                  % backwall error

for ic = 1:Nc2L

    c2L = c2L_sweep(ic);

    %-- Beamform nonzero-angle data
    
    BFdata = zeros(Nz,Nx);
    ia = round((1+Na)/2);
    Dts = IQdata(1:Nt,1:Ns,ia);     % IQ channel data frame to process
    z_newO = zv(ia);                % virtual source z-offset (-)
    x_newO = xv(ia);                % virtual source x-offset (+)
    tx_distance = sqrt((zOffset - z_newO).^2 + (s_axis - x_newO).^2);
    tx_distance_min = min(tx_distance);

    % frontwall boundary (as seen by DW)
    boundary1 = sqrt((frontDist + zOffset - z_newO).^2 + (x_axis - x_newO).^2);
    alpha_newO1 = asin((x_axis - x_newO)./boundary1);


    % estimate of incident angle (as seen by DW)
    theta1s = calculate_angles_from_delays(params, frontDist, zErr, xErr, ia);
    theta1 = mean(theta1s);
    if x_newO + xErr > 0
        theta1 = -theta1;
    end
    sinTheta1 = sin(theta1);
    % refraction angle for layer2
    sinTheta2 = c2L/c1*sinTheta1;
    if abs(sinTheta2) >= 1
        fprintf('\nERROR: Refraction angle is out of bounds...\n\n');
        return;
    end
    theta2 = asin(sinTheta2);
    % Calculate TX/RX delays
    delTX = zeros(Nz,Nx);
    delRX = zeros(Nz,Nx,Ns);
    for ix = xIndRange
        for iz = zIndMin(ix):zIndMax(ix)

            rTX = sqrt((z_axis(iz) + zOffset - z_newO)^2 + (x_axis(ix) - x_newO)^2);
            a = asin((x_axis(ix) - x_newO)/rTX);

            for ib = 1:Nx
                if alpha_newO1(ib) >= a
                    break;
                end
            end

            tx_path1 = boundary1(ib) - tx_distance_min;
            tx_path2 = (rTX - boundary1(ib))*cos(theta2-theta1);
            delTX(iz,ix) = tx_path1/c1 + tx_path2/c2L;

            for ifw = 1:Nx
                a_fw = alpha1(ifw);
                x_fw = x_axis(ifw);
                z_fw = frontDist(ifw);

                delta_z = z_axis(iz) - z_fw;
                delta_x = x_fw - x_axis(ix);
                rx_path2 = sqrt(delta_z^2 + delta_x^2);

                delta_a = abs(atan(delta_x/delta_z));
                if delta_x >= 0
                    rx_angle2 = delta_a + a_fw;
                    rx_angle1 = asin(c1/c2L*sin(rx_angle2)) - a_fw;
                    rx_path1 = z_fw/cos(rx_angle1);
                    x_s = x_fw + z_fw*tan(rx_angle1);
                    x_sInd = round((x_s-s_axis(1))/ds) + 1;
                    if x_sInd >= 1 && x_sInd <= Ns
                        ttRX = rx_path2/c2L + rx_path1/c1;
                        delRX(iz,ix,x_sInd) = ttRX;
                    end
                else
                    rx_angle2 = delta_a - a_fw;
                    rx_angle1 = asin(c1/c2L*sin(rx_angle2)) + a_fw;
                    rx_path1 = z_fw/cos(rx_angle1);
                    x_s = x_fw - z_fw*tan(rx_angle1);
                    x_sInd = round((x_s-s_axis(1))/ds) + 1;
                    if x_sInd >= 1 && x_sInd <= Ns
                        ttRX = rx_path2/c2L + rx_path1/c1;
                        delRX(iz,ix,x_sInd) = ttRX;
                    end
                end

            end

        end

    end

    % Apply delays and sum
    Dzx = zeros(Nz,Nx);
    for ix = xIndRange
        for iz = zIndMin(ix):zIndMax(ix)
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
    BFdata = BFdata + Dzx;
    env = abs(BFdata);
    env = (env-min(env(:)))/(max(env(:))-min(env(:)));

    % record L1 error for normalized envelope
    BWerr(ic) = sum(abs(env(:)-BWdata(:)));     
end


% interpolate L1 error curve
c2L_int = c2L_sweep(1):1:c2L_sweep(end);
BWerr_int = interp1(c2L_sweep,BWerr,c2L_int,'pchip');

% find SoS estimate
[~,loc] = min(BWerr_int);
c2L_est = c2L_int(loc);

figure; plot(c2L_int, BWerr_int); xlim([c2L_int(1) c2L_int(end)]); grid on;
set(gca,'FontName','Courier','FontSize',14,'FontWeight','demi');
xlabel('c2L (m/s)','FontName','Courier','FontSize',14,'FontWeight','bold');
title('L1 Error','FontName','Courier','FontSize',14,'FontWeight','bold');


end