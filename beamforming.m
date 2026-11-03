function [data] = beamforming(params, data, geom, delays, shift)
arguments
    params
    data
    geom
    delays
    shift = false;
end


Na = params.Na;
Nt = params.Nt;
Ns = params.Ns;
Nx = params.Nx;
Nz = params.Nz;
dt = params.dt;
t0 = params.t0;
fc = params.fc;
zv = params.zv;
xv = params.xv;
xIndRange = params.xIndRange;
zIndRange = params.zIndRange;
zOffset = params.zOffset;
layer2 = params.layer2;
s_axis = params.s_axis;
x_axis = params.x_axis;
z_axis = params.z_axis;
IQdata = params.IQdata;
apodWinRX = params.apodWinRX;

TXangles = data.TXangles;
c1 = data.c1;
c2L = data.c2L;
c2T = data.c2T;

frontDist = geom.frontDist;
backDist = geom.backDist;
pipeMask = geom.pipeMask;
zErr = geom.zErr;
xErr = geom.xErr;
pipe_IR = geom.pipe_IR;

delayRX_dL = delays.delayRX_dL;
delayRX_hL = delays.delayRX_hL;
delayRX_dT = delays.delayRX_dT;
delayRX_hT = delays.delayRX_hT;

NzBF = numel(zIndRange);
NxBF = numel(xIndRange);

BFdata_LdL = zeros(NzBF,NxBF);
BFdata_LLdLL = zeros(NzBF,NxBF);

BFdata_TdT_neg = zeros(NzBF,NxBF);
BFdata_TdT_pos = zeros(NzBF,NxBF);

BFdata_TTdT_neg = zeros(NzBF,NxBF);
BFdata_TTdT_pos = zeros(NzBF,NxBF);

BFdata_TdTT_neg = zeros(NzBF,NxBF);
BFdata_TdTT_pos = zeros(NzBF,NxBF);

BFdata_TTdTT_neg = zeros(NzBF,NxBF);
BFdata_TTdTT_pos = zeros(NzBF,NxBF);


for ia = 1:Na

    Dts = IQdata(1:Nt,1:Ns,ia);     % IQ channel data to process
    
    thetaTX = TXangles(ia);         % DW transmit angle

    z_newO = zv(ia);                % virtual source z-offset (-)
    x_newO = xv(ia);                % virtual source x-offset (+)

    % transmit delays estimate (to compare with 'TXdelays')
    tx_distance = sqrt((zOffset - z_newO).^2 + (s_axis - x_newO).^2);
    tx_distance_min = min(tx_distance);
    tx_delays_est = (tx_distance - tx_distance_min)/c1;
    % 
    % 
    % SoS for layer2
    if abs(thetaTX) < deg2rad(5)
        c2 = c2L;
        L_mode = true;
    else
        c2 = c2T;
        L_mode = false;
    end


    % frontwall/backwall boundary (as seen by DW)
    boundary1 = sqrt((frontDist + zOffset - z_newO).^2 + (x_axis - x_newO).^2);
    alpha_newO1 = asin((x_axis - x_newO)./boundary1);
     
    boundary2 = sqrt((backDist + zOffset - z_newO).^2 + (x_axis - x_newO).^2);
    alpha_newO2 = asin((x_axis - x_newO)./boundary2);


    % estimate of incident angle (as seen by DW)

    theta1s = calculate_angles_from_delays(params, frontDist, zErr, xErr, ia);
    theta1 = mean(theta1s);
    if x_newO + xErr > 0
        theta1 = -theta1;
    end
    sinTheta1 = sin(theta1);
    % refraction angle for layer2
    sinTheta2 = c2/c1*sinTheta1;
    if (abs(sinTheta2)) >= 1
        fprintf('\nERROR: Refraction angle is out of bounds...\n\n');
        return;
    end
    theta2 = asin(sinTheta2);

    % % record estimated angles
    % data.theta1(ia) = theta1;
    % data.theta2(ia) = theta2;
    % 
    %%% Calculate transmit/receive path delays
    fprintf('\nTX angle = %.1f deg\n',...
        rad2deg(TXangles(ia)));
    
    fprintf('Calculating TX delays...\n');

    delayTX_d = zeros(Nz,Nx);           % direct incoming (TX) delay
    delayTX_h = zeros(Nz,Nx);           % halfskip incoming (TX) delay

    layerMask = zeros(Nz,Nx);
    for iz = zIndRange
        for ix = xIndRange
            xDiff = s_axis - x_axis(ix);
            rTX = sqrt((z_axis(iz) + zOffset - z_newO)^2 + (x_axis(ix) - x_newO)^2);
            if z_axis(iz) <= frontDist(ix)
                % inside layer1
                layerMask(iz,ix) = 1;
                delayTX_d(iz,ix) = (rTX - tx_distance_min)/c1;
                delayTX_h(iz,ix) = delayTX_d(iz,ix);

            else
                % outside layer1
                layerMask(iz,ix) = 2;
                a = asin((x_axis(ix) - x_newO)/rTX);
                ta = tan(a);
                beta = x_axis(ix) - ta * (z_axis(iz) + zOffset - zErr);
                A = ta^2 + 1;
                B = 2 * (ta^2*xErr - beta);
                C = beta^2 + ta^2*(xErr^2 - pipe_IR^2);
                % solve for Ax^2 + Bx + C = 0
                discriminant = B^2 - 4*A*C;
                assert(discriminant >= 0, 'No real intersection.');                
                x1 = (-B + sqrt(discriminant))/(2*A);
                x2 = (-B - sqrt(discriminant))/(2*A);
                xCandidate = [x1, x2];  % Two real roots
                valid = xCandidate >= min(x_axis(ix), x_newO) & ...
                        xCandidate <= max(x_axis(ix), x_newO);            
                assert(nnz(valid) == 1, 'Expected one intersection between V and G.');
                xb = xCandidate(valid);
                zb = sqrt(pipe_IR^2 - (xb + xErr)^2) - zOffset + zErr;
                % for ib = 1:Nx
                %     if alpha_newO1(ib) >= a
                %         break;
                %     end
                % end
                boundaryDist = sqrt((zb + zOffset - z_newO).^2 + (xb - x_newO).^2);

                tx_path1 = boundaryDist - tx_distance_min;

                % direct TX path/delay
                tx_path2_d = (rTX - boundaryDist)*cos(theta2-theta1);
                delayTX_d(iz,ix) = tx_path1/c1 + tx_path2_d/c2;

                % halfskip TX path/delay
                hs_path = 2*layer2*cos(theta2);
                tx_path2_h = hs_path - (rTX - boundaryDist)*cos(theta2+theta1);
                delayTX_h(iz,ix) = tx_path1/c1 + tx_path2_h/c2;
    
            end
        end
    end


    fprintf('DAS beamforming...\n');

    %%% DAS beamforming: zero TX angle
    if ~(ia == 1 || ia == Na)

        % LdL mode
        Dzx = zeros(Nz,Nx);
        for iz = zIndRange
            for ix = xIndRange
                if layerMask(iz,ix) > 0
                    snapshot = zeros(1,Ns);
                    for is = 1:Ns                            
                        tau = delayTX_d(iz,ix) + delayRX_dL(iz,ix,is) + t0;
                        tInd = tau/dt + 1;
                        tInd_int = floor(tInd);
                        tInd_frac = tInd - tInd_int;
                        snapshot(is) = (Dts(tInd_int,is) + ...
                            (Dts(tInd_int+1,is) - Dts(tInd_int,is)) * tInd_frac);
                        if shift
                            snapshot(is) = snapshot(is) * exp(1i*2*pi*fc*tau);
                        end
                    end
                    % DAS output without coherence factor
                    Dzx(iz,ix) = sum(snapshot.*apodWinRX)/Ns;
                end
            end
        end
        % record beamformed data
        BFdata_LdL = Dzx(zIndRange,xIndRange);

        % LLdLL mode
        Dzx = zeros(Nz,Nx);
        for iz = zIndRange
            for ix = xIndRange
                if layerMask(iz,ix) > 0
                    snapshot = zeros(1,Ns);
                    for is = 1:Ns
                        tau = delayTX_h(iz,ix) + delayRX_hL(iz,ix,is) + t0;
                        tInd = tau/dt + 1;
                        tInd_int = floor(tInd);
                        tInd_frac = tInd - tInd_int;
                        snapshot(is) = Dts(tInd_int,is) + ...
                            (Dts(tInd_int+1,is) - Dts(tInd_int,is)) * tInd_frac;
                        if shift
                            snapshot(is) = snapshot(is) * exp(1i*2*pi*fc*tau);
                        end
                    end
                    % DAS output without coherence factor
                    Dzx(iz,ix) = sum(snapshot.*apodWinRX)/Ns;
                end
            end
        end
        % record beamformed data
        BFdata_LLdLL = Dzx(zIndRange,xIndRange);        
    end

    %%% DAS beamforming: non-zero TX angles
    if ia == 1 || ia == Na

        % TdT mode
        Dzx = zeros(Nz,Nx);
        for iz = zIndRange
            for ix = xIndRange    
                if layerMask(iz,ix) > 0    
                    snapshot = zeros(1,Ns);
                    for is = 1:Ns                            
                        tau = delayTX_d(iz,ix) + delayRX_dT(iz,ix,is) + t0;
                        tInd = tau/dt + 1;
                        tInd_int = floor(tInd);
                        tInd_frac = tInd - tInd_int;
                        snapshot(is) = (Dts(tInd_int,is) + ...
                            (Dts(tInd_int+1,is) - Dts(tInd_int,is)) * tInd_frac);
                        if shift
                            snapshot(is) = snapshot(is) * exp(1i*2*pi*fc*tau);
                        end
                    end
                    % DAS output without coherence factor
                    Dzx(iz,ix) = sum(snapshot.*apodWinRX)/Ns;
                end
            end
        end
        % record beamformed data
        if ia == 1, BFdata_TdT_neg = Dzx(zIndRange,xIndRange);
        elseif ia == Na, BFdata_TdT_pos = Dzx(zIndRange,xIndRange);
        end


        % TTdT mode
        Dzx = zeros(Nz,Nx);
        for iz = zIndRange
            for ix = xIndRange
                if layerMask(iz,ix) > 0    
                    snapshot = zeros(1,Ns);
                    for is = 1:Ns                            
                        tau = delayTX_h(iz,ix) + delayRX_dT(iz,ix,is) + t0;
                        tInd = tau/dt + 1;
                        tInd_int = floor(tInd);
                        tInd_frac = tInd - tInd_int;                           
                        snapshot(is) = Dts(tInd_int,is) + ...
                            (Dts(tInd_int+1,is) - Dts(tInd_int,is)) * tInd_frac;
                        if shift
                            snapshot(is) = snapshot(is) * exp(1i*2*pi*fc*tau);
                        end
                    end

                    % DAS output without coherence factor
                    Dzx(iz,ix) = sum(snapshot.*apodWinRX)/Ns;

                end

            end
        end
        % record beamformed data
        if ia == 1, BFdata_TTdT_neg = Dzx(zIndRange,xIndRange);
        elseif ia == Na, BFdata_TTdT_pos = Dzx(zIndRange,xIndRange);
        end

        % TdTT mode
        Dzx = zeros(Nz,Nx);
        for iz = zIndRange
            for ix = xIndRange
                if layerMask(iz,ix) > 0
                    snapshot = zeros(1,Ns);
                    for is = 1:Ns
                        tau = delayTX_d(iz,ix) + delayRX_hT(iz,ix,is) + t0;
                        tInd = tau/dt + 1;
                        tInd_int = floor(tInd);
                        tInd_frac = tInd - tInd_int;
                        snapshot(is) = Dts(tInd_int,is) + ...
                            (Dts(tInd_int+1,is) - Dts(tInd_int,is)) * tInd_frac;
                        if shift
                            snapshot(is) = snapshot(is) * exp(1i*2*pi*fc*tau);
                        end
                    end
                    % DAS output without coherence factor
                    Dzx(iz,ix) = sum(snapshot.*apodWinRX)/Ns;
                end
            end
        end
        % record beamformed data
        if ia == 1, BFdata_TdTT_neg = Dzx(zIndRange,xIndRange);
        elseif ia == Na, BFdata_TdTT_pos = Dzx(zIndRange,xIndRange);
        end

        % TTdTT mode
        Dzx = zeros(Nz,Nx);
        for iz = zIndRange
            for ix = xIndRange
                if layerMask(iz,ix) > 0
                    snapshot = zeros(1,Ns);
                    for is = 1:Ns
                        tau = delayTX_h(iz,ix) + delayRX_hT(iz,ix,is) + t0;
                        tInd = tau/dt + 1;
                        tInd_int = floor(tInd);
                        tInd_frac = tInd - tInd_int;
                        snapshot(is) = Dts(tInd_int,is) + ...
                            (Dts(tInd_int+1,is) - Dts(tInd_int,is)) * tInd_frac;
                        if shift
                            snapshot(is) = snapshot(is) * exp(1i*2*pi*fc*tau);
                        end
                    end
                    % DAS output without coherence factor
                    Dzx(iz,ix) = sum(snapshot.*apodWinRX)/Ns;
                end
            end
        end
        % record beamformed data
        if ia == 1, BFdata_TTdTT_neg = Dzx(zIndRange,xIndRange);
        elseif ia == Na, BFdata_TTdTT_pos = Dzx(zIndRange,xIndRange);
        end

    end   
end
data.LdL = BFdata_LdL;                  % zero-angle LdL data
data.LLdLL = BFdata_LLdLL;              % zero-angle LLdLL data
data.TdT_neg = BFdata_TdT_neg;          % negative-angle TdT data
data.TdT_pos = BFdata_TdT_pos;          % positive-angle TdT data
data.TdTT_neg = BFdata_TdTT_neg;        % negative-angle TdTT data
data.TdTT_pos = BFdata_TdTT_pos;        % positive-angle TdTT data
data.TTdT_neg = BFdata_TTdT_neg;        % negative-angle TTdT data
data.TTdT_pos = BFdata_TTdT_pos;        % positive-angle TTdT data
data.TTdTT_neg = BFdata_TTdTT_neg;      % negative-angle TTdTT data
data.TTdTT_pos = BFdata_TTdTT_pos;      % positive-angle TTdTT data
end