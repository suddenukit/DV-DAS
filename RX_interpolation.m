function [delayRX] = RX_interpolation(delayRX, pipeMask)
% RX_SMOOTH interpolates the zero points in RX delay
% The function does not interpolate zero points in the sensors on each side
% because in that case the sensors on the side do receive signal from the
% point. 
Nz = size(delayRX, 1);
Nx = size(delayRX, 2);
Ns = size(delayRX, 3);
zIndRange = 1 : Nz;
xIndRange = 1 : Nx;
sIndRange = 1 : Ns;

for iz = zIndRange
    for ix = xIndRange
        if pipeMask(iz,ix) > 0
            start_pixel = 1;
            end_pixel = Ns;
            % Pass side sensors with zero data
            while ~isfinite(delayRX(iz, ix, start_pixel))
                start_pixel = start_pixel + 1;
            end
            while ~isfinite(delayRX(iz, ix, end_pixel))
                end_pixel = end_pixel - 1;
            end
            % scan from the starting sensor to the end sensor
            current_pixel = start_pixel;
            while current_pixel <= end_pixel
                % find the segment where data are zeroes and interpolate
                if ~isfinite(delayRX(iz, ix, current_pixel))
                    p = current_pixel;
                    while ~isfinite(delayRX(iz, ix, p))
                        p = p + 1;
                    end
                    x = [current_pixel - 1, p];
                    v = [delayRX(iz, ix, current_pixel - 1), delayRX(iz, ix, p)];
                    delayRX(iz, ix, current_pixel:p-1) = interp1(x, v, current_pixel:p-1);
                    current_pixel = p;
                end
                current_pixel = current_pixel + 1;
            end
        end
    end
end

winLength = 2;
for iz = zIndRange
    for ix = xIndRange
        if pipeMask(iz,ix) > 0
            for is = sIndRange
                if ~isfinite(delayRX(iz, ix, is))
                    xIndRangeLoc = max(xIndRange(1),ix-winLength):min(ix+winLength,xIndRange(end));
                    zIndRangeLoc = max(zIndRange(1),iz-winLength):min(iz+winLength,zIndRange(end));
                    regionLoc = delayRX(zIndRangeLoc,xIndRangeLoc,is);
                    delayRX(iz,ix,is) = median(regionLoc(:));
                end
            end

        end
    end
end

for iz = zIndRange
    for ix = xIndRange
        for is = sIndRange
            if ~isfinite(delayRX(iz, ix, is))
                delayRX(iz, ix, is) = 0;
            end
        end
    end
end
