function [angles] = calculate_angles_from_delays(params, frontDist, zErr, xErr, ia)

pipe_IR = params.pipe_IR;               % pipe's inner radius [m]
TXdelays = params.TXdelays;
c1 = params.c1;
s_axis = params.s_axis;
x_axis = params.x_axis; 
zOffset = params.zOffset;
Nx = params.Nx;

skip = round(0.25*Nx);
xIndRange = ((1+skip):(Nx-skip))';
Nxx = numel(xIndRange);

angles = zeros(Nxx,1);
for ix = 1:Nxx
    dist = sqrt(frontDist(xIndRange(ix))^2 + (s_axis -  x_axis(xIndRange(ix))).^2);
    time = TXdelays(:,ia) + dist / c1;
    [~, loc] = min(time);
    s2center = sqrt((zOffset-zErr)^2 + (s_axis(loc)+xErr)^2);
    angles(ix) = acos((dist(loc).^2 + pipe_IR.^2 - s2center.^2) / (2 * dist(loc) * pipe_IR));
end


end