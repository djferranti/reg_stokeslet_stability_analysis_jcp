function [d2Xds, d4Xds] = d2XdsPeriodic(Xg, ds, L, M2, p, varargin)

bending = 0;
if nargin > 5
    bending = varargin{1};
end

d2Xgds = Xg(3 : M2, :) - 2 .* Xg(2 : M2 - 1, :) + Xg(1 : M2 - 2, :);
d2Xgds = d2Xgds ./ ( ds.^2 );

d2XgdsFirst = (Xg(2, :) - 2 .* Xg(1, :) + Xg(M2, :));
if p
    d2XgdsFirst = d2XgdsFirst - L;
end
d2XgdsFirst = d2XgdsFirst ./ (ds .^ 2);

d2XgdsLast = (Xg(1, :) - 2 .* Xg(M2, :) + Xg(M2 - 1, :));
if p
    d2XgdsLast = d2XgdsLast + L;
end
d2XgdsLast = d2XgdsLast ./ (ds .^ 2);

d2Xds = [d2XgdsFirst; d2Xgds; d2XgdsLast];

if bending %this is assuming periodicity
    Cg = [d2XgdsLast; d2Xds; d2XgdsFirst];
    Cg = Cg(3 : end, :) - 2 .* Cg(2 : end-1, :) + Cg(1 : end - 2, :) ;
    d4Xds = Cg ./ (ds.^2);
end
end