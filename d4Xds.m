function d4Xds = d4Xds(Xg, ds, M2) 

zeroRow = zeros(1, size(Xg, 2));
Cg = (Xg(3 : M2, :) - 2 .* Xg(2 : M2-1, :) + Xg(1 : M2 - 2, :)) ./ ds.^2 ;


d4Xds = [Cg; zeroRow; zeroRow] - 2 .* [zeroRow; Cg; zeroRow] ...
    + [zeroRow; zeroRow; Cg];
d4Xds = d4Xds ./ ds.^2;

end