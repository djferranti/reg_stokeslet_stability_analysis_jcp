function d2Xds = d2Xds(taug, ds, M2, normTau)

zeroRow = zeros(1, size(taug, 2));
% Ap12 = [Xg(2 : M2, :) - Xg(1 : M2 - 1, :) ; zeroRow - Xg(M2, :)]; 
% Am12 = [Xg(1, :) - zeroRow; Xg(2 : M2, :) - Xg(1 : M2 - 1, :)]

A12 = (normTau ./ ds - 1) .* (taug ./ normTau);
B = [A12; zeroRow] - [zeroRow; A12];

d2Xds = B ./ ds;
% d2Xds = (Ap12 - Am12) ./ (ds .^ 2);

end