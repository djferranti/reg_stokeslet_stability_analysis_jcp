function d4Xdqds = d4Xdqds(Xg, dq, ds, M1, M2) 

zeroColumn = zeros(size(Xg, 1), 1);
A12 = Xg(:, 2 : M1) - Xg(:, 1 : M1 -1);
d2Xdq = ( [A12, zeroColumn] - [zeroColumn, A12] ) ./ dq.^2;

zeroRow = zeros(1, size(Xg, 2));
A12 = d2Xdq(2 : M2, :) - d2Xdq(1 : M2 -1, :);

d4Xdqds = [A12; zeroRow] - [zeroRow; A12];
d4Xdqds = d4Xdqds ./ ds.^2;
end