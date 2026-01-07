function d4Xdq = d4Xdq(Xg, dq, M1)

zeroColumn = zeros(size(Xg,1), 1); 
Cg = (Xg(:, 3 : M1) - 2 .* Xg(:, 2 : M1 - 1) + Xg(:, 1 : M1 -2))./ dq.^2; 

d4Xdq = [zeroColumn, zeroColumn, Cg] - 2 .* [zeroColumn, Cg, zeroColumn] ... 
    + [Cg, zeroColumn, zeroColumn]; 

d4Xdq = d4Xdq ./ dq.^2;