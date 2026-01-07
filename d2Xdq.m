function d2Xdq = d2Xdq(taug, dq, M1, normTau)


%ala peskin lecture 3 immersed boundary notes
zeroColumn = zeros(size(taug,1), 1); 
% Ap12 = [Xg(:, 2 : M1) - Xg(:, 1 : M1-1), zeroColumn - Xg(:, M1)];
% Am12 = [Xg(:, 1) - zeroColumn, Xg(:, 2 : M1) - Xg(:, 1 : M1-1)];
A12 = (normTau ./ dq - 1) .* (taug ./ normTau);
B = [A12, zeroColumn] - [zeroColumn, A12];

d2Xdq = B ./ dq;