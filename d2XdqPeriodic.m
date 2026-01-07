function [d2Xdq, d4Xdq] = d2XdqPeriodic(Xg, dq, L, M1, p, varargin)

bending = 0; 
if nargin > 5
    bending = varargin{1};
end

d2Xgdq = Xg(:, 3 : M1) - 2 .* Xg(:, 2 : M1 - 1) + Xg(:, 1 : M1 - 2);
d2Xgdq = d2Xgdq ./ ( dq.^2 );

d2XgdqFirst = (Xg(:, 2) - 2 .* Xg(:, 1) + Xg(:, M1));
if p
    d2XgdqFirst = d2XgdqFirst - L;
end
d2XgdqFirst = d2XgdqFirst ./ (dq .^ 2);

d2XgdqLast = (Xg(:, 1) - 2 .* Xg(:, M1) + Xg(:, M1 - 1));
if p
    d2XgdqLast = d2XgdqLast + L;
end
d2XgdqLast = d2XgdqLast ./ (dq .^ 2);

d2Xdq = [d2XgdqFirst d2Xgdq d2XgdqLast];

if bending %this is assuming periodicity
    Cg = [d2XgdqLast d2Xdq d2XgdqFirst];
    Cg = Cg(:, 3 : end) - 2 .* Cg(:, 2 : end-1) + Cg(:, 1 : end - 2) ;
    d4Xdq = Cg ./ (dq.^2);
end

end