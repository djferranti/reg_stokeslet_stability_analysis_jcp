function [Fq, Fs, Fz, dS] = computetensilebendingforces(K, Xq, Xs, Xz, ...
    M1, M2, dq, ds, L, varargin)
%tension in first material direction
p = 1; np =0;
bending = 0;
if nargin > 9
    bending = varargin{1};
    Kb = varargin{2};
end

if ~bending
    d2Xqdq = d2XdqPeriodic(Xq, dq, L, M1, p);
    d2Xsdq = d2XdqPeriodic(Xs, dq, L, M1, np);
    d2Xzdq = d2XdqPeriodic(Xz, dq, L, M1, np);

    d2Xqds = d2XdsPeriodic(Xq, ds, L, M2, np);
    d2Xsds = d2XdsPeriodic(Xs, ds, L, M2, p);
    d2Xzds = d2XdsPeriodic(Xz, ds, L, M2, np);

    Fq = K .* (d2Xqdq + d2Xqds);
    Fs = K .* (d2Xsdq + d2Xsds);
    Fz = K .* (d2Xzdq + d2Xzds);

else
    [d2Xqdq, d4Xqdq] = d2XdqPeriodic(Xq, dq, L, M1, p,  bending);
    [d2Xsdq, d4Xsdq] = d2XdqPeriodic(Xs, dq, L, M1, np, bending);
    [d2Xzdq, d4Xzdq] = d2XdqPeriodic(Xz, dq, L, M1, np, bending);

    [d2Xqds, d4Xqds]=  d2XdsPeriodic(Xq, ds, L, M2, np, bending);
    [d2Xsds, d4Xsds] = d2XdsPeriodic(Xs, ds, L, M2, p,  bending);
    [d2Xzds, d4Xzds] = d2XdsPeriodic(Xz, ds, L, M2, np, bending); 

    %mixed partial 
    d4Xqdqds = d2XdsPeriodic(d2Xqdq, ds, L, M2, np);
    d4Xsdqds = d2XdsPeriodic(d2Xsdq, ds, L, M2, np);
    d4Xzdqds = d2XdsPeriodic(d2Xzdq, ds, L, M2, np);

    Fq = K .* (d2Xqdq + d2Xqds) - Kb .* (d4Xqdq + d4Xqds + 2 * d4Xqdqds);
    Fs = K .* (d2Xsdq + d2Xsds) - Kb .* (d4Xsdq + d4Xsds + 2 * d4Xsdqds);
    Fz = K .* (d2Xzdq + d2Xzds) - Kb .* (d4Xzdq + d4Xzds + 2 * d4Xzdqds);
end

dS = dq * ds .* ones( size(Xq) );
end