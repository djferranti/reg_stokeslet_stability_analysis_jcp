function [Fq, Fs, Fz, dS] = computetensilebendingforcesaperiodic(K, Xq, ... 
    Xs, Xz, M1, M2, dq, ds, Kb)

    tauqXq  = (Xq(:,2 : M2) - Xq(:,1 : M2 - 1));
    tauqXs  = (Xs(:,2 : M2) - Xs(:,1 : M2 - 1));
    tauqXz  = (Xz(:,2 : M2) - Xz(:,1 : M2 - 1)); 
    normTauq = sqrt(tauqXq.^2 + tauqXs.^2 + tauqXz.^2); 

    d2Xqdq = d2Xdq(tauqXq, dq, M1, normTauq);
    d2Xsdq = d2Xdq(tauqXs, dq, M1, normTauq);
    d2Xzdq = d2Xdq(tauqXz, dq, M1, normTauq);

    tausXq  = (Xq(2 : M1, :) - Xq(1 : M1 - 1, :));
    tausXs  = (Xs(2 : M1, :) - Xs(1 : M1 - 1, :));
    tausXz  = (Xz(2 : M1, :) - Xz(1 : M1 - 1, :)); 

    normTaus = sqrt(tausXq.^2 + tausXs.^2 + tausXz.^2); 

    d2Xqds = d2Xds(tausXq, ds, M2, normTaus);
    d2Xsds = d2Xds(tausXs, ds, M2, normTaus);
    d2Xzds = d2Xds(tausXz, ds, M2, normTaus);

    d4Xqdq = d4Xdq(Xq, dq, M1);
    d4Xsdq = d4Xdq(Xs, dq, M1);
    d4Xzdq = d4Xdq(Xz, dq, M1);

    d4Xqds = d4Xdq(Xq, ds, M2);
    d4Xsds = d4Xdq(Xs, ds, M2);
    d4Xzds = d4Xdq(Xz, ds, M2);

    d4Xqdqds = d4Xdqds(Xq, dq, ds, M1, M2);
    d4Xsdqds = d4Xdqds(Xs, dq, ds, M1, M2);
    d4Xzdqds = d4Xdqds(Xz, dq, ds, M1, M2);

    Fq = K .* (d2Xqdq + d2Xqds) - Kb .* (d4Xqdq + d4Xqds + 2 * d4Xqdqds);
    Fs = K .* (d2Xsdq + d2Xsds) - Kb .* (d4Xsdq + d4Xsds + 2 * d4Xsdqds);
    Fz = K .* (d2Xzdq + d2Xzds) - Kb .* (d4Xzdq + d4Xzds + 2 * d4Xzdqds);
  
    %if there is no bending forces, add external pulling forces to flatten
    %the surface
    if Kb == 0
        Kpull = 75; 
    else 
        Kpull = 0;
    end

    Fq(:,1) = Fq(:, 1) - Kpull.*ones(M2, 1); 
    Fq(:,M1) = Fq(:, M1) + Kpull.*ones(M2,1);
    Fs(1, :) = Fs(1, :) - Kpull.*ones(1, M1);
    Fs(M2, :) = Fs(M2, :) + Kpull.*ones(1, M1);

    dS = dq * ds .* ones( size(Xq) );
end