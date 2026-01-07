function [uLocal, vLocal, wLocal] = evaluatelocalvelocitynewgaussian(X, Y, ...
    Z, L, epsilon, xi, rCutoff, mu, x0, y0, z0, fx0, fy0, fz0)

x0 = mod(x0, L); y0 = mod(y0, L);
%for force in domain [0,L] x [0, L]

[U00, V00, W00] = localvelocitynewgaussian(X, Y, Z, epsilon, xi, ...
    rCutoff, x0, y0, z0, fx0, fy0, fz0);

%start at 12 o'clock and go in clockwise order

%for force in domain [0, L] x [L, 2L]
[U0L, V0L, W0L] = localvelocitynewgaussian(X, Y, Z, epsilon, xi, ...
    rCutoff, x0, y0 + L, z0, fx0, fy0, fz0);

%for force in domain [L, 2L] x [L, 2L]
[ULL, VLL, WLL] = localvelocitynewgaussian(X, Y, Z, epsilon, xi, ...
    rCutoff, x0 + L, y0 + L, z0, fx0, fy0, fz0);

%for force in domain [L, 2L] x [0, L]
[UL0, VL0, WL0] = localvelocitynewgaussian(X, Y, Z, epsilon, xi, ...
    rCutoff, x0 + L, y0, z0, fx0, fy0, fz0);

%for force in domain [L, 2L] x [-L, 0]
[ULmL, VLmL, WLmL] = localvelocitynewgaussian(X, Y, Z, epsilon, xi, ...
    rCutoff, x0 + L, y0 - L, z0, fx0, fy0, fz0);

%for force in domain [0, L] x [-L, 0]
[U0mL, V0mL, W0mL] = localvelocitynewgaussian(X, Y, Z, epsilon, xi, ...
    rCutoff, x0, y0 - L, z0, fx0, fy0, fz0);

%for force in domain [-L, 0] x [-L, 0]
[UmLmL, VmLmL, WmLmL] = localvelocitynewgaussian(X, Y, Z, epsilon, xi, ...
    rCutoff, x0 - L, y0 - L, z0, fx0, fy0, fz0);

%for force in domain [-L, 0] x [0, L]
[UmL0, VmL0, WmL0] = localvelocitynewgaussian(X, Y, Z, epsilon, xi, ...
    rCutoff, x0 - L, y0, z0, fx0, fy0, fz0);

%for force in domain [-L, 0] x [L, 2L]
[UmLL, VmLL, WmLL] = localvelocitynewgaussian(X, Y, Z, epsilon, xi, ...
    rCutoff, x0 - L, y0 + L, z0, fx0, fy0, fz0);

uLocal = 1 ./ (8 .* pi .* mu) .* ( U00 + U0L + ULL + UL0 + ULmL + ...
    U0mL + UmLmL + UmL0 + UmLL );
vLocal = 1 ./ (8 .* pi .* mu) .* ( V00 + V0L + VLL + VL0 + VLmL + ...
    V0mL + VmLmL + VmL0 + VmLL );
wLocal = 1 ./ (8 .* pi .* mu) .* ( W00 + W0L + WLL + WL0 + WLmL + ...
    W0mL + WmLmL + WmL0 + WmLL );


    function [U, V, W] = ...
            localvelocitynewgaussian(X, Y, Z, epsilon, xi, rCutoff, ...
            x0, y0, z0, fx0, fy0, fz0)
        U = zeros(size(X)); V = U; W = U;
        Rx = X - x0;
        Ry = Y - y0;
        Rz = Z - z0;

        R2 = Rx.^2 + Ry.^2 + Rz.^2;
        R = sqrt(R2);

        %find points within rCutoff - only need to compute
        % velocity field for these
        withinCutoff = R < rCutoff;
        rc = R(withinCutoff);

        %only evaluate if there are any points within cutoff radius
        if ~isempty(rc)

            rcm1 = rc .^ (-1); rcm3 = rc .^ (-3);
            rxc = Rx(withinCutoff); ryc = Ry(withinCutoff);
            rzc = Rz(withinCutoff);
            expXi = exp( - rc.^2 ./ xi.^2);
            expEps = exp( - rc.^2 ./ epsilon.^2);
            erfXi = erf( rc ./ xi);
            erfEps = erf( rc ./ epsilon);

            Cloc = rc / sqrt(pi) .* ( (-6 ./ xi + 4 .* rc.^2 ./ xi.^3  )...
                .* expXi + (26 ./ 3 ./ epsilon ...
                - 12 .* rc.^2 ./ epsilon.^3  ...
                + 8 ./ 3 .*rc.^4 ./ epsilon.^5)  .* expEps) ...
                - (erfXi - erfEps);
            Dloc = rc / sqrt(pi) .* ( (2 ./ xi - 4 .* rc.^2 ./ xi.^3  ) ...
                .* expXi  + (- 2 ./ epsilon  ...
                + 28 ./ 3 .* rc.^2 ./ epsilon.^3 ...
                - 8 ./ 3 .* rc.^4 ./ epsilon.^5 ).* expEps) ...
                - (erfXi - erfEps);


            FdotX = fx0 .* rxc + fy0 .* ryc + fz0 .* rzc;

            uCutoff = fx0 .* rcm1 .* Cloc + FdotX .* rxc .* rcm3 .* Dloc;
            vCutoff = fy0 .* rcm1 .* Cloc + FdotX .* ryc .* rcm3 .* Dloc;
            wCutoff = fz0 .* rcm1 .* Cloc + FdotX .* rzc .* rcm3 .* Dloc;


            %use taylor series for evaluation at r <epsilon/10 
            rClose = rc < epsilon/10;
            FdotXClose = FdotX(rClose);
            rcClose = rc(rClose);
            rxcClose = rxc(rClose);
            rycClose = ryc(rClose);
            rzcClose = rzc(rClose);

            C0 =  - 8 / 3 / sqrt(pi) / epsilon .* ...
                ( 3 * epsilon / xi - 4 ) .* ones(size(rcClose));
            D0 = - 16 / 3 /sqrt(pi) / epsilon^3 .* ...
                ( ( epsilon / xi )^3 - 2 ) .* ones(size(rcClose));
            C1 = 32 / 3 / sqrt(pi) ./ epsilon^3 .* ((epsilon / xi)^3 - 2 ) .* rcClose.^2;
            D1 = 8 / 5 / sqrt(pi) ./ epsilon^5 .* (3 * (epsilon / xi)^5 - 8 ) .* rcClose.^2;
            C2 = - 12 / 5 / sqrt(pi) ./ epsilon^5 .*(3 * (epsilon/xi)^5 - 8) .* rcClose.^4;
            D2 = - 16/7 / sqrt(pi) ./ epsilon^7 .*  ((epsilon/xi)^7 - 10/3) .* rcClose.^4;

            uCutoff(rClose) = fx0 .* (C0 + C1 + C2) + FdotXClose .* ...
                rxcClose .* (D0 + D1 + D2 ) ;
            vCutoff(rClose) = fy0 .* (C0 + C1 + C2) + FdotXClose .* ...
                rycClose .* (D0 + D1 + D2) ;
            wCutoff(rClose) = fz0 .* (C0 + C1 + C2) + FdotXClose .* ...
                rzcClose .* (D0 + D1 + D2 );



            U(withinCutoff) = uCutoff;
            V(withinCutoff) = vCutoff;
            W(withinCutoff) = wCutoff;

        end
    end
end