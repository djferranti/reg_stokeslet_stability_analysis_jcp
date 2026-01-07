function [U2dp, V2dp, W2dp] ...
    = dpregularizedstokeslethighordergauss(C, Kx, Ky, X, Y, Z, ...
    L, xyz0, forceK, epsilon, mu, NInterp)

%initialize the three dimensional arrays for the domain grid
U2dp = zeros(NInterp, NInterp);
V2dp = U2dp; 
W2dp = U2dp; 

uhatPad = zeros(NInterp, NInterp);
vhatPad = zeros(NInterp, NInterp);
whatPad = zeros(NInterp, NInterp);

x0 = xyz0(1); y0 = xyz0(2); z0 = xyz0(3);
numberZlevels = size(Z, 3); 

fx0 = forceK(1); fy0 = forceK(2); fz0 = forceK(3);

N = size(C, 1);

for j = 1 : numberZlevels 
    
    Xj = X(:,:,j); Yj = Y(:,:,j); Zj = Z(:,:,j);
    z = Zj - z0; z = z(1);

    % compute long range periodic forces
    Gfour = GHatGaussMoments(C, Kx, Ky, z, epsilon, x0, y0);
    Bfour = BHatGaussMoments(C, Kx, Ky, z, epsilon, x0, y0);
    Bprimefour = dBhatDzGaussMoments(C, Kx, Ky, z, epsilon, x0, y0);

    %then compute velocity for nonzero wavenumbers
    uhat = 1 / mu .* ( (-Kx .^ 2 .* Bfour - Gfour) .* fx0 + ...
        (-Kx .* Ky .* Bfour) .* fy0 + ...
        (sqrt(-1) .* Kx .* Bprimefour) .* fz0 );

    vhat = 1 / mu .* ( (-Kx .* Ky .* Bfour) .* fx0 + ...
        (-Ky .^ 2 .* Bfour - Gfour).* fy0 + ...
        (sqrt(-1) .* Ky .* Bprimefour) .* fz0 );

    what = 1 / mu .* ( (sqrt(-1) .* Kx .* Bprimefour) .* fx0 + ...
        (sqrt(-1) .* Ky.* Bprimefour) .* fy0 + ...
        (C .^ 2 .* Bfour) .* fz0);

    %then velocity for the wavenumber (0,0)
    G0 = (1/12).*epsilon.^(-5).*pi.^(-3/2).*(6.*exp(1).^((-1).*epsilon.^( ...
        -2).*z.^2).*epsilon.^4.*z.^2+(-2).*exp(1).^((-1).*epsilon.^(-2).* ...
        z.^2).*epsilon.^2.*z.^4+3.*epsilon.^5.*pi.^(1/2).*z.*erf( ...
        epsilon.^(-1).*z));

    uhat(1,1) = 1 / mu .* G0 .* (-fx0);

    vhat(1,1) = 1 / mu .* G0 .* (-fy0);

    what(1,1) = 0; 

    uhatPad(1 : N/2 + 1, 1 : N/2 + 1) =  uhat(1 : N/2 + 1, 1 : N/2 + 1);
    uhatPad(1 : N/2 + 1 , NInterp - N/2 + 2 : end) = ...
        uhat(1 : N/2 + 1, N/2 + 2 : end );
    uhatPad(NInterp - N/2 + 2: end, 1 : N/2 + 1) = ...
        uhat(N/2 + 2 : end, 1 : N/2 + 1);
    uhatPad(NInterp - N/2 + 2: end, NInterp - N/2 + 2:end) = ...
        uhat(N/2 + 2 : end, N/2 + 2 : end);

    vhatPad(1 : N/2 + 1, 1 : N/2 + 1) =  vhat(1 : N/2 + 1, 1 : N/2 + 1);
    vhatPad(1 : N/2 + 1 , NInterp - N/2 + 2 : end) = ...
        vhat(1 : N/2 + 1, N/2 + 2 : end );
    vhatPad(NInterp - N/2 + 2: end, 1 : N/2 + 1) = ...
        vhat(N/2 + 2 : end, 1 : N/2 + 1);
    vhatPad(NInterp - N/2 + 2: end, NInterp - N/2 + 2:end) = ...
        vhat(N/2 + 2 : end, N/2 + 2 : end);

    whatPad(1 : N/2 + 1, 1 : N/2 + 1) =  what(1 : N/2 + 1, 1 : N/2 + 1);
    whatPad(1 : N/2 + 1 , NInterp - N/2 + 2 : end) = ...
        what(1 : N/2 + 1, N/2 + 2 : end );
    whatPad(NInterp - N/2 + 2: end, 1 : N/2 + 1) = ...
        what(N/2 + 2 : end, 1 : N/2 + 1);
    whatPad(NInterp - N/2 + 2: end, NInterp - N/2 + 2 : end) = ...
        what(N/2 + 2 : end, N/2 + 2 : end);

    uLong=2*pi*(NInterp/L).^2.*real(ifft2(uhatPad));
    vLong=2*pi*(NInterp/L).^2.*real(ifft2(vhatPad));
    wLong=2*pi*(NInterp/L).^2.*real(ifft2(whatPad));

    U2dp(:,:,j) = uLong;
    V2dp(:,:,j) = vLong;
    W2dp(:,:,j) = wLong;

end
end