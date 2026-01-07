function [U2dp, V2dp, W2dp, ULocal, VLocal, WLocal] ...
    = dpregularizedstokesletewaldgauss1m(C, Kx, Ky, X, Y, Z, ...
    Xq, Xs, Xz, L, xyz0, forceK, epsilon, xi, rCutoff, mu, NInterp)

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

%evaluate local velocity at surface points only (no need to use grid) 
Xq = mod(Xq, L); Xs = mod(Xs, L);
[ULocal, VLocal, WLocal] = evaluatelocalvelocitygaussian(Xq, Xs, ...
    Xz, L, epsilon, xi, rCutoff, mu, x0, y0, z0, fx0, fy0, fz0);

for j = 1 : numberZlevels 
    
    Zj = Z(:,:,j);
    z = Zj - z0; z = z(1);

    % compute long range periodic forces
    Gfour = GhatXi(C, Kx, Ky, z, xi, x0, y0);
    Bfour = BhatXi(C, Kx, Ky, z, xi, x0, y0);
    Bprimefour = DBhatXiDz(C, Kx, Ky, z, xi, x0, y0);

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
    G0 = 1 / (2 .* pi) .* 1 / 2 .* ( z.^2 ./ (sqrt(pi) .* xi) .* ...
        exp( - z .^ 2 ./ xi .^ 2) + z .* erf( z./ xi) );

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