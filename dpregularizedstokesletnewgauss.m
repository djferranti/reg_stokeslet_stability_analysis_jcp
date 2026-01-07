function [U2dp, V2dp, W2dp] = dpregularizedstokesletnewgauss(C, Kx, Ky, ...
    X, Y, Z, L, xyz0, forceK, epsilon, mu, NInterp)

%initialize the three dimensional arrays for the domain grid
U2dp = zeros(size(X));
V2dp = U2dp; 
W2dp = U2dp;

x0 = xyz0(1); y0 = xyz0(2); z0 = xyz0(3);
numberZlevels = size(Z, 3); 

N = size(C,1);
uhatPad = zeros(NInterp, NInterp);
vhatPad = zeros(NInterp, NInterp); 
whatPad = zeros(NInterp, NInterp);

for j = 1 : numberZlevels 
    
    Xj = X(:,:,j); Yj = Y(:,:,j); Zj = Z(:,:,j);
    % XjShift = Xj + xshift; YjShift = Yj + yshift;
    z = Zj - z0; z = z(1);

    fx0 = forceK(1); fy0 = forceK(2); fz0 = forceK(3);

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
     uhat(1,1)=1/mu.*( exp(-z.^2./reg.^2) .* ...
        (-1 .* z.^4 + 3 .* z.^2 .*reg.^2) ./ 6 ./ pi^(3/2) ./ reg.^3 + ...
        z ./ 4 ./ pi .*erf(z./reg)).*(-fx0);
    vhat(1,1)= 1/mu.*( exp(-z.^2./reg.^2) .* ...
        (-1 .* z.^4 + 3 .* z.^2 .*reg.^2) ./ 6 ./ pi^(3/2) ./ reg.^3 + ...
        z ./ 4 ./ pi .*erf(z./reg)).*(-fy0);
    what(1,1)=0;  


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

    %physical space sum
    % [uLocal, vLocal, wLocal] = evaluatelocalvelocitynewalgebraic(Xj, Yj, ...
    %     Zj, L, epsilon, xi, rCutoff, mu, x0, y0, z0, fx0, fy0, fz0);

    % [uLocal, vLocal, wLocal] = evaluatelocalvelocityclassic(Xj, Yj, ...
    %         Zj, L, epsilon, xi, rCutoff, mu, x0, y0, z0, fx0, fy0, fz0);

    % U2dp(:,:,j) = uLong + uLocal;
    % V2dp(:,:,j) = vLong + vLocal;
    % W2dp(:,:,j) = wLong + wLocal;
end
end