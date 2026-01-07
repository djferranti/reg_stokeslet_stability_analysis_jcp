function dBdz=dBhat72dz(c,kx,ky,z,epsilon,x0,y0) 
    Z2 = z.^2 + epsilon.^2; Z = sqrt(Z2);
    dBdz = exp(-sqrt(-1).*(kx .* x0  + ky .* y0)) .* (-1) ./ (8 * pi) .* ...
        exp(-c .* Z) .* z ./ c;
end