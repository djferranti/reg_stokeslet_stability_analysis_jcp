function B=Bhat72(c,kx,ky,z,epsilon,x0,y0) 
    Z2 = z.^2 + epsilon.^2; Z = sqrt(Z2);
    B = exp(-sqrt(-1).*(kx .* x0  + ky .* y0)) .* (1) ./ (8 * pi) .* ...
        exp(-c .* Z) ./ (c.^3).* (1 + c .* Z);
end