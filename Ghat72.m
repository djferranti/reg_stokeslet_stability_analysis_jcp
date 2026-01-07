function G=Ghat72(c,kx,ky,z,epsilon,x0,y0) 
    Z2 = z.^2 + epsilon.^2; Z = sqrt(Z2);
    G = exp(-sqrt(-1).*(kx .* x0  + ky .* y0)) .* (-1) ./ (8 * pi) .* ...
        (2 .* exp(-c .* Z) ./ c + epsilon.^2 .* exp(-c .* Z ) ./ Z ) ;
end