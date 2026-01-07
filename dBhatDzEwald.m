function Bprime=dBhatDzEwald(c,kx,ky,z,xi,x0,y0)

Bprime=exp(-sqrt(-1).*(kx.*x0+ky.*y0)).*(-1./(16.*c.*pi).*z.*(exp(-c.*z).*erfc(-z./xi+c.*xi./2)...
    + exp(c.*z).*erfc(z./xi+c.*xi./2)) - (1/(8.*pi^(3/2))).*exp(1).^((-1/4).*c.^2.*xi.^2+(-1).*xi.^(-2).* ...
    z.^2).*xi.*z);

end