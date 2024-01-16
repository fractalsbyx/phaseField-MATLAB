Params
tex = zeros(Nx,Ny,Nz,'quaternion');
for i = 1:np
    D2=2;
    while D2>1
        a=2*rand()-1;
        b=2*rand()-1;
        c=2*rand()-1;
        d=2*rand()-1;
        D2= a*a+b*b+c*c+d*d;
    end
    tex(i) = quaternion(a,b,c,d);
end