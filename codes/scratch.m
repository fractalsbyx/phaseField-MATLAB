[x,y,z] = meshgrid([-3:0.25:3]); 
V = x.*exp(-x.^2 -y.^2 -z.^2);
s = isosurface(x,y,z,V,1e-4);
p = patch(s);
isonormals(x,y,z,V,p)
view(3);
set(p,'FaceColor',[0.5 1 0.5]);  
set(p,'EdgeColor','none');
%camlight;
lighting gouraud;