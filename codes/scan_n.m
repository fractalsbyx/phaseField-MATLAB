Params
slice = zeros(Nx,Ny);
fig = image(slice, 'CDataMapping', 'scaled');
drawnow
ID = id_from_ops(n);
for i = 1:Nz
    slice = ID(:,:,i,1);
    set(fig, 'Cdata', slice);
    drawnow
end