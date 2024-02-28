function C = qcolor(q)
    C = reshape(compact(q),[size(q),4]);
    C = C(:,:,2:4);
end