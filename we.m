function WE = we(x, nlv)
    nnlv = size(x,1);
    WE = zeros(nlv,nlv+nnlv);
    WE(:,1:nlv) = eye(nlv);
    for j = 1:nlv
        xindex1 = find(x(:,1) == j);
        c1 = x(xindex1,2);
        xindex2 = find(x(:,2) == j);
        c2 = x(xindex2,1);
        c = [c1;c2];
        xindex = [xindex1;xindex2];
        [c,order] = sort(c);
        xr = length(c);
        rowindex = xindex(order);
        rowindex = nlv + rowindex;
        if xr > 0
           WE(j,rowindex) = c;
        end
    end