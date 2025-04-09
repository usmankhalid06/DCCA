function xo = spatial(x,nV1,nV2)
    x = reshape(x,nV1,nV2);
    xr = padarray(x,[1 1],0,'both');
    TT = conv2(xr,1/9*[1 1 1;1 1 1;1 1 1],'valid');
    xo = reshape(TT,1,nV1*nV2);
end



