function L = computePointLaplacian(pts, type, rings, options)
%  º∆À„Laplacianæÿ’ÛL, L = -(Laplacian operator).

options.null = 0;
normalize = getoptions(options, 'normalize', 0);
symmetrize = getoptions(options, 'symmetrize', 1);

W = computePointWeight(pts, type, rings, options);
n = size(W,1);
    
if symmetrize==1 && normalize==0
    L = diag(sum(W,2)) - W;
elseif symmetrize==1 && normalize==1
    L = speye(n) - diag(sum(W,2).^(-1/2)) * W * diag(sum(W,2).^(-1/2));
elseif symmetrize==0 && normalize==1
    L = speye(n) - diag(sum(W,2).^(-1)) * W;
else
    error('Does not work with symmetrize=0 and normalize=0');    
end



