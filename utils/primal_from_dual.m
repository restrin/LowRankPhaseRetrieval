function [x,val,W] = primal_from_dual(prob, y, explicit)
%PRIMAL_FROM_DUAL   Compute primal solution given dual
%
% Computes most positive eigenvector of W = A'*(y.*A)
%
% Input:
%   prob        Problem data (struct)
%   y           Dual variable
%   explicit    Use eig (explicit=true) or eigs (explicit=false)
%
% Output:
%   x           Top eigenvector of A'*(y.*A)
%   val         Top eigenvalue of A'*(y.*A)
%   W           Dense matrix W=A'*(y.*A) (explicit=true), or operator
%               (explicit=false)

    opAopts = struct;
    opAopts.type = 'full';
    opAopts.symm = true;
    opAopts.explicit = explicit;

    [W,~] = opA(prob.A,y,true,opAopts);

    if explicit
        [V,D] = eig(W,'vector');
        [val, jmax] = max(D);
        x = V(:,jmax);
    else
        eigopts.issym = true;
        [x,val] = eigs(W,prob.n,1,'la',eigopts);
    end  
end