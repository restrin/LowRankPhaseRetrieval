function [x,val,W] = primal_from_dual(prob, y, explicit)

    samplestruct = struct;
    samplestruct.type = 'full';
    samplestruct.symm = true;

    [W,~] = opA(prob.A,y,true,explicit,false,samplestruct);

    if explicit
        [V,D] = eig(W,'vector');
        [val, jmax] = max(D);
        x = V(:,jmax);
    else
        eigopts.issym = true;
        [x,val] = eigs(W,prob.n,1,'la',eigopts);
    end  
end