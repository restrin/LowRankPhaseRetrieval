function [y,idx] = opA(A, x, adjoint, opts)
% OPA   Compute forward operator y = A*x or adjoint y = A'*x
%
% Input:
%   A         m-by-n operator
%   x         Input vector
%   adjoint   Compute adjoint operator?
%   opts      Options struct
%      type             sampling type for forward operator ('full', 'pos', 'topk')
%      total_samples    number of samples (for 'topk')
%      alpha            scaling factor for distribution (0 is uniform, bigger alpha more skewed)
%      explicit         if adjoint=true, return A'*x as matrix or operator?
%      symm             if adjoint=true and explicit=true, symmetrize matrix?
%
% Output:
%   y         if adjoint=false, Ax, else, A'*x
%   idx       indices of sampled entries (nonempty if adjoint=true)

    [m,n] = size(A);
    idx = [];

    if ~adjoint
        y = sum((A*x).^2,2);
    else
        if strcmpi(opts.type,'full')
            idx = 1:m;
            if opts.explicit
                y = A'*(x.*A);
            else
                y = @(v) A'*(x.*(A*v));
            end

        elseif strcmpi(opts.type,'pos')
            idx = x > 0;

%            scal = sum(x) / sum(x(idx));
            
            if opts.explicit
                y = A(idx,:)'*(x(idx).*A(idx,:));
%                y = scal*y;
            else
                Aidx = A(idx,:);
                y = @(v) (Aidx'*(x(idx).*(Aidx*v)));
            end

        elseif strcmpi(opts.type,'topk')
            % This is gonna be real slow
            idx_bank = find(x > 0);
            num_samples = min(opts.total_samples,length(idx_bank));
            idx = zeros(num_samples,1);
            
            for k = 1:num_samples
                if length(idx_bank) == 1
                    idx = idx_bank;
                    break
                end

                distr = x(idx_bank);
                distr = distr.^(opts.alpha);
                distr = distr / sum(distr);

                s = randsample(idx_bank,1,true,distr);
                idx(k) = s;
                idx_bank = setdiff(idx_bank,s);
            end

%            scal = sum(x) / sum(x(idx));
            
            if opts.explicit
                y = A(idx,:)'*(x(idx).*A(idx,:));
%                y = scal*y;
            else
                Aidx = A(idx,:);
                y = @(v) (Aidx'*(x(idx).*(Aidx*v)));
            end
        end

        % Even though y is symmetric, roundoff error makes MATLAB think its
        % nonsymmetric which makes eig slow
        % We can re-write y=A'*(x.*A) to
        %    sx = sqrt(x)
        %    sxA = sx.*A
        %    y = sxA'*sxA
        % MATLAB should think that the operator is symmetric then
        % Unclear which is more expensive
        if opts.explicit && isfield(opts,'symm') && opts.symm
            y = (y+y')/2;
        end
    end
end