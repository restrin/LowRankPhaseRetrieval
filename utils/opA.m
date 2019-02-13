function [y,idx] = opA(A, x, transpose, explicit, isfactor, opts)
% Options:
%   opts.type - 'full', 'pos', 'topk'
%   opts.symm - true or false

    [m,n] = size(A);
    idx = [];

    if ~transpose
        if isfactor
            y = sum((A*x).^2,2);
        else
            Ax  = A*x;
            y = sum(Ax.*A,2);
        end
    else
        if strcmpi(opts.type,'full')
            idx = 1:m;
            if explicit
                y = A'*(x.*A);
            else
                y = @(v) A'*(x.*(A*v));
            end

        elseif strcmpi(opts.type,'pos')
            idx = x > 0;

%            scal = sum(x) / sum(x(idx));
            
            if explicit
                y = A(idx,:)'*(x(idx).*A(idx,:));
%                y = scal*y;
            else
                Aidx = A(idx,:);
                y = @(v) (Aidx'*(x(idx).*(Aidx*v)));
            end

        elseif strcmpi(opts.type,'topk')
            idx_bank = find(x > 0);
            idx = [];
            
            for k = 1:min(opts.total_samples,length(idx_bank))
                if length(idx_bank) == 1
                    idx = [idx, idx_bank];
                    break
                end

                distr = x(idx_bank);
                distr = distr.^(opts.alpha);
                distr = distr / sum(distr);

                s = randsample(idx_bank,1,true,distr);
                idx = [idx, s];
                idx_bank = setdiff(idx_bank,s);
            end

%            scal = sum(x) / sum(x(idx));
            
            if explicit
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
        if explicit && isfield(opts,'symm') && opts.symm
            y = (y+y')/2;
        end
    end
end