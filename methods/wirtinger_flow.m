function [u,track] = wirtinger_flow(prob,opts, ustart)

    m = prob.m;
    b = prob.b;
    n = prob.n;
    A = prob.A;

    % Set parameters
    if isfield(opts, 'explicit')
        explicit = opts.explicit;
    else
        explicit = true;
    end
    if ~explicit
        eigopts.issym = true; 
    end
    opts.sampling_scheme.type = 'full';

    if isempty(ustart)

        [W,~] = opA(prob.A,b,true,explicit,false,opts.sampling_scheme);

        if explicit
            [V,D] = eig(W,'vector');
            [~, jmax] = max(D);
            u = V(:,jmax);
        else
            [u,~] = eigs(W,n,1,'la',eigopts);
        end
    else
        u = ustart;
    end

    stepsize = opts.stepsize;
    for iter = 1:opts.maxiter

        Au = A*u;
        y = Au.^2-b;
        track.obj(iter) = sum_square(y)/2/m;
        g = A'*(y.*Au)/m;
        stepsize = stepsize * opts.stepsize_decay;
        u = u - stepsize*g;

        if isfield(opts, 'callback')
            opts.callback([], iter, track.obj(iter), u, prob);
        end
    end
end