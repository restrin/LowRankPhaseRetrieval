function [u,track] = wirtinger_flow(prob, opts, ustart)
% WIRTINGER_FLOW    Run Wirtinger flow
%
% Input:
%   prob    Problem instance
%   ustart  Initial vector
%   opts    Options struct
%      maxiter          Max number of iterations (required)
%      stepsize         Initial stepsize (required)
%      stepsize_decay   stepsize *= stepsize_decay every iteration (required)
%
%      explicit         Compute initial point using dense eigendecomposition? (optional, default true)
%      callback         Callback function taking (y, iter, obj, vmax, prob, opts)
%                       returns stop = true if solver should exit (optional)
%      callbackopts     Options struct for callback
%
% Output:
%   u       Solution
%   track   Struct containing convergence history
%      obj      Objective value

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
    opA_opts.explicit = explicit;
    opA_opts.type = 'full';
    
    if isempty(ustart)
        [W,~] = opA(prob.A,b,true,opA_opts);

        if explicit
            [V,D] = eig(W,'vector');
            [~, jmax] = max(D);
            u = V(:,jmax);
        else
            [u,~] = eigs(W,n,1,'la',eigopts);
        end
        Au = (A*u).^2;
        u = u*sqrt(((Au'*b)/(norm(Au)^2)));
    else
        u = ustart;
    end

    if isfield(opts, 'callbackopts')
        callbackopts = opts.callbackopts;
    else
        callbackopts = struct();
    end
    
    stepsize = opts.stepsize;
    for iter = 1:opts.maxiter

        Au = A*u;
        y = Au.^2-b;
        track.obj(iter) = sum_square(y)/2/m;
        g = A'*(y.*Au)/m;
        stepsize = stepsize * opts.stepsize_decay;
        u = u - stepsize*g;

        stop = 0;
        if isfield(opts, 'callback')
            stop = opts.callback([], iter, track.obj(iter), u, prob, callbackopts);
        end
        if stop, return, end
    end
end