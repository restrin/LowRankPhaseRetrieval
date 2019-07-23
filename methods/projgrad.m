function [y,track] = projgrad(prob,opts)
% PROJGRAD    Run projected gradient
%
% Input:
%   prob    Problem instance
%   opts    Options struct
%      maxiter          Max number of iteraitons (required)
%      stepsize         Initial stepsize (required)
%      stepsize_decay   stepsize *= stepsize_decay every iteration (required)
%      opAopts          Struct for products with measurement operator (required)
%
%      y0               Initial y0 (optional, default 0)
%      explicit         Compute initial point using dense eigendecomposition? (optional, default true)
%      callback         Callback function taking (y, iter, obj, vmax, prob, opts)
%                       returns stop = true if solver should exit (optional)
%      callbackopts     Options struct for callback
%
% Output:
%   y       Solution
%   track   Struct containing convergence history
%      obj      Objective value

m = prob.m;
n = prob.n;
b = prob.b;
% Rescaling (should probably avoid for now)
%b = b / norm(prob.orig(:));
%prob.A = prob.A / m;
%b = b / m;
y = zeros(m,1);
if isfield(opts, 'y0')
    ybar = opts.y0;
else
    ybar = b/(b'*b);
end

if isfield(prob,'orig')
    [n1,n2] = size(prob.orig);
    prob.n1 = n1;
    prob.n2 = n2;
end

% Set parameters
if isfield(opts, 'explicit')
    explicit = opts.explicit;
else
    explicit = true;
end
if ~explicit
    eigopts.issym = true; 
end
if isfield(opts, 'callback')
    callback = opts.callback;
    if isfield(opts, 'callbackopts')
        callbackopts = opts.callbackopts;
    else
        callbackopts = struct();
    end
else
    callback = @(y, iter, vmax, prob, callbackopts) false;
end

track = [];

yfull = ybar;
for iter = 1:opts.maxiter
    
    [W,~] = opA(prob.A,yfull,true,opts.opAopts);
    
    if explicit
        [V,D] = eig(W,'vector');
        [objval, jmax] = max(D);
        vmax = V(:,jmax);
    else
        [vmax,objval] = eigs(W,n,1,'la',eigopts);
    end

    track.obj(iter) = objval;
    
    ygrad = (prob.A*vmax).^2;
    
    y = y - opts.stepsize*ygrad;
    y = y - ((y'*b)/(b'*b)*b);
    
    opts.stepsize = opts.stepsize*opts.stepsize_decay;
    
    yfull = y+ybar;
    stop = callback(yfull, iter, objval, vmax, prob, callbackopts);
    
    if stop
        break;
    end
    
end

% Return final dual iterate
y = yfull;

end