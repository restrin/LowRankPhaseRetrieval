function [y, track] = reducedgrad(prob,opts)
%REDUCEDGRAD    Run projected gradient
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

% Create sparse null-space matrix
% Assume b(i) != 0 for all i
% Z = [ b(2)/b(1) ...  b(m)/b(1)]
%     [-1         ...           ]
%     [           ...         -1]
[~, i] = max(abs(b)); i = i(1);
Z = -speye(m-1,m-1);
if i == m
   Z = [Z; b(1:end-1)'/b(end)];
else
   Z = [Z(1:i-1,:); b(1:i-1)'/b(i) b(i+1:end)'/b(i); Z(i:end,:)];
end
nZ = norm(Z,'fro');
Z = Z/nZ;

% Make ybar feasible
res = b'*ybar - 1;
ybar(i) = ybar(i) - res/b(i);

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

z = zeros(m-1,1);
y = ybar;
for iter = 1:opts.maxiter
    
    [W,~] = opA(prob.A,y,true,explicit,false,opts.sampling_scheme);
    
    if explicit
        [V,D] = eig(W,'vector');
        [objval, jmax] = max(D);
        vmax = V(:,jmax);
    else
        [vmax,objval] = eigs(W,n,1,'la',eigopts);
    end
    
    track.obj(iter) = objval;
    
    g = (prob.A*vmax).^2;
    
    z = z - opts.stepsize*Z'*g;
    
    opts.stepsize = opts.stepsize*opts.stepsize_decay;
    y = Z*z+ybar;
    stop = callback(y, iter, objval, vmax, prob, callbackopts);
    
    if stop
        break;
    end
    
end

end