function [y, track] = reduced_grad(prob,opts)
% Prob:
%    A, b, m, n, orig
% Opts:
%    Mandatory:
%       maxiter
%       sampling_scheme
%       stepsize
%       stepsize_decay
%       checkperiod
%       saveperiod
%    Optional:
%       y0       - default: 0
%       explicit - default: true
%       symm     - default: false
%       callback - stop = opts.callback(y, iter, obj, vmax, prob);

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
else
    callback = @(y, iter, vmax, prob) false;
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
    stop = callback(y, iter, objval, vmax, prob);
    
    if stop
        break;
    end
    
end

end