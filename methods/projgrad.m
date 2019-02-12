function [W,track] = projgrad(prob,opts)
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
%       callback - stop = opts.callback(y, iter, vmax, prob);

m = prob.m;
n = prob.n;
b = prob.b;
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
else
    callback = @(y, iter, vmax, prob) false;
end

for iter = 1:opts.maxiter
    
    [W,~] = opA(prob.A,y+ybar,true,explicit,false,opts.sampling_scheme);
    
    if explicit
        [V,D] = eig(W,'vector');
        [objval, jmax] = max(D);
        vmax = V(:,jmax);
    else
        [vmax,objval] = eigs(W,n,1,'la',eigopts);
    end
    
%    track.eiggap(iter,:) = D;
%    eiggap = D(1)-D(2);
    track.obj(iter) = objval;
    
    ygrad = (prob.A*vmax).^2;
    
    y = y - opts.stepsize*ygrad;
    y = y - ((y'*b)/(b'*b)*b);
    
    opts.stepsize = opts.stepsize*opts.stepsize_decay;
    stop = callback(y, iter, objval, vmax, prob);
    
    if stop
        break;
    end
    
end
end