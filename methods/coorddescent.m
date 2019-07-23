function [y,q,track] = coorddescent(prob,opts)
% COORDDESCENT    Run coordinate descent
%
% Input:
%   prob    Problem instance
%   opts    Options struct
%      maxiter          Max number of iteraitons (required)
%      stepsize         Initial stepsize (required)
%      stepsize_decay   stepsize *= stepsize_decay every iteration (required)
%      opAopts          Struct for products with measurement operator (required)
%      rank             Rank of approximation (required)
%      recycle          How often factorization is refreshed (required)
%
%      y0               Initial y0 (optional, default 0)
%      explicit         Compute initial point using dense eigendecomposition? (optional, default true)
%      callback         Callback function taking (y, iter, obj, vmax, prob, opts)
%                       returns stop = true if solver should exit (optional)
%      callbackopts     Options struct for callback
%      sample_strat     Strategy for coordinate sampling, 'unif' or 'greedy' (optional, default 'unif')
%      alpha            Weight for greedy sampling (optional, default 0.1)
%      Q0, T0           Initial low-rank factorization (optional)
%      blocklen         Size of block for update (optional, default 1)
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
[rows, ~, Zval] = find(Z);
rows = reshape(rows, 2, size(rows,1)/2)';
Zval = reshape(Zval, 2, size(Zval,1)/2)';

% Make ybar feasible
res = b'*ybar - 1;
ybar(i) = ybar(i) - res/b(i);

%% Set parameters
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
    callback = @(y, iter, vmax, prob, opts) false;
end

if ~isfield(opts, 'sample_strat')
    strat = 'unif';
else
    strat = opts.sample_strat;
end
if isfield(opts, 'alpha')
    alpha = opts.alpha;
else
    alpha = 0.1;
end
if ~isfield(opts, 'blocklen')
    blocklen = 1;
else
    blocklen = opts.blocklen;
end

%% Initialize stuff

z = zeros(m-1,1);
y = ybar;

samplestruct = struct;
samplestruct.type = 'full';
samplestruct.symm = true;
samplestruct.explicit = opts.explicit;

if isfield(opts, 'Q0') && isfield(opts, 'T0')
    Q = opts.Q0;
    T = opts.T0;
    
    [v,~] = eig(T);
    v = real(v);
    Q = Q*v;
else
    [W,~] = opA(prob.A,y,true,samplestruct);
    
    if explicit
        [Q,T] = eig(W,'vector');
        [d,ix] = sort(T, 'descend');
        T = diag(d(1:opts.rank));
        Q = Q(:,ix(1:opts.rank));
    else
        [Q,T] = eigs(W, n, opts.rank, 'la', eigopts);
    end
end

if strcmpi(strat,'greedy')
    nonunifsample('init', (ybar).^(alpha));
end

track.numrestarts = 0;
jmax = 1;
iter = 0;

%% Start solving
for e=1:opts.maxiter

    if strcmpi(strat,'unif')
        order = randperm(m-1);
    end
    
    for i=1:floor((m-1)/blocklen)

        iter = iter+1;

        if strcmpi(strat,'unif')
            k = order((i-1)*blocklen+1:i*blocklen);
        elseif strcmpi(strat, 'greedy')
            k = nonunifsample('sample', blocklen);
            k = max(k-1,1);
        end
        
        if(mod(iter, opts.recycle) == 0)
            [W,~] = opA(prob.A,Z*z+ybar,true,samplestruct);

            if explicit
                [Q,T] = eig(W,'vector');
                [d,ix] = sort(T, 'descend');
                T = diag(d(1:opts.rank));
                Q = Q(:,ix(1:opts.rank));
            else
                [Q,T] = eigs(W, n, opts.rank, 'la', eigopts);
            end

            jmax = 1;

            track.numrestarts = track.numrestarts + 1;
        end
        
        % Get approximate top eigenvector
        u = Q(:,jmax);
        
        i1 = rows(k,1);
        i2 = rows(k,2);
        gk = (prob.A(i1,:)*u).^2.*Zval(k,1) + (prob.A(i2,:)*u).^2.*Zval(k,2);
        
        beta = opts.stepsize*gk;
        z(k) = z(k) - opts.stepsize*gk;
        
        % Update sampling probabilities        
        if strcmpi(strat, 'greedy')
            w1 = z(k)*Zval(k,1);
            w2 = z(k)*Zval(k,2);
            for j=1:blocklen
                nonunifsample('update', w1(j)^(alpha), i1(j));
                nonunifsample('update', w2(j)^(alpha), i2(j));
            end
        end
        
        % Update low-rank approximation (will clean up later...)
        l = size(T,1);
        T = [T zeros(l,2*blocklen);
            zeros(blocklen,l) -diag(beta.*Zval(k,1)) zeros(blocklen);
            zeros(blocklen,l) zeros(blocklen) -diag(beta.*Zval(k,2))];
        [Q,R] = qr([Q prob.A(i1,:)' prob.A(i2,:)'],0);
        
        T = R*T*R';
        T = (T+T')/2;
        [v,d] = eig(T);
                    
        v = real(v); d = real(d);
        ix = abs(diag(d)) > 1e-8;
        v = v(:,ix); d = d(ix,ix);
        ll = min(sum(ix), opts.rank);
        [~,p] = sort(diag(d), 'descend');
        
        p = p(1:ll);
        T = d(p,p);
        Q = Q*v(:,p); Q = Q(:,1:ll);
        jmax = 1;
        
        
        track.obj(iter) = T(jmax);
    end
    
    opts.stepsize = opts.stepsize*opts.stepsize_decay;
        
    stop = callback(Z*z+ybar, e, 0, u, prob, callbackopts);
	
    if stop
        break;
    end 
end

% Return final y
y = Z*z + ybar;
q = Q(:,1);
end