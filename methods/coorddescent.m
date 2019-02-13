function [y,track] = coorddescent(prob,opts) 
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
%       rank
%       recycle  - how often to refresh factorization
%    Optional:
%       y0       - default: b/(b'*b)
%       explicit - default: true
%       symm     - default: false
%       callback - stop = opts.callback(y, iter, obj, vmax, prob);
%       sample_strat - unif (default), greedy
%       Q0,T0    - initial low-rank factorization

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
else
    callback = @(y, iter, vmax, prob) false;
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

%% Initialize stuff

z = zeros(m-1,1);
y = ybar;

samplestruct = struct;
samplestruct.type = 'full';
samplestruct.symm = true;

if isfield(opts, 'Q0') && isfield(opts, 'T0')
    Q = opts.Q0;
    T = opts.T0;

    [v,~] = eig(T);
    v = real(v);
    Q = Q*v;
else
    [W,~] = opA(prob.A,y,true,explicit,false,samplestruct);
    
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

%% Start solving

for e=1:opts.maxiter
    
    if strcmpi(strat,'unif')
        order = randperm(m-1);
    end
    
    for i=1:m-1      
        
        iter = (e-1)*(m-1)+i;
        
        if strcmpi(strat,'unif')
            k = order(i);
        elseif strcmpi(strat, 'greedy')
            r = rand(1);
            k = nonunifsample('sample', r);
            k = max(k-1,1);
        end
        
        if(mod(iter, opts.recycle) == 0)
            [W,~] = opA(prob.A,Z*z+ybar,true,explicit,false,samplestruct);

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
        gk = (prob.A(i1,:)*u)^2*Zval(k,1) + (prob.A(i2,:)*u)^2*Zval(k,2);

        beta = opts.stepsize*gk;
        z(k) = z(k) - opts.stepsize*gk;

        % Update sampling probabilities
        w1 = z(k)*Zval(k,1);
        w2 = z(k)*Zval(k,2);
        
        if strcmpi(strat, 'greedy')
            nonunifsample('update', w1^(alpha), i1);
            nonunifsample('update', w2^(alpha), i2);
        end

        opts.stepsize = opts.stepsize*opts.stepsize_decay;
        
        % Update low-rank approximation (will clean up later...)
        l = size(T,1);
        T = [T zeros(l,2);
            zeros(1,l) -beta*Zval(k,1) 0;
            zeros(1,l) 0 -beta*Zval(k,2)];
        [Q,R] = qr([Q prob.A(i1,:)' prob.A(i2,:)'],0);

        T = R*T*R';
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
    
    stop = callback(Z*z+ybar, e, 0, u, prob);

    if stop
        break;
    end    
end
   
% Return final y
y = Z*z + ybar;

end