% This script should be run from the root repo directory
clc
clear

%% Define problem
imagename = 'logo_ubc'; n = 5220; m = 15000;

%% Choose methods to use
% Available methods: projgrad, redgrad, coord, wf
method = 'coord';
% Available sample types for projgrad and redgrad:
%    full, pos, topk
sampletype = 'full';

%% Set options for initialization method (user defined)

% Use explicit matrix or matrix operator for eigenvalue computation?
% Aka: eig vs. eigs
opts.explicit = false;

callbackopts.checkperiod = 100;
%callbackopts.saveperiod = 1;
% Compute image during callback (aka compute top eigenvector)?
% Can be set to false for wf/projgrad/redgrad using full sampling
% Other methods/configurations need recomputation
callbackopts.recompute = false;
if strcmpi(method,'wf')
    callbackopts.recompute = false;
end

% Set the struct for applying measurement operator
opts.opAopts = struct;
opts.opAopts.type = sampletype;
opts.opAopts.total_samples = m/10;
opts.opAopts.alpha = 2;
opts.opAopts.symm = true;
opts.opAopts.explicit = opts.explicit;

% Total iterations (epochs for coord)
opts.maxiter = 10;
opts.stepsize = 1e1;%1e-2
opts.stepsize_decay = 1;

% Coordinate descent specific
opts.rank = 50;
opts.recycle = 100*m;
opts.sample_strat = 'unif';
opts.alpha = 0.2;
opts.blocklen = 10;

%% Set options for Wirtinger Flow
wfopts.maxiter = 3000;
wfopts.stepsize = 1e3;
wfopts.stepsize_decay = 1;

wfopts.explicit = opts.explicit;

%% Set some more stuff
eigopts.issym = true; 

callbackopts.explicit = opts.explicit;
callbackopts.eigopts = eigopts;
opts.callbackopts = callbackopts;

opts.callback = @(y, iter, obj, vmax, prob, opts) callback(y, iter, obj, vmax, prob, opts);

%% Load problem and start solving
load(sprintf('data/%s_prob_n%d_m%d.mat',imagename,n,m), 'prob')

figure;
subplot(2,1,1);
img = prob.orig;
imshow(img/max(max(abs(img))));

fig = subplot(2,1,2);
opts.callbackopts.fig = fig;

[~,imax] = max(prob.b);
prob.b = [prob.b(imax); prob.b(1:imax-1); prob.b(imax+1:end)];
prob.A = [prob.A(imax,:); prob.A(1:imax-1,:); prob.A(imax+1:end,:)];

prob.A = prob.A/sqrt(m);
prob.b = prob.b/m;

%% Create figure and start running
tic;
% projected gradient descent
if strcmpi(method,'projgrad')
    [y,track] = projgrad(prob,opts);
end

% Reduced gradient descent
if strcmpi(method,'redgrad')
    opts.sampling_scheme.type = sampletype;

    [y,track] = reduced_grad(prob,opts);
end

% Coordinate descent
if strcmpi(method,'coord')
    [y,track] = coorddescent(prob,opts);
end
toc;
fprintf('Finished initialization\n');

%% wirtinger flow
callbackopts.fig = fig;
callbackopts.checkperiod = 1;
callbackopts.recompute = false;
wfopts.callbackopts = opts.callbackopts;

opA_opts.explicit = opts.explicit;
opA_opts.type = 'full';

[W,~] = opA(prob.A,y,true,opA_opts);

if opts.explicit
    [V,D] = eig(W,'vector');
    [~, jmax] = max(D);
    u0 = V(:,jmax);
else
    [u0,lambda] = eigs(W,n,1,'la',eigopts);
end
Au = (prob.A*u0).^2;
u0 = u0*sqrt(((Au'*prob.b)/(norm(Au)^2)));

wfopts.callback = @(y, iter, obj, vmax, prob, opts) callback(y, iter, obj, vmax, prob, opts);
fprintf('Starting Wirtinger-Flow...\n');
tic;
[u,track] = wirtinger_flow(prob, wfopts, u0);
toc;

%% Final image
img = u;
img = sign(img(:)'*prob.orig(:))*img;
img = reshape(img, size(prob.orig));
figure
imshow(img/max(max(abs(img))));
drawnow

%% Define callback function
function stop = callback(y, iter, obj, vmax, prob, opts)
 
    if mod(iter,opts.checkperiod) == 0
        fprintf('%d\t%e\n', iter, obj)
    end
    stop = 0;
    return

    prob.obj = obj;
    prob.show = true;
    prob.fig = opts.fig;
    
    if opts.recompute
        [img,obj] = primal_from_dual(prob, y, opts.explicit);
        
        prob.obj = obj;
        %img = sign(img(:)'*prob.orig(:))*img;
    else
        img = vmax;
    end
    img = sign(img(:)'*prob.orig(:))*img;
    
    if iter == 1
        fprintf('iter\tobj\n');
    end
    if mod(iter,opts.checkperiod) == 0
        fprintf('%d\t%e\n', iter, obj)
    end
    
    if mod(iter, opts.checkperiod) == 0
        img = reshape(img, size(prob.orig));
        imshow(img/max(max(abs(img))), 'Parent', opts.fig);
        title(obj, 'Parent', opts.fig)
        drawnow
    end
    
    stop = 0;
end