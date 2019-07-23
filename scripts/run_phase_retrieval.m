% This script should be run from the root repo directory
clc
clear

global truetrack
truetrack.obj = [];

%% Define problem
imagename = 'logo_ubc'; n = 1276; m = 3828;

%% Choose methods to use
% Available methods: projgrad, redgrad, coord, wf
method = 'projgrad';
% Available sample types for projgrad and redgrad:
%    full, pos, topk
sampletype = 'full';

%% Set options (user defined)

% Use explicit matrix or matrix operator for eigenvalue computation?
% Aka: eig vs. eigs
opts.explicit = true;

callbackopts.checkperiod = 1;
%callbackopts.saveperiod = 1;
% Compute image during callback (aka compute top eigenvector)?
% Can be set to false for wf/projgrad/redgrad using full sampling
% Other methods/configurations need recomputation
callbackopts.recompute = true;
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
opts.maxiter = 300;
opts.stepsize = 5*1e-2;
opts.stepsize_decay = 1;

% Coordinate descent specific
opts.rank = 10;
opts.recycle = 10*m;
opts.sample_strat = 'greedy';
opts.alpha = 0.2;
opts.blocklen = 1;

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

%% Start running

% wirtinger flow
if strcmpi(method,'wf')
    [u,track] = wirtinger_flow(prob, opts, []);
end

% projected gradient descent
if strcmpi(method,'projgrad')
    [y,track] = projgrad(prob,opts);
end

% Reduced gradient descent
if strcmpi(method,'redgrad')
    opts.sampling_scheme.type = sampletype;

    [y,track] = reducedgrad(prob,opts);
end

% Coordinate descent
if strcmpi(method,'coord')
    [y,track] = coorddescent(prob,opts);
end

%% Define callback function
function stop = callback(y, iter, obj, vmax, prob, opts)

    global truetrack

    prob.obj = obj;
    prob.show = true;
    prob.fig = opts.fig;
    
    if opts.recompute
        [img,obj] = primal_from_dual(prob, y, opts.explicit);
        
        truetrack.obj = [truetrack.obj; obj];
        
        prob.obj = obj;
        img = sign(img(:)'*prob.orig(:))*img;
    else
        img = vmax;
    end
    
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