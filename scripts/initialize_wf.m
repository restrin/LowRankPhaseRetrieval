% This script should be run from the root repo directory
clc
clear

%% Define problem
imagename = 'logo_ubc'; n = 5220; m = 15000;

%% Choose methods to use
% Available methods: projgrad, redgrad, coord
methodvec = [{'redgrad'}];
% Available sample types for projgrad and redgrad:
%    full, pos, topk
sampletypevec = [{'pos'}];

%% Set options for initialization method (user defined)

callbackopts.checkperiod = 1;
callbackopts.saveperiod = 1;
% Compute image during callback (aka compute top eigenvector)?
% Can be set to false for wf/projgrad/redgrad using full sampling
% Other methods/configurations need recomputation
callbackopts.recompute = false;

% Set the sampling struct
opts.sampling_scheme = struct;
opts.sampling_scheme.total_samples = m/10;
opts.sampling_scheme.alpha = 2;
opts.sampling_scheme.symm = true;

% Use explicit matrix or matrix operator for eigenvalue computation?
% Aka: eig vs. eigs
opts.explicit = false;

% Total iterations (epochs for coord)
if strcmpi(methodvec(1), 'projgrad')
    opts.maxiter = 8;
    opts.stepsize = 5*1e-9;
elseif strcmpi(methodvec(1), 'redgrad')
    opts.maxiter = 8;
    opts.stepsize = 5*1e-5;
elseif strcmpi(methodvec(1), 'coord')
    opts.maxiter = 1;
    opts.stepsize = 1e-8;
end
opts.stepsize_decay = 1;

% Coordinate descent specific
opts.rank = 10;
opts.recycle = 5*m;
opts.sample_strat = 'greedy';
opts.alpha = 0.1;

%% Set options for Wirtinger Flow
wfopts.maxiter = 3000;
wfopts.stepsize = 5*1e-6;
wfopts.stepsize_decay = 1;

wfopts.explicit = opts.explicit;

%% Set some more stuff
eigopts.issym = true; 

callbackopts.explicit = opts.explicit;
callbackopts.eigopts = eigopts;

opts.callback = @(y, iter, obj, vmax, prob) callback(y, iter, obj, vmax, prob, callbackopts);

%% Load problem and start solving
load(sprintf('data/%s_prob_n%d_m%d.mat',imagename,n,m), 'prob')

[~,imax] = max(prob.b);
prob.b = [prob.b(imax); prob.b(1:imax-1); prob.b(imax+1:end)];
prob.A = [prob.A(imax,:); prob.A(1:imax-1,:); prob.A(imax+1:end,:)];

y = prob.b/(prob.b'*prob.b);

%% Create figure and start running
figure;

for mti = 1:length(methodvec)
    method = methodvec{mti};
    for sti = 1:length(sampletypevec)
        sampletype = sampletypevec{sti};
        
        %% projected gradient descent
        if strcmpi(method,'projgrad')
            opts.sampling_scheme.type = sampletype;
            
            [y,track] = projgrad(prob,opts);
        end
        
        if strcmpi(method,'redgrad')
            opts.sampling_scheme.type = sampletype;
            
            [y,track] = reduced_grad(prob,opts);
        end
        
        if strcmpi(method,'coord')
            [y,track] = coorddescent(prob,opts);
        end
        
    end
end
fprintf('Finished initialization\n');

%% wirtinger flow
[W,~] = opA(prob.A,y,true,opts.explicit,false,struct('type', 'full', 'symm', true));

if opts.explicit
    [V,D] = eig(W,'vector');
    [~, jmax] = max(D);
    ustart = V(:,jmax);
else
    [ustart,~] = eigs(W,prob.n,1,'la',eigopts);
end       

callbackopts.recompute = false;
wfopts.callback = @(y, iter, obj, vmax, prob) callback(y, iter, obj, vmax, prob, callbackopts);
fprintf('Starting Wirtinger-Flow...\n');
[u,track] = wirtinger_flow(prob, wfopts, ustart);

%% Define callback function
function stop = callback(y, iter, obj, vmax, prob, opts)

    prob.obj = obj;
    prob.show = true;
    
    if opts.recompute
        [img,obj] = primal_from_dual(prob, y, opts.explicit);
        prob.obj = obj;
    else
        img = vmax;
    end
    
    if iter == 1
        fprintf('iter\tobj\n');
    end
    if mod(iter,opts.checkperiod) == 0
        fprintf('%d\t%f\n', iter, obj)
    end
    
    if mod(iter, opts.checkperiod) == 0
        recover_image(img, prob);
    end
    
    stop = 0;
end