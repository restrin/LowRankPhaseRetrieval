% This script should be run from the root repo directory
clc
clear

%% Define problem
imagename = 'logo_ubc'; n = 1276; m = 3828;

%% Choose methods to use
% Available methods: projgrad, redgrad, coord, wf
methodvec = [{'coord'}];
% Available sample types for projgrad and redgrad:
%    full, pos, topk
sampletypevec = [{'full'}];

%% Set options (user defined)

callbackopts.checkperiod = 1;
callbackopts.saveperiod = 1;
% Compute image during callback (aka compute top eigenvector)?
% Can be set to false for wf/projgrad/redgrad using full sampling
% Other methods/configurations need recomputation
callbackopts.recompute = true;

% Set the sampling struct
opts.sampling_scheme = struct;
opts.sampling_scheme.total_samples = m/10;
opts.sampling_scheme.alpha = 2;
opts.sampling_scheme.symm = true;

% Use explicit matrix or matrix operator for eigenvalue computation?
% Aka: eig vs. eigs
opts.explicit = true;

% Total iterations (epochs for coord)
opts.maxiter = 50;
opts.stepsize = 1e6;
opts.stepsize_decay = 1;

% Coordinate descent specific
opts.rank = 5;
opts.recycle = inf;
opts.sample_strat = 'greedy';
opts.alpha = 0.4;
opts.blocklen = 1;

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

%% Create figure and start running
figure;

for mti = 1:length(methodvec)
    method = methodvec{mti};
    for sti = 1:length(sampletypevec)
        sampletype = sampletypevec{sti};

        %% wirtinger flow
        if strcmpi(method,'wf')
            [u,track] = wirtinger_flow(prob, opts, []);
        end
        
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
        fprintf('%d\t%e\n', iter, obj)
    end
    
    if mod(iter, opts.checkperiod) == 0
        recover_image(img, prob);
    end
    
    stop = 0;
end