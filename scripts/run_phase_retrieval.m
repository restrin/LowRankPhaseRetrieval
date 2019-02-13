% This script should be run from the root repo directory
clc
clear

%% Define problem
imagename = 'logo_ubc'; n = 5220; m = 15000;

%% Choose methods to use
methodvec = [{'projgrad'}];
sampletypevec = [{'pos'}];

%% Set options (user defined)

callbackopts.checkperiod = 1;
callbackopts.saveperiod = 1;
callbackopts.recompute = false;

opts.sampling_scheme = struct;
opts.sampling_scheme.total_samples = m/10;
opts.sampling_scheme.alpha = 2;
opts.sampling_scheme.symm = true;

opts.explicit = false;

opts.maxiter = 50;
opts.stepsize = 1e-9;
opts.stepsize_decay = 0.99;

% Coordinate descent specific
opts.rank = 10;
opts.recycle = 5*m;
opts.sample_strat = 'greedy';
opts.alpha = 0.1;

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
            opts.saveimage_period = 1;
            opts.checkperiod = 1;
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
        fprintf('%d\t%f\n', iter, obj)
    end
    
    if mod(iter, opts.checkperiod) == 0
        recover_image(img, prob);
    end
    
    stop = 0;
end