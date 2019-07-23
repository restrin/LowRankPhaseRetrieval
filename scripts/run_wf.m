% This script should be run from the root repo directory
clc
clear

%% Define problem
imagename = 'musicnote'; n = 121; m = 300;

%% Load problem and start solving
load(sprintf('data/%s_prob_n%d_m%d.mat',imagename,n,m), 'prob')

figure;
subplot(2,1,1);
img = prob.orig;
imshow(img/max(max(abs(img))));

fig = subplot(2,1,2);

[~,imax] = max(prob.b);
prob.b = [prob.b(imax); prob.b(1:imax-1); prob.b(imax+1:end)];
prob.A = [prob.A(imax,:); prob.A(1:imax-1,:); prob.A(imax+1:end,:)];

prob.A = prob.A/sqrt(m);
prob.b = prob.b/m;

y = prob.b/(prob.b'*prob.b);

%% Set options
wfopts.stepsize = 1e2;
wfopts.stepsize_decay = 1;
wfopts.explicit = false;
wfopts.maxiter = 500;

%%
callbackopts.fig = fig;
callbackopts.checkperiod = 1;
wfopts.callbackopts = callbackopts;

wfopts.callback = @(y, iter, obj, vmax, prob, opts) callback(y, iter, obj, vmax, prob, opts);
fprintf('Starting Wirtinger-Flow...\n');
[u,track] = wirtinger_flow(prob, wfopts, []);

%% Define callback function
function stop = callback(y, iter, obj, vmax, prob, opts)

    prob.obj = obj;
    prob.show = true;
    
    img = sign(vmax'*prob.orig(:))*vmax;
    
    if iter == 1
        fprintf('iter\tobj\n');
    end
    if mod(iter,opts.checkperiod) == 0
        fprintf('%d\t%f\n', iter, obj)
    end
    
    if mod(iter, opts.checkperiod) == 0
        img = reshape(img, size(prob.orig));
        imshow(img/max(max(abs(img))), 'Parent', opts.fig);
        title(obj, 'Parent', opts.fig)
        drawnow
    end
    
    stop = 0;
end