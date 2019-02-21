
clc
clear
addpath ../utils
addpath ../methods


% imagename = 'simple';
% n = 621;
% m = 2500;

% imagename = 'musicnote';
% n = 121;
% m = 500;



methodvec = [{'projgrad'},{'gradinnull'},{'coorddesc'}];



imagename = 'tree_sm';
n = 4824;
probtype = 'sdp';

opaopts.type = 'full';



% Set the sampling struct
opts.sampling_scheme = struct;
opts.sampling_scheme.alpha = 2;
opts.sampling_scheme.symm = true;

% Use explicit matrix or matrix operator for eigenvalue computation?
% Aka: eig vs. eigs
opts.explicit = false;

% Total iterations (epochs for coord)
opts.stepsize_decay = 0.99;

% Coordinate descent specific
opts.sample_strat = 'unif';
opts.alpha = 0.1;


mvec = [7500,10000,12500,15000,17500]



wfopts.stepsize_decay = 1;
wfopts.explicit = false;



wfopts.stepsize = 1e-6;
wfopts.maxiter =  10000;
wfopts.checkperiod = 1000;


figure(1)
clf

for mi = [1,2,3,5]%4%[4,5]%1:length(mvec)
    m = mvec(mi)
    load(sprintf('phaseretrieval/%s_prob_n%d_m%d.mat',imagename,n,m), 'prob')
    
    plotopts.orig = prob.orig;
    
    
    %% methods wirtinger random
    %     try
    %
    %         clear track
    %         load(sprintf('phaseretrieval/%s_prob_n%d_m%d_wirtinger_random.mat',imagename,n,m), 'u','track','ustart')
    %         track;
    %     catch
    %         tic
    %         ustart = randn(n,1);
    %         ustart = ustart / norm(ustart);
    %         overhead = toc;
    %         [u,track] = wirtinger_flow(prob,wfopts,ustart/10);
    %         track.overhead = overhead;
    %         save(sprintf('phaseretrieval/%s_prob_n%d_m%d_wirtinger_random.mat',imagename,n,m), 'u','track','wfopts','ustart')
    %
    %     end
    %         imrec = recover_image(u,plotopts);
    %         figure(m+2)
    %         clf
    %         subplot(2,1,1)
    %         imshow(imrec)
    %         title(m)
    %         subplot(2,1,2)
    %         semilogy(track.runtime,track.obj,'marker','.');
    %         drawnow
    
    %% methods wirtinger
    try
        sdfg
        clear track
        load(sprintf('phaseretrieval/%s_prob_n%d_m%d_wirtinger.mat',imagename,n,m), 'u','track','ustart')
        track;
        if length(track.obj) ~= 10000
            adsf
        end
    catch
        tic
        bA = repmat(prob.b,1,n);
        W = prob.A'*bA;
        W = (W+W')/2;
        optseig.maxiter = 1;
        [ustart,~] = eigs(W,1,'la',optseig);
        overhead = toc;
        [u,track] = wirtinger_flow(prob,wfopts,ustart/10);
        track.overhead = overhead;
        save(sprintf('phaseretrieval/%s_prob_n%d_m%d_wirtinger.mat',imagename,n,m), 'u','track','wfopts','ustart')
        
    end
%     imrec = recover_image(u,plotopts);
%     figure(mi*100)
%     clf
%     subplot(2,1,1)
%     imshow(imrec)
%     title(m)
%     subplot(2,1,2)
%     semilogy(track.runtime,track.obj,'marker','.');
%     drawnow
    
    
    %% methods fastwirtinger
    %     try
    %
    %         clear track
    %         load(sprintf('phaseretrieval/%s_prob_n%d_m%d_wirtinger_fast.mat',imagename,n,m), 'u','track','ustart')
    %         track;
    %     catch
    %         tic
    %         bA = repmat(prob.b,1,n);
    %         W = prob.A'*bA;
    %         W = (W+W')/2;
    %         ustart = randn(n,1); ustart = ustart/ norm(ustart);
    %         for iter = 1:5
    %             ustart = W*ustart; ustart = ustart/ norm(ustart);
    %         end
    %         overhead = toc;
    %         [u,track] = wirtinger_flow(prob,wfopts,ustart/2);
    %         track.overhead = overhead;
    %         save(sprintf('phaseretrieval/%s_prob_n%d_m%d_wirtinger_fast.mat',imagename,n,m), 'u','track','wfopts','ustart')
    %
    %     end
    %         imrec = recover_image(u,plotopts);
    %         figure(m+1)
    %         clf
    %         subplot(2,1,1)
    %         imshow(imrec)
    %         title(m)
    %         subplot(2,1,2)
    %         semilogy(track.runtime,track.obj,'marker','.');
    %         drawnow
    
    
    %% methods projgraddesc
    opts.decay = 1.1;
    sampletypevec = [{'full'},{'topk'}];
    for s = 2%1:length(sampletypevec)
        sampletype = sampletypevec{s}
        probstring = sprintf('phaseretrieval/%s_prob_n%d_m%d_%s_%s_projgrad.mat',imagename,n,m,probtype,sampletype);
        try
            asdf
            clear y track
            load(probstring,  'y','u','track')
            y;
        catch
            opts.sampling_scheme.total_samples = 1000;
            opts.sampling_scheme.type = sampletype;
            opts.sampling_scheme.alpha = 2;
            
            opts.stepsize = 1e-12;
            opts.maxiter = 300;
            
            opts.checkperiod = 100;
            [y,u,track] = projgrad(prob,opts);
            save(probstring, 'y','u','opts','track')
        end
        
%         imrec = recover_image(u,plotopts);
%         figure(100*mi)
%         clf
%         subplot(2,1,1)
%         imshow(imrec)
%         title(m)
%         subplot(2,1,2)
%         semilogy(track.runtime,track.obj,'marker','.');
%         drawnow
    end
    
    %% methods gradinnull
    
    opts.decay = 1;
    sampletypevec = [{'full'},{'topk'}];
    for s = 2%1%:length(sampletypevec)
        sampletype = sampletypevec{s}
        probstring = sprintf('phaseretrieval/%s_prob_n%d_m%d_%s_%s_gradinnull.mat',imagename,n,m,probtype,sampletype);
        
        
        try
            asdf
            clear y track
            load(probstring,  'y','u','track')
            y;
            
        catch
            
            opts.sampling_scheme.total_samples = 1000;
            opts.sampling_scheme.type = sampletype;
            opts.sampling_scheme.alpha = 2;
            
            opts.stepsize = 1e-7;
            opts.maxiter = 200;
            
            opts.checkperiod = 100;
            
            
            [y,u,track] = grad_in_null(prob,opts);
            save(probstring, 'y','u','opts','track')
        end
        
%         imrec = recover_image(u,plotopts);
%         figure(1000*mi)
%         clf
%         subplot(2,1,1)
%         imshow(imrec)
%         title(m)
%         subplot(2,1,2)
%         semilogy(track.runtime,track.obj,'marker','.');
%         drawnow
    end
    
    
    %% methods coorddescent
    
    
    opts.decay = 1.1;
    opts.recycle = inf;
    opts.blocklen = 5;
    opts.rank = 100;
    opts.checkperiod = 1000;
    opts.maxepoch = 2;
    sampletypevec = [{'full'},{'uniform'}];
    for s = 2%:length(sampletypevec)
        sampletype = sampletypevec{s}
        
        probstring = sprintf('phaseretrieval/%s_prob_n%d_m%d_%s_%s_coorddesc.mat',imagename,n,m,probtype,sampletype);
        try
            
            clear y track
            load(probstring,  'y','u','track')
            y;
            
        catch
            
            opts.stepsize = 1e-6;
            
            [y,u,track] = coorddescent(prob,opts);
            
            save(probstring, 'y','u','opts','track')
        end
        
        
        imrec = recover_image(u,plotopts);
        figure(10000*s+m)
        clf
        subplot(2,1,1)
        imshow(imrec)
        title(m)
        subplot(2,1,2)
        semilogy(track.runtime,track.obj,'marker','.');
        drawnow
    end
    
    
    %% wirtingerflow
    for methi = 1:2%3%1:length(methodvec)
        method = methodvec{methi};
        if strcmpi(method,'projgrad') ||  strcmpi(method,'gradinnull')
            sampletypevec = [{'full'},{'topk'}];
        end
        
        if strcmpi(method,'coorddesc')
            sampletypevec = [{'full'},{'uniform'}];
        end
        
        for s = 2
            sampletype= sampletypevec{s};
            
            wirtprobstring = sprintf('phaseretrieval/%s_prob_n%d_m%d_%s_%s_wirtingerinit_%s.mat',imagename,n,m,probtype,sampletype,method);
            try
                asdf
                clear u track
                load(wirtprobstring,  'ustart','u','track')
                u;
                
                if length(track.obj) ~= 10000
                    adsf
                end
            catch
                
                probstring = sprintf('phaseretrieval/%s_prob_n%d_m%d_%s_%s_%s.mat',imagename,n,m,probtype,sampletype, method);
                
                clear y track u
                load(probstring,  'y','u','track')
                tic
                opAopts.type = 'full';
                eigsopt.maxiter = 1;
                W = opA(prob.A, y, true, true, false, opAopts);
                [u,~] = eigs((W+W')/2,1,'la',eigsopt);
                overhead2 = toc;
                ustart = u;
                [u,track] = wirtinger_flow(prob,wfopts,ustart/10);
                track.overhead2 = overhead2;
                save(wirtprobstring,  'ustart','u','wfopts','track')
            end
            
            
            imrecstart = recover_image(ustart,plotopts);
            imrec = recover_image(u,plotopts);
            figure(methi*100+s*10)
            clf
            
            subplot(3,1,1)
            imshow(imrecstart)
            title(m)
            subplot(3,1,2)
            imshow(imrec)
            title(m)
            subplot(3,1,3)
            semilogy(track.runtime,track.obj,'marker','.');
            drawnow
        end
    end
    
end
