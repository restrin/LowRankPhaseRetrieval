
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
n = 4824; m = 10000;
probtype = 'sdp';


mvec = [10000,12500,15000,17500]




figure(1)
clf

for mi = 1:length(mvec)
    m = mvec(mi)
    
    load(sprintf('phaseretrieval/%s_prob_n%d_m%d.mat',imagename,n,m), 'prob')
    
    plotopts.orig = prob.orig;
    
    %% methods wirtinger
    %
    %     clear track
    %     load(sprintf('phaseretrieval/%s_prob_n%d_m%d_wirtinger_random.mat',imagename,n,m), 'u','track','ustart')
    %     imrec = recover_image(u,plotopts);
    %     subplot(6,4,1+4*(mi-1))
    %     imshow(imrec)
    %     title(track.overhead)
    %
    %
    %     clear track
    %     load(sprintf('phaseretrieval/%s_prob_n%d_m%d_wirtinger_fast.mat',imagename,n,m), 'u','track','ustart')
    %     imrec = recover_image(u,plotopts);
    %     subplot(6,6,2+6*(mi-1))
    %     imshow(imrec)
    %     title(track.overhead)
    %
    
    clear track
    load(sprintf('phaseretrieval/%s_prob_n%d_m%d_wirtinger.mat',imagename,n,m), 'u','track','ustart')
    imrec = recover_image(u,plotopts);
    subplot(length(mvec),4,4*(mi-1)+1)
    imshow(imrec)
    refobj = track.obj(end);
    
            title(sprintf('(%.2f,  %.2f)',track.obj(end),track.overhead),'fontsize',20)

    
    
    %% wirtingerflow
    for methi = 1:length(methodvec)
        method = methodvec{methi};
        if strcmpi(method,'projgrad') ||  strcmpi(method,'gradinnull')
            sampletype = 'topk';
        end
        
        if strcmpi(method,'coorddesc')
            sampletype = 'uniform';
        end
        
        
            try
        
        clear track
        probstring = sprintf('phaseretrieval/%s_prob_n%d_m%d_%s_%s_%s.mat',imagename,n,m,probtype,sampletype,method);
        load(probstring,  'y','u','track')
        %         figure(methi*100+mi)
        %         clf
        %         plot(track.runtime,track.obj)
        
        
            
            time = track.runtime(end) + track.overhead;
            
            wirtprobstring = sprintf('phaseretrieval/%s_prob_n%d_m%d_%s_%s_wirtingerinit_%s.mat',imagename,n,m,probtype,sampletype,method);
            clear u
            clear track
            load(wirtprobstring,  'ustart','u','wfopts','track')
            
            time = time + track.overhead2;
            imrec = recover_image(u,plotopts);
            subplot(length(mvec),4,4*(mi-1)+1+methi)
            imshow(imrec)      
            title(sprintf('(%.2f,  %.2f)',track.obj(end),time),'fontsize',20)

            end
        
        
    end
end


