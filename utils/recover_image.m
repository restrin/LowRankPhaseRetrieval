function im_rec = recover_image(u, opts)
% Options:
%   opts.orig   - original image
%   opts.fignum - figure number
%   opts.obj    - objective value
%   opts.n1, opts.n2 - image size

    if isfield(opts,'orig')
        [n1,n2] = size(opts.orig);
    else
        n1 = opts.n1;
        n2 = opts.n2;
    end
    im_rec = reshape(u,n1,n2);
        
    im_rec = im_rec - min(im_rec(:));
    im_rec = im_rec / max(im_rec(:));

    im_orig = opts.orig;
    im_orig = im_orig - min(im_orig(:));
    im_orig = im_orig /max(.0001, max(im_orig(:)));


    if norm(im_orig - im_rec,'fro') > norm(im_orig -(1- im_rec),'fro')
        im_rec = 1-im_rec;
    end
    
    if isfield(opts,'fignum')
        figure(opts.fignum)
        clf
        if isfield(opts,'orig')
            subplot(2,1,1)
            imshow(im_orig)
            subplot(2,1,2)
        end
        imshow(im_rec)
        title(opts.obj)
        drawnow
    end
end