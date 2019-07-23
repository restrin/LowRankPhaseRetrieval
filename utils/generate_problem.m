function prob = generate_problem(imagename, outdir, m, getopt)
%GENERATE_PROBLEM   Given image, create phase-recovery problem of the form
%
%   min_X  tr(X)
%   s.t.   A(X) = b,    A is m-by-n^2 operator
%
% Input:
%   imagename   path to PNG file
%   outdir      path to output file
%   m           number of measurements
%   getopt      compute (gauge-dual) optimal dual optimum?
%
% Output:
%   prob        .mat file containing:
%     A         measurement matrix
%     b         measurements
%     dual_opt  optimal dual value
%     ysol      optimal dual solution (if getopt=true)

im = mean(double(imread(imagename))/255,3)-.5;
n = length(im(:));

% Set problem dimensions and original image
prob.m = m;
prob.n = n;
prob.orig = im;

% Create hadamard measurements
A = hadamard(2^13);
idx = randperm(size(A,2));
A = A(:,idx(1:n));
idx = randperm(size(A,1));
A = A(idx(1:m),:);
b = (A*im(:)).^2;

prob.A = A;
prob.b = b;

if getopt
    cvx_precision default
    cvx_begin sdp
        variables y(m) t
        minimize( t )
        subject to:
            t*eye(n) - prob.A'*(y*((ones(1,m)*prob.A))) == semidefinite(n)
            prob.b'*y == 1
    cvx_end
    cvx_precision default
    
    prob.ysol = y;
end

prob.dual_opt = 1/norm(prob.orig(:))^2;

% Save image
[~,fname] = fileparts(imagename);
save(sprintf('%s/%s_prob_n%d_m%d.mat',outdir,fname,n,m), 'prob');

end