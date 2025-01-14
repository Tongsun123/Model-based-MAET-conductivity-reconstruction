function X = imtgvsmooth( Y, alpha, beta, nite )
%IMTGVSMOOTH performs smoothing by using TGV regularizer [4] and ADMM solver [2].
% The details of this implementation are described in [1].
%
% X = IMTGVSMOOTH( Y, alpha, beta, nite ) smooths a gray scale image 'Y', and
% outputs the smoothed image 'X'. The parameters 'alpha' and 'beta' are balancing
% weights for TGV regularization and they deal with the first and second order
% differentials respectively. 'nite' is the number of ADMM's iterations.
%
%
% REFERENCES
%
% [1] K. Shirai, M. Okuda,
%     "FFT based solution for multivariable l2 equations using KKT system
%      via FFT and efficient pixel-wise inverse calculation,"
%    in Proc. IEEE ICASSP, 2648-2652, 2014.
%
% [2] S.Boyd, N.Parikh, E.Chu, B.Peleato, J.Eckstein,
%     "Distributed Optimization and Statistical Learning via the Alternating
%      Direction Method of Multipliers,"
%      Foundations and Trends in Machine Learning, 3(1):1-122, 2011.
%
% [3] M.Tao, J. Yang,
%     "Alternating Direction Algorithms for Total Variation Deconvolution in
%      Image Reconstruction," Optimization Online.
%
% [4] K. Bredies, K. Kunisch, and T. Pock,
%     "Total generalized variation,"
%     SIAM J. Imaging Sci., 3(3), 492-526, 2010.
%
% [5] S. Ono and I. Yamada,
%     "Optimized JPEG image decompression with super-resolution interpolation
%      using multi-order total variation,"
%      in Proc. IEEE ICIP, 474-478, 2013.
%

% ----------------------------------
% parameters
% ----------------------------------
[sy,sx,sc] = size( Y );                 % Get the dimension of the image

global rho eta;
rho = 1;                                % Global parameters, Lagrange multipliers
eta = 1;                                % Global parameters, Lagrange multipliers

% ----------------------------------
% filter kernels
% ----------------------------------

% differential kernels
global FKh FKv FKh_ FKv_;
Kh = [1, -1, 0];                        % forward difference by convolution，difference in x direction
Kv = [1; -1; 0];                        % difference in y direction

% FFTed filters
FKh = psf2otf( Kh, [sy, sx] );          % Fourier transform Kh kernel function to specify the size of Fourier transform
FKv = psf2otf( Kv, [sy, sx] );          % Fourier transform Kv kernel function to specify the size of Fourier transform

FKh_ = conj(FKh);                       % for computation in the Fourier domain
FKv_ = conj(FKv);
FKl = FKh_.*FKh + FKv_.*FKv;            % Representation of Laplacian in frequency domain

global FY;
FY = fft2( Y );                         % Fourier transform the image to be denoised

%
% pixel-wise inverse by using the image-wise adjugate mathod
%
global iA11 iA12 iA13 iA21 iA22 iA23 iA31 iA32 iA33;
A11 = 1+rho*FKl;  A12 =     -rho*FKh_;  A13 =     -rho*FKv_;
A21 =  -rho*FKh;  A22 =   rho+eta*FKl;  A23 = eta*FKv.*FKh_;
A31 =  -rho*FKv;  A32 = eta*FKh.*FKv_;  A33 =   rho+eta*FKl;
%
detA = 1./ ((A11.*A22.*A33) + (A12.*A32.*A13) + (A31.*A12.*A23)...
	- (A11.*A32.*A23) - (A31.*A22.*A13) - (A21.*A12.*A33));
%
iA11 = (A22.*A33 - A23.*A32) .* detA;
iA12 = (A13.*A32 - A12.*A33) .* detA;
iA13 = (A12.*A23 - A13.*A22) .* detA;
iA21 = (A23.*A31 - A21.*A33) .* detA;
iA22 = (A11.*A33 - A13.*A31) .* detA;
iA23 = (A13.*A21 - A11.*A23) .* detA;
iA31 = (A21.*A32 - A22.*A31) .* detA;
iA32 = (A12.*A31 - A11.*A32) .* detA;
iA33 = (A11.*A22 - A12.*A21) .* detA;

clear A11 A12 A13 A21 A22 A23 A31 A32 A33;
clear FKl;


% ----------------------------------
% Initialization
% ----------------------------------
X = Y;

Z1h = diff_circ( Y, 2, 'forward' );         % 
Z1v = diff_circ( Y, 1, 'forward' );         % 
Z2h = zeros(sy,sx);                         % 
Z2d = zeros(sy,sx);                         % 
Z2v = zeros(sy,sx);                         % 

U1h = zeros(sy,sx);                         % 
U1v = zeros(sy,sx);
U2h = zeros(sy,sx);                         %
U2d = zeros(sy,sx);
U2v = zeros(sy,sx);

% ----------------------------------
% ADMM
% For simplicity, the stopping criteria of the ADMM is omitted.
% ----------------------------------
for t = 1:nite
	
	% ----------------------------------
	% solve J by ADMM and obtain X,v_h,v_v
	% ----------------------------------
	tic
	[X, Qh, Qv] = opt_X( Z1h, Z1v, U1h, U1v, Z2h, Z2d, Z2v, U2h, U2d, U2v );
	toc
	% 
	figure(100), imshow( X ,[]); colormap jet; drawnow;         % for DEBUG
	
	% ----------------------------------
	% solve Z1
	% ----------------------------------
	Xu = diff_circ( X, 2, 'forward' );              %
	Xv = diff_circ( X, 1, 'forward' );              %
	%
	T1h = (Xu - Qh);
	T1v = (Xv - Qv);
	
	Z1h = T1h + U1h;
	Z1v = T1v + U1v;

	[Z1h, Z1v] = shrinkage( alpha/rho , Z1h, Z1v );

	% ----------------------------------
	% solve Z2
	% ----------------------------------
	T2h = diff_circ( Qh, 2, 'backward' );
	T2d = diff_circ( Qh, 1, 'backward' ) + diff_circ( Qv, 2, 'backward' );
	T2v = diff_circ( Qv, 1, 'backward' );
	
	Z2h = T2h + U2h;
	Z2d = T2d + U2d;
	Z2v = T2v + U2v;

	[Z2h, Z2d, Z2v] = shrinkage( beta/eta, Z2h, Z2d, Z2v );

	% ----------------------------------
	% update U1
	% ----------------------------------
	U1h = U1h + (T1h - Z1h);
	U1v = U1v + (T1v - Z1v);
	
	% ----------------------------------
	% update U2
	% ----------------------------------
	U2h = U2h + (T2h - Z2h);
	U2d = U2d + (T2d - Z2d);
	U2v = U2v + (T2v - Z2v);
end

end

%
% -------------------------------------------------------------------------
%

function [X, Qh, Qv] = opt_X( Z1h, Z1v, U1h, U1v, Z2h, Z2d, Z2v, U2h, U2d, U2v )

global FKh FKv FKh_ FKv_;
global FY;
global iA11 iA12 iA13 iA21 iA22 iA23 iA31 iA32 iA33;
global rho eta;
% 
FZU1h = rho * fft2(Z1h - U1h);
FZU1v = rho * fft2(Z1v - U1v);
FZU2h = eta * fft2(Z2h - U2h);
FZU2d = eta * fft2(Z2d - U2d);
FZU2v = eta * fft2(Z2v - U2v);

B1 =     FY + FKh_.*FZU1h + FKv_.*FZU1v;
B2 = -FZU1h + FKh .*FZU2h + FKv .*FZU2d;
B3 = -FZU1v + FKh .*FZU2d + FKv .*FZU2v;
% 
FX  = iA11.*B1 + iA12.*B2 + iA13.*B3;
FQu = iA21.*B1 + iA22.*B2 + iA23.*B3;
FQv = iA31.*B1 + iA32.*B2 + iA33.*B3;

X  = real( ifft2(FX) );
Qh = real( ifft2(FQu) );
Qv = real( ifft2(FQv) );

end

%
% -------------------------------------------------------------------------
%

function varargout = shrinkage( t, varargin )
% 
% [X1,...,XN] = shrinkage( t, X1,...,XN );

varargout = cell( nargout, 1 );
n = length(varargin);

% pixel-wise norm
i = 1;
T = varargin{i};
N = T.*T;
for i = 2:n
	T = varargin{i};
	N = N + T.*T;
end
N = sqrt(N);
S = max(1 - t./N, 0);

for i = 1:n
	varargout{i} = S.*varargin{i};
end

end

%
% -------------------------------------------------------------------------
%

function J = diff_circ( I, dim, ori )
% Function：Circumferential difference

if dim == 1 % virtical
	switch ori
		case 'forward'
			J = diff( [I; I(1,:,:)], [], 1 );
		case 'backward'
			J = -diff( [I(end,:,:); I], [], 1 );
		otherwise
			J = diff( [I; I(1,:,:)], [], 1 );
	end
		
else % dim == 2 % horizontal
	switch ori
		case 'forward'
			J = diff( [I, I(:,1,:)], [], 2 );
		case 'backward'
			J = -diff( [I(:,end,:), I], [], 2 );
		otherwise
			J = diff( [I, I(:,1,:)], [], 2 );
	end
end

end