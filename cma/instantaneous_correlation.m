function sigma = instantaneous_correlation(x,y,filter_type,eta,delta)

% SIGMA = INSTANTANEOUS_CORRELATION(X,Y,FILTER_TYPE,ETA,DELTA) computes
% the instantaneous correlation coefficient SIGMA between signals X and
% Y. FILTER_TYPE specifies whether a causal (FILTER_TYPE = 'ss') or a
% non-causal (FILTER_TYPE = 'ds') filter should be used in computing the
% correlation (if not provided, it defaults to 'ds'). ETA specifies the
% filter sensitivity and must be a number between 0 and 1 (if not provided,
% it defaults to 0.5). DELTA specifies the range of offsets between the two
% signals; correlations are computed between x(k) and y(k+i) for all values
% of i in the range [-DELTA:+DELTA].
% 
% SIGMA is returned as a sparse matrix. SIGMA(I,J) is the correlation between
% "points" X(I) and Y(J). Thus, if signals X and Y have N samples each (they
% must have the same length), the full correlation will be a NxN matrix. This
% matrix can become rapidly unmanageable as the length of the input signals
% grow. However, we are usually not interested in computing the correlation
% between all possible pairs of points in the two signals. Rather, most times
% we are more interested in computing the correlation only between neighboring
% points in the two signals. DELTA is the parameter that controls the size of
% this vicinity.

% IMPORTANT: SIGMA(I,J) IS THE CORRELATION BETWEEN X(I) AND Y(J)

% References:
% 
% A. V. Barbosa, H. C. Yehia and E. Vatikiotis-Bateson, "Algorithm for
% Computing Spatiotemporal Coordination", Proceedings of the International
% Conference on Auditory-Visual Speech Processing (AVSP 2008), pp. 131-136,
% Sep. 2008.
% 
% R. M. Aarts, R. Irwan and A. J. E. M. Janssen, "Efficient Tracking
% of the Cross-Correlation Coefficient", IEEE Transactions on Speech
% and Audio Processing, vol. 10, no. 6, pp. 391-402, Sep. 2002.
% 
% Author: Adriano Vilela Barbosa (adriano.vilela@gmail.com)
% Copyright 2006-2014 Adriano Vilela Barbosa


% Non provided input arguments are assigned an empty matrix,
% which means they will be assigned default values later on
if (nargin<5), delta       = []; end;
if (nargin<4), eta         = []; end;
if (nargin<3), filter_type = []; end;

% Input arguments which are empty matrices are assigned default values
if isempty(delta), delta = 0; end;
if isempty(eta), eta = 0.5; end;
if isempty(filter_type), filter_type = 'ds'; end;

% The input signals must be numerical vectors
if ~((isvector(x) & isnumeric(x)) & (isvector(y) & isnumeric(y)))
	error('Input signals must be numerical vectors.');
end

% The input signals must have the same length
if ~(length(x)==length(y))
	error('The input signals must have the same length.');
end

% The filter type must be either 'ss' or 'ds'
if ~ismember(filter_type,{'ss' 'ds'})
	error('The filter type must be either ''ss'' or ''ds''.');
end

% The paremeter 'eta' must be a number between 0 and 1
if ((eta<0) | (eta>1))
	error('The paremeter ''eta'' must be a number between 0 and 1.');
end

% The signals' lengths (at this point we know both signals
% have the same length)
N = length(x);

% Set the appropriate value of DELTA if it has been passed
% as either 'diag' or 'full'
if isequal(delta,'diag'), delta = 0; end;
if isequal(delta,'full'), delta = N; end;

% If DELTA is not a numerical scalar at this point,
% it has been passed as an invalid input.
if (isnumeric(delta) & isscalar(delta))
	if ~((delta >= 0) & (delta <= N) & (round(delta) == delta))
		error('DELTA must be an integer between 0 and the signal length.');
	end
else
	error('DELTA must be either ''diag'', ''full'' or a numeric scalar.');
end

% The instantaneous covariances
s_xy = instantaneous_covariance(x,y,filter_type,eta,delta);
s_xx = instantaneous_covariance(x,x,filter_type,eta,0);
s_yy = instantaneous_covariance(y,y,filter_type,eta,0);

% For both s_xx and s_yy, we are only interested in the main diagonal
s_xx = diag(s_xx);
s_yy = diag(s_yy);

% The product s_xx*s_yy', computed only for the diagonals of interest
aux = product_diagonals(s_xx,s_yy,delta);

% We need to initialize matrix SIGMA as a sparse matrix with the
% same sparsity structure as S_XY. The easy way of doing this is
% by copying S_XY onto SIGMA
sigma = s_xy;

% The indices of the elements on the diagonals of interest
indices = find(sigma);

% The instantaneous correlation coefficient
% sigma = s_xy./sqrt(s_xx*s_yy');
sigma(indices) = s_xy(indices)./sqrt(aux(indices));

% If 'delta' is zero, then only the main diagonal of 'sigma'
% contains non-zero data and thus we can return it as a vector
if (delta==0), sigma = full(diag(sigma)); end;

% --------------------------------------------------------------------------------------------- %

function s_xy = instantaneous_covariance(x,y,filter_type,eta,delta)

% Computes the instantaneous covariance between every possible
% pair of samples in the signals x and y. The covariance between
% the samples x(i) and y(j) is computed in the following way
%
% s_xy(i,j) = sum_{k=-inf}^{k=+inf} c*exp(-eta*abs(k))*x0(i-k)*y0(j-k)
%
% that is, it is a weighted mean of the product between vectors x0
% and y0 centered around samples i and j, respectively. The weights
% have an exponential decay with a time constant tau = 1/eta.

% Signals x0 and y0 are zero mean versions of signals x and y,
% respectively, where the instantaneous means have been removed.
%
% x0(k) = x(k) - x_mean(k)
%
% y0(k) = y(k) - y_mean(k)
%
% In turn, the instantaneous mean signals are computed in the following way:
%
% x_mean(i) = sum_{k=-inf}^{k=+inf} c*exp(-eta*abs(k))*x(k-i)
%
% y_mean(i) = sum_{k=-inf}^{k=+inf} c*exp(-eta*abs(k))*y(k-i)
%
% Author: Adriano Vilela Barbosa


% The filter function (single side or double side exponential)
if isequal(filter_type,'ss'), exp_filter = @exp_filter_ss; end;
if isequal(filter_type,'ds'), exp_filter = @exp_filter_ds; end;

% Ensure signals are in column vector format
x = x(:);
y = y(:);

% The signals' lengths
Nx = length(x);
Ny = length(y);

% The instantaneous means
x_mean = exp_filter(x,eta);
y_mean = exp_filter(y,eta);

% % Removes the means
% x = x-x_mean;
% y = y-y_mean;

% Compute the diagonals of interest of the matrix x*y'. This is
% equivalent to doing [A,diags] = spdiags(x*y',[-delta:delta]),
% but it is much more memory efficient, especially if X and Y
% are large vectors and DELTA is small, as only the necessary
% diagonals of matrix x*y' are computed
[xy,diags] = product_diagonals(x,y,delta);

% Compute the diagonals of interest of the matrix x_mean*y_mean'
[xy_mean,diags] = product_diagonals(x_mean,y_mean,delta);

% Put the diagonals of 'xy' in the columns of 'A'
A = spdiags(xy,diags);

% Filter the columns of 'A'
B = exp_filter(A,eta);

% Put the columns of 'B' in the diagonals of 's_xy'
s_xy = spdiags(B,diags,Nx,Ny);

% Subtract the product of the means from 's_xy'
s_xy = s_xy - xy_mean;

% --------------------------------------------------------------------------------------------- %

function y = exp_filter_ss(x,eta)

% Single side (SS) exponential filter.

% The constant c
c = 1-exp(-eta);

% The first-order linear filter
a = exp(-eta);
num = 1;
den = [1 -a];

% Filter the input signal
y = c*filter(num,den,x);

% --------------------------------------------------------------------------------------------- %

function y = exp_filter_ds(x,eta)

% Double side (DS) exponential filter.

% The constant c
c = (1-exp(-eta))/(1+exp(-eta));

% The first-order linear filter
a = exp(-eta);
num = 1;
den = [1 -a];

% The causal filter
y1 = filter(num,den,x);

% The non-causal filter
y2 = flipud(filter(num,den,flipud(x)));

% Combine the causal and non-causal filters' outputs
y = y1 + y2 - x;

% Normalize the output
y = c*y;

% --------------------------------------------------------------------------------------------- %

function [xy,diags] = product_diagonals(x,y,delta)

% [XY,DIAGS] = PRODUCT_DIAGONALS(X,Y,DELTA), where X and Y are column vectors,
% computes the diagonals of the product matrix x*y' in the range [-DELTA:DELTA].
% The computed diagonals are returned in the columns of matrix XY, whereas their
% indices are returned in the vector DIAGS.
%
% Given two column vectors with the same length, this function computes some
% diagonals of the product matrix xy=x*y'. More specifically, the diagonals
% in the range [-delta:delta] are computed. This function is useful when we
% are computing the 2D instantaneous correlation between the signals x and y,
% but we are interested only in a specific range of deltas. The advantage of
% using this function is memory efficiency, as only the necessary diagonals
% are computed. Computing the product x*y' and then taking the diagonals of
% interest can be highly memory inefficient, especially if X and Y are long
% signals and the value of DELTA is small.
%
% Author: Adriano Vilela Barbosa

% The input signals must have the same length
if ~(length(x)==length(y))
	error('The input signals must have the same length.');
end

% Ensure signals are in column vector format
x = x(:);
y = y(:);

% The signals' lengths
N = length(x);

% Initialize matrix xy
xy = zeros([N 2*delta+1]);

% The center column of matrix xy
center_column = delta + 1;

% The main diagonal (center column)
xy(:,center_column) = x.*y;

% The remaining diagonals
for m=1:delta

	% The lower diagonals (columns to the left)
	xy(1:N-m,center_column-m) = x(m+1:N).*y(1:N-m);

	% The upper diagonals (columns to the right)
	xy(m+1:N,center_column+m) = y(m+1:N).*x(1:N-m);

end

% The indices of the computed diagonals, in Matlab notation
diags = -delta:delta;

% Put the columns of XY in the diagonals of XY
xy = spdiags(xy,diags,N,N);

% --------------------------------------------------------------------------------------------- %
