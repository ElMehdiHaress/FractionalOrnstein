function B = fbm1d(H, M, T)
%FBM1D  Fractional Brownian motion path via Davies–Harte.
%   B = FBM1D(H, M, T) returns a column vector of length M+1 representing
%   a fractional Brownian motion {B_t^H}_{t in [0,T]} sampled on a uniform
%   grid with M steps (dt = T/M). B(1)=0.
%
%   Inputs:
%     H : Hurst parameter, 0 < H < 1
%     M : number of time steps (positive integer)
%     T : time horizon > 0
%
%   Output:
%     B : (M+1) x 1 vector, B(1)=0
%
%   Notes:
%     - Uses the Davies–Harte circulant embedding to generate fractional
%       Gaussian noise (fGn), then cumulatively sums to fBm.
%     - For H = 0.5 this reduces to standard Brownian motion.
%
%   Example:
%     % match your call in eta_quantities:
%     % B = fbm1d(H, N*n, n*h);
%     B = fbm1d(0.4, 10000, 10);

    if nargin < 3, T = 1; end
    if ~(H > 0 && H < 1)
        error('fbm1d: H must be in (0,1).');
    end
    if ~(isscalar(M) && M == floor(M) && M > 0)
        error('fbm1d: M must be a positive integer.');
    end
    if ~(isscalar(T) && T > 0)
        error('fbm1d: T must be positive.');
    end

    % --- Autocovariance of fractional Gaussian noise (fGn)
    % gamma(k) = 0.5*(|k-1|^(2H) - 2|k|^(2H) + |k+1|^(2H)), k = 0..M-1
    k      = 0:(M-1);
    gamma0 = 0.5*((abs(k-1)).^(2*H) - 2*(abs(k)).^(2*H) + (abs(k+1)).^(2*H));

    % Circulant embedding vector (length 2M)
    % r = [gamma(0..M-1), 0, gamma(M-1..1)]
    r  = [gamma0, 0, gamma0(M:-1:2)];
    L  = 2*M;

    % Eigenvalues of the circulant covariance (should be >= 0)
    lam = real(fft(r));
    % Numerical floor for tiny negatives due to roundoff:
    tol = 1e-12 * max(lam);
    if any(lam < -tol)
        warning('fbm1d:DaviesHarte', ...
                'Negative eigenvalues encountered (min %.3g). Clamping to 0.', min(lam));
    end
    lam = max(lam, 0);

    % --- Generate complex Gaussian vector with required spectrum
    wk = zeros(L, 1) + 0i;

    % j = 0 term
    wk(1) = sqrt(lam(1)/(2*L)) * randn();

    % j = M term
    wk(M+1) = sqrt(lam(M+1)/(2*L)) * randn();

    % Pairs j = 1..M-1 with conjugate symmetry
    for j = 2:M
        a = randn(); b = randn();
        c = (a + 1i*b);
        scale = sqrt(lam(j)/(4*L));
        wk(j)         = scale * c;
        wk(L - j + 2) = conj(scale * c);
    end

    % Inverse FFT to get fGn (up to scaling)
    Z   = fft(wk);           % complex
    fGn = real(Z(1:M));      % length M real-valued fGn sample

    % Scale to fBm on [0,T]:
    % Standard Davies–Harte convention: fBm = T^H * M^{-H} * cumsum(fGn)
    fBm = (T^H) * (M^(-H)) * cumsum(fGn);

    % Prepend 0 to get B(1)=0, then B(2..M+1)
    B = [0; fBm(:)];
end
