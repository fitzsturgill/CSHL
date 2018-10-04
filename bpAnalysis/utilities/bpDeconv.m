function output = bpDeconv(data, kernel, eps, window, truncate)

% Fitz Sturgill 2018, Thomas Carey 2018 (see mathworks source at bottom)
if nargin < 5
    truncate = 0;
    if nargin < 4
        window = "none";
        if nargin < 3
            eps = 0.1;
        end
    end
end
% data must be a row vector or a series of row vectors to be averaged
% kernel must be a single row vector
% eps must be a scalar, or absent to be defaulted to .1
% returns a row vector normalized so that its max value is 1
n = size(data,1); % length of data vector
L = size(data,2); % number of data vectors to process
Lk = length(kernel);
if truncate
    Lx=L-Lk+1; 
else
    Lx = L; 
end
Lx2=pow2(nextpow2(Lx)); % closed power of two to deconvolved length, makes fft faster
kernel = kernel / trapz(kernel);
switch window
    case "none"
        W = ones(1, L);
        Wk = ones(1, Lk);
    case "blackman"
        W = blackman(L)';
        Wk = blackman(Lk)';        
    case "hann"
        W = hann(L)';
        Wk = hann(Lk)';
end

H= fft(kernel, Lx2); % transform kernel, kernel is constant so this only need be done once
out = zeros(n,Lx); % instantiate process vector
for counter = 1:n
    Y=fft(data(counter,:) .* W, Lx2); % transform data with window
    % deconvolve, regularizing with epsilon
    %  X = (Y.*conj(H))./(H.*conj(H)+eps*mean(H.*conj(H))); 
    % rearranging makes more sense to me
    X = Y ./ (H + eps * mean(H.*conj(H))./conj(H)); 
    x= real(ifft(X)); % transform back to original domain
    x=x(1:1:Lx); % take valid portion of deconvolved vector
    out(counter,:)=x;
    
%     figure; hold on;
%     plot(abs(H), 'k');
%     plot(abs(eps * mean(H.*conj(H))./conj(H)), 'b');
%     plot(abs(H +  (eps * mean(H.*conj(H))./conj(H))), 'g');
%     legend('raw', 'summand', 'raw + summand = regularized')
end


%output = nanmean(out);
%averages out vector ignoring nans (same as nanmean function)
output = zeros(1, Lx);
for i = 1:Lx
    avgs = out(:,i);
    avgs = avgs(~isnan(avgs));
    output(1,i) = sum(avgs) / length(avgs);
end
%output = output / max(output); % normalization
out2 = nan(n, 220);
for i = 1:Lx
    out2(:,i) = out(:,i);
end
output = out2;

%{ 
Reference- see this function on mathworks: 
    function[T,S]=deconvolution(uz1,uz2)
    %   DECONVOLUTION OF TWO DISCRETE TIME SIGNALS IN FREQUENCY DOMAIN
    %
    %   Deconvolution of two signals is defined in the frequency domain as
    %
    %                        S(w,z) = u(z1,w)/u(z2,w)
    %
    %   where u(z1,w) and u(z2,w) are the Fourier spectra of the signal
    %   recorded at z1 and z2, respectively.
    %
    %   The above equation may become ill-conditioned when the denominator
    %   approaches zero due to data corrupted by noise. In order to eliminate
    %   the instability, a regularized format following the Tikhonov
    %   deconvolution S_e(w) is used (Tikhonov and Arsenin, 1977; Petrovic and
    %   Parolai, 2016; Wen and Kalkan, 2017):
    %
    %   Syntax:
    %       [T,S] = deconvolution(uz1,uz2);
    %
    %   Input:
    %         uz1 = a time vector of signal recorded at location z1 
    %         uz2 = a time vector of signal recorded at location z2 
    %
    %   Output:
    %          T = Time axis of deconvolved waveform (1xn)
    %          S = Deconvolved waveform (1xn)
    %
    %          for plotting, use plot(T,W); 
    %
    %   Acknowledgement:
    %
    %   In preparing this function, I benefitted and inspired from MatLAB codes
    %   of Dr. Nori Nakata of University of Oklahoma, and Hasan Ulusoy.
    %
    %   References:
    %
    %   Tikhonov, A. N. and Arsenin, V. Y. (1977). Solution of Ill-Posed
    %   Problems, Wiston/Wiley, Washington, D.C.
    %   
    %   Petrovic, B. and Parolai, S. (2016). Joint Deconvolution of Building
    %   and Downhole Strong-Motion Recordings: Evidence for the Seismic
    %   Wavefield Being Radiated Back into the Shallow Geological Layers,
    %   Bull. Seismol. Soc. Am. 106(4), doi: 10.1785/0120150326.
    %
    %   Wen, W. and Kalkan, E. (2017). "System Identification Based on
    %   Deconvolution and Cross-Correlation?An Application To A Twenty-Story
    %   Instrumented Building In Anchorage, Alaska", Bulletin of Seismological
    %   Society of America. ol. 107(2), doi: 10.1785/0120160069.
    %
    %   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER "AS IS" AND ANY
    %   EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    %   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
    %   PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER BE LIABLE
    %   FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    %   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    %   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
    %   BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
    %   WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
    %   OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
    %   ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    %
    %   Written by Dr. Erol Kalkan, P.E. (ekalkan@usgs.gov)
    %   $Revision: 1.0 $  $Date: 2016/09/01 08:00:00 $

    % Regularization parameter, which is 10 percent of the average of the power
    % spectrum of uz2
    eps = 0.1; 

    L = numel(uz1);

    % Time vector
    T = [-(L/2+1):1:(L/2-2)];

    % L-point symmetric Hann window in the column vector W
    W = hann(L); 

    % Multiply input signals, uz1 and uz2, with Hann window and take FFT in
    % order to make sure that the ends of signal in frequency domain match up
    % while keeping everything reasonably smooth; this greatly diminish
    % spectral leakage
    uz1f = fft(W.*uz1,L); 
    uz2f = fft(W.*uz2,L);

    % Compute deconvolution
    Stmp = real(ifft((uz1f.*conj(uz2f))./(uz2f.*conj(uz2f)+eps*mean(uz2f.*conj(uz2f)))));
    S = [Stmp(L/2:L); Stmp(1:L/2-1)];
    return
%}
