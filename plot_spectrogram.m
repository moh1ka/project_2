load handel.mat
filename = 'enh_1_on_only.wma';
[y,Fs] = audioread(filename);

N = 1024;
n = 0:N-1;
s = spectrogram(x); %not working, x needs to be a array
%{
s = spectrogram(x) returns the short-time 
Fourier transform of the input signal, x.
Each column of s contains an estimate of 
the short-term, time-localized frequency content of x.
%}
spectrogram(y,'yaxis')
%{
Xtwz = zeros(N,nframes); % pre-allocate STFT output array
M = length(w);           % M = window length, N = FFT length
zp = zeros(N-M,1);       % zero padding (to be inserted)
xoff = 0;                % current offset in input signal x
Mo2 = (M-1)/2;           % Assume M odd for simplicity here
for m=1:nframes
  xt = x(xoff+1:xoff+M); % extract frame of input data
  xtw = w .* xt;         % apply window to current frame
  xtwz = [xtw(Mo2+1:M); zp; xtw(1:Mo2)]; % windowed, zero padded
  Xtwz(:,m) = fft(xtwz); % STFT for frame m
  xoff = xoff + R;       % advance in-pointer by hop-size R
end
%}