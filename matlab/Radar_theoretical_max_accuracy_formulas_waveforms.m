function Radar_theoretical_max_accuracy_formulas_waveforms(varargin)
% Toolbox-free hex planar array visualization + FMCW + Range-Doppler + 3D & FFT beampatterns
% Saves all charts to results/ and writes summary.txt + summary.pdf.
% NOW SUPPORTS:
%   - Named presets: Radar_theoretical_max_accuracy_formulas_waveforms('preset','high_res_range')
%   - Name-Value overrides: Radar_theoretical_max_accuracy_formulas_waveforms('B',1.5e9,'scanAz_deg',45)
%   - Or a single struct: Radar_theoretical_max_accuracy_formulas_waveforms(optsStruct)

%% ================== DEFAULTS ==================
opt = struct();
% Feature toggles
opt.doSpectrogram    = true;
opt.doRangeDoppler   = true;
opt.doBeampattern3D  = true;
opt.doFFTBeampattern = true;
opt.useTimestampSubfolder = false;

% RF & waveform
opt.fc        = 5e9;        % Hz
opt.B         = 1e9;        % Hz (RF bandwidth)
opt.PRF       = 5e3;        % Hz
opt.Npulses   = 128;        % coherent pulses per CPI
opt.SNRdB     = 20;         % context only
opt.c0        = 3e8;        % m/s

% Array geometry
opt.N_target  = 3333;       % desired elements (or subarrays)
opt.d         = 0.5*(opt.c0/opt.fc); % spacing (≈ λ/2)

% Steering
opt.scanAz_deg = 30;
opt.scanEl_deg = 0;

% Taper
opt.taper = struct('type','raisedcos','power',1); % 1 ~ Hann-like

% FMCW demo
opt.Tchirp   = 50e-6;       % s
opt.fs       = 20e6;        % Hz
opt.N_chirps = 16;          % for spectrogram & RD demo

% Grids
opt.uv_lim   = 0.9;         % visible region cap
opt.uv_N     = 401;
opt.azSpan   = [-80 80];
opt.elSpan   = [-80 80];
opt.NAzEl    = 201;
opt.apGridPx = 512;
opt.apOversamp = 1.2;
opt.cutElForHPBW_deg = 0;

%% ================== ARG PARSING ==================
presetName = '';
if nargin==1 && isstruct(varargin{1})
    opt = merge_opts(opt, varargin{1});
elseif mod(nargin,2)==0
    for k=1:2:nargin
        key = varargin{k}; val = varargin{k+1};
        if ischar(key) || isstring(key)
            key = char(key);
            if strcmpi(key,'preset'), presetName = char(val); continue; end
            opt = setfield(opt, key, val); %#ok<SFLD>
        end
    end
else
    error('Use either (optsStruct) or name/value pairs, e.g., (''preset'',''high_res_range'',''B'',1.5e9).');
end

% Apply preset overlays
if ~isempty(presetName)
    opt = apply_preset(opt, presetName);
end

%% ================== DERIVED ==================
lambda = opt.c0/opt.fc;

%% ================== OUTPUT FOLDER ==================
script_dir = fileparts(mfilename('fullpath')); if isempty(script_dir), script_dir = pwd; end
root_results = fullfile(script_dir,'results'); if ~exist(root_results,'dir'), mkdir(root_results); end
if opt.useTimestampSubfolder, outdir = fullfile(root_results, datestr(now,'yyyymmdd_HHMMSS')); else, outdir = root_results; end
if ~exist(outdir, 'dir'), mkdir(outdir); end
fprintf('[INFO] Results folder: %s\n', outdir);
% write-check
fid = fopen(fullfile(outdir,'writecheck.txt'),'w'); if fid>0, fprintf(fid,'ok %s', datestr(now)); fclose(fid); else, error('Cannot write to results/'); end

%% ================== BUILD HEX ARRAY ==================
K = ceil( ( -3 + sqrt(9 + 12*(opt.N_target-1)) )/6 ); % invert N=1+3K(K+1)
[q, r] = hex_axial_coords(K);
[x, y] = axial_to_xy(q, r, opt.d);
if numel(x) > opt.N_target
    R = sqrt(x.^2 + y.^2); [~, idx] = sort(R, 'ascend');
    keep = idx(1:opt.N_target); x = x(keep); y = y(keep);
end
N = numel(x);

%% ================== TAPER ==================
w_amp = ones(N,1);
if strcmpi(opt.taper.type,'raisedcos')
    R  = sqrt(x.^2 + y.^2); Rn = R / max(R + eps);
    w_amp  = (0.5*(1 + cos(pi*Rn))).^opt.taper.power;
end

%% ================== STEERING ==================
u0 = sind(opt.scanEl_deg)*cosd(opt.scanAz_deg);
v0 = sind(opt.scanEl_deg)*sind(opt.scanAz_deg);
k  = 2*pi/lambda;
w   = w_amp .* exp(1j*k*(x*u0 + y*v0));  % complex weights with steering

%% ================== PLOT: RADAR FACE ==================
figure('Name','Radar Face / Element Layout','Color','w','Visible','on');
scatter(x*1e3, y*1e3, 8, 'filled');
axis equal tight; grid on;
xlabel('x [mm]'); ylabel('y [mm]');
title(sprintf('Hex Planar Array (N=%d, d=%.2f mm, steer=[%g^\\circ,%g^\\circ])', ...
      N, opt.d*1e3, opt.scanAz_deg, opt.scanEl_deg));
savefig2(outdir, 'face');

%% ================== ARRAY FACTOR (2-D u–v) ==================
u = linspace(-opt.uv_lim, +opt.uv_lim, opt.uv_N);
v = linspace(-opt.uv_lim, +opt.uv_lim, opt.uv_N);
[UU,VV] = meshgrid(u,v);
AF = zeros(size(UU));
for n = 1:N
    AF = AF + w(n) .* exp(1j*k*(x(n).*UU + y(n).*VV));
end
AFdB = 20*log10( abs(AF)/max(abs(AF(:))) + eps );
AFdB(UU.^2 + VV.^2 > 1) = -Inf;
figure('Name','2-D Beam (Array Factor in u–v)','Color','w','Visible','on');
imagesc(u,v, AFdB); axis image; set(gca,'YDir','normal');
colormap(parula); cb = colorbar; cb.Label.String = 'AF (dB)';
xlabel('u = sin(el)cos(az)'); ylabel('v = sin(el)sin(az)');
title('Array Factor (normalized, dB)'); hold on;
th = linspace(0,2*pi,720); plot(cos(th), sin(th), 'k--'); hold off;
savefig2(outdir, 'uv_beam');

%% ================== BEAM CUT (AZ, fixed EL > 0) ==================
elCut = 5; azCut = linspace(-90,90,1001);
uCut  = sind(elCut).*cosd(azCut); vCut = sind(elCut).*sind(azCut);
AF_cut = zeros(size(azCut));
for n = 1:N
    AF_cut = AF_cut + w(n) .* exp(1j*k*(x(n).*uCut + y(n).*vCut));
end
AF_cut_dB = 20*log10( abs(AF_cut)/max(abs(AF_cut)) + eps );
figure('Name','Azimuth Beam Cut','Color','w','Visible','on');
plot(azCut, AF_cut_dB, 'LineWidth', 1.5); grid on;
xlabel(sprintf('Azimuth (deg) @ EL=%g^\\circ', elCut));
ylabel('Normalized AF (dB)'); ylim([-60 0]); xlim([min(azCut) max(azCut)]);
title('Azimuth Cut of Beam');
savefig2(outdir, 'az_cut');

%% ================== 3D BEAMPATTERN (sampled az–el) ==================
if opt.doBeampattern3D
    az3 = linspace(opt.azSpan(1), opt.azSpan(2), opt.NAzEl);
    el3 = linspace(opt.elSpan(1), opt.elSpan(2), opt.NAzEl);
    [AZ, EL] = meshgrid(az3, el3);
    U3 = sind(EL).*cosd(AZ); V3 = sind(EL).*sind(AZ);
    AF3 = zeros(size(AZ));
    for n = 1:N
        AF3 = AF3 + w(n) .* exp(1j*k*(x(n).*U3 + y(n).*V3));
    end
    AF3dB = 20*log10( abs(AF3)/max(abs(AF3(:))) + eps );
    AF3dB( (U3.^2 + V3.^2) > 1 ) = -Inf;
    figure('Name','3D Beampattern (Az–El)','Color','w','Visible','on');
    surf(AZ, EL, AF3dB, 'EdgeColor','none'); view(2); axis tight;
    xlabel('Azimuth (deg)'); ylabel('Elevation (deg)');
    title('Beampattern (normalized, dB)'); colormap(parula);
    cb2 = colorbar; cb2.Label.String = 'AF (dB)'; caxis([-60 0]); grid on;
    savefig2(outdir, 'beampattern3d');
end

%% ================== FMCW CHIRP (time) & DIY SPECTROGRAM ==================
S = opt.B/opt.Tchirp; Ns = round(opt.Tchirp*opt.fs);
t  = (0:Ns-1)/opt.fs; phi = 2*pi*( 0.5*S*t.^2 ); s1  = exp(1j*phi);
wave = repmat(s1, 1, opt.N_chirps);
figure('Name','FMCW Chirp','Color','w','Visible','on');
subplot(2,1,1); plot(t*1e6, real(s1)); grid on;
xlabel('Time within chirp [\mus]'); ylabel('Re\{s\}');
title(sprintf('One LFM chirp: B=%.1f MHz, T=%.1f \\mus, fs=%.1f MHz', opt.B/1e6, opt.Tchirp*1e6, opt.fs/1e6));
subplot(2,1,2); instf = S*t; plot(t*1e6, instf/1e6); grid on;
xlabel('Time within chirp [\mus]'); ylabel('Inst. freq [MHz]'); title('Instantaneous frequency');
savefig2(outdir, 'chirp');
if opt.doSpectrogram
    win = 512; hop = 128; nfft = 2048;
    [Sabs, Fm, TM] = simple_stft(wave, win, hop, nfft, opt.fs);
    figure('Name','Chirp Spectrogram (DIY)','Color','w','Visible','on');
    imagesc(TM*1e3, Fm/1e6, 20*log10(Sabs + 1e-12)); axis xy;
    xlabel('Time [ms]'); ylabel('Frequency [MHz]');
    colormap(parula); colorbar; title('Spectrogram (rect window)');
    savefig2(outdir, 'spectrogram');
end

%% ================== SUMMARY & METRICS ==================
TCPI     = opt.Npulses/opt.PRF;
rangeRes = opt.c0/(2*opt.B);
velRes   = lambda/(2*TCPI);
A_per   = sqrt(3)/2 * opt.d^2; A_hex = N * A_per; D_eq = 2*sqrt(A_hex/pi);
HPBW_deg = 70 * (lambda / max(D_eq,eps));
theta_scan = acosd(cosd(opt.scanAz_deg).*cosd(opt.scanEl_deg));
scanLoss_dB = -10*log10(cosd(theta_scan)+eps);
HPBW_scan_deg = HPBW_deg / max(cosd(theta_scan), 1e-6);
Runamb  = opt.c0/(2*opt.PRF); vunamb  = (opt.PRF*lambda)/4;
d_limit = lambda / (1 + abs(sind(theta_scan))); okGL = opt.d <= d_limit;

summary_text = sprintf([ ...
    '==== Summary ====\\n' ...
    'Timestamp: %s\\n' ...
    'Elements: N=%d (target %d), spacing d=%.3f m (%.2f \\lambda)\\n' ...
    'Steering: az=%g deg, el=%g deg\\n' ...
    'Range resolution                 ~ %.3f m\\n' ...
    'Velocity resolution (CPI=%.3f s) ~ %.4f m/s\\n' ...
    'HPBW (boresight)                ~ %.2f deg (D_{eq}=%.2f m)\\n' ...
    'HPBW at scan (%.1f deg)         ~ %.2f deg; scan loss ~ %.2f dB\\n' ...
    'Unambiguous range                ~ %.1f m\\n' ...
    'Unambiguous velocity             ~ %.2f m/s\\n' ...
    'Grating-lobe safe at this scan? %s (d=%.3f m, limit=%.3f m)\\n' ...
], datestr(now,'yyyy-mm-dd HH:MM:SS'), N, opt.N_target, opt.d, opt.d/lambda, ...
   opt.scanAz_deg, opt.scanEl_deg, rangeRes, TCPI, velRes, HPBW_deg, D_eq, ...
   theta_scan, HPBW_scan_deg, scanLoss_dB, Runamb, vunamb, string(okGL), opt.d, d_limit);
fprintf('%s\n', summary_text);
fid = fopen(fullfile(outdir,'summary.txt'), 'w'); if fid>0, fprintf(fid, '%s', summary_text); fclose(fid); end
fig = figure('Name','Run Summary','Color','w','Position',[100 100 900 700],'Visible','on');
annotation('textbox',[0.05 0.05 0.9 0.9], 'String', summary_text, ...
    'Interpreter','none','FontName','FixedWidth','FontSize',10,'EdgeColor','none');
savefig2(outdir, 'summary'); close(fig);

%% ================== RANGE–DOPPLER DEMO ==================
if opt.doRangeDoppler
    R0 = 120; v0 = 12; SNRdB_sim = 10;
    S = opt.B/opt.Tchirp;  Ns = round(opt.Tchirp*opt.fs);  t = (0:Ns-1)/opt.fs;
    f_beat = (2*S/opt.c0) * R0;   f_D = 2*v0/lambda;
    rx = zeros(Ns, opt.N_chirps);
    for m = 1:opt.N_chirps
        sig = exp(1j*2*pi*( (f_beat + f_D)*t ));
        sig = sig + 10^(-SNRdB_sim/20) * (randn(size(sig))+1j*randn(size(sig)))/sqrt(2);
        rx(:,m) = sig;
    end
    NrFFT = 4096;  NdFFT = 512;
    RD = fft(rx, NrFFT, 1); RD = RD(1:NrFFT/2, :);
    RD = fftshift(fft(RD, NdFFT, 2), 2);
    f_fast = (0:NrFFT/2-1) * (opt.fs/NrFFT);
    rng_axis = (opt.c0/(2*S)) * f_fast;
    f_dop = linspace(-opt.PRF/2, opt.PRF/2, NdFFT);
    dop_axis = (lambda/2) * f_dop;
    figure('Name','Range–Doppler Map','Color','w','Visible','on');
    imagesc(dop_axis, rng_axis, 20*log10(abs(RD)+1e-12)); axis xy;
    xlabel('Velocity [m/s]'); ylabel('Range [m]'); title(sprintf('Range–Doppler (R0=%.0f m, v0=%.1f m/s)', R0, v0));
    colormap(parula); colorbar; caxis([-80 0]);
    savefig2(outdir, 'range_doppler');
end

%% ================== APERTURE FFT BEAMPATTERN ==================
if opt.doFFTBeampattern
    xr = max(abs(x)); yr = max(abs(y));
    Rmax = opt.apOversamp * max(xr, yr);
    gx = linspace(-Rmax, +Rmax, opt.apGridPx);
    gy = linspace(-Rmax, +Rmax, opt.apGridPx);
    dxy = gx(2)-gx(1);
    A = zeros(opt.apGridPx, opt.apGridPx);
    ix = round( (x - gx(1))/dxy ) + 1;
    iy = round( (y - gy(1))/dxy ) + 1;
    valid = ix>=1 & ix<=opt.apGridPx & iy>=1 & iy<=opt.apGridPx;
    linIdx = sub2ind(size(A), iy(valid), ix(valid));
    for kidx = 1:numel(linIdx)
        A(linIdx(kidx)) = A(linIdx(kidx)) + w(valid(kidx));
    end
    pad = 2; Nx = pad*opt.apGridPx; Ny = pad*opt.apGridPx;
    F = fftshift(fft2(A, Ny, Nx)); P = abs(F).^2;
    PdB = 10*log10( P / max(P(:)) + 1e-15 );
    fx = (-Nx/2:Nx/2-1) / (Nx*dxy);
    fy = (-Ny/2:Ny/2-1) / (Ny*dxy);
    uF = lambda * fx; vF = lambda * fy; [UF, VF] = meshgrid(uF, vF);
    PdB( (UF.^2 + VF.^2) > 1 ) = -Inf;
    figure('Name','Aperture FFT Beampattern','Color','w','Visible','on');
    imagesc(uF, vF, PdB); axis image; set(gca,'YDir','normal');
    colormap(parula); colorbar; caxis([-60 0]); xlabel('u'); ylabel('v');
    title('Beampattern via 2D FFT of Aperture'); savefig2(outdir, 'fft_uv');
    el0 = opt.cutElForHPBW_deg; azList = linspace(-90,90,2049);
    uCut2 = sind(el0).*cosd(azList); vCut2 = sind(el0).*sind(azList);
    patCut_dB = interp2(UF, VF, PdB, uCut2, vCut2, 'linear', -Inf);
    patCut_dB = movmean(patCut_dB, 5);
    [~, idx0] = min(abs(azList - opt.scanAz_deg));
    [BW_deg, azL, azR] = hpbw_from_cut(azList, patCut_dB, idx0);
    figure('Name','FFT Azimuth Cut & HPBW','Color','w','Visible','on');
    plot(azList, patCut_dB, 'LineWidth', 1.5); grid on; ylim([-60 0]);
    xlabel(sprintf('Azimuth (deg) @ EL=%g^\\circ', el0)); ylabel('Normalized (dB)');
    title(sprintf('Az Cut (FFT beampattern). HPBW \\approx %.2f^\\circ', BW_deg));
    hold on; yline(-3,'k--'); if ~isnan(BW_deg), xline(azL,'r--'); xline(azR,'r--'); end; hold off;
    savefig2(outdir, 'fft_az_cut');
    fid = fopen(fullfile(outdir,'summary.txt'), 'a'); if fid>0
        fprintf(fid, 'Exact HPBW from FFT cut @ EL=%g deg: %.2f deg (−3 dB at %.2f° / %.2f°)\n', el0, BW_deg, azL, azR); fclose(fid);
    end
end

% List files
d = dir(outdir); names = {d(~[d.isdir]).name};
fprintf('[INFO] Saved files in results/: %s\n', strjoin(names, ', '));
end

%% ================== PRESETS / UTILS ==================
function opt = apply_preset(opt, name)
name = lower(string(name));
switch name
    case {"default","boresight_tight"}
        % Keep defaults; mild taper for sidelobe control.
    case "wide_beam_low_elements"
        opt.N_target = 721;             % fewer elements
        opt.taper = struct('type','none','power',1);
    case "high_res_range"
        opt.B = 1.5e9;                  % larger bandwidth -> finer range resolution
        opt.Tchirp = 60e-6;
    case "long_range_high_prf"
        opt.PRF = 10e3;                 % higher PRF
        opt.Npulses = 256;              % longer CPI for finer velRes
    case "fast_scan_az45"
        opt.scanAz_deg = 45; opt.scanEl_deg = 0;
    otherwise
        warning('Unknown preset "%s" — using defaults.', name);
end
end

function out = merge_opts(base, delta)
out = base;
fn = fieldnames(delta);
for i=1:numel(fn)
    f = fn{i};
    if isstruct(delta.(f)) && isfield(base,f)
        out.(f) = merge_opts(base.(f), delta.(f));
    else
        out.(f) = delta.(f);
    end
end
end

function [q, r] = hex_axial_coords(K)
q = 0; r = 0;
for k = 1:K
    cq =  k; cr =  0;
    dirs = [ -1  1; -1  0; 0 -1; 1 -1; 1  0; 0  1 ];
    for d = 1:6
        dq = dirs(d,1); dr = dirs(d,2);
        for step = 1:k
            q = [q; cq]; %#ok<AGROW>
            r = [r; cr]; %#ok<AGROW>
            cq = cq + dq; cr = cr + dr;
        end
    end
end
end

function [x, y] = axial_to_xy(q, r, d)
x = d*( q + 0.5*r );
y = d*( (sqrt(3)/2)*r );
end

function [Sabs, F, Tm] = simple_stft(x, win, hop, nfft, fs)
x = x(:).'; L = numel(x); nFrames = 1 + floor((L - win)/hop);
Sabs = zeros(nfft/2+1, nFrames); F = (0:nfft/2)*(fs/nfft);
Tm = ((0:nFrames-1)*hop + win/2)/fs; idx = 1;
for m = 1:nFrames
    seg = zeros(1,win); i2 = min(idx+win-1, L);
    seg(1:(i2-idx+1)) = x(idx:i2);
    X = fft(seg, nfft);
    Sabs(:,m) = abs(X(1:nfft/2+1)).';
    idx = idx + hop;
end
end

function [BW, xL, xR] = hpbw_from_cut(x, pat_dB, idxCenter)
BW = NaN; xL = NaN; xR = NaN;
if isempty(idxCenter) || idxCenter<2 || idxCenter>numel(x)-1, return; end
p0 = pat_dB(idxCenter); pat = pat_dB - p0;
i = idxCenter; while i>1 && pat(i) > -3, i = i-1; end
if i==1, return; end
x1 = interp1(pat(i:i+1), x(i:i+1), -3);
j = idxCenter; while j<numel(x) && pat(j) > -3, j = j+1; end
if j==numel(x), return; end
x2 = interp1(pat(j-1:j), x(j-1:j), -3);
BW = x2 - x1; xL = x1; xR = x2;
end

function savefig2(outdir, basename)
f = gcf; set(f,'Color','w'); drawnow;
png = fullfile(outdir, [basename '.png']);
pdf = fullfile(outdir, [basename '.pdf']);
fig = fullfile(outdir, [basename '.fig']);
ok = false;
try
    if exist('exportgraphics','file') == 2
        exportgraphics(f, png, 'Resolution', 200);
        exportgraphics(f, pdf);
        ok = true;
    end
end
if ~ok
    try, set(f,'PaperPositionMode','auto'); print(f, png, '-dpng', '-r200'); ok=true; catch, end
    try, print(f, pdf, '-dpdf'); catch, end
end
try, saveas(f, fig); catch, end
if exist(png,'file')==2, fprintf('[SAVE] %s\n', png); else, warning('[WARN] PNG not saved: %s', png); end
if exist(pdf,'file')==2, fprintf('[SAVE] %s\n', pdf); end
if exist(fig,'file')==2, fprintf('[SAVE] %s\n', fig); end
end
