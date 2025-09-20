function selfcheck()
% Self-check for CI (Octave-friendly, no graphics).
% Verifies basic math and array-factor routines.

% RF / geometry (reduced sizes to stay fast)
fc = 5e9; c0 = 3e8; lambda = c0/fc;
d  = 0.5*lambda;
N_target = 127;   % small array for speed
K = ceil( ( -3 + sqrt(9 + 12*(N_target-1)) )/6 );
[q, r] = hex_axial_coords(K);
[x, y] = axial_to_xy(q, r, d);
if numel(x) > N_target
  R = sqrt(x.^2 + y.^2);
  [~, idx] = sort(R, 'ascend');
  keep = idx(1:N_target);
  x = x(keep); y = y(keep);
end
N = numel(x);

% Steering 0 deg
scanAz_deg = 0; scanEl_deg = 0;
u0 = sind(scanEl_deg)*cosd(scanAz_deg);
v0 = sind(scanEl_deg)*sind(scanAz_deg);
k  = 2*pi/lambda;
w  = ones(N,1) .* exp(1j*k*(x*u0 + y*v0));

% Evaluate AF at boresight and off-boresight
u_test = [0, 0.2]; v_test = [0, 0];
AF = zeros(1,2);
for n=1:N
  AF = AF + w(n) .* exp(1j*k*(x(n).*u_test + y(n).*v_test));
end
AF = abs(AF);
if ~(AF(1) > AF(2))
  error('AF peak check failed: AF(boresight)=%.3g, AF(off)=%.3g', AF(1), AF(2));
end

% Range/velocity resolution sanity
B = 1e9; PRF = 5e3; Npulses = 64;
TCPI = Npulses/PRF;
rangeRes = c0/(2*B);
velRes   = lambda/(2*TCPI);
if ~(rangeRes > 0 && rangeRes < 1 && velRes > 0 && velRes < 10)
  error('Resolution sanity failed: rangeRes=%.3g, velRes=%.3g', rangeRes, velRes);
end

disp('selfcheck: OK');
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
x = d*( q + 0.5*r ); y = d*( (sqrt(3)/2)*r );
end
