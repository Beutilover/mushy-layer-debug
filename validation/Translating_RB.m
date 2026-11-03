%% Forced-throughflow Rayleigh-Benard benchmark for SOFTBALL
% Compare the dominant Fourier-Galerkin eigenmode with a SOFTBALL snapshot.
% The benchmark uses relative velocity and a uniform imposed throughflow.
% Requires MATLAB with local functions (R2016b or later) and SOFTBALL HDF5 output.
% NOTE: HDF5 layout and variable indices are specific to the supplied case.
clearvars;

%% Case parameters (match the SOFTBALL input file)
cfg.thetaLower = 2.1;
cfg.thetaUpper = 1.1;
cfg.throughflow = 10;
cfg.Pr = 1;
cfg.Ra = 21000;
cfg.horizontalWavenumber = 4.86;
cfg.nModes = 512;
cfg.nz = 256;
cfg.nx = 664;
cfg.nVariables = 15;
cfg.xIndex = 249;
cfg.outputDir = 'MFRBV10test';
cfg.temperatureIndex = 3;
cfg.verticalVelocityIndex = 11;

%% Growth-rate comparison (temperature and vertical velocity)
cfg.checkGrowthRate = true;
% Specify [tStart, tEnd] to fit a known linear-growth interval.
% Selected linear-growth interval, shared by both fields.
% Both fits use identical snapshots; no automatic window optimisation.
% Set [NaN NaN] to inspect all snapshots and use the middle 20%-70%
% of available times as a provisional window; review the diagnostic plot.
cfg.growthFitTimeRange = [0.15006152 0.55017871];
cfg.minGrowthSnapshots = 4;

%% Extract the final SOFTBALL snapshot
files = dir(fullfile(cfg.outputDir, 'plt*.2d.hdf5'));
assert(~isempty(files), 'No SOFTBALL plot files found in %s', cfg.outputDir);
indices = nan(numel(files), 1);
for j = 1:numel(files)
    token = regexp(files(j).name, '^plt(\d+)\.2d\.hdf5$', 'tokens', 'once');
    if ~isempty(token), indices(j) = str2double(token{1}); end
end
[~, latest] = max(indices);
snapshotFile = fullfile(cfg.outputDir, files(latest).name);

z = ((1:cfg.nz)' - 0.5) / cfg.nz;
V = cfg.throughflow;
% Base temperature profile for imposed vertical throughflow.
if abs(V) < 1e-12
    thetaBase = cfg.thetaLower + (cfg.thetaUpper - cfg.thetaLower) * z;
else
    thetaBase = (exp(V)*cfg.thetaLower - cfg.thetaUpper ...
        - (cfg.thetaLower - cfg.thetaUpper)*exp(V*z)) / (exp(V)-1);
end
T = readSnapshotVariable(snapshotFile, cfg.temperatureIndex, cfg);
w = readSnapshotVariable(snapshotFile, cfg.verticalVelocityIndex, cfg);
assert(cfg.xIndex >= 1 && cfg.xIndex <= cfg.nx, 'xIndex is out of bounds');
Tpert = T(:, cfg.xIndex) - thetaBase;
wpert = w(:, cfg.xIndex);

%% Independent Fourier-Galerkin linear stability calculation
[Tmode, wmode, growthRate] = getgrowth(cfg.nModes, ...
    cfg.horizontalWavenumber^2, V, cfg.Pr, cfg.Ra, ...
    cfg.thetaLower, cfg.thetaUpper);
Ttheory = Tmode(z);
wtheory = wmode(z);

% Determine ONE signed amplitude scale from the temperature profile.
% Apply the same scale to velocity to retain the relative mode amplitudes.
[~, iData] = max(abs(Tpert));
[~, iTheory] = max(abs(Ttheory));
assert(abs(Ttheory(iTheory)) > 0, 'Zero theoretical temperature eigenmode');
scale = Tpert(iData) / Ttheory(iTheory);
Ttheory = scale * Ttheory;
wtheory = scale * wtheory;

fprintf('Dominant LSA growth rate: %.9g\n', growthRate);
fprintf('SOFTBALL snapshot: %s\n', snapshotFile);

%% Plot both fields with consistent styling
fig = figure('Color','w', 'Position',[100 100 1000 430]);
tiledlayout(fig, 1, 2, 'TileSpacing','compact', 'Padding','compact');
nexttile;
plot(Ttheory, z, '-', 'LineWidth',1.8); hold on;
plot(Tpert, z, '--', 'LineWidth',1.6);
xlabel('Temperature perturbation'); ylabel('z');
title('(a) Temperature eigenmode'); grid on; box on;
legend('Fourier-Galerkin LSA','SOFTBALL','Location','best');
nexttile;
plot(wtheory, z, '-', 'LineWidth',1.8); hold on;
plot(wpert, z, '--', 'LineWidth',1.6);
xlabel('Relative vertical velocity'); ylabel('z');
title('(b) Vertical-velocity eigenmode'); grid on; box on;
legend('Fourier-Galerkin LSA','SOFTBALL','Location','best');
exportgraphics(fig, 'forced_throughflow_eigenmodes.png', 'Resolution', 200);

%% Compare temperature and vertical-velocity growth rates on the SAME time window
% Both amplitudes are projections onto the independent LSA eigenmodes.
% The fit interval is shared, so any discrepancy is not caused by selecting
% different linear-growth windows for the two variables.
if cfg.checkGrowthRate
    [times, ampsT, ampsW] = getBothModeAmplitudeSeries( ...
        files, indices, cfg, thetaBase, Tmode(z), wmode(z));
    good = isfinite(times) & isfinite(ampsT) & isfinite(ampsW) & ...
        abs(ampsT)>0 & abs(ampsW)>0;
    times = times(good); ampsT = ampsT(good); ampsW = ampsW(good);
    [times, order] = sort(times);
    ampsT = ampsT(order); ampsW = ampsW(order);

    if numel(times) < cfg.minGrowthSnapshots
        warning('Too few usable snapshots for a growth-rate fit.');
    else
        if all(isfinite(cfg.growthFitTimeRange))
            % Allow for rounded timestamps entered in the configuration.
            timeTol = 1e-7 * max(1, max(abs(times)));
            fitMask = times >= cfg.growthFitTimeRange(1)-timeTol & ...
                      times <= cfg.growthFitTimeRange(2)+timeTol;
        else
            t0 = times(1); t1 = times(end);
            fitMask = times >= t0+0.20*(t1-t0) & ...
                      times <= t0+0.70*(t1-t0);
            fprintf('Using provisional middle 20%%-70%% time window.\n');
        end
        if nnz(fitMask) < cfg.minGrowthSnapshots
            warning('Only %d snapshots in selected fit interval.', nnz(fitMask));
        else
            [sigmaT, r2T, fitT] = fitLogAmplitude(times, ampsT, fitMask);
            [sigmaW, r2W, fitW] = fitLogAmplitude(times, ampsW, fitMask);
            fprintf('LSA growth rate:                %.8f\n', growthRate);
            fprintf('SOFTBALL temperature growth:    %.8f (R^2 = %.6f)\n', sigmaT, r2T);
            fprintf('SOFTBALL vertical-velocity growth: %.8f (R^2 = %.6f)\n', sigmaW, r2W);
            fprintf('Temperature vs LSA difference: %.4f%%\n', ...
                100*abs(sigmaT-growthRate)/abs(growthRate));
            fprintf('Velocity vs LSA difference:    %.4f%%\n', ...
                100*abs(sigmaW-growthRate)/abs(growthRate));
            fprintf('Temperature vs velocity slope difference: %.8f\n', sigmaT-sigmaW);
            fprintf('Temperature vs LSA absolute difference: %.8f\n', abs(sigmaT-growthRate));
            fprintf('Velocity vs LSA absolute difference: %.8f\n', abs(sigmaW-growthRate));
            fprintf('Temperature vs velocity relative difference: %.4f%%\n', ...
                200*abs(sigmaT-sigmaW)/(abs(sigmaT)+abs(sigmaW)));
            fprintf('Fit interval: [%.8g, %.8g] (%d snapshots)\n', ...
                min(times(fitMask)), max(times(fitMask)), nnz(fitMask));

            figGrowth = figure('Color','w','Position',[150 150 1050 420]);
            tiledlayout(figGrowth,1,2,'TileSpacing','compact','Padding','compact');
            nexttile;
            plotGrowthPanel(times, ampsT, fitMask, fitT, sigmaT, r2T, ...
                'Temperature perturbation');
            nexttile;
            plotGrowthPanel(times, ampsW, fitMask, fitW, sigmaW, r2W, ...
                'Relative vertical velocity');
            exportgraphics(figGrowth,'forced_throughflow_growth_diagnostic.png', ...
                'Resolution',200);
        end
    end
end

function [times, amplitudesT, amplitudesW] = getBothModeAmplitudeSeries( ...
    files, indices, cfg, thetaBase, theoryT, theoryW)
    n = numel(files);
    times = nan(n,1); amplitudesT = nan(n,1); amplitudesW = nan(n,1);
    normT = sum(theoryT(:).^2);
    normW = sum(theoryW(:).^2);
    assert(normT > 0 && normW > 0, 'An LSA eigenmode has zero norm');
    for j = 1:n
        if ~isfinite(indices(j)), continue; end
        filePath = fullfile(cfg.outputDir,files(j).name);
        try
            times(j) = double(h5readatt(filePath,'/','time'));
            Tfield = readSnapshotVariable(filePath,cfg.temperatureIndex,cfg);
            wfield = readSnapshotVariable(filePath,cfg.verticalVelocityIndex,cfg);
            Tpert = Tfield(:,cfg.xIndex) - thetaBase;
            wpert = wfield(:,cfg.xIndex);
            amplitudesT(j) = sum(Tpert(:).*theoryT(:))/normT;
            amplitudesW(j) = sum(wpert(:).*theoryW(:))/normW;
        catch ME
            warning('Skipping %s: %s',filePath,ME.message);
        end
    end
end

function [slope, rSquared, fittedLogAmp] = fitLogAmplitude(times, amplitudes, fitMask)
    logAmp = log(abs(amplitudes));
    coefficients = polyfit(times(fitMask),logAmp(fitMask),1);
    slope = coefficients(1);
    fittedLogAmp = polyval(coefficients,times(fitMask));
    observed = logAmp(fitMask);
    ssTot = sum((observed-mean(observed)).^2);
    if ssTot > 0
        rSquared = 1-sum((observed-fittedLogAmp).^2)/ssTot;
    else
        rSquared = NaN;
    end
end

function plotGrowthPanel(times, amplitudes, fitMask, fitted, slope, r2, fieldName)
    plot(times,log(abs(amplitudes)),'o-','LineWidth',1.1, ...
        'DisplayName','SOFTBALL projection'); hold on;
    plot(times(fitMask),log(abs(amplitudes(fitMask))),'o','MarkerSize',7, ...
        'DisplayName','Fit samples');
    plot(times(fitMask),fitted,'--','LineWidth',1.8, ...
        'DisplayName','Linear fit');
    xlabel('Simulation time'); ylabel('log |projected mode amplitude|');
    title(sprintf('%s: slope %.5f, R^2 %.5f',fieldName,slope,r2));
    legend('Location','best'); grid on; box on;
end


function snapshot = readSnapshotVariable(fileName, variableIndex, cfg)
    % Case-specific SOFTBALL/Chombo HDF5 storage layout: one level, one box.
    raw = h5read(fileName, '/level_0/data:datatype=0');
    expected = cfg.nx * cfg.nz * cfg.nVariables;
    assert(numel(raw) == expected, ...
        'Unexpected HDF5 size: expected %d entries, got %d', expected, numel(raw));
    assert(variableIndex >= 1 && variableIndex <= cfg.nVariables);
    array = reshape(raw, [cfg.nx, cfg.nz, cfg.nVariables]);
    snapshot = array(:,:,variableIndex).'; % [nz, nx]
end

function [Tmod, wmod, gr] = getgrowth(N,k2,V,Pr,Ra,theta_d,theta_u)
L1 = getL1(N,k2);
R11 = getR11(N,V,k2,Pr);
R12 = getR12(N,k2,Pr,Ra);
L2 = getL2(N);
R21 = getR21(N,V,theta_d,theta_u);
R22 = getR22(N,V,k2);

M11 = L1 \ R11;
M12 = L1 \ R12;
M21 = L2 \ R21;
M22 = L2 \ R22;

M = [M11, M12;
     M21, M22];
[V, D] = eig(M);
lambda = diag(D);
[gr, idx] = max(real(lambda));
% Eigenvector associated with the largest real growth rate
v = V(:, idx);
    
% Ensure column vector
v = v(:);
vupper = v(1:N);
v = v(N+1:2*N); % Temperature expansion coefficients

wmod = @(z) composew_mode(z, vupper, N);
Tmod = @(z) composeT_mode(z, v, N);
end

function y = composeT_mode(z, v, N)
    % Preserve original size, work on row-vector for sine evaluation
    orig_size = size(z);
    z_vec = z(:).';               % 1 x M row vector
    y_vec = zeros(1, numel(z_vec));

    % Sum using for-loop to avoid dimension mismatch
    for n = 1:N
        y_vec = y_vec + v(n) * sin(n * pi * z_vec);
    end

    % Return real part and reshape to original input shape
    y = reshape(real(y_vec), orig_size);
end

function y = composew_mode(z, v, N)
    % Preserve original size, work on row-vector for sine evaluation
    orig_size = size(z);
    z_vec = z(:).';               % 1 x M row vector
    y_vec = zeros(1, numel(z_vec));

    % Sum using for-loop to avoid dimension mismatch
    for n = 1:N
        y_vec = y_vec + v(n) *z_vec.*(1-z_vec).*sin(n * pi * z_vec);
    end

    % Return real part and reshape to original input shape
    y = reshape(real(y_vec), orig_size);
end


% Assemble Fourier-Galerkin matrices
function L1 = getL1(N,k2)
L1 = zeros(N,N);
for m=1:N
    for n=1:N
        L1(m,n) = - getdelta(m,n) - (k2+n^2*pi^2)*obtain_A(m,n) + 2*n*pi*obtain_C(m,n);
    end    
end
end

function R11 = getR11(N,V,k2,Pr)
R11 = zeros(N,N);
for m=1:N
    for n=1:N
        R11(m,n) = V*((k2+3*n^2*pi^2)*obtain_B(m,n) + n*pi*(k2+n^2*pi^2)*obtain_D(m,n) + 6*n*pi*obtain_E(m,n)) + ...
            + Pr*((2*k2+6*n^2*pi^2)*getdelta(m,n) + (k2^2+2*k2*n^2*pi^2+n^4*pi^4)*obtain_A(m,n) - 4*n*pi*(k2+n^2*pi^2)*obtain_C(m,n));
    end    
end
end

function R12 = getR12(N,k2,Pr,Ra)
R12 = zeros(N,N);
for m=1:N
    for n=1:N
        R12(m,n) = -1/2*k2*Pr*Ra*getdelta(m,n);
    end    
end
end

function L2 = getL2(N)
L2 = zeros(N,N);
for m=1:N
    for n=1:N
        L2(m,n) = 1/2*getdelta(m,n);
    end
end
end


function R21 = getR21(N,V,theta_d,theta_u)
R21 = zeros(N,N);
if V==0
    for m=1:N
        for n=1:N
            R21(m,n) = obtain_A(m,n)*(theta_d-theta_u);
        end
    end
else
    for m=1:N
        for n=1:N
            Phi = @(k) real( (exp(V + 1i*k*pi).*(V + 1i*k*pi - 2) + V + 1i*k*pi + 2) ...
                 ./ (V + 1i*k*pi).^3 );

            I = 0.5*(Phi(abs(n-m)) - Phi(n+m));
            R21(m,n) = I*V*(theta_d-theta_u)/(exp(V)-1);
        end
    end
end
end

function R22 = getR22(N,V,k2)
R22 = zeros(N,N);
for m=1:N
    for n=1:N
        R22(m,n) = -1/2*getdelta(m,n)*(k2+n^2*pi^2) - n*pi*V*obtain_E(m,n);
    end
end
end


function delta_mn = getdelta(m,n)
if m==n
    delta_mn = 1;
else
    delta_mn = 0;
end
end


% Analytical Galerkin integral coefficients
function A_mn = obtain_A(m,n)
if m==n
    A_mn = 1/12 + 1/(4*n^2*pi^2);
else
    A_mn = 2*m*n*((-1)^(m+n+1)-1)/pi^2/(n^2-m^2)^2;
end
end

function B_mn = obtain_B(m,n)
if m==n
    B_mn = 0;
else
    B_mn = 4*m*n*(1-(-1)^(m+n))/pi^2/((n^2-m^2)^2);
end
end

function C_mn = obtain_C(m,n)
if m==n
    C_mn = 1/(2*n*pi);
else
    C_mn = m*(1+(-1)^(m+n))/pi/(m^2-n^2);
end
end

function D_mn = obtain_D(m,n)
if m==n
    D_mn = 0;
else
    D_mn = (1-(-1)^(m+n))/pi^3*(1/(m+n)^3 + 1/(m-n)^3);
end
end

function E_mn = obtain_E(m,n)
if m==n
    E_mn = 0;
else
    E_mn = m*(1-(-1)^(m+n))/pi/(m^2-n^2);
end
end














