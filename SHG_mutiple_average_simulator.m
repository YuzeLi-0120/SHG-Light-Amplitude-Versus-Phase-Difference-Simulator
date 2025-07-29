% Author: Yuze Li
% Date: July 2025
% Description: Simulation of average SHG amplitude.


%%% Global Parameters
c = 3e8;
fai_min = 0;
fai_max = 20;
N = 1000;
fai = linspace(fai_min * pi, fai_max * pi, N);

lambda = 296e-9;
lambda_2w = lambda / 2;

n_omega = 1.589;
n2w = 1.826; 

deff0 = 0.3e-12;

k_omega = 2*pi*n_omega / lambda;
k_2omega = 2*pi*n2w / lambda_2w;
deltak = k_2omega - 2 * k_omega;

mode_map = {'mismatch', 'QPM', 'APP'};

% local_params
% 1: errtype (1=independent, 2=dependent)
% 2: faberr
% 3: deltak
% 4: mode (1=mismatch, 2=QPM, 3=APP)
local_params = [
    2, 0.0, deltak, 3;
    2, 0.05, deltak, 3;
    2, 0.10, deltak, 3;
    2, 0.15, deltak, 3;
    2, 0.20, deltak, 3;
    2, 0.25, deltak, 3;
    2, 0.30, deltak, 3;
    2, 0.35, deltak, 3;
    2, 0.40, deltak, 3;
];

referencecurve = false;  % Decide whether to draw the reference line (perfectly phase matched)
colors = lines(size(local_params,1) + referencecurve);

if referencecurve
    deltak_ref = deltak;
    omega_ref = 2 * pi * c / lambda;
    Ew_ref = 1;
    z_ref = fai / deltak_ref;
    B_ref = (2 * omega_ref * deff0 * Ew_ref^2) / (c * n2w);
    E2w_ref = B_ref * z_ref; 
    
    plot(fai/pi, E2w_ref, '--k', 'LineWidth', 2, 'DisplayName', 'Reference (deltak=0 slope)');
end

num_trials = 100;
target_index = N;

mean_values = zeros(size(local_params,1), 1);
std_values = zeros(size(local_params,1), 1);
labels = strings(size(local_params,1), 1);

for i = 1:size(local_params,1)
    errtype_num = local_params(i,1);
    faberr = local_params(i,2);
    deltak = local_params(i,3);
    mode_num = local_params(i,4);

    errtype_map = {'independent', 'dependent'};
    errtype = errtype_map{errtype_num};

    mode = mode_map{mode_num};
    labels(i) = sprintf('%s, %.2f, %s', errtype, faberr, mode);

    omega = 2 * pi * c / lambda;
    Ew = 1;

    if deltak == 0
        z = zeros(size(fai));
    else
        z = fai / deltak;
    end

    E2w_end_values = zeros(num_trials, 1);

    for trial = 1:num_trials
        E2w_amp = simulate_SHG_local(z, omega, Ew, n2w, deff0, deltak, mode, faberr, errtype, c);
        E2w_end_values(trial) = E2w_amp(target_index);
    end

    mean_values(i) = mean(E2w_end_values);
    std_values(i) = std(E2w_end_values);
end

% Draw the bars.
figure;
b = bar(mean_values, 'FaceColor', 'flat');
hold on;

for k = 1:length(b.CData)
    b.CData(k,:) = colors(referencecurve + k, :);
end

% Error bars.
errorbar(1:length(mean_values), mean_values, std_values, ...
         'k.', 'LineWidth', 1.5, 'CapSize', 10);

xticks(1:length(labels));
xticklabels(labels);
xtickangle(30); 
ylabel('E_{2\omega} at \phi=20\pi');
title('SHG Output at \phi=20\pi with Varying Fab Errors');
grid on;

function E2w_amp = simulate_SHG_local(z, omega, Ew, n2w, deff0, deltak, mode, faberr, errtype, c)
    if deltak == 0
        B = (2 * omega * deff0 * Ew.^2) / (c * n2w);
        E2w_amp = B * z;
        return
    end

    switch mode
        case 'mismatch'
            deff = deff0 * ones(size(z));
            phase_factor = exp(1i * deltak * z);
        case 'QPM'
            Lambda = 2 * pi / deltak;
            N_flip = floor(max(z) / (Lambda / 2));
            deff = zeros(size(z));
            switch errtype
                case 'independent'
                    flip_pos = ((0:N_flip-1) * (Lambda / 2)) + (2*rand(1, N_flip)-1) * faberr * (Lambda/2);
                case 'dependent'
                    errors = (2*rand(1, N_flip) - 1) * faberr * (Lambda/2);
                    flip_pos = cumsum([0, repmat(Lambda/2, 1, N_flip-1)]) + cumsum(errors);
                otherwise
                    error('Invalid errtype.');
            end
            for i = 1:N_flip
                z_start = flip_pos(i);
                z_end = z_start + (Lambda / 2);
                polarity = (-1)^(i+1);
                idx = (z >= z_start) & (z < z_end);
                deff(idx) = polarity * deff0;
            end
            phase_factor = exp(1i * deltak * z);
        case 'APP'
            Lambda_app = 2 * pi / deltak;
            N_app = floor(max(z) / Lambda_app);
            deff = zeros(size(z));
            switch errtype
                case 'independent'
                    base_pos = ((0:N_app-1) * Lambda_app);
                    start_order = base_pos + (2*rand(1, N_app)-1) * faberr * Lambda_app;
                case 'dependent'
                    errors = (2*rand(1, N_app) - 1) * faberr * Lambda_app;
                    start_order = cumsum([0, repmat(Lambda_app, 1, N_app-1)]) + cumsum(errors);
                otherwise
                    error('Invalid errtype.');
            end
            for i = 1:N_app
                z_start = start_order(i);
                z_end = z_start + Lambda_app / 2;
                idx = (z >= z_start) & (z < z_end);
                deff(idx) = deff0;
            end
            phase_factor = exp(1i * deltak * z);
        otherwise
            error('Unsupported mode.');
    end

    B = (2 * omega .* deff .* Ew.^2) ./ (c * n2w);
    E2w = cumtrapz(z, B .* phase_factor);
    E2w_amp = abs(E2w);
end
