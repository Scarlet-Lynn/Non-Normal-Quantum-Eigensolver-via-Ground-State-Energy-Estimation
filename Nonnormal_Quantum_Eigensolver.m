function [mu, sample] = eigensolver(A, n, m, kappa, epsilon)
    % eigenvalue estimation for non-normal matrices

    I = eye(n);
    center = 0; % initial sampling center
    R_l = 1; % initial sampling radius
    sample = 0; % count the number of samples
    for l = 1:ceil(log2(1/epsilon))
        delta_l = 1/(kappa*(3*2^l)^m); % grid spacing
        s = ceil(R_l/delta_l);
        break_inner = false;
        for i = -s:s
            for j = -s:s
                sample = sample+1;
                mu = center + (i*delta_l+j*delta_l*1i);
                sv = svd(A-mu*I);
                sigma_0 = min(sv);
                if sigma_0 < delta_l
                    center = mu;
                    R_l = 3*(kappa*sigma_0)^(1/m);
                    break_inner = true;
                end
                if break_inner
                    break;
                end
            end
            if break_inner
                break;
            end
        end
    end
end

function mu = ext_eig(A, n, m, kappa, epsilon)
    % extreme eigenvalue estimation

    I = eye(n);
    sv = svd(A);
    sigma_0 = min(sv);
    R_1 = sigma_0;
    R_2 = min([1,3*(kappa*sigma_0)^(1/m)]);
    while R_2-R_1 > epsilon
        delta = 1/kappa * ((R_2-R_1)/(3*2))^m;
        M = ceil(pi/asin((sqrt(1-0.9^2)*delta)/R_1)); % set c = 0.9
        case_1 = false;
        case_2 = true;
        while case_2 
            for j = 0:M-1
                mu = R_1*exp(1i*j*2*pi/M);
                sv = svd(A-mu*I);
                sigma_0 = min(sv);
                if sigma_0 <= delta
                    case_1 = true;
                    case_2 = false;
                    break;
                end
            end
            if case_1
                R_2 = R_1 + 3*(kappa*sigma_0)^(1/m);
            else
                R_1 = R_1+0.9*delta;
            end
        end
    end
end

% ═════════════════════════════════════════════════════════════════════
% Asymptotic Performance
% ═════════════════════════════════════════════════════════════════════

% fix m=2 and epsilon = 0.01, vary kappa
n = 4;
m = 2;
kappa = [1, 5, 10];
epsilon = 0.01;
J = [1+1i 1 0 0
     0 1+1i 0 0
     0 0 -2 1
     0 0 0 -2]; % the Jordan matrix
eig_est = zeros(1, 3);
sample_count = zeros(1, 3);
for t = 1:3
    P = gallery('randsvd', n, kappa(t)); generate an invertible matrix so that A = P * J / P has Jordan condition number=kappa(i)
    A = P * J / P;
    [eig_est(t), sample_count(t)] = eigensolver(A, n, m, kappa(t), epsilon);
end

% result
% sample_count = [1.6526e+07, 4.0000e+08, 1.5274e+09];
% eig_est = [0.9974 + 1.0000i, 0.9980 + 1.0000i, 0.9988 + 1.0000i];

figure;
loglog(kappa, sample_count, 'o', 'MarkerFaceColor', 'k');
hold on;
x_mid = 10^(mean(log10(kappa)));
y_mid = 10^(mean(log10(sample_count)));
C = y_mid / (x_mid^2);
x_fit = logspace(log10(min(kappa)), log10(max(kappa)), 100);
y_fit = C * x_fit.^2;
loglog(x_fit, y_fit, 'b', 'LineWidth', 1);
xlim([min(x_fit)/1.1, max(x_fit)*1.1])
ylim([min(y_fit)/2, max(y_fit)*2])
xlabel('kappa (log scale)');
ylabel('sample count (log scale)');
hold off;

% fix m=2 and kappa=1, vary epsilon
n = 4;
m = 2;
kappa = 1;
epsilon = [0.1, 0.01, 0.001];
J = [1+1i 1 0 0
     0 1+1i 0 0
     0 0 -2 1
     0 0 0 -2];
P = gallery('randsvd', n, 1);
A = P * J / P;
eig_est = zeros(1, 3);
sample_count = zeros(1, 3);
for t = 1:3
    [eig_est(t), sample_count(t)] = eigensolver(A, n, m, kappa, epsilon(t));
end

% result
% eig_est = [0.9792 + 0.9996i, 0.9974 + 1.0000i, 0.9997 + 1.0000i];
% sample_count = [2.6001e+05, 1.6526e+07, 1.0570e+09];

figure;
loglog(epsilon, sample_count, 'o', 'MarkerFaceColor', 'k');
hold on;
x_mid = 10^(mean(log10(epsilon)));
y_mid = 10^(mean(log10(sample_count)));
C = y_mid / (x_mid^(-2));
x_fit = logspace(log10(min(epsilon)), log10(max(epsilon)), 100);
y_fit = C * x_fit.^(-2);
loglog(x_fit, y_fit, 'b', 'LineWidth', 1);
xlim([min(x_fit)/1.1, max(x_fit)*1.1]) 
xlabel('epsilon (log scale)');
ylabel('sample count (log scale)');
hold off;

% ═════════════════════════════════════════════════════════════════════
% Extreme Eigenvalue Estimation and Extreme Eigenvector Preparation
% ═════════════════════════════════════════════════════════════════════

% generating matrix instances
% instance 1
D = [1+1i 0 0 0
     0 -2-2i 0 0 
     0 0 8 0
     0 0 0 8];
n = 4;
m = 1;
% % instance 2
% D = [1 1 0 0
%      0 1 0 0
%      0 0 2 1
%      0 0 0 2];
% n = 4;
% m = 2;
% % instance 3
% D = [1+1i 1 0 0 0 0 0 0
%      0 1+1i 0 0 0 0 0 0
%      0 0 -2 1 0 0 0 0
%      0 0 0 -2 0 0 0 0
%      0 0 0 0 2i 1 0 0
%      0 0 0 0 0 2i 0 0
%      0 0 0 0 0 0 3 0
%      0 0 0 0 0 0 0 4];
% n = 8;
% m = 2;
D = D/norm(D);
disp(diag(D));
H = randn(n)+1i*randn(n);
H = H+H';
U = expm(1i*H);
kappa = 1;
epsilon = 0.01;
A = U*D*U^(-1);
mu = ext_eig(A, n, m, kappa, epsilon);
err = abs(mu-D(1,1));

% eigenstate preparation
H = sqrtm((A-mu*eye(n))'*(A-mu*eye(n)));
[V, eig_val] = eig(H);
eigenvalues = diag(eig_val);
[~, idx] = min(eigenvalues);
ground_state = V(:, idx);
eig_state_true = U(:,1);   % the true eigenstate corresponding to the extreme eigenvalue
fidelity = abs(eig_state_true' * ground_state)^2;
