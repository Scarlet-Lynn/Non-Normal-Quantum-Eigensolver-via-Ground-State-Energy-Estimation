function [mu, sample] = eigensolver_upper(A, n, m, kappa, epsilon)
    % search eigenvalues within the upper half plane

    I = eye(n);
    center = 0; % initial sampling center
    R_l = 1; % initial sampling radius
    sample = 0; % count the number of samples
    for l = 1:ceil(log2(1/epsilon))
        if l == 1
            delta_l = 1/(kappa*(3*2^l)^m); % grid spacing
            s = ceil(R_l/delta_l);
            break_inner = false;
            for i = -s:s
                for j = 1:s   % search eigenvalues within the upper half plane
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
        else
            delta_l = 1/(kappa*(3*2^l)^m); % grid spacing
            s = ceil(R_l/delta_l);
            break_inner = false;
            for i = -s:s
                for j = -s:s
                    sample = sample+1;
                    mu = center + (i*delta_l+j*delta_l*1i);
                    if imag(mu_candidate) <= 0
                        continue;
                    end
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
end

function [mu, sample] = eigensolver_lower(A, n, m, kappa, epsilon)
    % search eigenvalues within the lower half plane

    I = eye(n);
    center = 0; % initial sampling center
    R_l = 1; % initial sampling radius
    sample = 0; % count the number of samples
    for l = 1:ceil(log2(1/epsilon))
        if l == 1
            delta_l = 1/(kappa*(3*2^l)^m); % grid spacing
            s = ceil(R_l/delta_l);
            break_inner = false;
            for i = -s:s
                for j = -s:-1   % search eigenvalues within the lower half plane
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
        else
            delta_l = 1/(kappa*(3*2^l)^m); % grid spacing
            s = ceil(R_l/delta_l);
            break_inner = false;
            for i = -s:s
                for j = -s:s
                    sample = sample+1;
                    mu = center + (i*delta_l+j*delta_l*1i);
                    if imag(mu_candidate) >= 0
                        continue;
                    end                    
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
end

function [mu, sample] = eigensolver_real_inner(A, n, m, kappa, epsilon)
    % search eigenvalues on the real line
    % restrict the search region to the interval [-1/2,1/2]

    I = eye(n);
    center = 0; % initial sampling center
    R_l = 1; % initial sampling radius
    sample = 0; % count the number of samples
    for l = 1:ceil(log2(1/epsilon))
        if l == 1
            delta_l = 1/(kappa*(3*2^l)^m); % grid spacing
            s = ceil(R_l/delta_l);
            break_inner = false;
            for i = -s/2:s/2   % restrict the search region to the interval [-1/2,1/2]
                for j = -1:1   % search eigenvalues on the real line
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
        else
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
end

function [mu, sample] = eigensolver_real_outer(A, n, m, kappa, epsilon)
    % search eigenvalues on the real line
    % restrict the search region to intervals [-1,-1/2] and [1/2,1]

    I = eye(n);
    center = 0; % initial sampling center
    R_l = 1; % initial sampling radius
    sample = 0; % count the number of samples
    for l = 1:ceil(log2(1/epsilon))
        if l == 1
            delta_l = 1/(kappa*(3*2^l)^m); % grid spacing
            s = ceil(R_l/delta_l);
            break_inner = false;
            for i = [-s,-s/2, s/2,s]   % restrict the search region to intervals [-1,-1/2] and [1/2,1]
                for j = -1:1   % search eigenvalues on the real line
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
        else
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
end

% Liouvillian operator
gamma = 0.4;
omega = 0.25;
L = [-gamma 0 0 gamma
     0 -gamma-2*1j*omega 0 0
     0 0 -gamma+2*1j*omega 0
     gamma 0 0 -gamma];
[P,J] = jordan(L);
kappa = norm(P,2)*norm(inv(P),2);
eig_1 = eigensolver_upper(L, 4, 1, kappa, 0.001);
eig_2 = eigensolver_lower(L, 4, 1, kappa, 0.001);
eig_3 = eigensolver_real_inner(L, 4, 1, kappa, 0.001);
eig_4 = eigensolver_real_outer(L, 4, 1, kappa, 0.001);
