clear; clc; close all;

% loading variables form stress_log.m
stress_log;

% ---- choice of element and Gauss point ----
e0  = 1;
ig0 = 1;

mask = (e == e0) & (ig == ig0);

n_sel   = n(mask);
sigx_sel = sig_x(mask);
eps_sel = eps_x(mask);
sigy_sel = sig_y(mask);      % <-- must exist in stress_log.m
txy_sel  = tau_xy(mask);     % <-- or sig_xy, adapt if needed

svm_sel = sqrt( sigx_sel.^2 - sigx_sel.*sigy_sel + sigy_sel.^2 + 3*txy_sel.^2 );


[~, idx] = sort(n_sel);
n_sel   = n_sel(idx);
sigx_sel = sigx_sel(idx);
eps_sel = eps_sel(idx);
svm_sel = svm_sel(idx);

idx1 = (n_sel <= n1);
idx2 = (n_sel > n1) & (n_sel <= n1+n2); 
idx3 = (n_sel > n1+n2);               

% ---- Plots ----
figure(1); hold on; grid on;

plot(eps_sel(idx1), sigx_sel(idx1), '-o', 'DisplayName', 'Phase 1');
plot(eps_sel(idx2), sigx_sel(idx2), '-o', 'DisplayName', 'Phase 2');
plot(eps_sel(idx3), sigx_sel(idx3), '-o', 'DisplayName', 'Phase 3');

%plot(eps_sel(idx1), svm_sel(idx1), '-o', 'DisplayName', 'Phase 1');
%plot(eps_sel(idx2), svm_sel(idx2), '-o', 'DisplayName', 'Phase 2');
%plot(eps_sel(idx3), svm_sel(idx3), '-o', 'DisplayName', 'Phase 3');

xlabel('\epsilon');
ylabel('\sigma');
title('\bf stress-strain curve','Interpreter','tex');
legend('Location','best');





stress_log;

e0  = 1;
ig0 = 1;
mask = (e == e0) & (ig == ig0);

n_s    = n(mask);
epsp_s = epsp(mask);

[~, idx] = sort(n_s);
n_s    = n_s(idx);
epsp_s = epsp_s(idx);

figure; hold on; grid on;
plot(n_s, epsp_s, '-o');
xlabel('n'); ylabel('\epsilon_p');

xline(50,'--k');   % fin phase 1 si n1=50
xline(100,'--k');  % fin phase 2 si n2=50
title('\bf eps_p vs increment');