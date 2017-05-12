clear;
ndims = [300, 300, 50];
r = [5, 5, 20];
nmodes = length(r);

disp('Generate synthetic data...');
[data, ~] = synthetic(ndims, r,  'CORE');
data = sptensor(data);

err_type = 'AUC';
density = 0.05;

[data_train, P_train, data_test_post, P_test, data_mean, data_std] = tensplitspdata(data, density);
data_test = data .* P_test;

options.maxIter = 200;
options.maxK = 200 * ones(3, 1);
options.shrinkIter = 50;
tau = 100;
options.ex_dims = [];
options.eps = 5e-5;
options.rho = 0.4;
options.l_search = 1;
options.compute_test_per = 1;
options.accelerated = 0;
options.verbose = 1;
options.err_type = err_type;
options.AUC_size = 1e6;
[Xout, fw_obj, fw_info] = acctenfw(data_train, tau, data_test, P_test, options);
fw_runtime = fw_info.itime(1:fw_info.iter -1);
fw_test_err = fw_info.test_p(1:fw_info.iter-1);
plot(fw_runtime, fw_test_err , '-r', 'LineWidth', 1.5);
xlabel('\fontsize{14}time (sec)');
ylabel(strcat(['\fontsize{14}testing ',' ', err_type]));
grid on;
box on;
