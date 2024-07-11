clear;clc;

true_rts = zeros(4, 50);
max_rts = zeros(4, 50);
%test_rts = zeros(2, 50);

mu1 = 0;
mu2 = 2;
sigma1 = 1;
sigma2 = 1;
beta = 0.8;

count = 0;

for k = 0.01:0.01:0.5
    count = count + 1;
    
    %Rescaling of dielectric function
    w0 = -sigma1*k*1i;
    w1 = k*(mu2-mu1) - sigma2*k*1i;

    wdiff = w0-w1;

    %Let y = w - w0 and compute coefficients of polynomial in y
    coeffvec = [1, 2*wdiff, wdiff^2-1, -2*beta*wdiff, -beta*wdiff^2];
    
    %Solve for roots of quartic polynomial
    rts = roots(coeffvec);

%For Test - extraneous roots if beta = 1
    %abs(rts + wdiff)

%For Test - remove extraneous roots if beta = 1
    %rts_test = rts(find(abs(rts + wdiff) > 1e-4));
    %true_rts(:, count) = rts_test + w0 + k*mu1;
    
    % Transform roots of y back to roots of Omega + i*gamma
    true_rts(:, count) = rts + w0 + k*mu1;

    %Sort roots by greatest imaginary part
    temp = sort(1i*true_rts(:, count),'ComparisonMethod','real');
    true_rts(:, count) = -1i*temp;

%For Test - exact roots if beta = 1
    %test_rts(:, count) = [k*mu1 - 1 - sigma1*k*1i; k*mu1 + 1 - sigma1*k*1i];
end

%diff = norm(true_rts - test_rts, 2)

re_rts = real(true_rts);
im_rts = imag(true_rts);

kplot = 0.01:0.01:0.5;
fig = figure;
plot(kplot, re_rts(1,:)); % real part of root with maximal imaginary part
%plot(kplot, re_rts(1:4,:));
textInterp = 'latex';
set(fig, 'defaultTextInterpreter', textInterp);
xlabel('$k$')
ylabel('$\Omega(k)$')

fig = figure;
textInterp = 'latex';
set(fig, 'defaultTextInterpreter', textInterp);
plot(kplot, im_rts(1,:)); % imag part of root with maximal imaginary part
%plot(kplot, im_rts(1:4,:));
xlabel('$k$')
ylabel('$\gamma(k)$')


%For Test - plot real/imag parts of computed roots vs test roots (beta = 1)
%kplot = 0.01:0.01:0.5;
%figure;
%plot(kplot, real(true_rts(1, :)));
% hold on
% plot(kplot, real(test_rts(1, :)));
% figure;
% plot(kplot, real(true_rts(2, :)));
% hold on
% plot(kplot, real(test_rts(2, :)));
% figure;
% plot(kplot, imag(true_rts(1, :)));
% hold on
% plot(kplot, imag(test_rts(1, :)));
% figure;
% plot(kplot, imag(true_rts(2, :)));
% hold on
% plot(kplot, imag(test_rts(2, :)));