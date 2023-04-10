
clear 
close all

pi2 = 2 * pi;
domain = [0, pi2];
time = pi2/10;
a = 1;
lbound = 5;
ubound = 13;

exactFrequence = 2^ubound;
h = pi2 / exactFrequence;
k = h / 2;
n = single(time / k);
v0 = initial2(@testF2, h, k, domain, a);
x = (0:h:pi2)';
t = 0:k:time;
realvn = extSolution(@testF2, x, t, a);

% cal error convergence
loopVec = lbound:1:ubound;
error1 = zeros(size(loopVec));
for p = loopVec
clear h k n v0 t
localFrequence = 2^p;
h = pi2 / localFrequence;
k = h / 2;
n = single(time / k);
v0 = initial2(@testF2, h, k, domain, a);
vn = leapFrogPeriod(v0, h, k, n, a);
j = exactFrequence / localFrequence;
error1(p - lbound + 1) = L2(vn(:, end) - realvn(1:j:end, end));
end
% take k = 1.1h prove CFL conditino is necessary.
clear h k n v0 t
i = 6;
localFrequence = 2^loopVec(i);
h = pi2 / localFrequence;
k = 1.1 * h;
n = single(time / k);
v0 = initial2(@testF2, h, k, domain, a);
vn = leapFrogPeriod(v0, h, k, n, a);
j = exactFrequence / localFrequence;
CFLError1 = L2(vn(:, end) - realvn(1:j:end, end));
% plot(loopVec(i), CFLError, '*');


% initial function 2.
clear h k n v0 t realvn
exactFrequence = 2^ubound;
h = pi2 / exactFrequence;
k = h / 2;
n = single(time / k);
v0 = initial2(@testF3, h, k, domain, a);
x = (0:h:pi2)';
t = 0:k:time;
realvn = extSolution(@testF3, x, t, a);

% cal error convergence
loopVec = lbound:1:ubound;
error2 = zeros(size(loopVec));
for p = loopVec
clear h k n v0 t
localFrequence = 2^p;
h = pi2 / localFrequence;
k = h / 2;
n = single(time / k);
v0 = initial2(@testF3, h, k, domain, a);
vn = leapFrogPeriod(v0, h, k, n, a);
j = exactFrequence / localFrequence;
error2(p - lbound + 1) = L2(vn(:, end) - realvn(1:j:end, end));
end
% take k = 1.1h prove CFL conditino is necessary.
clear h k n v0 t
i = 6;
localFrequence = 2^loopVec(i);
h = pi2 / localFrequence;
k = 1.1 * h;
n = single(time / k);
v0 = initial2(@testF3, h, k, domain, a);
vn = leapFrogPeriod(v0, h, k, n, a);
j = exactFrequence / localFrequence;
CFLError2 = L2(vn(:, end) - realvn(1:j:end, end));
% plot(loopVec(i), CFLError, '*');


figure
yyaxis right
plot(loopVec, error2,'-'); hold on
ylabel('initial2 numerical error');
yyaxis left
plot(loopVec, error1,'-'); hold on
ylabel('initial1 numerical error');


figure
logError1 = log(error1);
logError2 = log(error2);
yyaxis right
plot(loopVec, logError2,'-');
ylabel('initial2 log(error)');
yyaxis left
plot(loopVec, logError1,'-');
ylabel('initial1 log(error)');
