% exact solution
excatFrequence = 5000;
jump = excatFrequence / 10;
h = pi2 / excatFrequence;
k = h / 2;
n = single(time / k);
v0 = initial(@testF1, h, domain);
realvn = laxWendroffPeriod(v0, h, k, n);

