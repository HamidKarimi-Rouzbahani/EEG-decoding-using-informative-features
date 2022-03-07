%% 15- CSP: this feature provides one two numbers for each pairs of trials of the two opposing categories
% % for all channels simultaneously
%
dataCat1 = randn(31,100);
dataCat2 = randn(31,100);


[W] = f_CSP(dataCat1,dataCat1);
J=W([1 end],:)*dataCat1;
CSP1=log([var(J(1,:)) var(J(2,:))])

% [PTranspose] = CSP(dataCat1,dataCat2);
% [Class1Class2] = spatFilt(dataCat1,PTranspose,2);
% CSP2=log([var(Class1Class2(1,:)) var(Class1Class2(2,:))])
