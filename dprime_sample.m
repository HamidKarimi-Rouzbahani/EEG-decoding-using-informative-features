clc;
clear all;
close all;
sample_data=randn(4,10); % rows*columns = neural space dimension * trials of the two conditions 
                         % here each condition has 5 trials 
% sample_data=[randn(4,5) 3*randn(4,5)+10]; % another example with variable
% mean and variance
                                                  
mult=50; % some constant which is used
         % when finding the connecting line between two clusters, it should be large
X_mean1=mean(sample_data(:,1:5),2)';
X_mean2=mean(sample_data(:,6:10),2)';
lent=(X_mean2-X_mean1)';
t=[-mult mult];
dim1=X_mean1'+t(1)*lent;
dim2=X_mean1'+t(2)*lent;
X_projected=zeros(size(sample_data,2),1);
ccc
for i=1:size(sample_data,2)
    X_projected(i,1)=pdist2(project_point_to_line_segment(dim1',dim2',sample_data(:,i)'),dim1');
end
mu2=nanmean(X_projected(1:5));
mu1=nanmean(X_projected(6:10));
var2=nanvar(X_projected(1:5));
var1=nanvar(X_projected(6:10));
dprime=(mu1-mu2)./sqrt(0.5*(var1+var2));
