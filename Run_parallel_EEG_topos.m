

clear
clc

allsubs = [1:10];

% open pool with right number of workers
nworkers = length(allsubs);
% open cbu parallel pool
P=cbupool(nworkers);
walltime_req = '50:00:00'; % good idea to check how long one iteration takes
%P.ResourceTemplate='-l nodes=^N^,mem=16GB,walltime=1:00:00'; %check whether increasing memory makes any difference to speed - it doesn't, and seems to be OK to prespecify the output matrix with just 4GB, so going with that

% try
%     P.ResourceTemplate=['-l nodes=^N^,mem=8GB,walltime=' walltime_req];
% catch %for SLURM (new cluster), it's a different field...
    P.SubmitArguments=['--ntasks=' num2str(nworkers) ...
        ' --mem-per-cpu=8G --time=' walltime_req];
% end
parpool(P,nworkers);

% send job for each sub to a different worker using parfor
tic();

parfor c = 1:length(allsubs)
    % CS remember to check if overwriting!! 
     [subject] = Revise_features_whole_trial(allsubs(c))
end
toc();

% close pool now
% matlabpool('close');
delete(gcp('nocreate'));


