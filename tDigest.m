% Use t-digest to combine two distributions and then calculate the quantiles.

% To make sure we get the same results as in paper we need to set the seed
% for the random number generator
rng(1)

%% This first part is just looking at the scaling function
% You can change delta, and change which scaling function is being use, to
% see how this places the quantiles of the t-Digest.
delta=10^4;

k0=@(q,delta) (delta/2)*q;
k1=@(q,delta) (delta/(2*pi))*asin(2*q-1);
k2=@(q,delta) (delta/(4*log(n/delta)+24))*log(q/(1-q));
k3=@(q,delta) (delta/(4*log(n/delta)+24))*(log(2*q)*(q<=0.5)-log(2*(1-q))*(q>0.5));
% Choose which k-scaling function to use
kfn=@(q,delta) k1(q,delta);

k0inv=@(k,delta) k*2/delta;
k1inv=@(k,delta) (sin(k/(delta/(2*pi)))+1)/2;
% k2inv=@(k,delta) 
% k3inv=@(k,delta) 
% Also need the inverse fn of k
kinvfn=@(k,delta) k1inv(k,delta);

qgrid=0:0.01:1;

kgrid=kfn(qgrid,delta);

qgridagain=kinvfn(kgrid,delta);

max(abs(qgrid-qgridagain)) % Check: This will be zero if I define the inverse functions correctly

% Note to self: following looks correct
figure(1)
plot(qgrid,kgrid)
ylabel('k')
xlabel('q')
title('Notice how more points are nearer the extreme of q')
% If you use other scaling functions k2 or k3 there will be even more points will be near extremes.

% Take a look at which percentiles will make up the digest
q=0;
qlimit=kinvfn(kfn(q,delta)+1,delta);
qlimitvec=[];
while qlimit<1
    qlimitvec=[qlimitvec,qlimit];
    qlimit=kinvfn(kfn(qlimit,delta)+1,delta);
end
qlimitvec
% delta=10:    length(qlimitvec) is 4
% delta=100:   length(qlimitvec) is 49
% delta=1000:  length(qlimitvec) is 499
% delta=10000: length(qlimitvec) is 4999
% Actually limit(qlimitvec) gets determined in practice by both delta and
% the dataset itself. It is often roughly this length, but can be slightly more.

%% Create first distribution: 10^6 observations of a Uniform[0,10] variable
% Distribution is a collection of grid points and their associated weights
evalgrid1=10*rand(10^6,1);
weights1=ones(length(evalgrid1),1)/length(evalgrid1);
weights1=weights1/sum(weights1); % Normalize mass to one

mean1=sum(evalgrid1.*weights1);

[sortevalgrid1,sortIndex]=sort(evalgrid1);
sortweights1=weights1(sortIndex);

mean2=sum(sortevalgrid1.*sortweights1);
% mean1 and mean2 return same, so sort is implemented correctly

% Calculate the median directly from the distribution
temp=cumsum(sortweights1);
[~,medianindex]=min(abs(temp-0.5));
temp(medianindex); % This should be roughly 0.5
medianvalue_1=sortevalgrid1(medianindex);
% Now the 99th percentile direct from the distribution
[~,index]=min(abs(temp-0.99));
temp(index); % This should be roughly 0.99
percentile99_1=sortevalgrid1(index);

%% Create second distribution: 10^7 observations of a Uniform[0,5] variable

% Distribution is a collection of grid points and their associated weights
evalgrid2=5*rand(10^7,1);
weights2=rand(length(evalgrid2),1);
weights2=weights2/sum(weights2); % Normalize mass to one

[sortevalgrid2,sortIndex2]=sort(evalgrid2);
sortweights2=weights2(sortIndex2);

% First, calculate the median directly from the distribution
temp=cumsum(sortweights2);
[~,medianindex]=min(abs(temp-0.5));
temp(medianindex); % This should be roughly 0.5
medianvalue_2=sortevalgrid2(medianindex);
% Now the 99th percentile direct from the distribution
[~,index]=min(abs(temp-0.99));
temp(index); % This should be roughly 0.99
percentile99_2=sortevalgrid2(index);

%% Do a full merge and calculate the percentiles of this

% To merge, we first need to give the relative weights of the two distributions
relativeweights=[0.6,0.4];

evalgridmerge=[evalgrid1;evalgrid2];
weightmerge=[relativeweights(1)*weights1;relativeweights(2)*weights2];

[sortevalgridmerge,sortIndexMerge]=sort(evalgridmerge);
sortweightsmerge=weightmerge(sortIndexMerge);

% First, calculate the median directly from the distribution
temp=cumsum(sortweightsmerge);
[~,medianindex]=min(abs(temp-0.5));
temp(medianindex); % This should be roughly 0.5
medianvalue_merge=sortevalgridmerge(medianindex);
% Now the 99th percentile direct from the distribution
[~,index]=min(abs(temp-0.99));
temp(index); % This should be roughly 0.99
percentile99_merge=sortevalgridmerge(index);


%% Now use t-Digests to estimate the median and 99th percentile of each of the two distributions, and of the merged distribution
delta=10000;

[C1,digestweights1,qlimitvec1]=createDigest(evalgrid1, weights1,delta);

[C2,digestweights2,qlimitvec2]=createDigest(evalgrid2, weights2,delta);

results1=interp1(qlimitvec1,C1,[0.5,0.99]);
results2=interp1(qlimitvec2,C2,[0.5,0.99]);


fprintf('For first distribution \n')
fprintf('Precise calculation gives median=%8.4f and 99th percentile=%8.4f \n',medianvalue_1,percentile99_1)
fprintf('Digest calculation gives median=%8.4f and 99th percentile=%8.4f \n',results1)

fprintf('For second distribution \n')
fprintf('Precise calculation gives median=%8.4f and 99th percentile=%8.4f \n',medianvalue_2,percentile99_2)
fprintf('Digest calculation gives median=%8.4f and 99th percentile=%8.4f \n',results2)


% relativeweights=[0.4,0.6]; % is defined above

Cmerge=[C1;C2];
digestweightsmerge=[relativeweights(1)*digestweights1;relativeweights(2)*digestweights2];

[C,digestweights,qlimitvec]=mergeDigest(Cmerge, digestweightsmerge, delta);

results_merge=interp1(qlimitvec,C,[0.5,0.99]);

quantiles=[0.5,0.99];

% How to do quantile cutoffs and quantile means.
quantilecutoffs=interp1(qlimitvec,C,quantiles);
quantilemeans=zeros(length(quantilecutoffs)+1,1);
Ctimesdisgestweights=C.*digestweights;
quantilemeans(1)=sum(Ctimesdisgestweights(qlimitvec<quantiles(1)))/sum(digestweights(qlimitvec<quantiles(1)));
for qq=2:length(quantilecutoffs)
    quantilemeans(qq)=sum(Ctimesdisgestweights(logical((qlimitvec>quantiles(qq-1)).*(qlimitvec<quantiles(qq)))))/sum(digestweights(logical((qlimitvec>quantiles(qq-1)).*(qlimitvec<quantiles(qq)))));
end
quantilemeans(end)=sum(Ctimesdisgestweights(qlimitvec>quantiles(end)))/sum(digestweights(qlimitvec>quantiles(end)));

fprintf('For merge of two distributions \n')
fprintf('Precise calculation gives median=%8.4f and 99th percentile=%8.4f \n',medianvalue_merge,percentile99_merge)
fprintf('Digest calculation gives median=%8.4f and 99th percentile=%8.4f \n',results_merge)

%% Add a third distribution
% Distribution is a collection of grid points and their associated weights
evalgrid3=3*randn(10^6,1);
weights3=ones(length(evalgrid3),1)/length(evalgrid3);
weights3=weights3/sum(weights3); % Normalize mass to one

[sortevalgrid3,sortIndex3]=sort(evalgrid3);
sortweights3=weights3(sortIndex3);

[C3,digestweights3,qlimitvec3]=createDigest(evalgrid3, weights3,delta);

results3=interp1(qlimitvec3,C3,[0.5,0.99]);

% First, calculate the median directly from the distribution
temp=cumsum(sortweights3);
[~,medianindex]=min(abs(temp-0.5));
temp(medianindex); % This should be roughly 0.5
medianvalue_3=sortevalgrid3(medianindex);
% Now the 99th percentile direct from the distribution
[~,index]=min(abs(temp-0.99));
temp(index); % This should be roughly 0.99
percentile99_3=sortevalgrid3(index);

fprintf('For third distribution \n')
fprintf('Precise calculation gives median=%8.4f and 99th percentile=%8.4f \n',medianvalue_3,percentile99_3)
fprintf('Digest calculation gives median=%8.4f and 99th percentile=%8.4f \n',results3)

%% Merge all three distribitons
relativeweights=[0.4,0.3,0.3];

Cmerge=[C1;C2;C3];
digestweightsmerge=[relativeweights(1)*digestweights1;relativeweights(2)*digestweights2;relativeweights(3)*digestweights3];

[C,digestweights,qlimitvec]=mergeDigest(Cmerge, digestweightsmerge, delta);

results=interp1(qlimitvec,C,[0.5,0.99]);

% Do a full merge and calculate the percentiles of this
% To merge, we first need to give the relative weights of the two distributions

evalgridmerge=[evalgrid1;evalgrid2;evalgrid3];
weightmerge=[relativeweights(1)*weights1;relativeweights(2)*weights2;relativeweights(3)*weights3];

[sortevalgridmerge,sortIndexMerge]=sort(evalgridmerge);
sortweightsmerge=weightmerge(sortIndexMerge);

% First, calculate the median directly from the distribution
temp=cumsum(sortweightsmerge);
[~,medianindex]=min(abs(temp-0.5));
temp(medianindex); % This should be roughly 0.5
medianvalue=sortevalgridmerge(medianindex);
% Now the 99th percentile direct from the distribution
[~,index]=min(abs(temp-0.99));
temp(index); % This should be roughly 0.99
percentile99=sortevalgridmerge(index);

fprintf('For merge of three distribitons \n')
fprintf('Precise calculation gives median=%8.4f and 99th percentile=%8.4f \n',medianvalue,percentile99)
fprintf('Digest calculation gives median=%8.4f and 99th percentile=%8.4f \n',results)





%% Just try out a very different kind of distribution: 0.3 with value 4, 0.2 with value 5, 0.5 with value 7
evalgrid4=[4*ones(3*10^5,1);5*ones(2*10^5,1);7*ones(5*10^5,1)];
% Distribution is a collection of grid points and their associated weights
weights4=ones(length(evalgrid4),1)/length(evalgrid4);
weights4=weights4/sum(weights4); % Normalize mass to one

mean4=sum(evalgrid4.*weights4); % Is correct: 0.3*4+0.2*5+0.5*7=5.7

[sortevalgrid4,sortIndex]=sort(evalgrid4); % Note the sort steps are redundant, but I will do them anyway
sortweights4=weights4(sortIndex);

mean4=sum(sortevalgrid4.*sortweights4);
% mean1 and mean2 return same, so sort is implemented correctly

% Calculate the 29.9 percentile directly from the distribution
temp=cumsum(sortweights4);
[~,p299index]=min(abs(temp-0.299));
temp(p299index); % This should be roughly 0.299
p299=sortevalgrid4(p299index);
% Calculate the 30.1 percentile directly from the distribution
temp=cumsum(sortweights4);
[~,p301index]=min(abs(temp-0.301));
temp(p301index); % This should be roughly 0.301
p301=sortevalgrid4(p301index);
% Calculate the 49.9 percentile directly from the distribution
temp=cumsum(sortweights4);
[~,p499index]=min(abs(temp-0.499));
temp(p499index); % This should be roughly 0.499
p499=sortevalgrid4(p499index);
% Calculate the 50.1 percentile directly from the distribution
temp=cumsum(sortweights4);
[~,p501index]=min(abs(temp-0.501));
temp(p501index); % This should be roughly 0.501
p501=sortevalgrid4(p501index);

% Now use the digest
[C4,digestweights4,qlimitvec4]=createDigest(evalgrid4, weights4,delta);

results4=interp1(qlimitvec4,C4,[0.299,0.301,0.499,0.501]);


fprintf('For fourth distribution \n')
fprintf('Precise calculation gives percentiles: 29.9th=%8.4f, 30.1th=%8.4f, 49.9th=%8.4f, 50.1th=%8.4f \n',p299,p301,p499,p501)
fprintf('Digest calculation gives percentiles: 29.9th=%8.4f, 30.1th=%8.4f, 49.9th=%8.4f, 50.1th=%8.4f \n',results4)


