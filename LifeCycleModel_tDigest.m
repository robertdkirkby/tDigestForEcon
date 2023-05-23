%% Life-Cycle Model 24: Using Permanent Type to model fixed-effects
% The exogenous process on labor efficiency units now uses an approach common in the literature:
% Labor efficiency units are a combination of four components:
% 1) kappa_j, a deterministic profile of age
% 2) z, a persistent markov shock
% 3) e, a transitory i.i.d. shock
% 4) alpha_i, a fixed effect
%
% All of these were already present in Life-Cycle Model 11, except the fixed-effect alpha_i
%
% We want to have five different possible values of alpha_i, and to do this we use the 'permanent type', PType, feature.
%
% We here use the easiest approach, we will create alpha_i as a vector with
% five values. And set N_i=5 (N_i is the number of permanent types). The
% codes will automatically realise that because alpha_i has N_i values it
% is different by permanent type.
%
% To compute the agent distribution we need to say how many agents are of
% each type. We denote this 'alphadist' and it is a vector of weights that
% sum to one (we need to put the name of this in PTypeDistParamNames, like
% we would for a discount factor in DiscountFactorParamNames).
%
% All the commands we run look slightly different as they are the PType
% versions of the commands until now. Permanent types are much more
% powerful than just fixed-effects, and later models will show more options.
%
% Model statistics, like the life-cycle profiles we calculate here, are
% reported both for each permanent type (that is to say, conditional on the
% permanent type), and 'grouped' across the permanent types. 

%% How does VFI Toolkit think about this?
%
% One decision variable: h, labour hours worked
% One endogenous state variable: a, assets (total household savings)
% Two stochastic exogenous state variables: z and e, persistent and transitory shocks to labor efficiency units, respectively
% Age: j
% Permanent types: i

%% Begin setting up to use VFI Toolkit to solve
% Lets model agents from age 20 to age 100, so 81 periods

Params.agejshifter=19; % Age 20 minus one. Makes keeping track of actual age easy in terms of model age
Params.J=100-Params.agejshifter; % =81, Number of period in life-cycle

% Grid sizes to use
n_d=51; % Endogenous labour choice (fraction of time worked)
n_a=201; % Endogenous asset holdings
n_z=21; % Exogenous labor productivity units shocks, persistent and transitiory
n_e=3;
N_i=5; % Permanent type of agents
N_j=Params.J; % Number of periods in finite horizon

%% The parameter that depends on the permanent type
% Fixed-effect (parameter that varies by permanent type)
Params.alpha_i=exp([0.5,0.3,0,-0.3,-0.5]); % Roughly: increase earnings by 50%, 30%, 0, -30%, -50%

% This is only seperate from the other parameters to make it easier to see
% what has changed in the codes.

%% How many of each permanent type are there
PTypeDistParamNames={'alphadist'};
Params.alphadist=[0.1,0.2,0.4,0.1,0.2]; % Must sum to one
% Note: this is not relevant to solving the value function, but is needed for
% stationary distribition. It then gets encoded into the StationaryDist and
% so is not needed for things like life-cycle profiles.

%% Parameters

% Discount rate
Params.beta = 0.96;
% Preferences
Params.sigma = 2; % Coeff of relative risk aversion (curvature of consumption)
Params.eta = 1.5; % Curvature of leisure (This will end up being 1/Frisch elasty)
Params.psi = 10; % Weight on leisure

% Prices
Params.w=1; % Wage
Params.r=0.05; % Interest rate (0.05 is 5%)

% Demographics
Params.agej=1:1:Params.J; % Is a vector of all the agej: 1,2,3,...,J
Params.Jr=46;

% Pensions
Params.pension=0.3;

% Age-dependent labor productivity units
Params.kappa_j=[linspace(0.5,2,Params.Jr-15),linspace(2,1,14),zeros(1,Params.J-Params.Jr+1)];
% persistent AR(1) process on idiosyncratic labor productivity units
Params.rho_z=0.9;
Params.sigma_epsilon_z=0.02;
% transitiory iid normal process on idiosyncratic labor productivity units
Params.sigma_epsilon_e=0.2; % Implictly, rho_e=0

% Conditional survival probabilities: sj is the probability of surviving to be age j+1, given alive at age j
% Most countries have calculations of these (as they are used by the government departments that oversee pensions)
% In fact I will here get data on the conditional death probabilities, and then survival is just 1-death.
% Here I just use them for the US, taken from "National Vital Statistics Report, volume 58, number 10, March 2010."
% I took them from first column (qx) of Table 1 (Total Population)
% Conditional death probabilities
Params.dj=[0.006879, 0.000463, 0.000307, 0.000220, 0.000184, 0.000172, 0.000160, 0.000149, 0.000133, 0.000114, 0.000100, 0.000105, 0.000143, 0.000221, 0.000329, 0.000449, 0.000563, 0.000667, 0.000753, 0.000823,...
    0.000894, 0.000962, 0.001005, 0.001016, 0.001003, 0.000983, 0.000967, 0.000960, 0.000970, 0.000994, 0.001027, 0.001065, 0.001115, 0.001154, 0.001209, 0.001271, 0.001351, 0.001460, 0.001603, 0.001769, 0.001943, 0.002120, 0.002311, 0.002520, 0.002747, 0.002989, 0.003242, 0.003512, 0.003803, 0.004118, 0.004464, 0.004837, 0.005217, 0.005591, 0.005963, 0.006346, 0.006768, 0.007261, 0.007866, 0.008596, 0.009473, 0.010450, 0.011456, 0.012407, 0.013320, 0.014299, 0.015323,...
    0.016558, 0.018029, 0.019723, 0.021607, 0.023723, 0.026143, 0.028892, 0.031988, 0.035476, 0.039238, 0.043382, 0.047941, 0.052953, 0.058457, 0.064494,...
    0.071107, 0.078342, 0.086244, 0.094861, 0.104242, 0.114432, 0.125479, 0.137427, 0.150317, 0.164187, 0.179066, 0.194979, 0.211941, 0.229957, 0.249020, 0.269112, 0.290198, 0.312231, 1.000000]; 
% dj covers Ages 0 to 100
Params.sj=1-Params.dj(21:101); % Conditional survival probabilities
Params.sj(end)=0; % In the present model the last period (j=J) value of sj is actually irrelevant

% Warm glow of bequest
Params.warmglow1=0.3; % (relative) importance of bequests
Params.warmglow2=3; % bliss point of bequests (essentially, the target amount)
Params.warmglow3=Params.sigma; % By using the same curvature as the utility of consumption it makes it much easier to guess appropraite parameter values for the warm glow

%% Grids
% The ^3 means that there are more points near 0 and near 10. We know from
% theory that the value function will be more 'curved' near zero assets,
% and putting more points near curvature (where the derivative changes the most) increases accuracy of results.
a_grid=10*(linspace(0,1,n_a).^3)'; % The ^3 means most points are near zero, which is where the derivative of the value fn changes most.

% First, the AR(1) process z
if Params.rho_z<0.99
    [z_grid,pi_z]=discretizeAR1_FarmerToda(0,Params.rho_z,Params.sigma_epsilon_z,n_z);
elseif Params.rho_z>=0.99 % Rouwenhourst performs better than Farmer-Toda when the autocorrelation is very high
    [z_grid,pi_z]=discretizeAR1_Rouwenhorst(0,Params.rho_z,Params.sigma_epsilon_z,n_z);
end
z_grid=exp(z_grid); % Take exponential of the grid
[mean_z,~,~,~]=MarkovChainMoments(z_grid,pi_z); % Calculate the mean of the grid so as can normalise it
z_grid=z_grid./mean_z; % Normalise the grid on z (so that the mean of z is 1)
% Now the iid normal process e
[e_grid,pi_e]=discretizeAR1_FarmerToda(0,0,Params.sigma_epsilon_e,n_e);
e_grid=exp(e_grid); % Take exponential of the grid
pi_e=pi_e(1,:)'; % Because it is iid, the distribution is just the first row (all rows are identical). We use pi_e as a column vector for VFI Toolkit to handle iid variables.
mean_e=pi_e'*e_grid; % Because it is iid, pi_e is the stationary distribution (you could just use MarkovChainMoments(), I just wanted to demonstate a handy trick)
e_grid=e_grid./mean_e; % Normalise the grid on z (so that the mean of e is 1)
% To use e variables we have to put them into the vfoptions and simoptions
vfoptions.n_e=n_e;
vfoptions.e_grid=e_grid;
vfoptions.pi_e=pi_e;
simoptions.n_e=vfoptions.n_e;
simoptions.e_grid=vfoptions.e_grid;
simoptions.pi_e=vfoptions.pi_e;


% Grid for labour choice
h_grid=linspace(0,1,n_d)'; % Notice that it is imposing the 0<=h<=1 condition implicitly
% Switch into toolkit notation
d_grid=h_grid;

%% Now, create the return function 
DiscountFactorParamNames={'beta','sj'};

% Notice: have added alpha_i to inputs (relative to Life-Cycle Model 11 which this extends)
ReturnFn=@(h,aprime,a,z,e,w,sigma,psi,eta,agej,Jr,pension,r,alpha_i,kappa_j,warmglow1,warmglow2,warmglow3,beta,sj) LifeCycleModel24_ReturnFn(h,aprime,a,z,e,w,sigma,psi,eta,agej,Jr,pension,r,alpha_i,kappa_j,warmglow1,warmglow2,warmglow3,beta,sj)

%% Now solve the value function iteration problem, just to check that things are working before we go to General Equilbrium
disp('Test ValueFnIter')
tic;
vfoptions.verbose=1; % Just so we can see feedback on progress
[V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
toc

% V is now a structure containing the value function for each permanent type
V
% So for example for permanent type 1, V_i(a,z,e,j)
size(V.ptype001)
% Policy likewise depends on type, e.g.,
size(Policy.ptype001)
% which is the same as
[length(n_d)+length(n_a),n_a,n_z,n_e,N_j]
% The other permanent types are numbered through N_i, so the last in our current case is
% V.ptype005
% The n_a,n_z,n_e,N_j represent the state on which the decisions/policys
% depend, and there is one decision for each decision variable 'd' and each
% endogenous state variable 'a', and one for the markov exogenous state variable
% 'z', and one for the markov exogenous state variable 'e'.
% The different permanent types, i, are essentially just seperate value function problems.

%% We won't plot the Value and Policy fn, but thinking out how you would might be a good way to check you understand the form of V and Policy

%% Now, we want to graph Life-Cycle Profiles

%% Initial distribution of agents at birth (j=1)
% Before we plot the life-cycle profiles we have to define how agents are
% at age j=1. We will give them all zero assets.
jequaloneDist=zeros([n_a,n_z,n_e],'gpuArray'); % Put no households anywhere on grid
jequaloneDist(1,floor((n_z+1)/2),floor((n_e+1)/2))=1; % All agents start with zero assets, and the median value of each shock

% Anything that is not made to depend on the permanent type is
% automatically assumed to be independent of the permanent type (that is,
% identical across permanent types). This includes things like the initial
% distribution, jequaloneDist

%% We now compute the 'stationary distribution' of households
% Start with a mass of one at initial age, use the conditional survival
% probabilities sj to calculate the mass of those who survive to next
% period, repeat. Once done for all ages, normalize to one
Params.mewj=ones(1,Params.J); % Marginal distribution of households over age
for jj=2:length(Params.mewj)
    Params.mewj(jj)=Params.sj(jj-1)*Params.mewj(jj-1);
end
Params.mewj=Params.mewj./sum(Params.mewj); % Normalize to one
AgeWeightsParamNames={'mewj'}; % So VFI Toolkit knows which parameter is the mass of agents of each age
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,N_i,pi_z,Params,simoptions);
% Again, we will explain in a later model what the stationary distribution
% is, it is not important for our current goal of graphing the life-cycle profile

%% FnsToEvaluate are how we say what we want to graph the life-cycles of
% Like with return function, we have to include (h,aprime,a,z) as first
% inputs, then just any relevant parameters.
FnsToEvaluate.fractiontimeworked=@(h,aprime,a,z,e) h; % h is fraction of time worked
FnsToEvaluate.earnings=@(h,aprime,a,z,e,w,kappa_j,alpha_i) w*h*kappa_j*alpha_i*z*e; % w*kappa_j*h*z*e is the labor earnings
FnsToEvaluate.assets=@(h,aprime,a,z,e) a; % a is the current asset holdings
FnsToEvaluate.alpha_i=@(h,aprime,a,z,e,alpha_i) alpha_i; % alpha_i is the fixed effect
FnsToEvaluate.agej=@(h,aprime,a,z,e,agej) agej; % alpha_i is the fixed effect
FnsToEvaluate.savingsrate=@(h,aprime,a,z,e,r,w,kappa_j,alpha_i) (aprime-a)/(r*a+w*h*kappa_j*alpha_i*z*e); % savings as a fraction of income

% notice that we have called these fractiontimeworked, earnings and assets
% Have added alpha_i so that we can see how this evaluates differently across the different permanent types of agents
% Note that alpha_i also appears in the function for earnings

%% Calculate the life-cycle profiles
% These are not actually used anywhere in the t-Digest paper
AgeConditionalStats=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate, Params,n_d,n_a,n_z,N_j,N_i,d_grid, a_grid, z_grid, simoptions);

rr=1;
% NOTE: To do the run time calculations mentioned in the paper just uncomment the following for loop over rr (and uncomment the 'end' around line 324)
% Repeat=100;
% for rr=1:Repeat
%     rr

%% Calculate some statistics of the agent distribution
tic;
AllStats=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate, Params,n_d,n_a,n_z,N_j,N_i,d_grid, a_grid, z_grid, simoptions);
timetdigest(rr)=toc
% Note: This command used t-Digests for the grouped statistics

%% Now calculate the statistics exactly (without using t-Digests, which is what VFI Toolkit does)
tic;
simoptions.parallel=1;
% This next example is only for comparing times of running on cpu to running on gpu. It is not actually used for anything.
ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate, Params,n_d,n_a,n_z,N_j,N_i,d_grid, a_grid, z_grid, simoptions);
timeexactcpu(rr)=toc

tic;
simoptions.parallel=2;
% First, evaluate the functions on the grid
ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate, Params,n_d,n_a,n_z,N_j,N_i,d_grid, a_grid, z_grid, simoptions);
timeexactgpu(rr)=toc


% Now calculate the median and 99th percentile directly from these
FullHoursValues=zeros(n_a,n_z,n_e,N_j,N_i);
FullHoursValues(:,:,:,:,1)=ValuesOnGrid.fractiontimeworked.ptype001;
FullHoursValues(:,:,:,:,2)=ValuesOnGrid.fractiontimeworked.ptype002;
FullHoursValues(:,:,:,:,3)=ValuesOnGrid.fractiontimeworked.ptype003;
FullHoursValues(:,:,:,:,4)=ValuesOnGrid.fractiontimeworked.ptype004;
FullHoursValues(:,:,:,:,5)=ValuesOnGrid.fractiontimeworked.ptype005;
FullHoursValues=reshape(FullHoursValues,[numel(FullHoursValues),1]);
% Now calculate the median and 99th percentile directly from these
FullEarningsValues=zeros(n_a,n_z,n_e,N_j,N_i);
FullEarningsValues(:,:,:,:,1)=ValuesOnGrid.earnings.ptype001;
FullEarningsValues(:,:,:,:,2)=ValuesOnGrid.earnings.ptype002;
FullEarningsValues(:,:,:,:,3)=ValuesOnGrid.earnings.ptype003;
FullEarningsValues(:,:,:,:,4)=ValuesOnGrid.earnings.ptype004;
FullEarningsValues(:,:,:,:,5)=ValuesOnGrid.earnings.ptype005;
FullEarningsValues=reshape(FullEarningsValues,[numel(FullEarningsValues),1]);
% Now calculate the median and 99th percentile directly from these
FullAssetsValues=zeros(n_a,n_z,n_e,N_j,N_i);
FullAssetsValues(:,:,:,:,1)=ValuesOnGrid.assets.ptype001;
FullAssetsValues(:,:,:,:,2)=ValuesOnGrid.assets.ptype002;
FullAssetsValues(:,:,:,:,3)=ValuesOnGrid.assets.ptype003;
FullAssetsValues(:,:,:,:,4)=ValuesOnGrid.assets.ptype004;
FullAssetsValues(:,:,:,:,5)=ValuesOnGrid.assets.ptype005;
FullAssetsValues=reshape(FullAssetsValues,[numel(FullAssetsValues),1]);
% Now calculate the median and 99th percentile directly from these
FullSavingsRateValues=zeros(n_a,n_z,n_e,N_j,N_i);
FullSavingsRateValues(:,:,:,:,1)=ValuesOnGrid.savingsrate.ptype001;
FullSavingsRateValues(:,:,:,:,2)=ValuesOnGrid.savingsrate.ptype002;
FullSavingsRateValues(:,:,:,:,3)=ValuesOnGrid.savingsrate.ptype003;
FullSavingsRateValues(:,:,:,:,4)=ValuesOnGrid.savingsrate.ptype004;
FullSavingsRateValues(:,:,:,:,5)=ValuesOnGrid.savingsrate.ptype005;
FullSavingsRateValues=reshape(FullSavingsRateValues,[numel(FullSavingsRateValues),1]);

%
FullAgentDist=zeros(n_a,n_z,n_e,N_j,N_i);
FullAgentDist(:,:,:,:,1)=StationaryDist.ptype001*StationaryDist.ptweights(1);
FullAgentDist(:,:,:,:,2)=StationaryDist.ptype002*StationaryDist.ptweights(2);
FullAgentDist(:,:,:,:,3)=StationaryDist.ptype003*StationaryDist.ptweights(3);
FullAgentDist(:,:,:,:,4)=StationaryDist.ptype004*StationaryDist.ptweights(4);
FullAgentDist(:,:,:,:,5)=StationaryDist.ptype005*StationaryDist.ptweights(5);
FullAgentDist=reshape(FullAgentDist,[numel(FullAgentDist),1]);

[SortedValues,SortedValues_index] = sort(FullHoursValues);
SortedStationaryDistVec=FullAgentDist(SortedValues_index);
median_index=find(cumsum(SortedStationaryDistVec)>=0.5,1,'first');
FractionTimeWorked_median=SortedValues(median_index);
p99_index=find(cumsum(SortedStationaryDistVec)>=0.99,1,'first');
FractionTimeWorked_99=SortedValues(p99_index);

[SortedValues,SortedValues_index] = sort(FullEarningsValues);
SortedStationaryDistVec=FullAgentDist(SortedValues_index);
median_index=find(cumsum(SortedStationaryDistVec)>=0.5,1,'first');
Earnings_median=SortedValues(median_index);
p99_index=find(cumsum(SortedStationaryDistVec)>=0.99,1,'first');
Earnings_99=SortedValues(p99_index);

[SortedValues,SortedValues_index] = sort(FullAssetsValues);
SortedStationaryDistVec=FullAgentDist(SortedValues_index);
median_index=find(cumsum(SortedStationaryDistVec)>=0.5,1,'first');
Assets_median=SortedValues(median_index);
p99_index=find(cumsum(SortedStationaryDistVec)>=0.99,1,'first');
Assets_99=SortedValues(p99_index);

[SortedValues,SortedValues_index] = sort(FullSavingsRateValues);
SortedStationaryDistVec=FullAgentDist(SortedValues_index);
median_index=find(cumsum(SortedStationaryDistVec)>=0.5,1,'first');
SavingsRate_median=SortedValues(median_index);
p99_index=find(cumsum(SortedStationaryDistVec)>=0.99,1,'first');
SavingsRate_99=SortedValues(p99_index);


timeexactextra(rr)=toc

% end

fprintf('Median of fraction of time worked is %8.4f with exact, %8.4f with tDigest \n',FractionTimeWorked_median,AllStats.fractiontimeworked.Median)
fprintf('99th percentile of fraction of time worked is %8.4f with exact, %8.4f with tDigest \n',FractionTimeWorked_99,AllStats.fractiontimeworked.Percentile99th)

fprintf('Median of earnings is %8.4f with exact, %8.4f with tDigest \n',Earnings_median,AllStats.earnings.Median)
fprintf('99th percentile of earnings is %8.4f with exact, %8.4f with tDigest \n',Earnings_99,AllStats.earnings.Percentile99th)

fprintf('Median of assets holdings is %8.4f with exact, %8.4f with tDigest \n',Assets_median,AllStats.assets.Median)
fprintf('99th percentile of assets holdings is %8.4f with exact, %8.4f with tDigest \n',Assets_99,AllStats.assets.Percentile99th)

fprintf('Median of savings rate is %8.4f with exact, %8.4f with tDigest \n',SavingsRate_median,AllStats.savingsrate.Median)
fprintf('99th percentile of savings rate is %8.4f with exact, %8.4f with tDigest \n',SavingsRate_99,AllStats.savingsrate.Percentile99th)

%% A couple of further checks of results, these do not involve t-Digests but are just about ensuring everything is working as expected
% We can also check the conditional on ptype stats (note that these do not involve t-Digest)
AssetDist_ptype001=ValuesOnGrid.assets.ptype001;
AssetDist_ptype001=reshape(AssetDist_ptype001,[numel(AssetDist_ptype001),1]);
AgentDist_ptype001=StationaryDist.ptype001;
[SortedValues_ptype001,SortedValues_index] = sort(AssetDist_ptype001);
SortedStationaryDistVec_ptype001=AgentDist_ptype001(SortedValues_index);
median_index=find(cumsum(SortedStationaryDistVec_ptype001)>=0.5,1,'first');
Assets_ptype001_median=SortedValues_ptype001(median_index);
p99_index=find(cumsum(SortedStationaryDistVec_ptype001)>=0.99,1,'first');
Assets_ptype001_99=SortedValues_ptype001(p99_index);

fprintf('Median of assets holdings of first permanent type is %8.4f with handcrafted, %8.4f with toolkit \n',Assets_ptype001_median,AllStats.assets.ptype001.Median)
fprintf('99th percentile of assets holdings of first permanent type is %8.4f with handcrafted, %8.4f with toolkit \n',Assets_ptype001_99,AllStats.assets.ptype001.Percentile99th)

% Check the means as this doesn't involve t-Digest, just as an extra test
% that underlying aspects are all correct.
% For assets
AssetMean=sum(SortedStationaryDistVec.*SortedValues)
AllStats.assets.Mean
% These are identical :)

% Check the means conditional on ptype
AssetMean_ptype001=sum(SortedStationaryDistVec_ptype001.*SortedValues_ptype001)
AllStats.assets.ptype001.Mean
% These are identical :)


%
CumSumSortedWeights=cumsum(SortedStationaryDistVec_ptype001);
cutoffvalues=prctile(CumSumSortedWeights,[50,90,95,99]); % Note: This could be faster done off the digest
cutoffvalues(1)


%% Create table of parameter values for paper

%% The following is unrelated to the discretization. It is simply to create a Table containing model parameter values.

% Table of life-cycle model parameters (except the earnings process)
FID = fopen('./SavedOutput/LatexInputs/Table_LifeCycleModelParameters.tex', 'w');
fprintf(FID, 'Parameters of the Life-Cycle Model \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lcc} \n \\hline \\hline \n');
fprintf(FID, ' \\multicolumn{3}{l}{\\textbf{Age structure}} \\\\ \n');
fprintf(FID, ' Total periods & $J$ & %i \\\\ \n', Params.J);
fprintf(FID, ' Retirement period & $Jr$ & %i \\\\ \n', Params.Jr);
fprintf(FID, ' Age in years at j=0 & $Jr$ & %i \\\\ \n', Params.agejshifter);
fprintf(FID, ' \\multicolumn{3}{l}{\\textbf{Preferences}} \\\\ \n');
fprintf(FID, ' Discount factor & $\\beta$ & %8.2f \\\\ \n', Params.beta);
fprintf(FID, ' CRRA & $\\sigma$ & %8.2f \\\\ \n', Params.sigma);
fprintf(FID, ' Curv. of leisure & $\\eta$ & %8.2f \\\\ \n', Params.eta);
fprintf(FID, ' Weight on leisure & $\\psi$ & %8.2f \\\\ \n', Params.psi);
fprintf(FID, ' \\multicolumn{3}{l}{\\textbf{Prices and policies}} \\\\ \n');
fprintf(FID, ' Wage & $w$ & %8.2f \\\\ \n', Params.w);
fprintf(FID, ' Interest rate & $r$ & %8.2f \\\\ \n', Params.r);
fprintf(FID, ' Pension & $pension$ & %8.2f \\\\ \n', Params.pension);
fprintf(FID, ' \\multicolumn{3}{l}{\\textbf{Warm glow of bequests}} \\\\ \n');
fprintf(FID, ' relative importance of bequests & $warmglow1$ & %8.1f \\\\ \n', Params.warmglow1);
fprintf(FID, ' bliss point of bequests & $warmglow2$ & %8.1f \\\\ \n', Params.warmglow2);
fprintf(FID, ' curvature for bequests & $warmglow3$ & %8.1f \\\\ \n', Params.warmglow3);
fprintf(FID, ' \\multicolumn{3}{l}{\\textbf{Permanent Types}} \\\\ \n');
fprintf(FID, ' $\\alpha_i$ & \\multicolumn{2}{l}{[%8.1f,%8.1f,%8.1f,%8.1f,%8.1f]} \\\\ \n', Params.alpha_i(1), Params.alpha_i(2), Params.alpha_i(3), Params.alpha_i(4), Params.alpha_i(5));
fprintf(FID, ' $alphadist$ & \\multicolumn{2}{l}{[%8.1f,%8.1f,%8.1f,%8.1f,%8.1f]} \\\\ \n', Params.alphadist(1), Params.alphadist(2), Params.alphadist(3), Params.alphadist(4), Params.alphadist(5));
fprintf(FID, ' \\multicolumn{3}{l}{\\textbf{Exogenous processes}} \\\\ \n');
fprintf(FID, ' Autocorrelation of $z$ & $\\rho_z$ & %8.1f \\\\ \n', Params.rho_z);
fprintf(FID, ' Std. dev. of innovs. to (log) $z$ & $\\sigma_{\\epsilon,z}$ & %8.2f \\\\ \n', Params.sigma_epsilon_z);
fprintf(FID, ' Std. dev. of (log) $e$ & $\\sigma_{\\epsilon,e}$ & %8.2f \\\\ \n', Params.sigma_epsilon_e);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
% fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
% fprintf(FID, 'Note: notice. \n');
% fprintf(FID, '}} \\end{minipage}');
fclose(FID);


jstr='';
istr='';
cstr='';
for jj=1:10
    jstr=[jstr,'& %8.6f '];
    istr=[istr,'& %i '];
    cstr=[cstr,'c'];
end
jstrend='';
istrend='';
for jj=81:Params.J
    jstrend=[jstrend,'& %8.6f '];
    istrend=[istrend,'& %i '];
end
for jj=Params.J+1:90
    jstrend=[jstrend,'& '];
    istrend=[istrend,'& '];
end

temp_agej=1:1:Params.J;
temp_ageyear=Params.agejshifter+(1:1:Params.J);
temp_dj=Params.dj(Params.agejshifter+1+1:101);
temp_sj=Params.sj;
temp_cumsj=cumprod(Params.sj);

% Table of age conditional death and survival probabilities
FID = fopen('./SavedOutput/LatexInputs/Table_LifeCycleModelSurvivalProbs.tex', 'w');
fprintf(FID, 'Age-Conditional Survival/Death Probabilities for the Life-Cycle Model \\\\ \n');
fprintf(FID, ['\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lc',cstr,'} \n \\hline \\hline \n']);
fprintf(FID, ['Model age & $j$ ',istr,' \\\\ \n'], temp_agej(1:10));
fprintf(FID, ['Age in Years & ',istr,' \\\\ \n'], temp_ageyear(1:10));
fprintf(FID, [' Death & $d_j$ ',jstr,' \\\\ \n'], temp_dj(1:10));
fprintf(FID, [' Survival & $s_j=1-d_j$ ',jstr,' \\\\ \n'], temp_sj(1:10));
fprintf(FID, [' Cumulative survival & $\\Pi_{j=1}^J s_j$ ',jstr,' \\\\ \\hline \n'], temp_cumsj(1:10));
fprintf(FID, ['Model age & $j$ ',istr,' \\\\ \n'], temp_agej(11:20));
fprintf(FID, ['Age in Years & ',istr,' \\\\ \n'], temp_ageyear(11:20));
fprintf(FID, [' Death & $d_j$ ',jstr,' \\\\ \n'], temp_dj(11:20));
fprintf(FID, [' Survival & $s_j=1-d_j$ ',jstr,' \\\\ \n'], temp_sj(11:20));
fprintf(FID, [' Cumulative survival & $\\Pi_{j=1}^J s_j$ ',jstr,' \\\\ \\hline \n'], temp_cumsj(11:20));
fprintf(FID, ['Model age & $j$ ',istr,' \\\\ \n'], temp_agej(21:30));
fprintf(FID, ['Age in Years & ',istr,' \\\\ \n'], temp_ageyear(21:30));
fprintf(FID, [' Death & $d_j$ ',jstr,' \\\\ \n'], temp_dj(21:30));
fprintf(FID, [' Survival & $s_j=1-d_j$ ',jstr,' \\\\ \n'], temp_sj(21:30));
fprintf(FID, [' Cumulative survival & $\\Pi_{j=1}^J s_j$ ',jstr,' \\\\ \\hline \n'], temp_cumsj(21:30));
fprintf(FID, ['Model age & $j$ ',istr,' \\\\ \n'], temp_agej(31:40));
fprintf(FID, ['Age in Years & ',istr,' \\\\ \n'], temp_ageyear(31:40));
fprintf(FID, [' Death & $d_j$ ',jstr,' \\\\ \n'], temp_dj(31:40));
fprintf(FID, [' Survival & $s_j=1-d_j$ ',jstr,' \\\\ \n'], temp_sj(31:40));
fprintf(FID, [' Cumulative survival & $\\Pi_{j=1}^J s_j$ ',jstr,' \\\\ \\hline \n'], temp_cumsj(31:40));
fprintf(FID, ['Model age & $j$ ',istr,' \\\\ \n'], temp_agej(41:50));
fprintf(FID, ['Age in Years & ',istr,' \\\\ \n'], temp_ageyear(41:50));
fprintf(FID, [' Death & $d_j$ ',jstr,' \\\\ \n'], temp_dj(41:50));
fprintf(FID, [' Survival & $s_j=1-d_j$ ',jstr,' \\\\ \n'], temp_sj(41:50));
fprintf(FID, [' Cumulative survival & $\\Pi_{j=1}^J s_j$ ',jstr,' \\\\ \\hline \n'], temp_cumsj(41:50));
fprintf(FID, ['Model age & $j$ ',istr,' \\\\ \n'], temp_agej(51:60));
fprintf(FID, ['Age in Years & ',istr,' \\\\ \n'], temp_ageyear(51:60));
fprintf(FID, [' Death & $d_j$ ',jstr,' \\\\ \n'], temp_dj(51:60));
fprintf(FID, [' Survival & $s_j=1-d_j$ ',jstr,' \\\\ \n'], temp_sj(51:60));
fprintf(FID, [' Cumulative survival & $\\Pi_{j=1}^J s_j$ ',jstr,' \\\\ \\hline \n'], temp_cumsj(51:60));
fprintf(FID, ['Model age & $j$ ',istr,' \\\\ \n'], temp_agej(61:70));
fprintf(FID, ['Age in Years & ',istr,' \\\\ \n'], temp_ageyear(61:70));
fprintf(FID, [' Death & $d_j$ ',jstr,' \\\\ \n'], temp_dj(61:70));
fprintf(FID, [' Survival & $s_j=1-d_j$ ',jstr,' \\\\ \n'], temp_sj(61:70));
fprintf(FID, [' Cumulative survival & $\\Pi_{j=1}^J s_j$ ',jstr,' \\\\ \\hline \n'], temp_cumsj(61:70));
fprintf(FID, ['Model age & $j$ ',istr,' \\\\ \n'], temp_agej(71:80));
fprintf(FID, ['Age in Years & ',istr,' \\\\ \n'], temp_ageyear(71:80));
fprintf(FID, [' Death & $d_j$ ',jstr,' \\\\ \n'], temp_dj(71:80));
fprintf(FID, [' Survival & $s_j=1-d_j$ ',jstr,' \\\\ \n'], temp_sj(71:80));
fprintf(FID, [' Cumulative survival & $\\Pi_{j=1}^J s_j$ ',jstr,' \\\\ \\hline \n'], temp_cumsj(71:80));
fprintf(FID, ['Model age & $j$ ',istrend,' \\\\ \n'], temp_agej(81:end));
fprintf(FID, ['Age in Years & ',istrend,' \\\\ \n'], temp_ageyear(81:end));
fprintf(FID, [' Death & $d_j$ ',jstrend,' \\\\ \n'], temp_dj(81:end));
fprintf(FID, [' Survival & $s_j=1-d_j$ ',jstrend,' \\\\ \n'], temp_sj(81:end));
fprintf(FID, [' Cumulative survival & $\\Pi_{j=1}^J s_j$ ',jstrend,' \\\\ \\hline \n'], temp_cumsj(81:end));
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Source: $s_j$ is the probability of surviving to be age $j+1$, given alive at age $j$. Age-conditional death probabilites for the US are taken from "National Vital Statistics Report, volume 58, number 10, March 2010.", first column (qx) of Table 1 (Total Population). Conditional survival probabilities are calculated as then $s_j=1-d_j$. \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

temp_kappaj=Params.kappa_j;
% Table of deterministic labor productivity units
FID = fopen('./SavedOutput/LatexInputs/Table_LifeCycleModelAgePrody.tex', 'w');
fprintf(FID, 'Deterministic Age-dependent labor productivity units for the Life-Cycle Model \\\\ \n');
fprintf(FID, ['\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lc',cstr,'} \n \\hline \\hline \n']);
fprintf(FID, ['Model age & $j$ ',istr,' \\\\ \n'], temp_agej(1:10));
fprintf(FID, ['Age in Years & ',istr,' \\\\ \n'], temp_ageyear(1:10));
fprintf(FID, [' $\\kappa_j$ & ',jstr,' \\\\ \n'], temp_kappaj(1:10));
fprintf(FID, ['Model age & $j$ ',istr,' \\\\ \n'], temp_agej(11:20));
fprintf(FID, ['Age in Years & ',istr,' \\\\ \n'], temp_ageyear(11:20));
fprintf(FID, [' $\\kappa_j$ & ',jstr,' \\\\ \n'], temp_kappaj(11:20));
fprintf(FID, ['Model age & $j$ ',istr,' \\\\ \n'], temp_agej(21:30));
fprintf(FID, ['Age in Years & ',istr,' \\\\ \n'], temp_ageyear(21:30));
fprintf(FID, [' $\\kappa_j$ & ',jstr,' \\\\ \n'], temp_kappaj(21:30));
fprintf(FID, ['Model age & $j$ ',istr,' \\\\ \n'], temp_agej(31:40));
fprintf(FID, ['Age in Years & ',istr,' \\\\ \n'], temp_ageyear(31:40));
fprintf(FID, [' $\\kappa_j$ & ',jstr,' \\\\ \n'], temp_kappaj(31:40));
fprintf(FID, ['Model age & $j$ ',istr,' \\\\ \n'], temp_agej(41:50));
fprintf(FID, ['Age in Years & ',istr,' \\\\ \n'], temp_ageyear(41:50));
fprintf(FID, [' $\\kappa_j$ & ',jstr,' \\\\ \n'], temp_kappaj(41:50));
fprintf(FID, ['Model age & $j$ ',istr,' \\\\ \n'], temp_agej(51:60));
fprintf(FID, ['Age in Years & ',istr,' \\\\ \n'], temp_ageyear(51:60));
fprintf(FID, [' $\\kappa_j$ & ',jstr,' \\\\ \n'], temp_kappaj(51:60));
fprintf(FID, ['Model age & $j$ ',istr,' \\\\ \n'], temp_agej(61:70));
fprintf(FID, ['Age in Years & ',istr,' \\\\ \n'], temp_ageyear(61:70));
fprintf(FID, [' $\\kappa_j$ & ',jstr,' \\\\ \n'], temp_kappaj(61:70));
fprintf(FID, ['Model age & $j$ ',istr,' \\\\ \n'], temp_agej(71:80));
fprintf(FID, ['Age in Years & ',istr,' \\\\ \n'], temp_ageyear(71:80));
fprintf(FID, [' $\\kappa_j$ & ',jstr,' \\\\ \n'], temp_kappaj(71:80));
fprintf(FID, ['Model age & $j$ ',istrend,' \\\\ \n'], temp_agej(81:end));
fprintf(FID, ['Age in Years & ',istrend,' \\\\ \n'], temp_ageyear(81:end));
fprintf(FID, [' $\\kappa_j$ & ',jstrend,' \\\\ \n'], temp_kappaj(81:end));
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Source: $\\kappa_j$ is the labor productivity units as a deterministic function age. \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);



