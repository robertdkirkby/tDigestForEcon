%% These are just two snippets of code that I used to decide two of my minor implementation
% details for the createDigest() and mergeDigest() commands.
% The first was about qlimit getting stupidly close to one.
% The second was about how big should the preallocated matrices be for things like C.
% Note that a third, checking that the new qlimit is creater than q, when
% updating qlimit, came from resolving a problem with an actual (simulated) dataset.

%% Based on this code I decided to implement the if statement around 1-qlimit<10^(-7)
% Otherwise qlimit just goes ridiculously close to 1

qlimit=0.999999

for ii=1:10
    qlimit=kinvfn(kfn(qlimit,delta)+1,delta);
    if qlimit<0.5
        qlimit
    else
        1-qlimit
    end
end

qlimit=0

% I should stop once 1-qlimit<10^(-7)



%% Try get not Nq=4999
% The remaining part of this code was used to determine the likely
% maximum number of quantiles created by a t-Digest to allow the
% createDigest and mergeDigest functions to preallocated memory.


for N=[10^3,10^7] % I toyed with different values for N on grounds that might matter for how many quantiles end up recorded
    % Distribution is a collection of grid points and their associated weights
    evalgrid=5*rand(N,1);
    weights=rand(length(evalgrid),1);
    weights=weights/sum(weights); % Normalize mass to one

    [sortevalgrid,sortIndex]=sort(evalgrid);
    sortweights=weights(sortIndex);

    [C,digestweights,qlimitvec]=createDigest(evalgrid, weights,delta);

    size(C)
    size(digestweights)
    size(qlimitvec)

end

% Note: It seemed like I cannot exceed 4999 (with delta=10000), but it
% turned out in practice that actually for very specific underlying
% datasets it is possible. So in the end I use 5100, and I know the last
% few entries will be left as zeros so then I just trim them.
