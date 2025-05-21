function xyb_est = SimulateReceiver(xyb,M,N,receiver,varargin)
% Description:
% An interface for choosing which receiver to run.

% valid receiver names
validReceivers = {'DirectImaging','StaticSPADE'};

% input validators
validateScene = @(xyb) isnumeric(xyb) & all(xyb(:,3)>0) & (abs(sum(xyb(:,3))-1)<1e-10);
validatePositiveInteger = @(x) isnumeric(x) & (rem(x,1)==0) & (x>0);
validatePositiveReal = @(x) isnumeric(x) & (x>0);
validateSplittingRatio = @(x) isnumeric(x) & (0<x) & (x<=1);
validateReceiver = @(x) any(validatestring(x,validReceivers));
validateHGIndices = @(nm) all(nm(:) == int8(nm(:))) & all(nm(:)>=0);

% default inputs
default_sigma = 1;
default_visualize_flag = 1;

% input parser
p = inputParser;
addRequired(p,'xyb',validateScene)
addRequired(p,'M',validatePositiveInteger)
addRequired(p,'N',validatePositiveInteger)
addRequired(p,'receiver',validateReceiver)
addOptional(p,'sigma',default_sigma,validatePositiveReal)
addOptional(p,'visualize_flag',default_visualize_flag)

switch receiver
    case 'DirectImaging'
        parse(p,xyb,M,N,receiver,varargin{:})
        xyb_est = Simulate_DirectImaging_Receiver(xyb,M,N,p.Results.sigma,p.Results.visualize_flag);

    case 'StaticSPADE'
        addRequired(p,'splitting_ratio',validateSplittingRatio)
        addRequired(p,'nm',validateHGIndices)
        parse(p,xyb,M,N,receiver,varargin{:})
        xyb_est = Simulate_StaticSPADE_Receiver(xyb,M,N,p.Results.sigma,p.Results.splitting_ratio,p.Results.nm,p.Results.visualize_flag);
end
end
