clc;
addpath('./mcm');
clear variables;
%close all;
seed = 10;
s = RandStream('mt19937ar','Seed',seed);
rand('seed',seed);
randn('seed',seed);    
numNodes = 50;
numLayers = 1;
carSpeed = 50 * 1000/60; 
timeSteps = 10000;
seed = ceil(rand()*100);
tau = 1;
%generate a directed weighted network
[A1] = triu(generateRandomGraphGivenEdgeProb(numNodes,0.1));
[A2] = tril(generateRandomGraphGivenEdgeProb(numNodes,0.1));
A = sparse(A1 + A2);    
A(find(A(:)==1)) = rand(length(find(A(:)==1)),1)*2000;
ATime = A;
nonZeroValues = find(ATime(:)>0);
ATime(nonZeroValues) = ATime(nonZeroValues)*(1/carSpeed);      
ATime = ATime';
A = ATime;    
%Structures that define the routing policy
[SPModelData.spBW,...
        SPModelData.numExternalPacketRoutings,...
        SPModelData.numLocalPacketRoutings,...
        SPModelData.numStartsAtNode,...
        SPModelData.numEndsAtNode]=...
        SPEdgeNodeBetweennessC_BrandesWeightedTestEdgeBW_local_fheap(A,numNodes,1);          
SPModelData.spBW = sum(SPModelData.spBW,2);        
%compute critical injection rate theorical
%only accurate if the genration rate is the same for all node pairs
rho_cs = tau*((numNodes-1)./(SPModelData.spBW + 2*(numNodes-1)));
rho_c = min(rho_cs);    
[SPModelData.DistRG,SPModelData.NumPathsRG,SPModelData.PredRG,SPModelData.PredProbRG,...
     SPModelData.srcMultRG,SPModelData.dstMultRG,SPModelData.probMultRG,...           
     ] = SPDataPathDegeneration_staticMemory(A,numNodes,1);  
%possGenRates = linspace(0.01,3,100);
possGenRates = logspace(-2,log10(3),100);

congestedNodes = [];
h = waitbar(0,'Initializing waitbar...');    
for rhosIndex = 1:length(possGenRates)
    waitbar(rhosIndex/length(possGenRates),h,'Halfway there...')
    genRatePerMinuteGates = zeros(numNodes,1);
    genRatePerMinuteGates(:) = possGenRates(rhosIndex);            
    %genRatePerMinuteGates(1) = possGenRates(rhosIndex);   
    %genRatePerMinuteGates(2) = possGenRates(rhosIndex);            
    %genRatePerMinuteGates(1:3) = possGenRates(rhosIndex);
    %define the "OD" matrix somehow, this is all to all but with the
    %given generationrate
    [SPModelData.possDest,SPModelData.possDestCumProb,SPModelData.possDestProb] =...
        matrixOfAccesibleDestinations_AllToAll(...
            find(genRatePerMinuteGates>0),...
            find(genRatePerMinuteGates==0),...
            genRatePerMinuteGates,...
            SPModelData.DistRG);                
    genRatePerMinute = genRatePerMinuteGates;        
    processingRatePerMinute = ones(numNodes,1);                        
    % monte carlo simulations
    [numPackets,numPacketOutUnitTime,B_i_ext,e_i_ext,B_i,e_i,s_i,PExp,sigma,DExp,...
    extractedPackets,packetsThatContinue]...
        = cSPCongestion_statMem_dir_weighted_local(...
            genRatePerMinute,...
            processingRatePerMinute,...
            [],...%SPModelData.DistRG,...
            [],...%SPModelData.NumPathsRG,...
            SPModelData.PredRG,...
            SPModelData.PredProbRG,...
            SPModelData.srcMultRG,...
            SPModelData.dstMultRG,...
            SPModelData.probMultRG,...
            SPModelData.possDest,...
            SPModelData.possDestCumProb,...
            numNodes,1,timeSteps,seed); 
    % compute microscopic vars obtained by the montecarlo
    normFactor = 1./sum(PExp);
    normFactor(normFactor == Inf) = 0;
    PExp = PExp*diag(normFactor);    
    pExp = packetsThatContinue./extractedPackets;
    DExp = mean(DExp,2);
    % compute queue increments obtained by the montecarlo
    for j=1:numNodes
        p = polyfit(1:size(numPackets(j,:),2),numPackets(j,:),1);
        DeltaNExp(j,rhosIndex) = (p(1));                
    end                    
    % compute theoretical predictions for the microscopic and
    % macroscopic vars
    [DeltaNTeor(:,rhosIndex),rateNumPacketsDelivered,...
    B_i_obs, e_i_obs, EB_Ext_obs, EB_Int_obs,...
    DTeor, pTeor, PTeor,congestedNodes] = SP_computeEtaGivenRho_dir_weigh_nonHomStart(...
            A,numNodes,numLayers,...
            genRatePerMinute,processingRatePerMinute,...
            SPModelData.possDest,SPModelData.possDestProb,SPModelData,congestedNodes);  
    packetsGenPerMinute(rhosIndex) = sum(genRatePerMinute);
end
close(h);
%% Plots for the eta
figure ('Position', [50 50 650 500]);
semilogx(possGenRates,sum(DeltaNTeor) ./ (packetsGenPerMinute),'b-','LineWidth',2); hold on;
semilogx(possGenRates,sum(DeltaNExp) ./ (packetsGenPerMinute),'ro','LineWidth',2); hold on;
scatter(rho_c,0,100,'MarkerFaceColor', 'g')

xlabel('Generation rate $\rho$','interpreter','latex');
ylabel('Order parameter $\eta$','interpreter','latex');
legend('Theoretical prediction',...
       'Monte carlo simulation',...
       'Critical injection rate $\rho_c$',...
       'Location','northwest','interpreter','latex');

axis([min(possGenRates) max(possGenRates) 0 1]);
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis
scale = 1.2;
makePlotNice;
