function makeTimeSeriesHistograms(timeFrame,W,threshold,times,binSize,folderName)
% Plots tweet count vs time histograms for the time series for the entire 
% data set with timestamps represented in nx1 vector times, and for 
% each topic represented by columns nxk matrix W, an output of NMF.
%
% To determine which tweets belong to which topic normalizes rows of 
% W to create a distribution of topics on tweets. Tweets with values for
% topics greater than specified topicTreshold are included in that topic's 
% time series.
%
% Inputs:
% timeFrame      = 1x2 vector of date-time objects specifying date range 
%                  to consider [minTime,maxTime]
% W              = nxk representation of topic distribution on tweets
% threshold      = the desired threshold for a Tweet to be part of a topic
% times          = column vector of datetime tweet timestamps with indices 
%                  corresponding to rows of W (i.e. tweets)
% binSize        = number of days worth of tweets to include in a bin
% folderName     = folder in which to save figures
%
timeFrame = datenum(timeFrame);
times = datenum(times);
binSize = datenum(0,0,binSize,0,0,0);

[Number_Documents,~] = size(W);  
sums = sum(W,2);
invSums = sums.^(-1);
invSumsMat = sparse(1:Number_Documents,1:Number_Documents,invSums);
normedW = invSumsMat*W;

newFolderName = sprintf('timeHistograms_binsize%ddays_topicThreshold%.2f',binSize,topicThreshold);
mkdir(folderName,newFolderName)

meetMin = find(times>timeFrame(1));
meetMax = find(times<timeFrame(2));
allTweetsTimesRange = times(meetMin(1):meetMax(end));
numBins = ceil((allTweetsTimesRange(end)-allTweetsTimesRange(1))/binSize);
[allTweetsBins,binEdges] = histcounts(allTweetsTimesRange,numBins);
binCenters = binEdges(2:end) - (binSize/2) * ones(size(binEdges(2:end)));

figure(1);
plot(binCenters,allTweetsBins)
title('All Tweets')
minTime = times(1);
maxTime = times(length(times));
xDivision = floor(maxTime-minTime)/6;
xVals = minTime:xDivision:maxTime;
xLabels = datestr(xVals,'mm/dd/yy');
set(gca, 'XTick', xVals);
set(gca, 'XTickLabel', xLabels);
xlabel('Time')
ylabel('Number of Tweets')
savefig(sprintf('%s/%s/All_Tweets',folderName,newFolderName))
saveas(gcf,sprintf('%s/%s/All_Tweets',folderName,newFolderName),'tiff')

for ii = 1:numTopics
    topicTweetIndices = find(normedW(:,ii)>threshold);
    topicTimes = times(topicTweetIndices);
    meetMin = find(topicTimes>timeFrame(1));
    meetMax = find(topicTimes<timeFrame(2));
    topicTimesRange = topicTimes(meetMin(1):meetMax(end));
    topicBins = histcounts(topicTimesRange,binEdges);
    weightedTopicBins = topicBins./allTweetsBins;
    
    figure(ii+1)
    plot(binCenters,weightedTopicBins)
    title(sprintf('Topic %d',ii))
    minTime = topicTimes(1);
    maxTime = topicTimes(length(topicTimes));
    xDivision = floor(maxTime-minTime)/6;
    xVals = minTime:xDivision:maxTime;
    xLabels = datestr(xVals,'mm/dd/yy');
    set(gca, 'XTick', xVals);
    set(gca, 'XTickLabel', xLabels);
    xlabel(sprintf('Time (Bin = %d Days)',binSize))
    ylabel('Percentage of Total Tweets')
    savefig(sprintf('%s/%s/Topic_%003d',folderName,newFolderName,ii))
    saveas(gcf,sprintf('%s/%s/Topic_%003d',folderName,newFolderName,ii),'tiff')
end

end




