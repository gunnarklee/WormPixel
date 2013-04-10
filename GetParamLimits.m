% get parameter limits from a list of a property
function [upper, lower, stats] = GetParamLimits(DataList, PRCbracket)
stats.mean=mean(DataList);
stats.range=mean(DataList);
stats.upper=stats.mean+stats.mean*PRCbracket;upper=stats.upper;
stats.lower=stats.mean-stats.mean*PRCbracket;lower=stats.lower;
stats.std=std(DataList);

end