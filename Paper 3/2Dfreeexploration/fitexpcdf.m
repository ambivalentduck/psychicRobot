function [cost,fit]=fitexpcdf(lambda)

global fitmecounts fitmebins

ecdf=expcdf(fitmebins,lambda);
fit=ecdf(2:end)-ecdf(1:end-1);
cost=sum((fitmecounts-fit).^2);