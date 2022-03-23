function [dates, discounts] = bootstrap(datesSet, ratesSet)
% bootstrap for a given set of dates and quoted rates
%
% INPUT
% datesSet: struct where datesSet.settlement -> settlement date
%                        datesSet.depos      -> expiry dates for depos                      
%                        datesSet.futures    -> [settlement dates of futures, expiry dates of futures] 
%                        datesSet.swaps      -> expiry dates for swaps                               
%
% ratesSet: struct where ratesSet.depos   -> [BID, ASK] for depos
%                        ratesSet.futures -> [BID, ASK] for futures
%                        ratesSet.swaps   -> [BID, ASK] for swaps
%                     
% all the rates are already divided by 100 (converted from percentage)                           
% futures BID-ASK are already transformed with (1-BID/ASK)
%
% OUTPUT
% dates:        expiry dates of the most liquid contracts 
%               settlement date                    (first date 19-feb2008 )
%               depos until the first month        (until 19-march2008)
%               expiries of first seven futures    (until 16-dec2009)                                           
%               expiries of the swaps              (until 19-feb2058)                                          
% discounts:    rates computed via bootstrap at above dates


%% Day conventions
DepoDayCount = 2;  %yearfrac Act/360
IBDayCount   = 3;  %yearfrac Act/365
SwapDayCount = 6;  %yearfrac 30/360 European



%% depos

depo_rates  = 0.5*(ratesSet.depos(1:3,1)+ ratesSet.depos(1:3,2)); %needed depo rates
sett_date = datesSet.settlement;    %settlement date
% depo_dates = datesSet.depos(:,1);
fractions_of_year = yearfrac(sett_date, datesSet.depos(1:3,1), DepoDayCount);
discounts_depos   = 1./(1 + fractions_of_year .* depo_rates);
zero_rates_depos  = -log(discounts_depos)./yearfrac(sett_date, datesSet.depos(1:3,1), IBDayCount);

%% futures

futures_dates = datesSet.futures(1:7,:);   %[setdate, expdate]
futures_rates = 0.5*(ratesSet.futures(1:7,1)+ ratesSet.futures(1:7,2));  %rates for futures
fwd_discounts = (1+yearfrac(datesSet.futures(1:7,1), datesSet.futures(1:7,2), DepoDayCount).*futures_rates).^(-1);  %forward discounts

%first future discount and first two zero rates 'by hands'
zero_rates           =  zeros(14,1);
zero_rates(1)        = -log(discounts_depos(end,1))/yearfrac(sett_date, datesSet.depos(3,1), IBDayCount); %zero rate at first setdate
futures_discounts    = zeros(length(futures_dates),1);
futures_discounts(1) = discounts_depos(end, 1)*fwd_discounts(1);
zero_rates(2)        = -log(futures_discounts(1))/yearfrac(sett_date, datesSet.futures(1,2), IBDayCount); %zero rate at first expdate

futures_dates_new = zeros(14,1);
for i=1:length(futures_dates)-1
    futures_dates_new(2*i-1) = futures_dates(i, 1);
    futures_dates_new(2*i)   = futures_dates(i, 2);
    [~, ind1] = unique(zero_rates(1:2*i), 'stable');
    [~, ind2] = unique(futures_dates_new(1:2*i), 'stable');
    %extrap_zero_rate = zero_rates(2*i-1) + ((futures_dates(i+1,1)-futures_dates(i,1))/(futures_dates(i,2)-futures_dates(i,1)))...
    %   *(zero_rates(2*i)-zero_rates(2*i-1));
    extrap_zero_rate  = interp1(futures_dates_new(ind1), zero_rates(ind2), futures_dates(i+1, 1), 'linear', 'extrap');
    zero_rates(2*i+1) = extrap_zero_rate;
    disc_factor = exp(-extrap_zero_rate*yearfrac(sett_date, futures_dates(i+1,1), IBDayCount));
    futures_discounts(i+1) = disc_factor * fwd_discounts(i+1);
    zero_rates(2*i+2) = -log(futures_discounts(i+1))/yearfrac(sett_date, futures_dates(i+1,2), IBDayCount);
end
futures_dates_new(end-1) = futures_dates(end,1);
futures_dates_new(end) = futures_dates(end,2);


% just a first check
% discounts = [discounts_depos; futures_discounts];
% dates=[datesSet.depos(1:3); futures_dates(:,2)];
% plot(dates, discounts, '-o')
% ax = gca;
% ax.XTick = dates;
% datetick('x', 1, 'keepticks')
% figure
% tot_zero_rates = [zero_rates_depos; zero_rates(2:2:end)];
% plot(dates, tot_zero_rates, '-')
% ax = gca;
% ax.XTick = dates;
% datetick('x', 1, 'keepticks')

%% swaps
% interpolation for the first swap discount (not returned)
swaps_dates      = datesSet.swaps(1:end);
swaps_zero_rates = zeros(length(swaps_dates),1);
swaps_discounts  = zeros(length(swaps_dates),1);
[~, ind1] = unique(zero_rates, 'stable');
[~, ind2] = unique(futures_dates_new, 'stable');
zero_swap_interp    = interp1(futures_dates_new(ind2), zero_rates(ind1), swaps_dates(1));
swaps_zero_rates(1) = zero_swap_interp;
swaps_discounts(1)  = exp(-swaps_zero_rates(1)*yearfrac(sett_date, swaps_dates(1), IBDayCount));

S_IR = 0.5*(ratesSet.swaps(:,1) + ratesSet.swaps(:,2)); 
% swaps_dates1 = swaps_dates(2:end)-swaps_dates(1:end-1);
% fractions_of_year = zeros(length(swaps_dates),1);
fractions_of_year = yearfrac(swaps_dates(1:end-1), swaps_dates(2:end), SwapDayCount);
for i=2:length(swaps_dates)
    swaps_discounts(i) = (1 - S_IR(i)*swaps_discounts(1:i-1)'*fractions_of_year(1:i-1))/(1+fractions_of_year(i-1)*S_IR(i));
end

%% outputs
discounts = [1; discounts_depos; futures_discounts; swaps_discounts(2:end)];
dates = [datesSet.settlement; datesSet.depos(1:3); futures_dates(:,2); swaps_dates(2:end)];



end