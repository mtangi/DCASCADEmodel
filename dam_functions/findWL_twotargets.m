function WL_Target = findWL_twotargets(ORparameters_dam, date)
%findSLtarger_MinMax find the target water level of a reservoir managed
%with an operating rule rule curve depends on four parameters, namely the
%minimum  and maximum water levels that a reservoir should reach within a
%year and  the time at which the two levels should be reached.
%Between these dates, the target water level is given by linearly
%interpolating the two targets according to the current date

%%

target_WL = ORparameters_dam.WL_twotargets(:,1);
target_dates = ORparameters_dam.WL_twotargets(:,2:3); %target month in the first column, day in the second   

%the first date should be earlier in the year than the
%second. If it's not, invert the dates and targets
if days(datetime(date(1),target_dates(2,1),target_dates(2,2)) - datetime(date(1),target_dates(1,1),target_dates(1,2)))<0
    target_dates = flip(target_dates);
    target_WL = flip(target_WL);
end

%number of days between the two MinMax dates
dist_1to2 = days(datetime(date(1),target_dates(2,1),target_dates(2,2)) - datetime(date(1),target_dates(1,1),target_dates(1,2)));
dist_2to1 = days(datetime(date(1)+1,target_dates(1,1),target_dates(1,2)) - datetime(date(1),target_dates(2,1),target_dates(2,2)));

%distance in days between the current timestep and the first target (given current year)
days_date(1) = days(datetime(date(1),date(2),date(3)) - datetime(date(1),target_dates(1,1),target_dates(1,2)));
%distance in days between the current timestep and the second target (given current year)
days_date(2) = days(datetime(date(1),date(2),date(3)) - datetime(date(1),target_dates(2,1),target_dates(2,2)));

diff = [target_WL(2)-target_WL(1);target_WL(1)-target_WL(2)]; %support matrix: difference between the SL of the two targets
if sign(days_date(1)) == sign(days_date(2)) % if the distances have the same sign, date is between target 2 and target 1 of the next year
    [~,I] = min(abs([days_date(1),days_date(2)]));
    WL_Target = target_WL(I) + diff(I)*abs(days_date(I))/dist_2to1;
else % if the distances have different signa, date is between target 1 and target 2 of the same year
    [~,I] = min(abs([days_date(1),days_date(2)])); 
    WL_Target = target_WL(I) + diff(I)*abs(days_date(I))/dist_1to2;
end

end

