% File: RetroLowess_Time_Seres_Analysis  November 2, 2019
% calculates smoothed data based on
% a weighted regression of past data (an asymmetrical lowess)
% copes with NaNs in the data set
% analyzes the distributions of:
%       duration of extreme events
%       intervals between events
%       extreme temperatures
%       deviations from "expected" temperatures
%       degree-days above threshold
% 

clear 

percent = 0.95; % cut off for "extreme" (as a fraction; e.g., 0.95)
pcnt = num2str(percent);

% load the data file

filename = input('Data File Name (with suffix) ','s'); % the time-series data
fname = input('Data File Name (without suffix) ','s'); % for figure labels
data = load(filename);
max = size(data,1);

interval = input('Data interval (days)'); % e.g., 1/24 for hourly data
interval_count = 0;

for d1 = 4:2:25 % 365 % the increment in which to loop the smoothing interval(days)
    
    interval_count = interval_count + 1
    
    interval_xaxis(interval_count) = d1; % for figure labels
    
    d = floor(d1/interval); % the first point where lowess can be calculated
    smth = num2str(d1);

    smoothedy = zeros(max,1);
    delta = zeros(max,1);
    
    for i = 1:d-1 % take care of initial points where delta is not calculated
        delta(i) = -9999;
    end
    
    x = zeros(max,1);
    y = zeros(max,1);
    rawdata = zeros(max,1);
    dev = zeros(max,1);
    flag = zeros(max,1);

    goodrawcount = 0;
    
    for i= 1:max % get rid of the NaNs by substituting a large negative number

        x(i) = i; % a vector of time, in # of measurement intervals

        if data(i) > -50 
            goodrawcount = goodrawcount + 1;
            y(i) = data(i); 
        else
            y(i) = -9999;
            flag(i) = 1; % flagged as bad data
        end

    end

    goodrawdata = NaN(goodrawcount,1); % avoids zeros at high numbers

    goodrawcount = 0;
    sum = 0;
    sum2 = 0;
    sum3 = 0;
    for i = 1:max
        if y(i) > -9999 
            goodrawcount = goodrawcount + 1;
            goodrawdata(goodrawcount) = y(i); % just the reliable data
            sum = sum + y(i);
            sum2 = sum2 + y(i)^2;
            sum3 = sum3 + y(i)^3;
        end
    end

    if interval_count == 1
        Mean_Raw_Data = sum/goodrawcount
        var_rawdata = (sum2 - sum^2/goodrawcount)/(goodrawcount - 1);
        SD_Raw_Data = sqrt(var_rawdata) % the standard deviation
        Gamma_Raw_Data = sum3/(goodrawcount * SD_Raw_Data^3) % the skew
    end

    % -----------------------------------------------------------------------
    % smooth the data using the retrolowess algorithm
    % uses vector y because it is the time series

    for i = d:max % start at the first point where a full retro is possible

        start =i-d+1; % where to start the running average

        sumx = 0;
        sumy = 0;
        sumx2 = 0;
        sumxy = 0;
        z = -1;
        count = 0; % keep track of how many non NaNs in the calculation

        for k = d-1:-1:0 % calculate the weighted data

            z = z+1;

            w = (1-(abs(k/d))^3)^3; % the lowess weighting function

            if y(start+z) > -9999 % only include reliable data in the calc.              
                count = count + 1;
                sumx = sumx + w*x(start+z); % weight the x values
                sumy = sumy + y(start+z); % don't weight the y's
                sumx2 = sumx2 + w*x(start+z)^2; % weighted sum of squares
                sumxy = sumxy + w*x(start+z)*y(start+z); % weighted xy-prod
            end

        end

        if count > 1 

            meanx = sumx/(count);
            meany = sumy/(count);
            ssx = sumx2 - sumx^2/(count); % machine formula for x sum of squares
            scp = sumxy - (sumx*sumy)/(count); % machine form. for sum of x-product

            b = scp/ssx; % best estimate of the slope

            a = meany -b*meanx; % best estimate of the y intercept

            smoothedy(i) = a + b*x(i); % smoothed estimate

        else

            smoothedy(i) = -9999; % bad data points

        end


    end

    % -----------------------------------------------------------------------
    % calculate the deviations from the smoothed data

    badcount = 0; % number of bad data points
    goodcount = 0; %number of good data points

    sum_deviation = 0;
    sum2_deviation = 0;
    sum3_deviation = 0;
    for i = d:max % subtract smoothed from raw

        if smoothedy(i) > -9999 & y(i) > -9999
            goodcount = goodcount + 1;
            delta(i) = y(i) - smoothedy(i); % a vector of the deviations
            sum_deviation = sum_deviation + delta(i);
            sum2_deviation = sum2_deviation + delta(i)^2;
            sum3_deviation = sum3_deviation + delta(i)^3;
            rawdata(i) = y(i); % the raw data corresponding to smoothed points
        else
            badcount = badcount + 1;
            delta(i) = -9999;
        end

    end

    Stats_Dev(interval_count,1) = sum_deviation/goodcount;
    var_deviation = (sum2_deviation - sum_deviation^2/goodcount)/(goodcount - 1);
    Stats_Dev(interval_count,2) = sqrt(var_deviation);
    Stats_Dev(interval_count,3) = sum3_deviation/(goodcount * Stats_Dev(interval_count,2)^3);

    good_ratio = goodcount/max; % the fraction of data that we use in calcs

    ranked = sort(delta); % an ascending sort of deviations

    rawranked = sort(rawdata); % an ascending sort of rawdata

    q = floor(percent * goodcount); % cut off for "extreme," by percentage

    coff = ranked(q + badcount); % the cut-off deviation defining "extreme"
    
    Stats_Raw(interval_count,1) = coff;

    rawcoff = rawranked(q + badcount); % cut-off for the rawdata
    
    Stats_Raw(interval_count,2) = rawcoff;

    for i = 1:max % set everything but the extreme deviations to -9999
        if delta(i) >= coff
            delta(i) = delta(i);
        else
            delta(i) = -9999;
        end
    end

    dev = NaN(goodcount,1); % avoids zeros at high numbers

    deviationcount = 0;
    sum_extdev = 0;
    sum2_extdev = 0;
    sum3_extdev = 0;
    for i = 1:max % vector of extreme deviations for histogramming
        if delta(i) >= coff
            deviationcount = deviationcount + 1;
            dev(deviationcount) = delta(i);
            sum_extdev = sum_extdev + delta(i);
            sum2_extdev = sum2_extdev + delta(i)^2;
            sum3_extdev = sum3_extdev + delta(i)^3;
        end
    end

    Stats_Ex_Dev(interval_count,1) = sum_extdev/deviationcount;
    var_extdev = (sum2_extdev - sum_extdev^2/deviationcount)/(deviationcount - 1);
    Stats_Ex_Dev(interval_count,2) = sqrt(var_extdev);
    Stats_Ex_Dev(interval_count,3) = sum3_extdev/(deviationcount * Stats_Ex_Dev(interval_count,2)^3);

    crawdata= NaN(max-q-badcount,1);
    rawcount = 0;
    sum_craw = 0;
    sum2_craw = 0;
    sum3_craw = 0;
    for i = 1:max % vector of extreme rawdata for histogramming
        if delta(i) >= coff
            rawcount = rawcount + 1;
            crawdata(rawcount) = rawdata(i);
            sum_craw = sum_craw + rawdata(i);
            sum2_craw = sum2_craw + rawdata(i)^2;
            sum3_craw = sum3_craw + rawdata(i)^3;
        end
    end

    Stats_Craw(interval_count,1) = sum_craw/rawcount;
    var_craw = (sum2_craw - sum_craw^2/rawcount)/(rawcount - 1);
    Stats_Craw(interval_count,2) = sqrt(var_craw);
    Stats_Craw(interval_count,3) = sum3_craw/(rawcount * Stats_Craw(interval_count,2)^3);

   
    % -----------------------------------------------------------------------
    % calculate the duration of extreme events

    newcount = 0;
    duratn = 0;
    degree_intervals = 0;
    sum_dur = 0;
    sum2_dur = 0;
    sum3_dur = 0;
    sum_intensity = 0;
    sum2_intensity = 0;
    sum3_intensity = 0;
    distribution = zeros(max,1);
    intensity = zeros(max,1);
    
    for i = 1:max-1 % detect and record intervals of extreme conditions

        if delta(i) == -9999 & delta(i+1) > -9999 & flag(i) == 0  
            duratn = 1;
            degree_intervals = delta(i+1)* interval/2; % note the half interval
        end

        if delta(i) > -9999 & delta(i+1) > -9999 
            duratn = duratn + 1;
            degree_intervals = degree_intervals + delta(i+1) * interval;
        end

        if delta(i) > -9999 & delta(i+1) == -9999
            newcount = newcount + 1;
            distribution(newcount) = duratn * interval; % duration in days
            intensity(newcount) = degree_intervals + delta(i) * interval/2; % degree days
        end

    end    

    Stats_Raw(interval_count,3) = newcount;
    
    for i = 1:newcount % because the last interval is ongoing
        sum_dur = sum_dur + distribution(i);
        sum2_dur = sum2_dur + distribution(i)^2;
        sum3_dur = sum3_dur + distribution(i)^3;
        sum_intensity = sum_intensity + intensity(i);
        sum2_intensity = sum2_intensity + intensity(i)^2;
        sum3_intensity = sum3_intensity + intensity(i)^3;
    end
    
    Stats_Dur(interval_count,1) = sum_dur/newcount;
    var_Dur = (sum2_dur - sum_dur^2/newcount)/(newcount - 1);
    Stats_Dur(interval_count,2) = sqrt(var_Dur);
    Stats_Dur(interval_count,3) = sum3_dur/(newcount * Stats_Dur(interval_count,2)^3);

    Stats_Intensity(interval_count,1) = sum_intensity/newcount;
    var_Intensity = (sum2_craw - sum_intensity^2/newcount)/(newcount - 1);
    Stats_Intensity(interval_count,2) = sqrt(var_Intensity);
    Stats_Intensity(interval_count,3) = sum3_intensity/(newcount * Stats_Intensity(interval_count,2)^3);

  
    % ---------------------------------------------------------------------
    % Calculate the distribution of inter-event intervals

    intercount = 0;
    interduratn = 0;
    sum_inter = 0;
    sum2_inter = 0;
    sum3_inter = 0;
    
    for i = 1:max-1 % detect and record intervals

        if delta(i) > -9999 & flag(i) == 0 & delta(i+1) == -9999 & flag(i+1) == 0    
            interduratn = 1/2;
        end

        if delta(i) == -9999 & flag(i) == 0 & delta(i+1) == -9999 & flag(i+1) == 0
            interduratn  = interduratn  + 1;
        end

        if delta(i) == -9999 & flag(i) == 0 & delta(i+1) > -9999
            intercount = intercount + 1;
            if intercount > 0
                interdist(intercount) = (interduratn + 1/2) * interval; % duration in days
            end
        end

    end    

    for i = 1:intercount-1 % because the last interval is ongoing
        sum_inter = sum_inter + interdist(i);
        sum2_inter = sum2_inter + interdist(i)^2;
        sum3_inter = sum3_inter + interdist(i)^3;
    end
    
    Stats_Inter(interval_count,1) = sum_inter/intercount;
    var_inter = (sum2_inter - sum_inter^2/intercount)/(intercount - 1);
    Stats_Inter(interval_count,2) = sqrt(var_inter);
    Stats_Inter(interval_count,3) = sum3_inter/(intercount * Stats_Inter(interval_count,2)^3);
    
end

figure
plot(interval_xaxis(:),Stats_Dev(:,1))
xlabel('Smoothing Interval (hrs) ')
ylabel('Mean Deviation (oC)')
phrase = [filename];
title(phrase);


figure
plot(interval_xaxis(:),Stats_Dev(:,2))
xlabel('Smoothing Interval (hrs) ')
ylabel('SD of Deviations (oC)')
phrase = [filename];
title(phrase);

figure
plot(interval_xaxis(:),Stats_Dev(:,3))
xlabel('Smoothing Interval (hrs) ')
ylabel('Gamma of Deviations (oC)')
phrase = [filename];
title(phrase);

figure
plot(interval_xaxis(:),Stats_Ex_Dev(:,1))
xlabel('Smoothing Interval (hrs) ')
ylabel('Mean Extr. Deviation (oC)')
phrase = [filename];
title(phrase);

%---------------------------------------------------------------
% output data to a text file

for i = 1:interval_count
    output_table(i,1) = interval_xaxis(i);
    output_table(i,2) = Stats_Ex_Dev(i,1);
end

fname = [fname,'_RetroLow_1yr_table.txt']
dlmwrite(fname,output_table);
%----------------------------------------------------------------

figure
plot(interval_xaxis(:),Stats_Ex_Dev(:,2))
xlabel('Smoothing Interval (hrs) ')
ylabel('SD of Extr. Deviations (oC)')
phrase = [filename];
title(phrase);

figure
plot(interval_xaxis(:),Stats_Ex_Dev(:,3))
xlabel('Smoothing Interval (hrs) ')
ylabel('Gamma of Extr. Deviations (oC)')
phrase = [filename];
title(phrase);

figure
plot(interval_xaxis(:),Stats_Craw(:,1))
xlabel('Smoothing Interval (hrs) ')
ylabel('Mean Raw Extremes(oC)')
phrase = [filename];
title(phrase);

figure
plot(interval_xaxis(:),Stats_Craw(:,2))
xlabel('Smoothing Interval (hrs) ')
ylabel('SD Raw Extremes(oC)')
phrase = [filename];
title(phrase);

figure
plot(interval_xaxis(:),Stats_Craw(:,3))
xlabel('Smoothing Interval (hrs) ')
ylabel('Gamma of Raw Extremes(oC)')
phrase = [filename];
title(phrase);

figure
plot(interval_xaxis(:),Stats_Inter(:,1))
xlabel('Smoothing Interval (hrs) ')
ylabel('Mean Interevent Interval(hrs)')
phrase = [filename];
title(phrase);

figure
plot(interval_xaxis(:),Stats_Inter(:,2))
xlabel('Smoothing Interval (hrs) ')
ylabel('SD Interevent Interval(hrs)')
phrase = [filename];
title(phrase);

figure
plot(interval_xaxis(:),Stats_Inter(:,3))
xlabel('Smoothing Interval (hrs) ')
ylabel('Gamma of Interevent Interval(hrs)')
phrase = [filename];
title(phrase);

figure
plot(interval_xaxis(:),Stats_Dur(:,1))
xlabel('Smoothing Interval (hrs) ')
ylabel('Mean Event Duration(hrs)')
phrase = [filename];
title(phrase);

figure
plot(interval_xaxis(:),Stats_Dur(:,2))
xlabel('Smoothing Interval (hrs) ')
ylabel('SD of Event Duration(hrs)')
phrase = [filename];
title(phrase);

figure
plot(interval_xaxis(:),Stats_Dur(:,3))
xlabel('Smoothing Interval (hrs) ')
ylabel('Gamma of Event Duration(hrs)')
phrase = [filename];
title(phrase);

figure
plot(interval_xaxis(:),Stats_Intensity(:,1))
xlabel('Smoothing Interval (hrs) ')
ylabel('Mean Intensity (degree hrs)')
phrase = [filename];
title(phrase);

figure
plot(interval_xaxis(:),Stats_Intensity(:,2))
xlabel('Smoothing Interval (hrs) ')
ylabel('SD of Event Intensity(degree hrs)')
phrase = [filename];
title(phrase);

figure
plot(interval_xaxis(:),Stats_Intensity(:,3))
xlabel('Smoothing Interval (hrs) ')
ylabel('Gamma of Event Intensity(degree hrs)')
phrase = [filename];
title(phrase);

figure
plot(interval_xaxis(:),Stats_Raw(:,1))
xlabel('Smoothing Interval (hrs) ')
ylabel('Cutoff (oC)')
phrase = [filename];
title(phrase);

figure
plot(interval_xaxis(:),Stats_Raw(:,2))
xlabel('Smoothing Interval (hrs) ')
ylabel('Raw Cutoff (oC)')
phrase = [filename];
title(phrase);

figure
plot(interval_xaxis(:),Stats_Raw(:,3))
xlabel('Smoothing Interval (hrs) ')
ylabel('Number of Extreme Events')
phrase = [filename];
title(phrase);


