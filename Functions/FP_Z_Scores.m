function [output_table,individual_trials] = FP_Z_Scores(blockpath, EXPERIMENT, SUBJECT)
%  FP_Z_SCORES UNTITLED2 Summary of this function goes here
% 
% Detailed explanation goes here
    data = TDTbin2mat(blockpath);
    
    % Potential for pulling into input variables
    Pre = 10;
    Post = 10;
    BL_Width = 5;
    BL_Start = 5;
    bin = 100;

    
    % Colors for plotting
    GCAMP_color = [0.8500, 0.3250, 0.0980];
    ISOS_color = [0.4660, 0.6740, 0.1880];
    cyan = [0.3010, 0.7450, 0.9330];
    gray1 = [.7 .7 .7];
    gray2 = [.8 .8 .8];

    % Stream Names:
    GCAMP = 'x465A';
    ISOS = 'x405A';
    REF_CHANNEL = EXPERIMENT;
    

    %% Setting up event timestamps
    % Timestamp refers to the moment the pellet is removed from the feeder
    data.epocs.(EXPERIMENT).onset = data.epocs.Note.onset;

    data.epocs.(EXPERIMENT).offSet = data.epocs.(EXPERIMENT).onset + .1;
    data.epocs.(EXPERIMENT).name = EXPERIMENT; 

    data.epocs.(EXPERIMENT).data = ones(1, length(data.epocs.(EXPERIMENT).onset));



    FS = data.streams.(GCAMP).fs;
    StartBuffer = 2000;
    EndBuffer = 2000;
    Ts=((1:numel(data.streams.(GCAMP).data(1,:))) /FS);
    TS = Ts(StartBuffer:(length(data.streams.(GCAMP).data) - EndBuffer));
    EventTS = data.epocs.(EXPERIMENT).onset;
    CH_GCAMP = data.streams.(GCAMP).data;
    CH_ISOS = data.streams.(ISOS).data;

    % Set Up BaseLine Window Length
    BaselineWind=round(BL_Start*FS); 
    BaselineWind2=round(BL_Width*FS); 
    
    
    %% ADD DEFAULT BIN!~
    % Create Window Size and eliminate events that overlap within a given window
        PreWind=round(Pre*FS);
        PostWind=round(Post*FS);
                
    % Eliminate events with overlapping windows.
    %NEED TO REMOVE -1 IF POSSIBLE!!!

    tmp=[];CurrEvent=EventTS(1);
    %REMOVE -1 HERE
    for i=1:length(EventTS)-1
        if EventTS(i+1)-CurrEvent>Pre+Post
            tmp(end+1,1)=CurrEvent;
            CurrEvent=EventTS(i+1);
        else
        end
    end

    tmp(end+1,1)=EventTS(length(EventTS),1);
    EventTS=tmp;

    %Find Time=0 for the event within the photometry data
    CSidx=[];
    %Ts=Ts'; %FIX FOR WRONG ARRAY ORIENTATION-REMOVED AND FIXED BELOW
    for i=1:length(EventTS)
       [MinVal, CSidx(i,1)]=min(abs(TS(:,:)-EventTS(i)));
    end

    %% Obtain the DeltaF/F for each event window
    CS_ISOS=[];CSTS=[];CS_GCAMP=[] ;
    CH_GCAMP=CH_GCAMP';
    CH_ISOS=CH_ISOS';
    counter=1;
    for i=1:length(CSidx)
        %%NEED TO CHECK IF INDEX IS <0 or GREATER THAN LENGTH OF DELTA490, NOT
        %%IF THE VALUE IS LESS THAN 0 OR GREATER THAN THE MAX VALUE OF DELTA490
        if CSidx(i)-(BaselineWind+BaselineWind2)<=0 || CSidx(i)+PostWind > length(Ts)
        else
            %Obtain data within baseline and event window for 405
            %and 490 channels.
            CSTS=(-PreWind:PostWind)./FS;
            CSBL(1,:)=CH_GCAMP((CSidx(i)-(BaselineWind+BaselineWind2)):(CSidx(i)-(BaselineWind))); 
            CSBL2(1,:)=CH_ISOS((CSidx(i)-(BaselineWind+BaselineWind2)):(CSidx(i)-(BaselineWind)));
            
            if i>length(CSidx)
                break
            elseif i<=length(CSidx) 
            CS_ISOS(1,:)=CH_ISOS((CSidx(i)-PreWind):(CSidx(i)+PostWind)); 
            CS_GCAMP(1,:)=CH_GCAMP((CSidx(i)-PreWind):(CSidx(i)+PostWind)); 
            end

            %Smooth to eliminate high frequency noise.
            F_GCAMP=smooth(CS_GCAMP,0.002,'lowess');  %Was 0.002- DJB
            F_ISOS=smooth(CS_ISOS,0.002,'lowess');
            F_GCAMPCSBL=smooth(CSBL,0.002,'lowess');  %Was 0.002- DJB
            F_ISOSCSBL=smooth(CSBL2,0.002,'lowess');

            %Scale and fit data
            bls=polyfit(F_ISOS(1:end),F_GCAMP(1:end),1);
            blsbase=polyfit(F_ISOSCSBL(1:end),F_GCAMPCSBL(1:end),1);
            Y_Fit=bls(1).*F_ISOS+bls(2);
            Y_Fit_base=blsbase(1).*F_ISOSCSBL+blsbase(2);
 
            %Center data and generate Delta F/F (DF_F) by dividing
            %event window by median of baseline window.
            DF_Event(:,i)=(F_GCAMP(:,1)-Y_Fit(:,1));
            DF_F(:,i)=DF_Event(:,i)./(Y_Fit); %Delta F/F
            DF_F(:,i)=DF_F(:,i)-DF_F(1,i);
             
            DF_Event(:,i)=(F_GCAMP(:,1)-Y_Fit(:,1));
            DF_Event(:,i)=DF_Event(:,i)-DF_Event(1,i); %zero first point
            DF_Base(:,i)=(F_GCAMPCSBL(:,1)-Y_Fit_base(:,1));
            DF_Base(:,i)=DF_Base(:,i)-DF_Base(1,i); %zero to first point
            DF_ZScore(:,counter)=(DF_Event(:,i)-median(DF_Base(:,i)))./mad(DF_Base(:,i)); %Z-Score
            counter=counter+1;
            %Center data and generate Delta F/F (DF_F) by dividing
            %event window by median of baseline window.
            DF_Event(:,i)=(F_GCAMP(:,1)-Y_Fit(:,1));
            DF_F(:,i)=DF_Event(:,i)./(Y_Fit); %Delta F/F
            DF_F(:,i)=DF_F(:,i)-DF_F(1,i);
            
            DF_Event(:,i)=(F_GCAMP(:,1)-Y_Fit(:,1));
            DF_Event(:,i)=DF_Event(:,i)-DF_Event(1,i); %zero first point
            DF_Base(:,i)=(F_GCAMPCSBL(:,1)-Y_Fit_base(:,1));
            DF_Base(:,i)=DF_Base(:,i)-DF_Base(1,i); %zero to first point
            DF_ZScore(:,counter)=(DF_Event(:,i)-median(DF_Base(:,i)))./mad(DF_Base(:,i)); %Z-Score
            counter=counter+1;
      
            clear CS_ISOS CS_GCAMP F_GCAMP F_ISOS bls Y_Fit

         end
     end

    %% BINNING OF ALL DATA
    DF_ZScore=DF_ZScore';
    tmp=[];tmp2=[];
    for i=1:bin:length(CSTS)
        if i+bin>length(CSTS)
        tmp(1,end+1)=median(CSTS(i:end));
        tmp2(:,end+1)=median(DF_ZScore(:,i:end),2);
        else
        tmp(1,end+1)=median(CSTS(i:i+bin));
        tmp2(:,end+1)=median(DF_ZScore(:,i:i+bin),2);
        end
    end
    tsGCAMP=tmp;
    DF_ZScore=tmp2;
    zStdError = std(DF_ZScore)/sqrt(size(DF_ZScore,1));

    %% Find mean z-score
    zGCAMP_Mean=(mean(DF_ZScore,1));

    %% Exporting data formated as a table
    output_table = table( ...
        tsGCAMP', zGCAMP_Mean', zStdError', ...
        'VariableNames', {'Time', 'Z Score', 'Std Dev'} ...
    );
    time_series = num2cell(tsGCAMP');
    dims = size(DF_ZScore');
    trial_data = mat2cell(DF_ZScore', ones(1,dims(1)), ones(1,dims(2)));

    individual_trials = horzcat(time_series,trial_data);
    
end