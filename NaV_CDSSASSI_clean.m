%   For NaV protocol :)

%   Variables
T2P_upperLim = 1100 ; % set time to peak limit
T2P_lowerLim = 500 ; % set time to peak starting point (set to 500us to avoid capacitance errors at start of depolarisation)
pA_LowerLim = -0.0000000004 ; %400pA
pA_UpperLim = -0.000000002 ; %2nA

SSA_timestart = 7000 ; %250000us is the start of depol from -120mV (from voltage protocol)
SSA_timeend = 7200 ; % cursor end of analysis sector
SSI_timestart = 17000 ; % protocol goes for 1s, therefore ends at 12500
SSI_timeend = 17060 ; % cursor end of analysis sector
LK_timestart = 991 ; %350000us is the start of depol from -120mV
LK_timeend = 2500 ;
voltages = [-140:5:60] ; % What voltages are tested?

%   Paths - import
Parent_dir = '\\Data\Vandenberg-Lab\Syncropatch\Functional_Genomics_Syncropatch\SCN5A\2.Control Mutants\' ;
Exp_dir = '20221128_JM1\' ; % update with exp date
Protocol_dir = 'NaV_ssAct_InAct_500ms_JM_10.39.45\' ; % update with exp protocol


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Should not need to change below

mkdir(fullfile(Parent_dir,Exp_dir,Protocol_dir,'PostQC'))
PostQC_dir = 'PostQC\' ; 
QCdWells_dir = 'QCd_Wells\' ;
QCd_wells = dir(fullfile(Parent_dir,Exp_dir,Protocol_dir,PostQC_dir,QCdWells_dir,'*.csv')) ;


%   Paths - export1 Peaks&Times
folderdir = fullfile(Parent_dir,Exp_dir,Protocol_dir,PostQC_dir) ;
mkdir(folderdir,'PeaksAndTimes')
PnT_dir = 'PeaksAndTimes\' ;
folderPTdir = fullfile(Parent_dir,Exp_dir,Protocol_dir,PostQC_dir,PnT_dir) ;
mkdir(folderPTdir,'Raw_SSA')
Raw_SSA_dir = 'Raw_SSA\' ;
mkdir(folderPTdir,'Raw_SSI')
Raw_SSI_dir = 'Raw_SSI\' ;

PeaksAtMV = 'SSA_PeaksAtMV.csv' ; % Enter name of file wanted
PeaksAtMV2 = 'SSI_PeaksAtMV.csv' ; % Enter name of file wanted

%   Paths - export2 Current Density
mkdir(folderdir,'CurrentDensity')
CD_dir = 'CurrentDensity' ;

CD_perSweep = 'CD_persweep.csv' ; % Enter name of file wanted
CD_summary = 'CD_perVariant.csv' ; % Enter name of file wanted
CD_calc = 'CD_calculated.csv' ; % Enter name of file wanted

% Paths - export3 SSA Boltz

mkdir(folderdir,'SteadyStateActivations_plots')
SSAplot_dir = 'SteadyStateActivations_plots\' ;
mkdir(folderdir,'SteadyStateActivation_V50')
SSA_dir = 'SteadyStateActivation_V50\' ;

Peaks_Norm_SSA = 'SSA.csv' ; % Enter name of file wanted
SSAsummary = 'SSAsummary.csv' ; % Enter name of file wanted

% Paths - export4 SSI Boltz
mkdir(folderdir,'SteadyStateInactivations_V50')
SSI_dir = 'SteadyStateInactivations_V50\' ;
mkdir(folderdir,'SteadyStateInactivations_plots')
SSIplot_dir = 'SteadyStateInactivations_plots\' ;

Peaks_Norm_SSI = 'SSI.csv' ; % Enter name of file wanted
SSIsummary = 'SSIsummary.csv' ; % Enter name of file wanted



%TimeAtPeaks = 'SSA_TimeAtPeaks.csv' ; % Enter name of file wanted
%TimeAtPeaks2 = 'SSI_TimeAtPeaks.csv' ; % Enter name of file wanted

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Times_us_SSA = [SSA_timestart*50-350000:50:SSA_timeend*50-350000]' ;
Times_us_SSI = [SSI_timestart*50-850000:50:SSI_timeend*50-850000]' ;
voltagesSSA = voltages(9:end)
voltagesSSI = voltages(1:21)

%   Import data
opts = delimitedTextImportOptions("NumVariables", 43);
opts.DataLines = [4, Inf];
opts.Delimiter = "\t";
opts.VariableNames = ["Var1", "Var2", "RawDataA", "RawDataA1", "RawDataA2", "RawDataA3", "RawDataA4", "RawDataA5", "RawDataA6", "RawDataA7", "RawDataA8", "RawDataA9", "RawDataA10", "RawDataA11", "RawDataA12", "RawDataA13", "RawDataA14", "RawDataA15", "RawDataA16", "RawDataA17", "RawDataA18", "RawDataA19", "RawDataA20", "RawDataA21", "RawDataA22", "RawDataA23", "RawDataA24", "RawDataA25", "RawDataA26", "RawDataA27", "RawDataA28", "RawDataA29", "RawDataA30", "RawDataA31", "RawDataA32", "RawDataA33", "RawDataA34", "RawDataA35", "RawDataA36", "RawDataA37", "RawDataA38", "RawDataA39", "RawDataA40"];
opts.SelectedVariableNames = ["RawDataA", "RawDataA1", "RawDataA2", "RawDataA3", "RawDataA4", "RawDataA5", "RawDataA6", "RawDataA7", "RawDataA8", "RawDataA9", "RawDataA10", "RawDataA11", "RawDataA12", "RawDataA13", "RawDataA14", "RawDataA15", "RawDataA16", "RawDataA17", "RawDataA18", "RawDataA19", "RawDataA20", "RawDataA21", "RawDataA22", "RawDataA23", "RawDataA24", "RawDataA25", "RawDataA26", "RawDataA27", "RawDataA28", "RawDataA29", "RawDataA30", "RawDataA31", "RawDataA32", "RawDataA33", "RawDataA34", "RawDataA35", "RawDataA36", "RawDataA37", "RawDataA38", "RawDataA39", "RawDataA40"];
opts.VariableTypes = ["string", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, ["Var1", "Var2"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Var2"], "EmptyFieldRule", "auto");

%   Import SSA/SSI data 
n_wellsAll = size(QCd_wells);
for wells = 1:n_wellsAll(1)
    fid = fullfile(Parent_dir,Exp_dir,Protocol_dir,PostQC_dir,QCdWells_dir,QCd_wells(wells).name);
    fidopen = fopen(fid);
    QC_rawWells_full = readtable(fid,opts);
    QC_rawWells_data{wells,1} = QCd_wells(wells).name ;
    QC_rawWells_data{wells,2} = QC_rawWells_full(SSA_timestart+1:SSA_timeend+1,9:end); % only take voltages -100mV onwards
    QC_rawWells_data{wells,3} = QC_rawWells_full(SSI_timestart+1:SSI_timeend+1,1:21); % only take voltages -140mV to -40mV
    QC_rawWells_data{wells,4} = QC_rawWells_full(LK_timestart+1:LK_timeend+1,1:21); % only take voltages -140mV to -40mV
    fclose("all");
end
clear opts


%   Extraction peaks and times indexes per sweeps - SSA
eRev_complete = [] ;
wellID_full = [] ;
titleID_full = [] ;
peaks_sweep = nan(n_wellsAll(1),33) ;
times_sweep = nan(n_wellsAll(1),33) ;

for wells = 1:n_wellsAll(1)    
    temp = QC_rawWells_data{wells,2} ;
    titleWell3 = char(string(QC_rawWells_data{wells,1})) ;
    
    dataplot = table2array(temp(:,:)) ; % to plot SSA raw sweep data
    plot(Times_us_SSA,dataplot,"-k","LineWidth",0.5) ;
    hold on
    plot(Times_us_SSA(17),dataplot(:,17),"-r","LineWidth",2) ;
    xlabel("Time (ms)") ;
    ylabel("Amplitude (pA)") ;
    yline(0,'--','r') ;
    wellID = string(titleWell3(39:41)) ;
    titleWell = strcat("Current per sweep ",wellID,".jpg")  ;
    title(titleWell) ;
    hold off
    
    pid = fullfile(Parent_dir,Exp_dir,Protocol_dir,PostQC_dir,PnT_dir,Raw_SSA_dir,titleWell) ;
    exportgraphics(gca,pid) ;

    wellID_full = [wellID_full;wellID] ;
    titleID_full = [titleID_full;titleWell3] ;
    
    datasweeps_fromSSAstart = table2array(temp(11:end,:)) ; %ind_start:end is to remove data from first 500ms from analysis as they consist of artefacts
    time_SSA = Times_us_SSA(11:end) ;
    ind_start = 1 ;
    
    for sweeps = 21:-1:1
        [minval,minind] = min(datasweeps_fromSSAstart(ind_start:end,sweeps)) ; % for peaks
        peaktime = time_SSA(minind+ind_start-1) ;
        times_sweep(wells,sweeps) = peaktime ;
        peaks_sweep(wells,sweeps) = minval   ; 
        if sweeps == 21, ind_start = minind ;
        else 
            ind_start = ind_start-1+minind ; % if second peak at same location as first peak then second idx=1 (relative to first) but 1+4 (e.g. if 4 was first index) = 5, hence minus 1
        end   
    end

    for sweeps = 22:33   
        [maxval,maxind] = max(datasweeps_fromSSAstart(:,sweeps)) ;
        [minval,minind] = min(datasweeps_fromSSAstart(:,sweeps)) ;
        if abs(maxval) > abs(minval)
            peak = maxval ;
            peaktime = time_SSA(maxind) ;
        else peak = minval ;
            peaktime = time_SSA(minind) ;
        end
        times_sweep(wells,sweeps) = peaktime ;
        peaks_sweep(wells,sweeps) = peak   ;     
    end

    try 
    eRev = interp1(peaks_sweep(wells,21:33),voltagesSSA(21:33),0,"linear","extrap");
    catch
        eRev = string("eRev calculation failed") ;
    end 
    eRev_complete = [eRev_complete;eRev] ;
      
    plot(voltages(9:end),peaks_sweep(wells,:),"ko-","MarkerFaceColor","auto", "LineWidth",1) ;
    hold on
    try 
    plot(eRev,0,"ko","MarkerFaceColor","b","LineWidth",1) ;
    catch
    end
    xlabel("Voltages (mV)") ;
    ylabel("pA") ;
    legend("Activation","Reversal Potential","location","southwest") ;
    titleWell3 = strcat("Activation by Voltage ",wellID,".jpg") ; %%%% can set to pdf etc
    title(titleWell3)    ;
    hold off
    
    pid = fullfile(Parent_dir,Exp_dir,Protocol_dir,PostQC_dir,PnT_dir,Raw_SSA_dir,titleWell3) ;
    exportgraphics(gca,pid)     ;
end

eRev_completetable = array2table(eRev_complete) ;
eRev_completetable.Properties.RowNames = wellID_full ;
eRev_completetable.titleID_full = titleID_full
writetable(eRev_completetable,fullfile(Parent_dir,Exp_dir,Protocol_dir,PostQC_dir,PnT_dir,"ActivationeRevAll.csv"),"WriteRowNames",true)


SSA_table = array2table(peaks_sweep)
SSA_table.Properties.RowNames = wellID_full ;
SSA_table.titleID_full = titleID_full

SSAtimes_table = array2table(times_sweep)
SSAtimes_table.Properties.RowNames = wellID_full ;
SSAtimes_table.titleID_full = titleID_full


%	Extraction peaks and times indexes per sweeps - SSI
    peaks_sweep = nan(n_wellsAll(1),21) ;
    times_sweep = nan(n_wellsAll(1),21) ;

for wells = 1:n_wellsAll(1)    
    temp = QC_rawWells_data{wells,3} ;
    titleWell3 = char(string(QC_rawWells_data{wells,1})) ;
    dataplot = table2array(temp(:,:)) ; % to plot SSI raw sweep data
    
    plot(Times_us_SSI,dataplot) 
    hold on
    xlabel("Time (s)")
    ylabel("pA")
    yline(0,'--')
    wellID = string(titleWell3(33:35)) ;
    titleWell = strcat("Current per sweep ",wellID,".jpg")  ;
    title(titleWell)
    hold off
    
    pid = fullfile(Parent_dir,Exp_dir,Protocol_dir,PostQC_dir,PnT_dir,Raw_SSI_dir,titleWell) ;
    exportgraphics(gca,pid)   
    
    datasweeps_fromSSIstart = table2array(temp(11:end,:)) ; %ind_start:end is to remove data from first 500ms from analysis as they consist of artefacts
    time_SSI = Times_us_SSI(11:end) ;
        
    for sweeps = 1:21
        [minval,minind] = min(datasweeps_fromSSIstart(:,sweeps)) ; % for peaks
        peaktime = time_SSI(minind) ;
        times_sweep(wells,sweeps) = peaktime ;
        peaks_sweep(wells,sweeps) = minval   ;       
    end

    plot(voltages,peaks_sweep(wells,:),"ko-")
    hold on
    xlabel("Voltages (mV)")
    ylabel("pA")
    titleWell3 = strcat("Inactivation by Voltage ",wellID,".jpg") ;
    title(titleWell3)
    hold off    
    pid = fullfile(Parent_dir,Exp_dir,Protocol_dir,PostQC_dir,PnT_dir,Raw_SSI_dir,titleWell3) ;
    exportgraphics(gca,pid)                 
end 

SSI_table = array2table(peaks_sweep)
SSI_table.Properties.RowNames = wellID_full ;
SSI_table.titleID_full = titleID_full


%   import variable names/QC
opts = delimitedTextImportOptions("NumVariables", 125);
opts.DataLines = [3, Inf];
opts.Delimiter = "\t";
opts.VariableNames = ["WellID", "CellType", "SealResistance", "Capacitance", "SeriesResistance", "SealResistance1", "Capacitance1", "SeriesResistance1", "SealResistance2", "Capacitance2", "SeriesResistance2", "SealResistance3", "Capacitance3", "SeriesResistance3", "SealResistance4", "Capacitance4", "SeriesResistance4", "SealResistance5", "Capacitance5", "SeriesResistance5", "SealResistance6", "Capacitance6", "SeriesResistance6", "SealResistance7", "Capacitance7", "SeriesResistance7", "SealResistance8", "Capacitance8", "SeriesResistance8", "SealResistance9", "Capacitance9", "SeriesResistance9", "SealResistance10", "Capacitance10", "SeriesResistance10", "SealResistance11", "Capacitance11", "SeriesResistance11", "SealResistance12", "Capacitance12", "SeriesResistance12", "SealResistance13", "Capacitance13", "SeriesResistance13", "SealResistance14", "Capacitance14", "SeriesResistance14", "SealResistance15", "Capacitance15", "SeriesResistance15", "SealResistance16", "Capacitance16", "SeriesResistance16", "SealResistance17", "Capacitance17", "SeriesResistance17", "SealResistance18", "Capacitance18", "SeriesResistance18", "SealResistance19", "Capacitance19", "SeriesResistance19", "SealResistance20", "Capacitance20", "SeriesResistance20", "SealResistance21", "Capacitance21", "SeriesResistance21", "SealResistance22", "Capacitance22", "SeriesResistance22", "SealResistance23", "Capacitance23", "SeriesResistance23", "SealResistance24", "Capacitance24", "SeriesResistance24", "SealResistance25", "Capacitance25", "SeriesResistance25", "SealResistance26", "Capacitance26", "SeriesResistance26", "SealResistance27", "Capacitance27", "SeriesResistance27", "SealResistance28", "Capacitance28", "SeriesResistance28", "SealResistance29", "Capacitance29", "SeriesResistance29", "SealResistance30", "Capacitance30", "SeriesResistance30", "SealResistance31", "Capacitance31", "SeriesResistance31", "SealResistance32", "Capacitance32", "SeriesResistance32", "SealResistance33", "Capacitance33", "SeriesResistance33", "SealResistance34", "Capacitance34", "SeriesResistance34", "SealResistance35", "Capacitance35", "SeriesResistance35", "SealResistance36", "Capacitance36", "SeriesResistance36", "SealResistance37", "Capacitance37", "SeriesResistance37", "SealResistance38", "Capacitance38", "SeriesResistance38", "SealResistance39", "Capacitance39", "SeriesResistance39", "SealResistance40", "Capacitance40", "SeriesResistance40"];
opts.VariableTypes = ["string", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, "WellID", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["WellID", "CellType"], "EmptyFieldRule", "auto");

QC_dir = 'QC\' ;
QCimport = readtable(fullfile(Parent_dir,Exp_dir,Protocol_dir,QC_dir,'QC'),opts,"ReadRowNames",true);
VarNames = QCimport(:,2);
clear opts


%   Time to peak
VarNames = VarNames(contains(VarNames.Properties.RowNames,SSA_table.Properties.RowNames),:)  ;
VarNames = sortrows(VarNames,"Row","ascend")
SSA_table.Variant = VarNames.CellType

toKeepNAVar = SSA_table.peaks_sweep21 > -0.00000000005 ;% all cells below 50pa passes (regardless of Brs or bad cell) and all above 50pA has to go through t2p calculations
toKeepT2Peak = SSA_table.peaks_sweep21 <= -0.00000000005 & SSAtimes_table.times_sweep21 <= T2P_upperLim & SSAtimes_table.times_sweep21 > T2P_lowerLim 
toKeep_SSAPeaks = SSA_table(logical(toKeepT2Peak+toKeepNAVar),:) 
toKeep_SSIPeaks = SSI_table(contains(SSI_table.Properties.RowNames,toKeep_SSAPeaks.Properties.RowNames),:) % Match SSI wellls to SSA wells
toKeep_SSIPeaks.Variant = toKeep_SSAPeaks.Variant

WellNames = char(toKeep_SSAPeaks.Properties.RowNames) 
WellNames_alphabet = WellNames(:,1)
WellNames_numeric = WellNames(:,2:3)

toKeep_SSAPeaks.WellAlph = WellNames_alphabet ;
toKeep_SSAPeaks.WellNum = WellNames_numeric ;
toKeep_SSIPeaks.WellAlph = WellNames_alphabet ;
toKeep_SSIPeaks.WellNum = WellNames_numeric ;
toKeep_SSAPeaks = sortrows(toKeep_SSAPeaks,["WellNum","WellAlph"],{'ascend','ascend'}) ;
toKeep_SSIPeaks = sortrows(toKeep_SSIPeaks,["WellNum","WellAlph"],{'ascend','ascend'}) ;
toKeep_SSAPeaks.WellAlph = [] ;
toKeep_SSAPeaks.WellNum = [] ;
toKeep_SSIPeaks.WellAlph = [] ;
toKeep_SSIPeaks.WellNum = [] ;


%   Exporting Peaks
writetable(toKeep_SSAPeaks,fullfile(Parent_dir,Exp_dir,Protocol_dir,PostQC_dir,PnT_dir,PeaksAtMV),"WriteRowNames",true)
writetable(toKeep_SSIPeaks,fullfile(Parent_dir,Exp_dir,Protocol_dir,PostQC_dir,PnT_dir,PeaksAtMV2),"WriteRowNames",true)
 
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%   Current Density % CDsqrt
QC_comp = QCimport
QC_comp = movevars(QC_comp,"CellType","Before","WellID")
QC_Cslow = QC_comp(:,1:3:end)  % capacitances

QC_Cslow_filtered = QC_Cslow(contains(QC_Cslow.Properties.RowNames,toKeep_SSAPeaks.Properties.RowNames),:)  ;
Peaks_filtered = toKeep_SSAPeaks(contains(toKeep_SSAPeaks.Properties.RowNames,QC_Cslow_filtered.Properties.RowNames),:)  ;

Peaks_Cslow = array2table(table2array(Peaks_filtered(:,1:end-2))./table2array(QC_Cslow_filtered(:,2:34)))
Peaks_Cslow.Properties.RowNames = Peaks_filtered.Properties.RowNames
Peaks_Cslow.Properties.VariableNames = Peaks_filtered.Properties.VariableNames(1:33)
Peaks_Cslow.Variant = Peaks_filtered.Variant

Peaks_Cslow = movevars(Peaks_Cslow,"Variant","Before","peaks_sweep1") ; % Per sweep dataset
Peaks_CslowSummary = groupsummary(Peaks_Cslow,"Variant","mean")

CurrentDensitysqrt17 = table(Peaks_Cslow.Variant,Peaks_Cslow.peaks_sweep17,'VariableNames',["Variant","CD_sweep17"])
CurrentDensitysqrt17.Properties.RowNames = Peaks_Cslow.Properties.RowNames
CurrentDensitysqrt17.ABSSQRT17 = sqrt(abs(Peaks_Cslow.peaks_sweep17)) % Calculated dataset

WTmeanSQRT4 = groupsummary(CurrentDensitysqrt17,"Variant","mean")
WTmeanSQRT4.Properties.RowNames = string(Peaks_CslowSummary.Variant)
WTmeanSQRT2 = WTmeanSQRT4{'SCN5A_DelQ1077_WT', 'mean_ABSSQRT17'}
WTmean2 = WTmeanSQRT4{'SCN5A_DelQ1077_WT', 'mean_CD_sweep17'}

CurrentDensitysqrt17.normWTmean = CurrentDensitysqrt17.CD_sweep17 / WTmean2
CurrentDensitysqrt17.normSQRTWTmean = CurrentDensitysqrt17.ABSSQRT17 / WTmeanSQRT2

%   Exporting CD, CD summary per variant, calculated CD for sweep 17 (-20mV)
writetable(Peaks_Cslow,fullfile(Parent_dir,Exp_dir,Protocol_dir,PostQC_dir,CD_dir,CD_perSweep),"WriteRowNames",true,"WriteVariableNames",true)
writetable(Peaks_CslowSummary,fullfile(Parent_dir,Exp_dir,Protocol_dir,PostQC_dir,CD_dir,CD_summary),"WriteRowNames",true,"WriteVariableNames",true)
writetable(CurrentDensitysqrt17,fullfile(Parent_dir,Exp_dir,Protocol_dir,PostQC_dir,CD_dir,CD_calc),"WriteRowNames",true,"WriteVariableNames",true)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Steady State Activations

voltagesSSA = voltages(9:29)

PeaksAtMVcomp = toKeep_SSAPeaks(:,1:21)
PeaksAtMV_filterLowertemp = PeaksAtMVcomp.peaks_sweep17 < pA_LowerLim & PeaksAtMVcomp.peaks_sweep17 > pA_UpperLim;
PeaksAtMV_filterLow_Upp = PeaksAtMVcomp(PeaksAtMV_filterLowertemp,:) ;

Erev = eRev_completetable
Erev = Erev(contains(Erev.Properties.RowNames,PeaksAtMV_filterLow_Upp.Properties.RowNames),:) ;

WellID = char(string(Erev.Properties.RowNames))
WellNames_alphabet = WellID(:,1)
WellNames_numeric = WellID(:,2:3)
Erev.WellAlph = WellNames_alphabet ;
Erev.WellNum = WellNames_numeric ;
Erev = sortrows(Erev,["WellNum","WellAlph"],{'ascend','ascend'}) ;
Erev.WellAlph = [] ;
Erev.WellNum = [] 
test = str2double(table2array(Erev(:,1)))
test2 = table(test,'RowNames',Erev.Properties.RowNames)
test2(any(ismissing(test2),2), :) = [];   %if the syncropatch cannot calculate Erev, it will be NaN

Erev_arr = table2array(test2)
PeaksAtMV_filterLow_Upp = PeaksAtMV_filterLow_Upp(contains(PeaksAtMV_filterLow_Upp.Properties.RowNames,test2.Properties.RowNames),:) ;


Peaks_arr = table2array(PeaksAtMV_filterLow_Upp) ;
n_Peaks = size(PeaksAtMV_filterLow_Upp) ;

 for iWells = 1:n_Peaks(1)
    for iSweeps = 1:21 
        Peaks_eRev(iWells,iSweeps) = Peaks_arr(iWells,iSweeps)/(voltagesSSA(iSweeps)-Erev_arr(iWells,1)) ;
    end 
    IMax(iWells) = max(Peaks_eRev(iWells,:)) ;
    Peaks_Norm(iWells,:) = Peaks_eRev(iWells,:)./IMax(iWells) ;
 end
 
 
Peaks_Norm2 = table(Peaks_Norm) ;
Peaks_Norm2.Properties.RowNames = PeaksAtMV_filterLow_Upp.Properties.RowNames
Peaks_Norm2 = Peaks_Norm2(contains(Peaks_Norm2.Properties.RowNames,PeaksAtMV_filterLow_Upp.Properties.RowNames),:) ;
Peaks_Norm2 = table2array(Peaks_Norm2) ;
voltagesSSA2 = [-100:5:0];

BoltzmannEq = 'BoltzMin + ((BoltzMax-BoltzMin)./(1+exp((V50-x)./k)))';
Bfitdata_V50 = [] ;
Bfitdata_Slope = []  ;
Bfitdata_R2 = [] ;
Bfitdata_Bmin = [] ;
Bfitdata_BMax = [] ;


for iWells = 1:n_Peaks(1)
    BoltzFit = Peaks_Norm(iWells,:)' ;
    [fittedCurve,fittedData] = fit(voltagesSSA2',BoltzFit,BoltzmannEq,'StartPoint',[0,1,-40,5]);  

    Bfitdata_V50 = [Bfitdata_V50;fittedCurve.V50] ;
    Bfitdata_Slope = [Bfitdata_Slope;fittedCurve.k] ;
    Bfitdata_R2 = [Bfitdata_R2;fittedData.rsquare] ;
    Bfitdata_Bmin = [Bfitdata_Bmin;fittedCurve.BoltzMin] ;
    Bfitdata_BMax = [Bfitdata_BMax;fittedCurve.BoltzMax] ;

    clf
    plot(voltagesSSA2,Peaks_Norm2(iWells,1:21)',"k:*") % change well as required
    hold on
    plot(fittedCurve,'c-') % change well as required
    xlabel("Voltage (mV)")
    ylabel("Normalised Current")
    yline(0.5,'--')
    ylim([-0.1,1.1])
    yticks(0:0.5:1)
    titleWell = string(PeaksAtMV_filterLow_Upp.Properties.RowNames(iWells)) ;
    titleWell = strcat("Current by Voltage ",titleWell,".jpg") ; %change between pdf and jpg if preferred
    title(titleWell);
    xtickangle(45)
    legend("Data", "Fitted Boltzmann","Location","northwest")
    legend("boxoff")
    hold off

    pid = fullfile(Parent_dir,Exp_dir,Protocol_dir,PostQC_dir,SSAplot_dir,titleWell) ;
    exportgraphics(gca,pid)
end


BfitdataSA = table(Bfitdata_V50,Bfitdata_Slope,Bfitdata_R2,Bfitdata_Bmin,Bfitdata_BMax,Peaks_Norm2) ;
BfitdataSA.Properties.RowNames = PeaksAtMV_filterLow_Upp.Properties.RowNames

BfitdataR2 = BfitdataSA.Bfitdata_R2 >0.99 ;
BfitdataSSA = BfitdataSA(BfitdataR2,:) ;
VarNames = QCimport(:,2)
VarNames = VarNames(contains(VarNames.Properties.RowNames,BfitdataSSA.Properties.RowNames),:)
BfitdataSSA.Variant = VarNames.CellType

WTmean1 = groupsummary(BfitdataSSA,"Variant","mean")
WTmean1.Properties.RowNames = string(WTmean1.Variant)
WTmean2 = WTmean1{'SCN5A_DelQ1077_WT', 'mean_Bfitdata_V50'}
BfitdataSSA.V50subWTmean = BfitdataSSA.Bfitdata_V50-WTmean2
Bfitdata_SSA_reorg = movevars(BfitdataSSA,"Variant","Before","Bfitdata_V50") ;
Bfitdata_SSA_reorg = movevars(Bfitdata_SSA_reorg,"V50subWTmean","After","Bfitdata_V50") 

writetable(Bfitdata_SSA_reorg,fullfile(Parent_dir,Exp_dir,Protocol_dir,PostQC_dir,SSA_dir,Peaks_Norm_SSA),"WriteRowNames",true,"WriteVariableNames",true) ;
writetable(WTmean1,fullfile(Parent_dir,Exp_dir,Protocol_dir,PostQC_dir,SSA_dir,SSAsummary),"WriteRowNames",true,"WriteVariableNames",true) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Steady State Inactivation

PeaksAtMVcomp = toKeep_SSIPeaks
PeaksAtMV_filterLowertemp = PeaksAtMVcomp.peaks_sweep1 < pA_LowerLim & PeaksAtMVcomp.peaks_sweep1 > pA_UpperLim ;
PeaksAtMV_filterLow_Upp = PeaksAtMVcomp(PeaksAtMV_filterLowertemp,:) ;

n_Peaks = size(PeaksAtMV_filterLow_Upp) ;
PeaksAtMV_filtered_MAT = table2array(PeaksAtMV_filterLow_Upp(:,1:21))
clear Peaks_Norm

 for iWells = 1:n_Peaks(1)
    IMax(iWells,1) = max(PeaksAtMV_filtered_MAT(iWells,:)) ;
    Peaks_e2(iWells,:) = PeaksAtMV_filtered_MAT(iWells,:)-IMax(iWells,1) ;
    IMin(iWells,1) = min(Peaks_e2(iWells,:)) ;
    Peaks_Norm(iWells,:) = Peaks_e2(iWells,:)./IMin(iWells,1) ; 
 end

Peaks_Norm2 = table(Peaks_Norm) ;
Peaks_Norm2.Properties.RowNames = PeaksAtMV_filterLow_Upp.Properties.RowNames
Peaks_Norm2 = Peaks_Norm2(contains(Peaks_Norm2.Properties.RowNames,PeaksAtMV_filterLow_Upp.Properties.RowNames),:) ;
Peaks_Norm2 = table2array(Peaks_Norm2) ;

BoltzmannEq = 'BoltzMin + ((BoltzMax-BoltzMin)./(1+exp((V50-x)./k)))';
Bfitdata_V50 = [] ;
Bfitdata_Slope = []  ;
Bfitdata_R2 = [] ;
Bfitdata_Bmin = [] ;
Bfitdata_BMax = [] ;

for iWells = 1:n_Peaks(1)
    BoltzFit = Peaks_Norm2(iWells,:)' ;
    [fittedCurve,fittedData] = fit(voltagesSSI',BoltzFit,BoltzmannEq,'StartPoint',[0,1,-80,-9]);  
    Bfitdata_V50 = [Bfitdata_V50;fittedCurve.V50] ;
    Bfitdata_Slope = [Bfitdata_Slope;fittedCurve.k] ;
    Bfitdata_R2 = [Bfitdata_R2;fittedData.rsquare] ;
    Bfitdata_Bmin = [Bfitdata_Bmin;fittedCurve.BoltzMin] ;
    Bfitdata_BMax = [Bfitdata_BMax;fittedCurve.BoltzMax] ;

    clf
    plot(voltagesSSI,Peaks_Norm2(iWells,:)',"k:*") % change well as required
    hold on
    plot(fittedCurve,'c-') % change well as required
    xlabel("Voltage (mV)")
    ylabel("Normalised Current")
    yline(0.5,'--')
    ylim([-0.1,1.1])
    yticks(0:0.5:1)
    titleWell = string(PeaksAtMV_filterLow_Upp.Properties.RowNames(iWells)) ;
    titleWell = strcat("Current by Voltage ",titleWell,".jpg") ;
    title(titleWell);
    xtickangle(45)
    legend("Data", "Fitted Boltzmann","Location","northeast")
    legend("boxoff")
    hold off

    pid = fullfile(Parent_dir,Exp_dir,Protocol_dir,PostQC_dir,SSIplot_dir,titleWell) ;
    exportgraphics(gca,pid)
end

BfitdataSI = table(Bfitdata_V50,Bfitdata_Slope,Bfitdata_R2,Bfitdata_Bmin,Bfitdata_BMax,Peaks_Norm2) ;
BfitdataSI.Properties.RowNames = PeaksAtMV_filterLow_Upp.Properties.RowNames
BfitdataR2 = BfitdataSI.Bfitdata_R2 >0.99 ;
BfitdataSSI = BfitdataSI(BfitdataR2,:) ;
VarNames = QCimport(:,2)
VarNames = VarNames(contains(VarNames.Properties.RowNames,BfitdataSSI.Properties.RowNames),:)
BfitdataSSI.Variant = VarNames.CellType

WTmeanSSI = groupsummary(BfitdataSSI,"Variant","mean")
WTmeanSSI.Properties.RowNames = string(WTmeanSSI.Variant)
WTmeanSQRTSSI2 = WTmeanSSI{'SCN5A_DelQ1077_WT', 'mean_Bfitdata_V50'}
BfitdataSSI.V50subWTmean = BfitdataSSI.Bfitdata_V50-WTmeanSQRTSSI2 ;
Bfitdata_SSI_reorg = movevars(BfitdataSSI,"Variant","Before","Bfitdata_V50") ;
Bfitdata_SSI_reorg = movevars(Bfitdata_SSI_reorg,"V50subWTmean","After","Bfitdata_V50")

writetable(Bfitdata_SSI_reorg,fullfile(Parent_dir,Exp_dir,Protocol_dir,PostQC_dir,SSI_dir,Peaks_Norm_SSI),"WriteRowNames",true,"WriteVariableNames",true) ;
writetable(WTmeanSSI,fullfile(Parent_dir,Exp_dir,Protocol_dir,PostQC_dir,SSI_dir,SSIsummary),"WriteRowNames",true,"WriteVariableNames",true) ;

clear
