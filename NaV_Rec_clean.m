% Recoveries
% Export all wells and a QC file (Rseal, Cslow, Rseries) of all 4 protocols.
% Export raw uncorrected leak file as well as raw file.
% QCs in this file - Rseal, Cslow, corrected leak, uncorrected leak, peak
% cut off, mean recovery>1 cut off, inverse slope cut off??


% Variables
LeakCheck = 95000 
LeakCheck2 = LeakCheck/50+1

Rec_control1_start = 1720500 ;
Rec_1ms_start = 2291500 ;
Rec_2ms_start = 2863500 ;
Rec_3ms_start = 3436500 ;
Rec_4ms_start = 4010500 ;
Rec_5ms_start = 4585500 ;
Rec_6ms_start = 5160500 ;
Rec_7ms_start = 5738500 ;

Rec_control2_start = 1720500 ;
Rec_8ms_start = 2875500 ;
Rec_9ms_start = 3454500 ;
Rec_10ms_start = 4034500 ;
Rec_15ms_start = 4619500 ;
Rec_20ms_start = 5209500 ;
Rec_25ms_start = 5804500 ;
Rec_30ms_start = 6404500 ;

Rec_control3_start = 1720500 ;
Rec_35ms_start = 2925500 ;
Rec_40ms_start = 3535500 ;
Rec_45ms_start = 4150500 ;
Rec_50ms_start = 4770500 ;
Rec_75ms_start = 5415500 ;
Rec_100ms_start = 6085500 ; 
Rec_0ms_start = 0 ; %EXCLUDE
%Rec_200ms_start = 5755500 %6265500 %5755500 ;

Rec_control4_start = 1720500 ;
Rec_200ms_start = 3160500 ;
Rec_500ms_start = 4230500 ;
Rec_1000ms_start = 5800500 ;
Rec_0ms_start = 0 ;
Rec_0ms_start = 0 ;
Rec_0ms_start = 0 ; 
Rec_0ms_start = 0 ; %EXCLUDE

Rec_table = table(Rec_control1_start,Rec_1ms_start,Rec_2ms_start,Rec_3ms_start,Rec_4ms_start,Rec_5ms_start,Rec_6ms_start,Rec_7ms_start,Rec_control2_start,Rec_8ms_start,Rec_9ms_start,Rec_10ms_start,Rec_15ms_start,Rec_20ms_start,Rec_25ms_start,Rec_30ms_start,Rec_control3_start,Rec_35ms_start,Rec_40ms_start,Rec_45ms_start,Rec_50ms_start,Rec_75ms_start,Rec_100ms_start,Rec_0ms_start,Rec_control4_start,Rec_200ms_start,Rec_500ms_start,Rec_1000ms_start,Rec_0ms_start,Rec_0ms_start,Rec_0ms_start,Rec_0ms_start) ;
%Rec_table.Properties.VariableNames = ["Rec_control1","Rec_1ms","Rec_2ms","Rec_3ms","Rec_4ms","Rec_5ms","Rec_6ms","Rec_7ms","Rec_control2","Rec_8ms","Rec_9ms","Rec_10ms","Rec_15ms","Rec_20ms","Rec_25ms","Rec_30ms","Rec_control3","Rec_35ms","Rec_40ms","Rec_45ms","Rec_50ms","Rec_75ms","Rec_100ms","Rec_200ms"]

%  QC
Rseal_lowerLim = 500000000 ;
Rseal_upperLim = Inf ;
Cslow_lowerLim = 2e-12 ;
Cslow_upperLim = 30e-12 ;
QCpA_uLim = 40e-12 ;
QCpA_lLim = -40e-12 ;
ControlPeakpA_Lim = -400e-12 ;
Leak_start = 1001 ; % 45ms after step to -120 from -80
Leak_end = 1900 ; % 5ms before end of leak correction step (i.e. mean over 50ms used to QC leak correction - 95 to 145ms)
QCUncorpA_uLim = 0 ;

% Import paths
Parent_dir = '\\Data\Vandenberg-Lab\Syncropatch\Functional_Genomics_Syncropatch\SCN5A\2.Control Mutants\' ;
Exp_dir = '20221128_JM1\' ;
Raw_dir = 'QCd_Wells\' ;
QC_dir = 'QC\' ;
Protocol_dir = ["NaRecInac1_10.47.31\","NaRecInac2_10.48.12\","NaRecInac3_10.48.53\","NaRecInac4_10.49.33\"] ;
%UncorrectedRaw_dir = 'Raw_NoLeakCorrection\' ;

% Export paths
folderdir = fullfile(Parent_dir,Exp_dir,Protocol_dir(1)) ;
mkdir(folderdir,'Recoveries')
recoveries_dir = 'Recoveries\' ;

folderdir = fullfile(Parent_dir,Exp_dir,Protocol_dir(1),recoveries_dir) ;
mkdir(folderdir,'RecoveryPlots')


% Import QC first, if it doesnt pass QC there is no point importing file

opts = delimitedTextImportOptions("NumVariables", 8);

% Specify range and delimiter
opts.DataLines = [3, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["WellID", "CellType", "SealResistance", "Capacitance", "SeriesResistance", "SealResistance1", "Capacitance1", "SeriesResistance1"];
opts.VariableTypes = ["string", "categorical", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "WellID", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["WellID", "CellType"], "EmptyFieldRule", "auto");

% Import the data
Raw_wells_4comb = [] ;
QC_wells_4comb = [] ;
for i = 1:4
    QC_wells_temp = dir(fullfile(Parent_dir,Exp_dir,Protocol_dir(i),QC_dir,'QC.csv')) ;
%    QC_wells_temp2 = QC_wells_temp(:,6:8)
    QC_wells_4comb = [QC_wells_4comb;QC_wells_temp] ;
    Raw_wells_temp = dir(fullfile(Parent_dir,Exp_dir,Protocol_dir(i),Raw_dir,'*.csv')) ;
    Raw_wells_4comb = [Raw_wells_4comb;Raw_wells_temp] ;
end


%   Import of QC sweeps data
for i = 1:4
    fid = fullfile(Parent_dir,Exp_dir,Protocol_dir(i),QC_dir,QC_wells_4comb(i).name);
    fidopen = fopen(fid);
    QCWells_data{1,i} = QC_wells_4comb(i).name ;
    QCWells_data{2,i} = readtable(fid,opts) ;
    fclose("all");
end  


% Clear variables
clear opts



% Filtering for wells which pass QC across ALL 4 protocols

for i = 1:4
    QCWells_temp = QCWells_data{2,i} ;
    QCWells_tempR = QCWells_temp.SealResistance > Rseal_lowerLim ;
    QCWells_temp = QCWells_temp(QCWells_tempR,:) ;
    QCWells_tempC = QCWells_temp.Capacitance > Cslow_lowerLim & QCWells_temp.Capacitance < Cslow_upperLim ;
    QCWells_temp = QCWells_temp(QCWells_tempC,:) ;
    QCWells_temp.Properties.RowNames = QCWells_temp.WellID ;
    QCWells_data{3,i} = QCWells_temp ;
end

QCWells_temp1 = QCWells_data{3,1} ;
QCWells_temp2 = QCWells_data{3,2} ;
QCWells_temp3 = QCWells_data{3,3} ;
QCWells_temp4 = QCWells_data{3,4} ;
QCWells_filter1_2 = QCWells_temp1(contains(QCWells_temp1.Properties.RowNames,QCWells_temp2.Properties.RowNames),:) ;
QCWells_filter1_2_3 = QCWells_filter1_2(contains(QCWells_filter1_2.Properties.RowNames,QCWells_temp3.Properties.RowNames),:) ;
QCWells_filter1_2_3_4 = QCWells_filter1_2_3(contains(QCWells_filter1_2_3.Properties.RowNames,QCWells_temp4.Properties.RowNames),:) 


%Boltz_filtered = Boltz_allWells(contains(Boltz_allWells.Properties.RowNames,PeaksAtMV_filterLow_Upp.Properties.RowNames),:) ;


% Extracting wellIDs

WellID_RawID1 = char(Raw_wells_4comb.name) ;
WellID_Raw = string(zeros(384*4,1)) ;
%20221117 n_wells = size(WellID_Raw) ;
n_wells = size(Raw_wells_4comb) ;

for i = 1:n_wells(1) 
    WellID_RawID2(i) = string(WellID_RawID1(i,25:27)) ; %25:27
end

WellID_RawID3 = WellID_RawID2' ;
%WellID_RawID2 = table(WellID_RawID2)

Raw_wells_4comb = struct2table(Raw_wells_4comb) ;
Raw_wells_4comb.WellIDs = WellID_RawID3 ;
Raw_wells_4_filtered = Raw_wells_4comb(contains(Raw_wells_4comb.WellIDs,QCWells_filter1_2_3_4.Properties.RowNames),:) ;
Raw_wells_4_filtered = table2struct(Raw_wells_4_filtered) ;



opts = delimitedTextImportOptions("NumVariables", 5);

% Specify range and delimiter
opts.DataLines = [4, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["SamplePoint", "SampleTimeus", "Var4", "RawDataA1"];
opts.SelectedVariableNames = ["SamplePoint", "SampleTimeus", "RawDataA1"];
opts.VariableTypes = ["double", "double", "string", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Var4", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Var4", "EmptyFieldRule", "auto");


% Import the data
n_pass = size(Raw_wells_4_filtered) ;

%   Import of Raw sweeps data
for irec = 1:4
    wellStart = n_pass/4*irec-n_pass/4+1 ;
    wellEnd = wellStart+n_pass/4-1 ;
    LabelRow = (irec*irec)-((irec-1)*(irec-1)) ;
    DataRow = LabelRow+1 ;
    for i = wellStart(1):wellEnd(1)
        fid = fullfile(Parent_dir,Exp_dir,Protocol_dir(irec),Raw_dir,Raw_wells_4_filtered(i).name);
        fidopen = fopen(fid);  
            RawWells_data{1,i} = Raw_wells_4_filtered(i).name ;    
            RawWells_data{2,i} = readtable(fid,opts) ;                 
        fclose("all");         
    end 
end  

clear opts


% Organising
% RawWells_data2 = cell2table(RawWells_data)
% RawWells_data3 = [RawWells_data2{1,1:n_pass(1)/3-1};RawWells_data2{2,1:n_pass(1)/3-1};RawWells_data2{1,n_pass(1)/3:n_pass(1)/3*2-1};RawWells_data2{2,n_pass(1)/3:n_pass(1)/3*2-1};RawWells_data2{1,n_pass(1)/3*2:n_pass(1)-1};RawWells_data2{2,n_pass(1)/3*2:n_pass(1)-1}]
n = n_pass/4 ;
RawWells_data2 = reshape(RawWells_data',n(1),[]) ;

% QC Leak

LeakCheck2 = LeakCheck/50+1
n_pass = size(RawWells_data2)
LeakQC = []

for iii = 1:4 % PLATE
    for i = 1:n_pass(1) %WELL
    RawWells_peaks(i,1) = Raw_wells_4_filtered(i).WellIDs ;
    RawWell_leakDataTemp = RawWells_data2{i,iii+4} ;
    RawWell_leakDataTemp2 = table2array(RawWell_leakDataTemp(LeakCheck2:LeakCheck2+100,3)) ;

    [mean_temp] = mean(RawWell_leakDataTemp2) ;
    LeakQC = [LeakQC;mean_temp] ;
    end
end

LeakQC = reshape(LeakQC,n_pass(1),[]) ;
LeakQC_tab = array2table(LeakQC) ;
LeakQC_tab.Properties.RowNames = RawWells_peaks(:,1)  ;

LeakQC_tabFilter1 = LeakQC_tab.LeakQC1 < QCpA_uLim & LeakQC_tab.LeakQC1 > QCpA_lLim ;
LeakQC_tabFilter1 = LeakQC_tab(LeakQC_tabFilter1,:) ;
LeakQC_tabFilter2 = LeakQC_tabFilter1.LeakQC2 < QCpA_uLim & LeakQC_tabFilter1.LeakQC2 > QCpA_lLim ;
LeakQC_tabFilter2 = LeakQC_tabFilter1(LeakQC_tabFilter2,:) ;
LeakQC_tabFilter3 = LeakQC_tabFilter2.LeakQC3 < QCpA_uLim & LeakQC_tabFilter2.LeakQC3 > QCpA_lLim ;
LeakQC_tabFilter3 = LeakQC_tabFilter2(LeakQC_tabFilter3,:) ;
LeakQC_tabFilter4 = LeakQC_tabFilter3.LeakQC4 < QCpA_uLim & LeakQC_tabFilter3.LeakQC4 > QCpA_lLim ;
LeakQC_tabFilter4 = LeakQC_tabFilter3(LeakQC_tabFilter4,:) ;

% Extracting values at times of interest per recovery
% Final layout of cell - variables: wellID, control1, 1ms, 2ms,
% 3ms, ....7ms
% Final layout of cell - rows = wells that have passed QC

Rec_table2 = reshape(table2cell(Rec_table),8,4) ;
n_rec = size(Rec_table2) ;
n_pass = size(RawWells_data2) ;

for iii = 1:4 % PLATE
    for i = 1:n_pass(1) %WELL
    RawWells_peaks(i,1) = Raw_wells_4_filtered(i).WellIDs ;
    RawWell_temp = RawWells_data2{i,iii+4} ;
    for ii = 1:n_rec(1) % RECOVERY 
        R_start = cell2mat(Rec_table2(ii,iii)) ;
        R_start_idx = R_start/50+1 ;
        R_end_idx = R_start_idx+200 ;
        Rec_range = RawWell_temp(R_start_idx:R_end_idx,3) ;
        Rec_min = min(table2array(Rec_range)) ;
        RawWells_peaks(i,ii+8*(iii-1)+1) = Rec_min ;        
    end
    end 
end


% Normalising to control step and conversion to table

s = size(RawWells_peaks)
%for i = 1:8
RawWells_peaksTab = array2table(str2double(RawWells_peaks(:,2:end))) ;
RawWells_peaksTab.Properties.RowNames = RawWells_peaks(:,1) ;

RawWells_peaksTab_PeaksQC1 = RawWells_peaksTab.Var1 < ControlPeakpA_Lim  ;
RawWells_peaksTab_PeaksQC1 = RawWells_peaksTab(RawWells_peaksTab_PeaksQC1,:) ;
RawWells_peaksTab_PeaksQC2 = RawWells_peaksTab_PeaksQC1.Var9 < ControlPeakpA_Lim  ;
RawWells_peaksTab_PeaksQC2 = RawWells_peaksTab_PeaksQC1(RawWells_peaksTab_PeaksQC2,:) ;
RawWells_peaksTab_PeaksQC3 = RawWells_peaksTab_PeaksQC2.Var17 < ControlPeakpA_Lim  ;
RawWells_peaksTab_PeaksQC3 = RawWells_peaksTab_PeaksQC2(RawWells_peaksTab_PeaksQC3,:) ;
RawWells_peaksTab_PeaksQC4 = RawWells_peaksTab_PeaksQC3.Var25 < ControlPeakpA_Lim  ;
RawWells_peaksTab_PeaksQC4 = RawWells_peaksTab_PeaksQC3(RawWells_peaksTab_PeaksQC4,:) ;

RawWells_peaksArr = table2array(RawWells_peaksTab_PeaksQC4) ;
Peaks_P1_Norm = RawWells_peaksArr(:,2:8)./RawWells_peaksArr(:,1) ;
Peaks_P2_Norm = RawWells_peaksArr(:,10:16)./RawWells_peaksArr(:,9) ;
Peaks_P3_Norm = RawWells_peaksArr(:,18:24)./RawWells_peaksArr(:,17) ;
Peaks_P4_Norm = RawWells_peaksArr(:,26:32)./RawWells_peaksArr(:,25) ;
Peaks_Pall_Norm = [Peaks_P1_Norm,Peaks_P2_Norm,Peaks_P3_Norm,Peaks_P4_Norm] ;
Peaks_Pall_NormTab = array2table(Peaks_Pall_Norm) ;
Peaks_Pall_NormTab.Properties.RowNames = RawWells_peaksTab_PeaksQC4.Properties.RowNames ;


Peaks_Pall_NormTab.Peaks_Pall_Norm21 = [] ;
Peaks_Pall_NormTab.Peaks_Pall_Norm25 = [] ;
Peaks_Pall_NormTab.Peaks_Pall_Norm26 = [] ;
Peaks_Pall_NormTab.Peaks_Pall_Norm27 = [] ;
Peaks_Pall_NormTab.Peaks_Pall_Norm28 = [] ;


% Filtering bad wells - when each recovery peak is the ~same size as
% controls, the norm data would show ~1. This step is to remove those bad data.

% UPDATED: also filter out the wells not passing leak correction check

Peaks_Pall_mean = mean(Peaks_Pall_Norm')' ;
Peaks_Pall_mean = array2table(Peaks_Pall_mean) ;
Peaks_Pall_mean.Properties.RowNames = RawWells_peaksTab_PeaksQC4.Properties.RowNames 
Peaks_Pall_meanF = Peaks_Pall_mean.Peaks_Pall_mean < 1 ;
Peaks_Pall_meanF = Peaks_Pall_mean(Peaks_Pall_meanF,:) ;

Peaks_Pall_filterM = Peaks_Pall_NormTab(contains(Peaks_Pall_NormTab.Properties.RowNames,Peaks_Pall_meanF.Properties.RowNames),:) ;
Peaks_Pall_filterM2 = Peaks_Pall_filterM(contains(Peaks_Pall_filterM.Properties.RowNames,LeakQC_tabFilter3.Properties.RowNames),:) ;


% cant use slope to filter in case there is a one-off rec that is out of
% place



time = [1,2,3,4,5,6,7,8,9,10,15,20,25,30,35,40,45,50,75,100,200,500,1000]
time = time'
%SpanFast = '(Plateau-Y0)*PercentFast*.01'
%SpanSlow = '(Plateau-Y0)*(100-PercentFast)*.01'
%Y = 'Y0+ SpanFast*(1-exp(-KFast*X)) + SpanSlow*(1-exp(-KSlow*X))'
%singleExpY=Y0 + (Plateau-Y0)*(1-exp(-K*x))
  

DXfitdata_KFast = [] ;
DXfitdata_KSlow = [] ;
DXfitdata_TauFast = [] ;
DXfitdata_TauSlow = [] ;
DXfitdata_HalfTimes = []
DXfitdata_adjR2 = [] ;

SXfitdata_K = [] ;
SXfitdata_Tau = [] ;
SXfitdata_HalfTimes = []
SXfitdata_adjR2 = [] ;

n_wells = size(Peaks_Pall_filterM2)
for wells = 1:1:n_wells(1)
    %try
    y1 = table2array(Peaks_Pall_filterM2(wells,:)) ;
    y1 = y1' 
    
    Y = fittype('Y0+ ((Plateau-Y0)*(100-PercentFast)*.01)*(1-exp(-KSlow*x)) + ((Plateau-Y0)*PercentFast*.01)*(1-exp(-KFast*x))') ;
    [DXfittedcurve,DXfitdata] = fit(time,y1,Y,'StartPoint',[0.5,0.05,1,0,0],'lower',[0,0,0,0,0], 'upper',[inf,inf,inf,inf,0]) 

    DXci = confint(DXfittedcurve,0.95) % some cells have unstable calculations (esp BrS variants with extemely long recoveries that cannot be plotted accurately. Matlab still gives a value but the CI will show as inf. This is added so that if inf is given anyywhere, then a single exp is used isntead.
    DXcii = find(DXci ==inf)
    
    singleY = fittype('Y0 + (Plateau-Y0)*(1-exp(-K*x))') 
    [SXfittedcurve,SXfitdata] = fit(time,y1,singleY,'StartPoint',[0.05,1,0],'lower',[0,0,0], 'upper',[inf,inf,0]) 

%    if DXcii >=1
    SXfitdata_K = [SXfitdata_K; SXfittedcurve.K] ;
    SXfitdata_adjR2 = [SXfitdata_adjR2; SXfitdata.adjrsquare] ;
    Tau = 1/SXfittedcurve.K ;
    SXfitdata_Tau = [SXfitdata_Tau; Tau] ;
    SXHalftime = (log(2))/(SXfittedcurve.K) ;
    SXfitdata_HalfTimes = [SXfitdata_HalfTimes;SXHalftime] ;

%    else 
    DXfitdata_KFast = [DXfitdata_KFast; DXfittedcurve.KFast] ;
    DXfitdata_KSlow = [DXfitdata_KSlow; DXfittedcurve.KSlow] ;
    DXfitdata_adjR2 = [DXfitdata_adjR2; DXfitdata.adjrsquare] ;
    TauFast = 1/DXfittedcurve.KFast ;
    TauSlow = 1/DXfittedcurve.KSlow ;
    DXfitdata_TauFast = [DXfitdata_TauFast; TauFast] ;
    DXfitdata_TauSlow = [DXfitdata_TauSlow; TauSlow] ;
    DXHalftime = (log(2)/(DXfittedcurve.KFast * (DXfittedcurve.PercentFast/100))) + (log(2)/(DXfittedcurve.KSlow * (1-DXfittedcurve.PercentFast/100))) ;
    DXfitdata_HalfTimes = [DXfitdata_HalfTimes;DXHalftime] ;
%    end




    clf
    plot(DXfittedcurve,'c',time,y1,'*k') 
    hold on
    plot(SXfittedcurve,'m') 
    
    set(gca,'Xscale','log') ;
    xlabel("Time (ms)") ;
    ylabel("Normalised Recovery from Inactivation") ;
    yline(0.5,'--') ;
    ylim([0,1.1]) ;
    yticks(0:0.5:1) ;
    titleWell = string(Peaks_Pall_filterM2.Properties.RowNames(wells)) ; 
    titleWell = strcat("Recovery from inactivation ",titleWell,".jpg") ;
    title(titleWell);
    
    legend("Recovery Data", "Double Exponential","Single Exponential","Location","southeast") ;
    legend("boxoff") ;
    hold off

    pid = fullfile(Parent_dir,Exp_dir,Protocol_dir(1),recoveries_dir,'RecoveryPlots\',titleWell) ;
    exportgraphics(gca,pid)



 % Try and catch eventually removed. if cannot plot on double exp, it
 % should plot on single.
  %  catch % try and catch is employed so that if a dataset cannot be fitted (bad cell resulting in poor values, e.g. poor initial values), the data is replaced with 0 so to not disrupt #wells and removed later
  %  DXfitdata_KFast = [DXfitdata_KFast; 0] ;
  %  DXfitdata_KSlow = [DXfitdata_KSlow; 0] ;
  %  DXfitdata_adjR2 = [DXfitdata_adjR2; 0] ;

  %  DXfitdata_TauFast = [DXfitdata_TauFast; 0] ;
  %  DXfitdata_TauSlow = [DXfitdata_TauSlow; 0] ;

  %  DXfitdata_HalfTimes = [DXfitdata_HalfTimes;0] ;

  %  end

end


   


%Turning into tables for tracking row names and variant IDs
alpcomp = [] ;
numcomp = [] ;
for i = 1:n_wells(1)
wellIDs = char(Peaks_Pall_filterM2.Properties.RowNames(i)) ;
alp = wellIDs(1) ;
num = wellIDs(2:3) ;
alpcomp = [alpcomp;alp] ;
numcomp = [numcomp;num] ;
end

DXfitfull = table(DXfitdata_KFast) ;
DXfitfull.dXKSlow = DXfitdata_KSlow ;
DXfitfull.dXTauFast = DXfitdata_TauFast ;
DXfitfull.dXTauSlow = DXfitdata_TauSlow ;
DXfitfull.dXHalfTime = DXfitdata_HalfTimes ;
DXfitfull.dXR2 = DXfitdata_adjR2 ;

DXfitfull.sXK = SXfitdata_K ;
DXfitfull.sXTau = SXfitdata_Tau ;
DXfitfull.sXHalfTime = SXfitdata_HalfTimes ;
DXfitfull.sXR2 = SXfitdata_adjR2 ;


DXfitfull.Properties.RowNames = Peaks_Pall_filterM2.Properties.RowNames
DXfitfull.WellID1 = alpcomp ;
DXfitfull.WellID2 = numcomp ;

DXfitfull = sortrows(DXfitfull,[12,11],{'ascend','ascend'})
DXfitfull(:,11:12) = []

toremove = find(DXfitfull.dXR2 <= 0.9) % remove failed/bad fit wells
DXfitfull(toremove,:) = []

% Add Variant ID

QCWells_filter1_2_3_4 = QCWells_filter1_2_3_4(contains(QCWells_filter1_2_3_4.Properties.RowNames,DXfitfull.Properties.RowNames),:)
DXfitfull.Variant = QCWells_filter1_2_3_4.CellType




DXfitfull_summ = groupsummary(DXfitfull,"Variant","mean")
DXfitfull_summ.Properties.RowNames = string(DXfitfull_summ.Variant)
HalfTimeNORM = DXfitfull_summ{'SCN5A_DelQ1077_WT', 'mean_sXHalfTime'}

DXfitfull.sXHTsubWTmean = DXfitfull.sXHalfTime-HalfTimeNORM



%Exporting
DXfitfull = movevars(DXfitfull,"Variant","Before","DXfitdata_KFast")
DXfitfull = movevars(DXfitfull,"sXHalfTime","After","Variant")
DXfitfull = movevars(DXfitfull,"sXHTsubWTmean","After","sXHalfTime")


writetable(Peaks_Pall_filterM2,fullfile(folderdir,'NormalisedRecovery.csv'),"WriteRowNames",true,"WriteVariableNames",true) ;
writetable(DXfitfull,fullfile(folderdir,'HT_KF_KS_TF_TS.csv'),"WriteRowNames",true,"WriteVariableNames",true) ;
writetable(DXfitfull_summ,fullfile(folderdir,'Summary.csv'),"WriteRowNames",true,"WriteVariableNames",true) ;

clear

