%function[]=mk_tseries_nc(); 
% (c) dmitry.aleynik@sams.ac.uk 2013.12.06; 2014.12.21 2016.03.10         & 
%  __o_O__¬                                                               &
%  \_____/ ~~~~~~<@})))< ~ ~ ~~~~~ ~~~ SAMS MIDAS CCZ; ETIVE27 s25,s22    &
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
close all; clear all;
  co2='n'; % co2='y';
this_dir=pwd;

DRV=this_dir(1:1);
addpath ([DRV ':\HECTOR\matlab\matlab\']);
addpath ([DRV ':\HECTOR\matlab\']);
addpath( ('M:\Mar_Phys\matlab\seawater'));
addpath( (['M:\Mar_Phys\matlab\m_map\']));
pause off
warning off
        %mk_mat='N';  elim1=[-1.9 1.9] ;
         mk_mat='Y';  elim1=[-1.9 1.9] ;

  mkgif     ='Y'; gi=0; gim=0; % mkgif='N'; gi=0;
          slayer='tanh';
 %%-----------------------------------------------------------
% set directories default values
 %            caseN ='26'; 
              caseN ='27'; 
  casename=['etive' caseN];  turb = 'y';co2='n';
 files_FVCOM =[casename '_' '0001.nc'];   turb = 'y';co2='n'; crun='v0/';     
 [sigma_levs,zl,meshdir,sms_model,FVCOM_data_dir,FVCOM_mat_dir,FVCOM_plot_dir] ...
       = set_vars_sms(caseN) ;
%FVCOM_data_dir =[ 'D:/ARCHER/etive27/OUT_2013_s25/']; crun='v1/' ;        
 FVCOM_data_dir =[ 'D:/ARCHER/etive27/OUT_2013_s22/']; crun='v2/' ;        
 matlab_FVCOM =[ DRV ':/HECTOR/matlab/matlab/']; 

% params_opts={'time','data_dir','files','nzopt','trnsotp','varnames','trns_idx','nz_idx'};
 time_offset = 678942; % from FVCOM time to matlab time

[IM_range_lat ,IM_range_lon] = deal([56.430  56.568],[ -5.529   -5.05]); cm=''   ; %all
 x_vec=-5-14/60; y_vec=56+27.25/60;  
 base_year = datenum(2013,1,1,0,0,0);  
 date_range={'01/01/2009 00:00:00','31/12/2037 23:00:00'}; % 1 if all available data wanted
 date_day  =datestr(date_range(1),'yymmdd');
 fig_name = [date_day,'_tide_',casename];
                                  pathout=[casename '_' crun ];             
  FVCOM_plot_dirg=[FVCOM_plot_dir pathout 'Probe_g' '/' ];      %  
 [~,~,messid] = mkdir (FVCOM_plot_dirg) ;
  try load([FVCOM_plot_dir pathout 'mesh.mat'],'-mat'); end;
% try load([FVCOM_plot_dir pathout 'CVM_090101.mat'],'-mat'); end
           
fts_dr=dir( [FVCOM_data_dir [casename '*station_timeseries.nc']]);
fnam1=[ fts_dr(1,1).name];
[finfo outstrct]=read_nc_file_struct([FVCOM_data_dir fnam1]);
FV_TS = outstrct; 

      clear outstrct sigma_levs finfo params_opts ;
FV_TS.name=FV_TS.name_station';
FV_TS.mtime=double(FV_TS.time)+time_offset;

xtm=double(FV_TS.mtime); mtime=xtm;

path_fig=FVCOM_plot_dirg;
 cvars={'x'; 'y';  'lon';  'lat';  'siglay';  'siglev';  'h';  'time';  'iint'; ...
'u'; 'v'; 'ww'; 'ua'; 'va'; 'temp'; 'salinity'; 'zeta'; ...
'uwind_speed'; 'vwind_speed'; 'name'; 'mtime';};

if 0,
% plot(FV_TS.lon,FV_TS.lat,'s');hold on
vrm= double(squeeze(FV_TS.(char('zeta'))(:,:)));
xtm=double(FV_TS.mtime); mtime=xtm;
plot(xtm, vrm ); hold on;
  set(gca,'YMinorTick','on','XMinorTick','on');
  set(gca,'TickDir','out');     
                  xlmt=[min(xtm) max(xtm)];
  set(gca,'xlim',[xlmt]);
  datetick('x','keeplimits');
  ylabel('zeta, m ')
end

        [ni,nk]=size(FV_TS.x);
         niss=[29, 20,          1,    25];
        % Ganavan, Dunstafnage, Hypox,RE5
          % t, s, el,h, u, v,ww,ua,va;
    cvi  = [15,16,17,7,10,11,12,13,15];
    cunits={'^oC';'psu';'m';'m';'ms^{-1}';'ms^{-1}';'ms^{-1}';'ms^{-1}';'ms^{-1}';};
for vai=1:2,...
 %  vai=1;%t
%     vai=2;%s
    civ  =cvi(vai); 
    cunit=char(cunits(vai));
    
 for iis=1:length(niss),...
 clear nis vrm xlmt db de ylmt stnm dv  yearsm yearsYS yearsc ylmdif
       figure (100+iis); clf
   nis = niss(iis); % 
   cvar= char(cvars(civ));
     
vrm = double(squeeze(FV_TS.(char(cvar))(nis,:,:)));

plot(xtm, vrm ); hold on;
  set(gca,'YMinorTick','on','XMinorTick','on');
  set(gca,'TickDir','out'); %in    
                  xlmt=[min(xtm) max(xtm)];
  set(gca,'xlim',[xlmt]);
  datetick('x','keeplimits');
  ylabel( [cvar ',' cunit ]);
   db=datestr(mtime(1)   , 'yyyymmdd');  
   de=datestr(mtime(end) , 'yyyymmdd'); 
          stnm =strtrim(FV_TS.name(nis,:));
   title([stnm ' ' cvar ' : ' db ' ' de ]);
   
 ylmt=get(gca,'YLim'); ylmdif=diff(ylmt)/10;
     dv=datevec(xlmt); years=[dv(1):1:dv(2)]; nya=length(years);
 yearsm=       (datenum(num2str(years(:)),'YYYY'));
 yearsm(1)=xlmt(1);
 yearsc=datestr(datenum(num2str(years(:)),'YYYY') ,'yyyy');
        yearsYS=[zeros(length(yearsm),1)]+(ylmt(1)+ylmdif*0.25)*1.0;
 text(yearsm+1, yearsYS' ,yearsc ); hold on;
 
namepng =[  path_fig ['tmser_' cvar '_' stnm '_'  db '-' de ] ];    
print(['-f'],'-dpng','-loose','-r500',[namepng  ,'.png']) ;    

 end
end
 save([FVCOM_plot_dir pathout 'FV_TS_'  db '-' de '.mat' ],'FV_TS','-v7.3');