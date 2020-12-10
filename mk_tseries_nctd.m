%function[]=mk_tseries_nctd(); 
% (c) dmitry.aleynik@sams.ac.uk 2013.12.06; 2014.12.21 2016.03.10         & 
%  __o_O__¬                                                               &
%  \_____/ ~~~~~~<@})))< ~ ~ ~~~~~ ~~~ SAMS MIDAS CCZ; ETIVE27 s25,s22    &
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
close all; clear all;
  co2='n'; % co2='y';
this_dir=pwd;

DRV=this_dir(1:1);   DRVd='D';
if strcmpi(DRV,'c'); DRVd='C';end
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
%FVCOM_data_dir=[DRVd ':/ARCHER/etive27/OUT_2013_s25/']; crun='v1/' ;        
%FVCOM_data_dir=[DRVd ':/ARCHER/etive27/OUT_2013_s22/']; crun='v2/' ;        
%FVCOM_data_dir=[DRVd ':/ARCHER/etive27/OUT_2013_s10/']; crun='v3/' ;% NO mean flow
%FVCOM_data_dir=[DRV ':/samhanach/sa01da/work/etive27/OUT_2013_s20mf/']; crun='v20/' ;        
%FVCOM_data_dir=[DRV ':/samhanach/sa01da/work/etive27/OUT_2013_s10mf/']; crun='v4/' ;        
%   FVCOM_data_dir =[ DRV ':/samhanach/sa01da/work/'...
%         'etive27/OUT_2013_s10mf/']; crun='v5/';  turb = 'y';co2='n';

% FVCOM_data_dir=[DRVd ':/ARCHER/etive27/OUT_2014_a/']; crun='2014/' ;% NO mean flow
% FVCOM_data_dir=[DRVd ':/ARCHER/etive27/OUT_2014_b/']; crun='2014_b/' ;% NO mean flow
% FVCOM_data_dir=[DRVd ':/ARCHER/etive27/OUT_2014_c/']; crun='2014_c/' ;% NO mean flow
% FVCOM_data_dir=[DRVd ':/ARCHER/etive27/OUT_2015_as10/']; crun='2015_as10/'   ;%s10

% FVCOM_data_dir=[DRVd ':/ARCHER/etive27/OUT_2015_a/']; crun='2015_a/' ;%s20
% FVCOM_data_dir=[DRVd ':/ARCHER/etive27/OUT_2015_b/']; crun='2015_b/' ;%s10
% FVCOM_data_dir=[DRVd ':/ARCHER/etive27/OUT_2015_c/']; crun='2015_c/' ;%s10
%   FVCOM_data_dir=[DRVd ':/ARCHER/etive27/OUT_2015_d/']; crun='2015_d/' ;%s20

 FVCOM_data_dir=[DRVd ':/ARCHER/etive27/OUT_2014_e/']; crun='2014_e/' ;%s20
%   FVCOM_data_dir=[DRVd ':/ARCHER/etive27/OUT_2016_c/']; crun='2016_c/' ;%s20

  
%FVCOM_data_dir=[DRV ':/samhanach/sa01da/work/etive27/OUT_recent/'];     crun='v4/' ;        
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
  path_fig=FVCOM_plot_dirg;

if 1,...           
 fts_dr=dir( [FVCOM_data_dir [casename '*station_timeserie*.nc']]);
 [fn,fi]=size(fts_dr);
 for fi=1:fn,...
clear FV_TS fnam1
fnam1=[fts_dr(fi,1).name];
[finfo outstrct]=read_nc_file_struct([FVCOM_data_dir fnam1]);

FV_TS = outstrct; 

      clear outstrct sigma_levs finfo params_opts  mtime xtm Date;
FV_TS.name=FV_TS.name_station';
     tmstep=0.15*10; %sec
     FV_TS.iint = double(FV_TS.iint);
             mtime1=datenum(double(fix(FV_TS.time(1)))+time_offset);  
  FV_TS.Date=mtime1+(FV_TS.iint-FV_TS.iint(1))*tmstep/(24*3600); 
% iter number=2016000*tmstep/(24*3600)-> day from year start
      datestr(   FV_TS.Date(1:3),31)      
%--------------
FV_TS.mtime = FV_TS.Date;
FV_TS.mtime1= double(FV_TS.time)+time_offset;
        xtm=double(FV_TS.mtime); 
mtime=xtm;
 cvars={'x'; 'y';  'lon';  'lat';  'siglay';  'siglev';  'h';  'time';  'iint'; ...
'u'; 'v'; 'ww'; 'ua'; 'va'; 'temp'; 'salinity'; 'zeta'; ...
'uwind_speed'; 'vwind_speed'; 'name'; 'mtime';};
FV_TS.cvars=cvars;
cvarsU=cvars(10:19);
   
if fi==1,
   FV_TS0=FV_TS;
          clear FV_TS;
   mtime0=FV_TS0.mtime;
else
    go=[]; go=find(mtime>mtime0(end));
    FV_TS0.mtime=cat(1,FV_TS0.mtime,mtime(go));
    for ii=1:length(cvarsU),...
            cva=char(cvarsU(ii));
    sa=[]; a=[]; a=FV_TS.([cva]);  sa=size(a);
     if length(sa)<=2,
      FV_TS0.([cva])=cat(2,FV_TS0.([cva]),a(:,go));
     end
     
     if length(sa)==3,
      FV_TS0.([cva])=cat(3,FV_TS0.([cva]),a(:,:,go));
     end
     
     if length(sa)==4,
      FV_TS0.([cva])=cat(4,FV_TS0.([cva]),a(:,:,:,:,go));
     end
    
    end   
end

end
FV_TS0.cvarsU=cvarsU;
FV_TS=FV_TS0;
      xtm=double(FV_TS.mtime); 
   db=datestr(FV_TS.mtime(1)   , 'yyyymmdd');  
   de=datestr(FV_TS.mtime(end) , 'yyyymmdd'); 
save([FVCOM_plot_dir pathout 'FV_TS_'  db '-' de '.mat' ],'FV_TS','-v7.3');
else
   dir_ma=dir([FVCOM_plot_dir pathout 'FV_TS_*' '.mat' ]);
         load([FVCOM_plot_dir pathout dir_ma(end,1).name ]);

end
 cvars=FV_TS.cvars;
 xtm=double(FV_TS.mtime); 
 mtime=xtm; 

%% mk_figures ================
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
%
datejb=FV_TS.mtime(1);
dateje=FV_TS.mtime(end);
datej= 0.5*(datejb+dateje);

check_Srivs;

check_CTDs;
    %  varn={'z';'t';'s';'r';'o';'pr'};
%        CTD_Prof.(char(vr));   vr=varn(iv); 

        [ni,nk]=size(FV_TS.x);
         niss=[29, 20,          1,    25];
        % Ganavan, Dunstafnage, Hypox,RE5
          % t, s, el,h, u, v,ww,ua,va;
    cvi  = [15,16,17,7,10,11,12,13,15];
    cunits={'^oC';'psu';'m';'m';'ms^{-1}';'ms^{-1}';'ms^{-1}';'ms^{-1}';'ms^{-1}';};
for vai=1:2,...
 %  vai=1;%t
%   vai=2;%s
    civ  =cvi(vai); 
    cunit=char(cunits(vai));
        CTDmt=CTD_Prof.mtime;
         CTDz=CTD_Prof.z;
         CTDv=CTD_Prof.(char(varn(vai+1)));
          [ctz,cts,ctm]=size(CTDv);
          ctz10=10; %cz=23 m
 for iis=1:length(niss),...
 clear nis vrm xlmt db de ylmt stnm dv  yearsm yearsYS yearsc ylmdif
       figure (100+iis); clf
   nis = niss(iis); % 
   cvar= char(cvars(civ));
     
vrm = double(squeeze(FV_TS.(char(cvar))(nis,:,:)));

plot(xtm, vrm ); hold on;
  set(gca,'YMinorTick','on','XMinorTick','on');
  set(gca,'TickDir','out'); %in    
                  xlmt=[min(xtm)-0.0 max(xtm)];
  set(gca,'xlim',[xlmt]);                

  ylabel( [cvar ',' cunit ]);
   db=datestr(mtime(1)   , 'yyyymmdd');  
   de=datestr(mtime(end) , 'yyyymmdd'); 
          stnm =strtrim(FV_TS.name(nis,:));
   slash=1; slash=strfind(FVCOM_data_dir,'/');
   title([stnm ' ' cvar ' : ' db ' ' de ',' FVCOM_data_dir(slash(end-1):end-1) ] );
   
 ylmt=get(gca,'YLim'); ylmdif=diff(ylmt)/10;
     dv=datevec(xlmt); years=[dv(1):1:dv(2)]; nya=length(years);
 yearsm=       (datenum(num2str(years(:)),'YYYY'));
 yearsm(1)=xlmt(1);
 yearsc=datestr(datenum(num2str(years(:)),'YYYY') ,'yyyy');
        yearsYS=[zeros(length(yearsm),1)]+(ylmt(1)+ylmdif*0.25)*1.0;
 text(yearsm+1, yearsYS' ,yearsc ); hold on;
 
         for itc=1:ctm,
          if (CTDmt(itc) >= min(xlmt) & CTDmt(itc) <= max(xlmt) ),
            ctms= repmat(CTDmt(itc),ctz-ctz10+1,cts-1);
            ctv = squeeze(CTDv(ctz10:ctz,1:cts-1,itc));
          plot(ctms,ctv,'s'); hold on;
          end
         end
if vai==2, 
    gmt=[];gmt=find(Srivs.mt >= min(xlmt) & Srivs.mt <= max(xlmt) );
    gmtd=[];gmtd=find(Srivs.mtmd >= min(xlmt) & Srivs.mtmd <= max(xlmt) );
    plot(Srivs.mt(gmt)   ,Srivs.S(gmt)+Srivs.dSre5,'k-');hold on; 
    plot(Srivs.mtmd(gmtd),Srivs.Sid(gmtd)+Srivs.dSre5,'k--');hold on; 
end
% datetick('x','m','keeplimits');
  datetick('x','m');
  set(gca,'xlim',[xlmt]);
if abs(diff([xlmt]))< 30,...
  datetick('x','keeplimits');
end

set(gcf,'color','w'); % set(gcf,'color','none')
set(gcf,'PaperPositionMode','auto')

namepng =[  path_fig ['tmser_' cvar '_' stnm '_'  db '-' de ] ];    
print(['-f'],'-dpng','-loose','-r500',[namepng  ,'.png']) ;    

 end
end

if 1, 
    
%     iye=datevec(max(xlmt));iye=iye(1)+0;
% set(gca,'xlim',[xlmt(1) datenum(num2str(iye), 'yyyy') ]);
set(gca,'xlim',[xlmt  ]);
xlm1=get(gca,'xlim'); 
   if vai==2, set(gca,'ylim',[24 31]); end;
  datetick('x','m');
% datetick('x','keeplimits');
  set(gca,'xlim',[xlm1]);

set(gcf,'PaperPositionMode','auto')

namepngx =[  path_fig ['tmser_' cvar '_' stnm '_'  db '-' datestr(xlm1(2),'yymmdd') ] ];    
print(['-f'],'-dpng','-loose','-r500',[namepngx  'x','.png']) ;    
end
