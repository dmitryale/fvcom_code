%function[]=mk_vect
% subplot(2,1,2);
             scaleu=2.0;
             cvl=['cm\cdots^{-1}'];
clear uu3 vv3 xx3 yy3   uu2 vv2 xx2 ntx2 yy2  xx yy vv uu ntx scit hu1 hu2 xlimx;
% xx=mto;   uu=uo; vv=vo;   %obs
 xx=mtdm;   uu=udo; vv=vdo; %obs resid
ntx=length(xx);yy=zeros(ntx,1);      
%% add the scale 
                   scit= -fix(ntx*1/4)+1; %-42;
                   if scit>0, scit=0; end;
xx(ntx+1)=xx(ntx + scit);         
xx(ntx+2)=xx(ntx + scit);
yy(ntx+1)=yy(ntx + scit)+scaleu*2; yy(ntx+2)=yy(ntx + scit) + scaleu*2;
uu(ntx+1)=scaleu;     vv(ntx+1)=0; uu(ntx+2)=0; vv(ntx+2)=scaleu;

 %xx2=mtm;  uu2=um;  vv2=vm ; %model
 xx2=mtdm; uu2=udm; vv2=vdm; %model resid
 ntx2=length(xx2);
 yy2=zeros(ntx2,1);  % yy3=zeros(ntx,1);         % yy4=zeros(length(go4),1);
 xx2(ntx2+1)=xx2(ntx2 + scit);          xx2(ntx2+2)=xx2(ntx2 + scit);
 yy2(ntx2+1)=yy2(ntx2 + scit)+scaleu*2; yy2(ntx2+2)=yy2(ntx2 + scit) + scaleu*2; 
 uu2(ntx2+1)=scaleu;   vv2(ntx2+1)=0;   uu2(ntx2+2)=0; vv2(ntx2+2)=scaleu;
 
     rlt=mean(lat); vp=vv*cosd(rlt);
  %xx3=[xx,xx2]; yy3=[yy,yy2]; uu3=[uu,uu2]; vv3=[vv,vv2];
 if 0,...
       hu1=quiver(xx ,yy-0        ,uu ,vv ,1,'b','filled');   hold on;
       hu2=quiver(xx2,yy2-scaleu*0,uu2,vv2,1,'g','filled');   hold on;
  else
     clear dx dy offset bp bp2 dx2 dy2 offset2 %hp hp2 hp1 
     scale = scaleu; dts=1;
     bp= [1:length(uu)]';    dx=dts*uu';   dy=dts*vv';   offset=[dx;dy]';
     bp2=[1:length(uu2)]';  dx2=dts*uu2'; dy2=dts*vv2'; offset2=[dx2;dy2]';    
     
     [hp ]=arrow([xx,yy],([xx,yy]+offset(bp,:)),...
          'length',10,'tipangle',15,'width',1.9,'BaseAngle',30); hold on;
       set(hp,'FaceColor','b'); set(hp,'EdgeColor','k');
  
  [hp2 ]=arrow([xx2,yy2],([xx2,yy2]+offset2(bp2,:)),...
          'length',10,'tipangle',15,'width',2.0,'BaseAngle',30);     hold on;
    set(hp2,'FaceColor','g');set(hp2,'EdgeColor','k');
    set(hp ,'FaceColor','b');set(hp ,'EdgeColor','k');

       end 
                   xlmtx=[min(xx)-1 max(xx)+ scaleu*3]; 
    set(gca,'xlim',xlmtx);        
   datetick('x','keeplimits');
%  datetick('x','m','keeplimits');
  set(gca,'YMinorTick','on','XMinorTick','on');
  set(gca,'TickDir','out'); %in
                                 fcb=['\color{blue} ' ]; % fc=['\color[rgb]{.7 .7 .7} \sl' ]
                                 fcg=['\color{green} '];
                                 fcr=['\color{red} '  ];
  dbe=[   datestr(xx(1),'yyyymmdd') ...
      '-' datestr(xx(end,1),'yyyymmdd') ];
    cid(1,:) = ['o' num2str(jp)];
    cid(2,:) = ['m' num2str(jp)];
text(xx(ntx+1)+1,yy(ntx+1)+scaleu*1.25,...
    {[''];[ ['\color{black} ' num2str(scaleu) ' ' cvl ...
           ':' fcb cid(1,:) ',' fcg cid(2,:) '.' ] ... ',' fcr  cid(3,:) '. '] ...
       ] });  
    

   set(gca,'dataaspectratio',[1 1 1]);
   set(gca,'position',[0.10 0.07  0.84  0.41]);
   set(gca,'yTicklabel','');
   ylmt=get(gca,'ylim');      
        set(gca,'ylim',ylmt*0.75);   
        