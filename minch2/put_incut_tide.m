% function[]=put_incut_tide(}
%(c)dmitry.aleynik@sams.ac.uk, 2013.05.15)
         clear kt tml tb te jw mt;
%  elim1=[-1.99 1.99] ;    
   elim1=[-2.75 2.75] ;    
   npob = np2;%oban
%  FVCOM=TE;
%% ----------------------- %set incut:
  posii1=get(gca,'Position'); 
  gca1=get(gca);
      ax1=gca;
  set(ax1,'XColor','k','YColor','k');     hold on;  
  posii2=get(gca,'Position');   %get(ax2,'Position');
%       posii2 =[ 0.21 0.21 0.26 0.072]; %exch2
        posii2 =[ 0.38 0.21 0.33 0.078];
  ax2 = axes('Position',posii2, 'XAxisLocation','top', 'YAxisLocation','right',...
             'Color','none',   'XColor','g','YColor','g');
        el_min=min(min(FVCOM.el(:,:)));    
        nl= length(FVCOM.el(:,1));
       ntd=  length(Date);
%        if it1 >0 & ntimes > 0, 
%            tb=it1; te=ntimes;
%        else
               tb=it-24;   if tb<1, tb=1;end
               te=it+24*1; 
%        end;     
       
       if te > ntd ,te=ntd; end;
       
         tml=Date(tb:te);%-base_year;-base_year
          kt=find(tml>=Date(it),1,'first');
          jw=[tb:te]';
         if tml(1) <  1000; tml=tml+base_year_date; end
                               nmx0=7*24;
     if length(jw) <  nmx0, nmx=length(tml); mt=[1 : nmx ]';end;         
     if length(jw) >= nmx0, nmx=length(tml); mt=[jw-nmx0 : nmx ]';end
%          jwi=find(tml >=  Date(it),1,'first');

 h2(1)= plot(tml(mt),      FVCOM.el(jw,   npob ),'-b.');   hold on;
 h2(3)= plot(tml(1:kt),    FVCOM.el(tb:it,np1 )   ,'-g.');   hold on;
 h2(4)= plot(tml(kt)   ,min(FVCOM.el(:,   np1 ) ),'vk','MarkerSize',10); hold on;            %timer
 h2(2)= plot(tml(kt)   ,    FVCOM.el(it,  np2 )  ,'or','MarkerSize' ,8);hold on;
 
% lh= legend (h2(1:2),'W  ', 'E  '); 
  clear lht1 lht ncl;     
  lht1=['\zeta ( ' num2str(1    ,'%4.4d') ' )']; ncl=[1:length(lht1)]; %14
      xlimd= [tml(mt(1))-0.1 tml(mt(end))+0.1];

  set(gca,'xlim',xlimd);  % xlabel('Days');
%                                           ntk=10;
%                  if (tml(end)-tml(1))<31, ntk=5;end
%                  if (tml(end)-tml(1))<10, ntk=1;end
%        set(gca,'xtick',[fix(tml(1)):ntk:ceil(tml(end))]);
    set(gca,'YMinorTick','on','XMinorTick','on');  
  datetick('x','keeplimits','keepticks');    
  set(gca,'ylim',elim1);   ylabel('\zeta, m');
hold on;
%% -------- 