function [Date2,Val2] = aver_3hc(Date1, Val1, Navh)
%(c) dmitry.aleynik@sams.ac.uk 2018.05.05
 clearvars Date2 Val2 ; %date1 val1 nex dth dt 
%  ====== extend data it both way by Navh/dth extra point 
   dt=mean(diff(Date1(1:3))); %in days
  dth=dt*24;                  %in hours
  nex=round(0.5*Navh/dth); a=1:nex;
  date1(1:nex,1)= Date1(1); 
  val1=[];
  for  i=1:nex, 
      date1(i)  =date1(i)-dt*(nex-i+1); 
       val1(i,1)=Val1(1);
  end
        date1=cat(1,date1,Date1); 
        val1 =cat(1,val1 ,Val1);
         ie=length(date1); j=0;
  for  i=ie+1:ie+nex, 
       j=j+1;
       date1(i)=date1(ie)+dt*j;
        val1(i)=Val1(end);
  end  
%    ==================
       T1=Date1(1)  ;
       T2=Date1(end);
       Ndays=ceil(T2)-fix(T1);
       ni=ceil(Ndays*(24/Navh));
  Val2=zeros(ni,1)*NaN;
 for i=1:ni,...              
   Tb2=fix(T1) + (i-1)*Navh/24;
   Te2=Tb2 + Navh/24;  
    indav=[]; 
    indav = find(date1 >= Tb2 & date1 < Te2); 
%   Date2(i,1)=(Tb2 + Te2 )*0.5;     
    Date2(i,1)=(Tb2);     % centerd at target (rounded) times 1:00 h
       if isempty(indav), continue; end
       if (Date2(i) <= T1-dt || Date2(i) > T2) , continue; end
    Val2( i,1)=nanmean(val1(indav,1));   
 end
end
