function [Date2,Val2] = aver_3h(Date1, Val1, Navh)
clear Date2 Val2;
       T1=Date1(1)  ;
       T2=Date1(end);
        Ndays=ceil(T2)-fix(T1);
for i=1:Ndays*(24/Navh),...              
   Tb2=fix(T1) + (i-1)*Navh/24;
   Te2=Tb2 + Navh/24;  
    indav=[]; 
    indav = find(Date1 >= Tb2 & Date1 < Te2);   
    Date2(i,1)=(Tb2 +Te2)*0.5;
    Val2(i,1)=mean(nanmean(Val1(indav,1)));   
end
end
