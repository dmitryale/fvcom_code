 function[B] = sinter(X,A,Y,M1,N1)
%based on sinter.f: get B(N1) at Y(N1) with call SINTER(DPTHSL,TA,ZM,TI,KSL,KBM1)
%  EXTRAPOLATION:

   for I=1:N1,...
       if Y(I) > X(1 ), B(I) = A(1) + ((A(1)-A(2))/(X(1)-X(2))) * (Y(I)-X(1)); end
       if Y(I) > X(M1), B(I) = A(M1); end       
   end
%  INTERPOLATION
     NM = M1 - 1;
   for I=1:N1,...
     for J=1:NM,...
        if Y(I) >= X(J) && Y(I) <= X(J+1), ...
        B(I) = A(J) - (A(J)- A(J+1)) * (X(J)-Y(I)) / (X(J)-X(J+1)) ;
        end
     end;   
   end;
 end   
