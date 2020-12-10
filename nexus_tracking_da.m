close all
clear

load FVsubsampled

starttime = FVsub.mt(1); % start time in matlab days

dtsec = 60;   % particle timestep in mins
dtday = dtsec/(60*60*24); % convert to days

metrestolon = 90/1e7/cosd(56);
metrestolat = 90/1e7;

xpart(1,:) = -5.43+rand(1,1000)*.01;
ypart(1,:) = 56.47+rand(1,1000)*.01;
time(1) = starttime;
% time stepping loop
it = 1;  % time step
itmax = 650; %1000;
load etive6_0coast

%da -----------
 time_end = time(1)+(itmax-1)*dtday;
 np=length(xpart(1,:)); % number of particles
 mkgif     ='Y'; gi=0;  
 path_fig='../fig_plTR/';
 [status,message,messageid]=mkdir(path_fig);   clear status message messageid
  load('../code/A_Grid.mat')
  G.h=-G.Vq; 
    hlv=[-250:25:-100 -90:10:-20 -15:5:-5.01 -5 -3 0]; %;-3 -1 0 ];
%da -----------

while(it<itmax)
it = it+1;
time(it) = time(it-1)+dtday;

% need to decide which is the model snapshot to use
fvstep = find(FVsub.mt<time(it), 1, 'last' );

% find u and v
uthis = interp2(FVsub.x,FVsub.y,squeeze(FVsub.u(fvstep,1,:,:)),xpart(it-1,:),ypart(it-1,:));
vthis = interp2(FVsub.x,FVsub.y,squeeze(FVsub.v(fvstep,1,:,:)),xpart(it-1,:),ypart(it-1,:));

xpart(it,:) = xpart(it-1,:)+uthis*dtsec*metrestolon;
ypart(it,:) = ypart(it-1,:)+vthis*dtsec*metrestolat;

% plot particles
clf
    contour(G.xq,G.yq,G.h,hlv); hold on
    cbar(1)=colorbar('location','east');
    set(cbar(1),'position',[ 0.84 0.20 0.013 0.65]);  % [ 0.053 0.18 0.02  0.56]);
    set(get(cbar(1),'title'),'string',[ 'h,m' ]);    
          
plot(xpart(it,:),ypart(it,:),'r.');hold on
% draw coastline
plot(ncst(:,1),ncst(:,2),'k')
plot(mean(xpart(1,:)),mean(ypart(1,:)),'go','markersize',4); hold on
                
set(gca,'xlim',[-5.4667   -5.3750] ,'ylim',[56.445 56.496667 ]);% da
% set(gca,'xlim',[-5.5 -5.35],'ylim',[56.43 56.5])
set(gca,'dataaspectratio',[1 cosd(mean(ylim)) 1]);

% text(-5.45,56.44,datestr(time(it))); % add a timestamp
  text(-5.43,56.447,datestr(time(it))); % add a timestamp %da

   set(gcf,'color','w'); % set(gcf,'color','none')
   set(gcf,'PaperPositionMode','auto');
 drawnow  % needed to make each frame appear as a movie

 %% %%%%%%%%%%%making gif vect:
  if  mkgif=='Y', ....
          
   if gi==0,...
   date_daye= [datestr(time(1),'yyyymmdd_HHMM') '-' datestr(time_end,'yyyymmdd_HHMM')];        
   namemtr=strcat([ path_fig 'Tracking_' num2str(np)  '_' date_daye   ] );
                   gifname=[namemtr '.gif'];               
   II2 = getframe(gcf);
   II2 = frame2im(II2);                
   [XX2, map] = rgb2ind(II2, 256); % need  Image Processing Toolbox
   imwrite(XX2, map, gifname, 'GIF', 'WriteMode', 'overwrite', 'DelayTime', 0, 'LoopCount', Inf);
   gi=gi+1;
   else 
   II2 = getframe(gcf);
   II2 = frame2im(II2);
   [XX2, map] = rgb2ind(II2, 128);%
   imwrite(XX2, map, gifname, 'GIF', 'WriteMode', 'append', 'DelayTime', 0);
   gi=gi+1;
   
   end; %gi
   clear II2 XX2
  end; %mkgif
 %% %%%%%%%%%%%%%%%%%%
 
end




