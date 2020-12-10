close all
clear

% read in dima's model states (subsampled)
load FVsubsampled

% set start time for particle tracking in matlab days
starttime = FVsub.mt(1);

% set particle tracking timestep
dtsec = 60;   % minutes
dtday = dtsec/(60*60*24); % convert to days

% conversion factors between degrees of longitude/latitude
metrestolon = 90/1e7/cosd(56);
metrestolat = 90/1e7;

% intialise particle postions
% **your code here!**

% time stepping loop
it = 1;  % initial time step
itmax = 2000; % number of timesteps for this simulation

while(it<itmax)
it = it+1;
time(it) = time(it-1)+dtday;

% need to decide which is the model snapshot to use
% no interpolaton in time used
% **your code here!**

% find u and v at the current time and particle positions
% **your code here!**

% step particle positions forward in time
% **your code here!**

% plot particles
clf
% **your code here!**

% draw coastline
load etive6_0coast
plot(ncst(:,1),ncst(:,2),'k')
set(gca,'xlim',[-5.5 -5.35],'ylim',[56.43 56.5])
set(gca,'dataaspectratio',[1 cosd(mean(ylim)) 1])
% add a timestamp
text(-5.45,56.44,datestr(time(it)))
drawnow  % needed to make each frame appear as a movie

end




