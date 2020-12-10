% the difficulties with FVCOM is that our netcdf (3 & 4)  files are in Ugrid convention (unstructured)  https://ugrid-conventions.github.io/ugrid-conventions/ 
% The best way would be to use native OPeNDAP option within Matlab scripts for fields extraction from remote files:
% a single-time snapshot (the last small file) of 2D or 3D temperature field: it takes 2.6 & 35.0 seconds: 

tic
url_fw='https://thredds.sams.ac.uk/thredds/dodsC/scoats-westcoms2/Archive_forecast/netcdf_2020F/R20200318/westcoms2_20200324_R20200318_0007.nc';
FV1.info=ncinfo(url_fw);
FV1.Times=ncread(url_fw,'Times')';
FV1.temp=ncread(url_fw,'temp');  % , startLoc,count,stride);
toc

tic
url_fw='https://thredds.sams.ac.uk/thredds/dodsC/scoats-westcoms2/Archive_forecast/netcdf_2020F/R20200318/westcoms2_20200323_R20200318_0006.nc';

FV.info=ncinfo(url_fw);
FV.Times=ncread(url_fw,'Times')';
FV.temp=ncread(url_fw,'temp');  % , startLoc,count,stride);
toc