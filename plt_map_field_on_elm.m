%  function[]=plt_map_field_on_elm (FVCOM1,Mesh)
%  addpath([homeM 'Mar_Phys/matlab/general/']);
var_plot=['UV']; 
var_plot_un =['cm\cdots^{-1}']; 
IZ=1; % surface
it=1; % first time slot

vrmfe=[]; 
u=double(squeeze(FVCOM1.u(it,IZ,:)));
v=double(squeeze(FVCOM1.v(it,IZ,:)));
[dir_deg,uv] = cart2compass(u,v);
   vrmfe=uv; %speed  

%get cells id to be used
 N=Mesh.Nelems; 
 M=Mesh.Nverts; % nodes
     celn=zeros(M,1); %==Mesh.ntve
     cell=zeros(M,8); %==Mesh.nbve(:,1:8);  - slightly reordered 
 for ii=1:M 
     no=ii; %node
     cels=[];ce =[];
  for jj=1:3
        cel=find(Mesh.trinodes(:,jj)==no);
    ce(jj,1:length(cel))=cel;
  end  
      ce=ce(ce>0);
     ceu=unique(ce);         
    cels=cat(1,cels,ceu);
      nc=length(cels);
    celn(ii,1)=nc;
    cell(ii,1:nc)=cels;
 end  
%% ===================
% do extraction of the field from 'elems' to the 'nodes'
  vrmfn=zeros(M,1);    

 for ii=1:M  %nodes
           clear cels nc;
      nc=celn(ii,1);
    cels=cell(ii,1:nc);    
   vrmfn(ii,1)=mean(vrmfe(cels));
 end 
  vrm=vrmfn;  
 
%
      HPatch= patch('Vertices',[Mesh.geog(:,1),Mesh.geog(:,2)],'Faces',Mesh.trinodes,...
         'Cdata',vrm ,'edgecolor','interp','facecolor','interp');
     hold on;
      colormap(flipud(lbmap(24,'Step5Seq')))
     cbr=colorbar('location','west');      % cbt=['T ^o C'];
                     cbt={[ var_plot ', '], [var_plot_un ]};
     set(get(cbr,'title'),'string',[cbt],'HorizontalAlignment','left');
%    set(cbr,'position',[ 0.82 0.47 0.02 0.35]); 
     set(cbr,'position',[  0.91  0.343 0.02 0.35]);
   
%    end
     