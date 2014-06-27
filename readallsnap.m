function A=readallsnap

% reads all snap*.out files into a matlab struct, 
% from which it can be easily plotted.
%

dat=dir('./snap*.out');

if length(dat)==0
    error('Error: no snap*.out found in current directory');
end

tmp=importdata(dat(1).name);

A.nspecies=tmp(1);
A.nfield=tmp(2);
A.nvgrid=tmp(3);
A.mpsi=tmp(4);
A.mtgrid=tmp(5);
A.mtoroidal=tmp(6);
A.tmax=tmp(7);

nspecies=A.nspecies;
nfield=A.nfield;
nvgrid=A.nvgrid;
mpsi=A.mpsi;
mtgrid=A.mtgrid;
mtoroidal=A.mtoroidal;
tmax=A.tmax;

rho = dir('./rho.out');
if length(rho)==0
    A.rho = 0:mpsi-1;
else    
    A.rho=importdata('./rho.out');
end    

A.profile=zeros(length(dat),mpsi,6,nspecies);
A.pdf=zeros(length(dat),nvgrid,4,nspecies);
A.poloidata=zeros(length(dat),mtgrid,mpsi,nfield);
A.fluxdata=zeros(length(dat),mtgrid,mtoroidal,nfield);
A.x=zeros(mtgrid,mpsi);
A.z=zeros(mtgrid,mpsi);
A.t=zeros(length(dat),1);
A.fieldrms=zeros(length(dat),mpsi,nfield);



A.totalnumber=length(dat);


for i=1:length(dat)
    fprintf('Now processing %s\n',dat(i).name);
    if (i>1) 
        tmp=importdata(dat(i).name);
    end
    A.profile(i,:,:,:)=reshape(tmp(8:A.mpsi*6*A.nspecies+7),[A.mpsi,6,A.nspecies]);
    index=A.mpsi*6*A.nspecies+8;
    A.pdf(i,:,:,:)=reshape(tmp(index:index+A.nvgrid*4*A.nspecies-1),[A.nvgrid,4,A.nspecies]);
    index=index+A.nvgrid*4*A.nspecies;
    A.poloidata(i,:,:,:)=reshape(tmp(index:index+A.mtgrid*A.mpsi*A.nfield-1),[A.mtgrid,A.mpsi,A.nfield]);
    index=index+A.mtgrid*A.mpsi*A.nfield;
    if(i==1)
        A.x(:,:)=reshape(tmp(index:index+A.mtgrid*A.mpsi-1),[A.mtgrid,A.mpsi]);
        index=index+A.mtgrid*A.mpsi;
        A.z(:,:)=reshape(tmp(index:index+A.mtgrid*A.mpsi-1),[A.mtgrid,A.mpsi]);
        index=index+A.mtgrid*A.mpsi;
    else
        index=index+A.mtgrid*A.mpsi*2;
    end
    A.fluxdata(i,:,:,:)=reshape(tmp(index:index+A.mtgrid*A.mtoroidal*A.nfield-1),[A.mtgrid,A.mtoroidal,A.nfield]);   
    A.t(i)=str2num(dat(i).name(5:9));
    for j=1:mpsi
        for k=1:nfield
            A.fieldrms(i,j,k)=A.poloidata(i,:,j,k)*transpose(A.poloidata(i,:,j,k))/mtgrid;
        end
    end
end


