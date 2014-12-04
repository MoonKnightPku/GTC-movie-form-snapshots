function A=readallsnap(target_dir, output)

% reads all snap*.out files into a matlab struct, 
% from which numerical analysis can be easily produced.
%
% target_dir: is the target directory you want to read
% 'a' or 'a/' work the same as reading all snap*.out in the a/ directory
% if not specified, the default value is './'
%
% output is a switch (0 or 1) whether you want output during the reading
% Sometimes the reading would take several minuites and you may want to
% monitor the process to make sure it is running but not freezing. 
% if not specified, the default value is output=0
%

%% initialization: handle optional argument and deal with potential format issue
if ~exist('target_dir','var') || isempty(target_dir)
    fprintf('Target directory not set, reading the current directory.\n')
    % default target directory
    filename = './';
end

if ~exist('output','var') || isempty(output)
    % default target directory
    output=0;
end


if target_dir(length(target_dir)) ~= '/'
    % if target_dir not ends with '/', add '/' to the end of target_dir
    target_dir = strcat(target_dir , '/');
end


%% read in the first file, and set up the dimensions
dat=dir(strcat(target_dir,'snap*0*.out'));

tmp=importdata(strcat(target_dir, dat(1).name));

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

rhofile = strcat(target_dir, 'rho.out');
if length(dir(rhofile))==0
    A.rho = 0:mpsi-1;
else
    A.rho=importdata(rhofile);
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

%% read in the rest data file
for i=1:length(dat)
    if (output>0)
        fprintf('Now processing %s\n',strcat(target_dir, dat(i).name));
    end
    if (i>1) 
        tmp=importdata(strcat(target_dir, dat(i).name));
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
            A.fieldrms(i,j,k)=rms(A.poloidata(i,:,j,k));
        end
    end
end


