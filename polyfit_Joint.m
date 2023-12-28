function [ph_de, ph_trend]=polyfit_Joint(ph_in,win_ksize)       
  a=ph_in;
[ma,na]=size(a);
b=2*ones(ma,na);

a1=find(ones(size(a)));
[az,range]=find(ones(size(a)));
obs=[a1,az,range,a(:),b(:)];

coord1=obs(:,2:3);
nlink=4;
NPS=size(coord1,1);
Pcls=knnsearch(coord1,coord1,'k',8);  %每个点周围的nlink个临近点 
Pcls(:,2:5)=[];

n_i=ma/win_ksize;
n_j=na/win_ksize;
n_total=n_i*n_j;

B_idx=zeros(n_total,n_total);
count=0;
idx_edge=[];
win0=zeros(win_ksize,win_ksize);
for i=1:win_ksize
    win0(i,:)=(i-1)*win_ksize+1:i*win_ksize;
end

for i=1:n_i
    for j=1:n_j
        wink=(i-1)*win_ksize+(j-1)*win_ksize+win0;
        idx_win=(i-1)*n_j+j;
        if mod(idx_win,n_j)~=0
            idx_edge=[idx_edge;wink(:,end-1:end)];
            count=count+1;
            B_idx(idx_win,idx_win+1)=count;
            B_idx(idx_win+1,idx_win)=-count;
        end
        if i<n_i
            idx_edge=[idx_edge;wink(end-1:end,:)'];
            count=count+1;
            B_idx(idx_win,idx_win+n_j)=count;
            B_idx(idx_win+n_j,idx_win)=-count;
        end
    end
end

idx_edge=unique(idx_edge);
idx_inte=setdiff([1:ma*na]',idx_edge);

Pcls1=Pcls(idx_inte,:);
Pcls2=Pcls(idx_edge,:);

Pcls1(:,3:4)=[];

NPS1=size(idx_inte,1);
sp=(repmat(Pcls1(:,1),1,nlink-1-2))';   
from = reshape(sp,(nlink-1-2)*NPS1,1); 
to   = reshape(Pcls1(:,2:end)',(nlink-1-2)*NPS1,1); 
Arctem1=[from,to]; %基于临近点构网

NPS2=size(idx_edge,1);
sp=(repmat(Pcls2(:,1),1,nlink-1))';   
from = reshape(sp,(nlink-1)*NPS2,1); 
to   = reshape(Pcls2(:,2:end)',(nlink-1)*NPS2,1); 
Arctem3=[from,to]; %基于临近点构网

 
DT=delaunayTriangulation(coord1);
    Arctem2= edges(DT);
   
    Arc=[Arctem1;Arctem2;Arctem3];
%     Arc=[Arctem1];
    IDX_to=Arc(:,2);
    IDX_from=Arc(:,1);
    


% nestpoint=knnsearch(coord1,coord1,'k',9);
% figure(2)
% triplot(DT);
pkeep=size(IDX_from,1);
ranges=obs(:,3);
azis=obs(:,2);
ph=obs(:,4);
dph=wrap(ph(IDX_to,:)-ph(IDX_from,:));

% for i=1:pkeep
%     X1=[azis(IDX_from(i)),azis(IDX_to(i))];
%     X2=[ranges(IDX_from(i)),ranges(IDX_to(i))];
%     line(X1,X2);
%     hold on;
% end
% hold off;
% amp=obs(:,5);

% win_ksize=5;

% A=zeros(pkeep,n_total*3+count-1);
% tmp=[];
% A=zeros(pkeep,n_total*6);


%estimate a0
a0_mat=zeros(n_i,n_j);
D=a0_mat;
for i=1:n_i
    for j=1:n_j
        i1=(i-1)*win_ksize+1;
        i2=i*win_ksize;
        j1=(j-1)*win_ksize+1;
        j2=j*win_ksize;
        [a0_mat(i,j),D(i,j)]=arcpoly(a(i1:i2,j1:j2));
        
    end
end
[m,n]=find(D==min(min(D)));

a0=a0_mat(m,n);
t=(m-1)*n_j+n;
nArc=size(Arc,1);
% arc_del=[];
si=[];sj=[];sv=[];
for i=1:nArc
     point1_x=mod(azis(IDX_from(i))-1,win_ksize)+1;
    point1_y=mod(ranges(IDX_from(i))-1,win_ksize)+1;
    point2_x=mod(azis(IDX_to(i))-1,win_ksize)+1;
    point2_y=mod(ranges(IDX_to(i))-1,win_ksize)+1;
    
    point1_x=(point1_x-1)/(win_ksize-1);
    point1_y=(point1_y-1)/(win_ksize-1);
    point2_x=(point2_x-1)/(win_ksize-1);
    point2_y=(point2_y-1)/(win_ksize-1);
    
    
    n_x1=ceil(azis(IDX_from(i))/win_ksize);
    n_y1=ceil(ranges(IDX_from(i))/win_ksize);
    n_idx1=(n_x1-1)*n_j+n_y1;
    n_x2=ceil(azis(IDX_to(i))/win_ksize);
    n_y2=ceil(ranges(IDX_to(i))/win_ksize);
    n_idx2=(n_x2-1)*n_j+n_y2;
%     idx_a0=B_idx(n_idx1,n_idx2);
%     tmp=[tmp;point1_x,point1_y,point2_x,point2_y,n_idx1,n_idx2];
    
    if n_idx1==n_idx2
          si=[si,i*ones(1,5)];
          sj=[sj,(n_idx1-1)*6+2:n_idx1*6];
          sv=[sv,point2_x-point1_x ,point2_y-point1_y,point2_x*point2_y-point1_x*point1_y,point2_x*point2_x-point1_x*point1_x,point2_y*point2_y-point1_y*point1_y];
    elseif n_idx1==t
          dph(i)=dph(i)+a0;
          si=[si,i*ones(1,12)];
          sj=[sj,(n_idx1-1)*6+1:n_idx1*6,(n_idx2-1)*6+1:n_idx2*6];
          sv=[sv,-1,-point1_x,-point1_y,-point1_x*point1_y,-point1_x^2,-point1_y^2,1,point2_x,point2_y,point2_x*point2_y,point2_x^2,point2_y^2];

    elseif n_idx2==t
          dph(i)=dph(i)-a0;
          si=[si,i*ones(1,12)];
          sj=[sj,(n_idx1-1)*6+1:n_idx1*6,(n_idx2-1)*6+1:n_idx2*6];
          sv=[sv,-1,-point1_x,-point1_y,-point1_x*point1_y,-point1_x^2,-point1_y^2,1,point2_x,point2_y,point2_x*point2_y,point2_x^2,point2_y^2];

    else
          si=[si,i*ones(1,12)];
          sj=[sj,(n_idx1-1)*6+1:n_idx1*6,(n_idx2-1)*6+1:n_idx2*6];
          sv=[sv,-1,-point1_x,-point1_y,-point1_x*point1_y,-point1_x^2,-point1_y^2,1,point2_x,point2_y,point2_x*point2_y,point2_x^2,point2_y^2];

    end
     
end

% A(arc_del,:)=[];
% S=sparse(A);
S=sparse(si,sj,sv);
xxx=sum(abs(S'));
xxxd=find(xxx~=0);
tempS=S(xxxd,:); 
tempS=[tempS(:,1:6*(t-1)),tempS(:,6*(t-1)+2:end)];
% dph(arc_del,:)=[];

% par=(S'*S)^-1*(S'*dph);
% par=lsmr(tempS,dph);
[par,res]=IRLS(tempS,dph);
par=[par(1:6*(t-1));a0;par(6*(t-1)+1:end)];


ph_est=zeros(ma,na);
for i=1:ma
    for j=1:na
        i1=mod(i-1,win_ksize)+1;
        j1=mod(j-1,win_ksize)+1;
        i1=(i1-1)/(win_ksize-1);
        j1=(j1-1)/(win_ksize-1);
        
        idx_win=(ceil(i/win_ksize)-1)*n_j+ceil(j/win_ksize);
%         ph_est(i,j)=[1,i1,j1,i1*j1]*[a0(idx_win);par(3*idx_win-2:3*idx_win)];
        ph_est(i,j)=[1,i1,j1,i1*j1,i1^2,j1^2]*par(6*idx_win-5:6*idx_win);
    end
end

% dph0=a-wrap(ph_est);
% dph0=mean(dph0(:));
% dph0=a(1,1)-ph_est(1,1);
% ph_trend=dph0+ph_est;
% ph_de=a-wrap(ph_trend);
ph_trend=ph_est;

ph_de=a-wrap(ph_trend);

  
end

% estimation method
function [x,res]=IRLS(A,b)

x0=lsmr(A,b);
r=b-A*x0;
cvg=1;
scale=median(abs(r-median(r)))/0.6745;
% scale=1;
count=0;
n=length(b);
r=r/scale;
idx=r>1;
W=ones(size(b));
W(idx)=1./r(idx);

for count=1:30
%     r=r/scale;
%     x=x0;
%     W=(abs(r)<1).*(1-r.^2).^2;
   
%     P=diag(W);
    x=x0;
    P=sparse(1:n,1:n,W);
    x0=(P*A)\(P*b);
    r=abs(b-A*x0);
    
    scale=median(abs(r-median(r)))/0.6745;
    r=r/scale;
    
    idx=r>1;
    W=ones(size(b));
    W(idx)=1./r(idx);
    
    cvg=max(abs(x-x0));
     disp(cvg);
    if cvg<1e-3
        break;
    end
   
end
    disp(count);
    res=r;
end


function [a0,D]=arcpoly(A)
[ma,na]=size(A);

id=find(ones(size(A)));
[az,range]=find(ones(size(A)));
obs_win0=[id,range,az,A(:)];


% range_win=pcluster_range(i,:);
% azi_win=pcluster_azi(i,:);


coord_win0=obs_win0(:,2:3);
obs_ph=obs_win0(:,4);
nlink=6;
NPS=size(coord_win0,1);
Pcls=knnsearch(coord_win0,coord_win0,'k',10);
Pcls(:,2:5)=[];

sp=(repmat(Pcls(:,1),1,nlink-1))';   
from = reshape(sp,(nlink-1)*NPS,1); 
to   = reshape(Pcls(:,2:end)',(nlink-1)*NPS,1); 
Arc1_win0=[from,to]; %基于临近点构网


DT=delaunayTriangulation(coord_win0);
Arc2_win0= edges(DT);
Arc_win0=[Arc1_win0;Arc2_win0];

    m = min(range);
    M = max(range);
    ranges_norm = (range-m)/(M-m) ;
    m = min(az); M = max(az);
    azis_norm= (az-m)/(M-m);
%     
    from=Arc_win0(:,1);
    to=Arc_win0(:,2);
   
    
    dr=ranges_norm(to)-ranges_norm(from);
    daz=azis_norm(to)-azis_norm(from);
    drdaz=ranges_norm(to).*azis_norm(to)-ranges_norm(from).*azis_norm(from);
%     dr=range(to)-range(from);
%     daz=az(to)-az(from);
%     drdaz=range(to).*az(to)-range(from).*az(from);
    arc_ploy=[dr,daz,drdaz];
      y=wrap(obs_ph(to)-obs_ph(from));
%     
    [b2,state2]=robustfit(arc_ploy,y,'bisquare',4.685,'off');
    
    obse2=[ranges_norm,azis_norm,ranges_norm.*azis_norm]*b2;
    ph_tre=reshape(obse2,ma,na);
    dif=A-ph_tre;
    a0=mean(dif(:));
    D=var(dif(:));
%     res2=state2.resid;
%     a0=1;
    
end

