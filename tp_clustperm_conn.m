oldp = p;

%%

dat = rand(90,90);
dat = dat<0.05;
t = 0.05;

clear idx

% dat = p<t;

for i = 1 : 90
  if sum(dat(:,i))<4
    idx(i)=1;
  end
end

dat(:,find(idx))=0;
% dat = tril(dat,-1);
k = 1 : 90;

%%
 i= 0;
while ~isempty(k)
  i = i + 1;
  [d]=bfs(dat,i);
  
  clust{i}=find(d>0);
  a(i) = sum(d>0);
  
  
  
end


%%
%
clear clust


j       = 0;
nclust  = 0;
conns   = [];
chain   = [];
cnt     = 0;
mat     = zeros(90,90);

done = [];
%
while length(done)<90 
  
    chain = [];
    j = j + 1;
     
    if ismember(j,done)
      continue
    end
    
    chain = [find(dat(j,:)) find(dat(:,j))'];
    done  = [done j];

    if ~isempty(chain)
      nclust = nclust + 1;
      clust{nclust} = [repmat(j,length(chain),1) chain'];
      fprintf('Cluster found!\n')
    else
      continue
    end
    
    idx = [chain' repmat(j,length(chain),1) ];
    
    for iidx = 1 : size(idx,1)
      mat(idx(iidx,1),idx(iidx,2)) = nclust;
    end
    
    % find all cluster components here    
    while ~isempty(chain) 
      
      curr  = chain(1);
      chain(1) = [];
      
      if ismember(curr,done)
        continue
      end
      
      done  = unique([done curr]);

      new_conns1 = find(dat(curr,:));
      new_conns2 = find(dat(:,curr));
    
      idx1 = [repmat(curr,length(new_conns1),1) new_conns1'];
      idx2 = [new_conns2 repmat(curr,length(new_conns2),1) ];
      
      for iidx = 1 : size(idx1,1)
        mat(idx1(iidx,1),idx1(iidx,2)) = nclust;
      end
      for iidx = 1 : size(idx2,1)
        mat(idx2(iidx,1),idx2(iidx,2)) = nclust;
      end
  
      new_conns1(ismember(new_conns1,done))=[];
      new_conns2(ismember(new_conns2,done))=[];
      
      chain = [chain new_conns1 new_conns2'];
            
    end
end

%%

dat = zeros(90,90);

clust1 = [  5 1; ...
  10 1
  20 1
  50 20
  40 20
  85 18
  20 18];
  
clust2 =  ...
[  3 2
  12 2
  70 60
  70 15
  15 2
  
  % three
  6 4
  28 4
  89 28
];

[a k]=find(dat);


for i = 1 : length(clust1)
  dat(clust1(i,1),clust1(i,2))=1;
end
for i = 1 : length(clust2)
  dat(clust2(i,1),clust2(i,2))=1;
end

% for i = 1 : length(

% COMBINE CLUSTERS NOW
% k = 1:size(clust,2);
% for i = 1 : size(clust,2)
%   
%   i = i + 1
%   clear comb
%   for j = k(k~=i)
%     
%     comb(j) = any(ismember(clust{i},clust{j}));
%     
%   end
%   
%   sum(comb)
% end
% 
%   if ~any(comb)
%     continue
%   else
%     
  
  
  

    
    
  

%%