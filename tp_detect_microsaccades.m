function y = tp_detect_microsaccades(x,fs)

% based on engbert et al. (2003) and fieldtrip code

kernel = [1 1 0 -1 -1].*fs/6;
y = [];
thresh = 6;


x(isnan(x(:,1)),:)=[];
n = size(x,1);
x=x-repmat(mean(x),[length(x) 1]);

if size(x,2)>1
if sum(x(:,2)==0) 
  x(:,2) = [];
end
end

if n < 100
  % signals are rarely shorter than 100 samples...
else
  vel = convn(x',kernel,'same');
end

medianstd = sqrt( median(vel.^2,2) - median(vel,2).^2 );
radius = thresh*medianstd;
test = sum((vel./radius(:,ones(1,n))).^2,1);
sacsmp = find(test>1);

% first find eye movements of n-consecutive time points
j = find(diff(sacsmp)==1);
j1 = [j; j+1];
com = intersect(j,j+1);
cut = ~ismember(j1,com);
sacidx = reshape(j1(cut),2,[]);

for k=1:size(sacidx,2);
  duration = sacidx(1,k):sacidx(2,k);
  if size(duration,2) >= 12;
    % finding peak velocity by Pitagoras
    begtrl = sacsmp(duration(1,1));
    endtrl = sacsmp(duration(1,end));
    
    [peakvel smptrl] = max(sqrt(sum(vel(:,begtrl:endtrl).^2,1)));
    veltrl = sacsmp(duration(1,smptrl));% peak velocity microsaccade sample -> important for spike conversion
    
    trlsmp = 1:length(x);
    begsample = trlsmp(1, begtrl); % begining microsaccade sample
    endsample = trlsmp(1, endtrl); % end microsaccade sample
    velsample = trlsmp(1, veltrl); % velocity peak microsaccade sample
    y(end+1,:) = [begsample endsample velsample];
  end
end