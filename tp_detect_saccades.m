% detect eye movements
function sacc = tp_detect_saccades(x)

fsample = 1000;

velocity2D.kernel = [1 1 0 -1 -1].*(fsample/6); % this is equivalent to Engbert et al (2003) Vis Res, eqn. (1)
velocity2D.mindur  =  12;  % minimum microsaccade duration in samples
velocity2D.velthres = 6;    

n = size(velocity2D.kernel,2);
pad = ceil(n/2);

ndatsample = size(x,2);

prepad    = min(pad, floor(size(x,2)/2));
edgeleft  = mean(x(:,1:prepad),2);
postpad    = min(pad, floor(size(x,2)/2));
edgeright = mean(x(:,1+end-postpad:end),2);
x       = [edgeleft*ones(1,pad) x edgeright*ones(1,pad)];
vel = convn(x,   velocity2D.kernel,   'same');
vel = vel(:, pad+1:end-pad);

medianstd = sqrt( nanmedian(vel.^2,2) - (nanmedian(vel,2)).^2 );

% Engbert et al (2003) Vis Res, eqn. (3)
radius = velocity2D.velthres*medianstd;

% compute test criterion: ellipse equation
test = sum((vel./radius(:,ones(1,ndatsample))).^2,1);
sacsmp = find(test>1);% microsaccade's indexing

% first find eye movements of n-consecutive time points
j = find(diff(sacsmp)==1);
j1 = [j; j+1];
com = intersect(j,j+1);
cut = ~ismember(j1,com);
sacidx = reshape(j1(cut),2,[]);
movement = []
for k=1:size(sacidx,2);
  duration = sacidx(1,k):sacidx(2,k);
  if size(duration,2) >= velocity2D.mindur;
    % finding peak velocity by Pitagoras
    begtrl = sacsmp(duration(1,1));
    endtrl = sacsmp(duration(1,end));
    
    [peakvel smptrl] = max(sqrt(sum(vel(:,begtrl:endtrl).^2,1)));
    veltrl = sacsmp(duration(1,smptrl));% peak velocity microsaccade sample -> important for spike conversion
    
    trlsmp =  1:size(x,2);
    begsample = trlsmp(1, begtrl); % begining microsaccade sample
    endsample = trlsmp(1, endtrl); % end microsaccade sample
    velsample = trlsmp(1, veltrl); % velocity peak microsaccade sample
    movement(end+1,:) = [begsample endsample velsample];
  end
end