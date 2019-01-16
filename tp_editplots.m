function k = tp_editplots(varargin)

if ~isempty(varargin)
  set(varargin{1},'linewidth',0.5,'ticklength',[0.02 0.02],'tickdir','out');
  set(varargin{1},'fontsize',7); %axis square
else
  set(gca,'linewidth',0.5,'ticklength',[0.02 0.02],'tickdir','out');
  set(gca,'fontsize',7); %axis square
end
