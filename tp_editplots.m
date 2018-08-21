function k = tp_editplots(varargin)

if ~isempty(varargin)
  set(varargin{1},'linewidth',1,'ticklength',[0.03 0.03],'tickdir','out');
  set(varargin{1},'fontsize',10); %axis square
else
  set(gca,'linewidth',1,'ticklength',[0.03 0.03],'tickdir','out');
  set(gca,'fontsize',10); %axis square
end
