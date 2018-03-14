function s = pconn_getbehavior(para)
% takes following input:
% para.cond (cnt/bttn/bth)
% para.subj (SUBJLIST)


addpath ~/pcbi/
addpath ~/pconn/matlab/

cond     = para.cond;
SUBJLIST =  para.subj;

ord = pconn_randomization;

if ~strcmp(cond,'bth')
  switch cond
    case 'cnt'

      cnt=pcbi_cnt(1:34);

      for isubj = 1:34
        for m = 1 : 3

          im = find(ord(isubj,:)==m);

          cnt_all(isubj,m) = nanmean(cnt(isubj,im*2-1:im*2));

        end
      end

      cnt_all = cnt_all(SUBJLIST,:);

      d_behav1 = cnt_all;

    case 'bttn'

      v_bttn = 4;

      load ~/pconn_bttn/proc/pconn_bttn_hist_v4.mat

      for isubj =  SUBJLIST
        for m = 1 : 3
          for iblock = 1 : 2

            try
              load(['~/pconn_bttn/proc/' sprintf('pconn_bttn_dur_s%d_b%d_m%d_v%d.mat',isubj,iblock,m,v_bttn)]) 
            catch 
              cnt(iblock,m,isubj) = nan;
                med_dur(iblock,m,isubj) = nan; 
                continue
            end

              if ~isempty(par.dur)
                cnt(iblock,m,isubj) = length(par.dur); 
                med_dur(iblock,m,isubj) = mean(par.dur); 
              else
                cnt(iblock,m,isubj) = nan;
                med_dur(iblock,m,isubj) = nan; 
              end
              clear par


          end    
        end
      end

      dur=squeeze(nanmean(med_dur(:,:,SUBJLIST)));
      cnt=squeeze(nanmean(cnt(:,:,SUBJLIST),1));

      d_behav1 = cnt(1,:)';

  end
  
else
    clear cnt_all
   cnt=pcbi_cnt(1:34);

    for isubj = SUBJLIST
      for m = 1 : 3

        im = find(ord(isubj,:)==m);

        cnt_all(isubj,m) = nanmean(cnt(isubj,im*2-1:im*2));

      end
    end

    cnt_all = cnt_all(SUBJLIST,:);

    d_behav = nanmean(cnt_all(:,1),2);
        
    v_bttn = 4; clear cnt_all cnt med_dur 

    load ~/pconn_bttn/proc/pconn_bttn_hist_v4.mat

    for isubj =  SUBJLIST
      for m = 1 : 3
        for iblock = 1 : 2

          try
            load(['~/pconn_bttn/proc/' sprintf('pconn_bttn_dur_s%d_b%d_m%d_v%d.mat',isubj,iblock,m,v_bttn)]) 
          catch 
            cnt(iblock,m,isubj) = nan;
              med_dur(iblock,m,isubj) = nan; 
              continue
          end

            if ~isempty(par.dur)
              cnt(iblock,m,isubj) = length(par.dur); 
              med_dur(iblock,m,isubj) = mean(par.dur); 
            else
              cnt(iblock,m,isubj) = nan;
              med_dur(iblock,m,isubj) = nan; 
            end
            clear par


        end    
      end
    end

    dur=squeeze(nanmean(med_dur(:,:,SUBJLIST)));
    cnt=squeeze(nanmean(cnt(:,:,SUBJLIST),1));

    d_behav1 = cnt(1,:)';
    
    d_behav1 = nanmean([d_behav1 d_behav],2);
    
end

s = d_behav1;