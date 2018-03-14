function y = tp_compute_degree(x,)


clear z1 z2

% reference (1) drug (1-3), (2) rest/task (1/2)
ref_cond = [1 1];

for ifoi = 1:13
  j = 1 : size(s_fc,1);

  m1 = squeeze(nanmean(s_fc(j,:,:,ref_cond(1),ref_cond(2),:),2));
  s1 = squeeze(nanstd(s_fc(j,:,:,ref_cond(1),ref_cond(2),:),[],2));
  m2 = squeeze(nanmean(s_fc(:,j,:,ref_cond(1),ref_cond(2),:),1));
  s2 = squeeze(nanstd(s_fc(:,j,:,ref_cond(1),ref_cond(2),:),[],1));

    for il = 1 : size(s_fc,1)
       fprintf('Computing location %d ...\n',il);

      jj = j(j~=il);


      for jl = 1:90


        x = squeeze(s_fc(il,jl,:,1,1,ifoi));

        z1 = (x' - squeeze(m1(jl,:,ifoi)))./squeeze(s1(jl,:,ifoi));
        z2 = (x' - squeeze(m2(jl,:,ifoi)))./squeeze(s2(jl,:,ifoi));

        [~,p1] = ttest(z1');
        [~,p2] = ttest(z2');

        th1(jl,:) = p1 < 0.01/2;
        th2(jl,:) = p2 < 0.01/2;     

      end

      th1 = th1(jj,:);
      th2 = th2(jj,:);

      % any connection significant?
      th(il,jj,:) = (th1 + th2) > 0; clear th1 th2

    end
        c(ifoi)= sum(th(:))/(size(th,1)*size(th,1))
        deg(:,ifoi) = mean(th);

end