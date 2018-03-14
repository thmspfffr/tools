c = rand(500,1)

dat = rand(500,1) + c;
ref = rand(500,1) + c;

X = [ones(length(dat),1) ref];
Y = [ones(length(dat),1) dat];

reg = X\Y;
dat_reg = dat - reg(2,2).*ref + reg(1,2);

subplot(1,4,1)
scatter(dat,ref,'.'); lsline; axis square

r1 = corr(dat,ref);
disp(norm(ref))
ref1 = ref/norm(ref);
% ref = zscore(ref);

cleandat = dat-(dat'*ref1)*ref1;

r2 = corr(cleandat,ref1);

subplot(1,4,2)
scatter(cleandat,ref,'.'); lsline; axis square

dat = zscore(dat);
ref1 = zscore(ref);
disp(norm(ref1))

ref1 = ref1/norm(ref1);

cleandat = (dat)-(zscore(dat)'*ref1)*ref1;

subplot(1,4,3)
scatter(cleandat,ref1,'.'); lsline; axis square

corrcoef(cleandat,ref1)


subplot(1,4,4)

scatter(dat_reg,ref,'.'); lsline; axis square

corrcoef(dat_reg,ref)


