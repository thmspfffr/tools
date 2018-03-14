function ExportGridOnMRI(pnt, mri)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

	for k = 68 : 5 : 188
		fid = figure('Position', [500 500 720 720]);
		fname = ['.\grid2Topo', num2str(k, '%03d')];
		im = squeeze(mri(k, :, :))';
		imagesc(im);
		colormap(gray);
		hold on;
		tmp = abs(pnt(:, 1) - k) < 5;
		plot(pnt(tmp, 2), pnt(tmp, 3), 'g.');
		axis equal;
		print(fid, '-dpng', fname);
		clear tmp fname fid;
	end
	
	fid = figure('Position', [500 500 720 720]);
	fname = '.\grid2TopoAll.png';
	im = squeeze(mri(128, :, :))';
	imagesc(im);
	colormap(gray);
	hold on;
	plot(pnt(:, 2), pnt(:, 3), 'g.');
	axis equal;
	print(fid, '-dpng', fname);
	clear fname fid;
	close all;
end

