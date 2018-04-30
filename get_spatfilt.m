function A1 = get_spatfilt(para)

if ~isfield(para,'cs')
  fprintf('CS not found. Computing CS ...')
end

if ndims(para.cs)==3
  % LCMV WITH CODE FROM GUIDO
  if strcmp(para.filt,'lcmv')
    if strcmp(para.grid,'xcoarse')
      load(para.sa);
      [~,A1] = mkfilt_lcmv(sa.L_xcoarse,nanmean(para.cs(:,:,para.foi(1,1):para.foi(1,2)),3));
    elseif strcmp(para.grid,'coarse')
      load(para.sa);
      [~,A1] = mkfilt_lcmv(sa.L_coarse,nanmean(para.cs(:,:,para.foi(1,1):para.foi(1,2)),3));
    elseif strcmp(para.grid,'cortex')
      load(para.sa);
      [~,A1] = mkfilt_lcmv(sa.L_coarse,nanmean(para.cs(:,:,round(para.foi(1,1)):round(para.foi(1,2))),3));
      A1 = getdipdir(nanmean(para.cs(:,:,round(para.foi(1,1)):round(para.foi(1,2))),3),A1);
    elseif strcmp(para.grid,'medium')
      load(para.sa);
      [~,A1] = mkfilt_lcmv(sa.L_medium,nanmean(para.cs(:,:,round(para.foi(1,1)):round(para.foi(1,2))),3));
      A1 = getdipdir(nanmean(para.cs(:,:,round(para.foi(1,1)):round(para.foi(1,2))),3),A1);
    end

  % ELORETA WITH CODE FROM GUIDO
  elseif strcmp(para.filt,'eloreta')
    if isstr(para.sa)
      load(para.sa);
    else
      sa = para.sa;
    end
    if strcmp(para.grid,'xcoarse')
      % v2
      A1 = mkfilt_eloreta_v2(sa.L_xcoarse);
      A1 = getdipdir(nanmean(para.cs(:,:,round(para.foi(1,1)):round(para.foi(1,2))),3),A1);
    elseif strcmp(para.grid,'coarse')
      % v1
      A1 = mkfilt_eloreta_v2(sa.L_coarse);
      A1 = getdipdir(nanmean(para.cs(:,:,round(para.foi(1,1)):round(para.foi(1,2))),3),A1);
    elseif strcmp(para.grid,'cortex')
      % v3
      % watch out! L_coarse is correct even for cortex grid
      A1 = mkfilt_eloreta_v2(sa.L_coarse);
      A1 = getdipdir(nanmean(para.cs(:,:,round(para.foi(1,1)):round(para.foi(1,2))),3),A1);
    elseif strcmp(para.grid,'aal')
      A1 = mkfilt_eloreta_v2(sa.L_aal);
      A1 = getdipdir(nanmean(para.cs(:,:,round(para.foi(1,1)):round(para.foi(1,2))),3),A1);
    elseif strcmp(para.grid,'medium')
      A1 = mkfilt_eloreta_v2(sa.L_medium);
      A1 = getdipdir(nanmean(para.cs(:,:,round(para.foi(1,1)):round(para.foi(1,2))),3),A1);
    end

    % LCMV WITH CODE FROM JOERG
  elseif strcmp(para.filt,'jh_lcmv')
    load(para.sa);
    if strcmp(para.grid,'xcoarse')
      A1 = pconn_beamformer(para.cs(:,:,round(para.foi(1,1)):round(para.foi(1,2))),sa.L_xcoarse);
    elseif strcmp(para.grid,'coarse')
      A1 = pconn_beamformer(para.cs(:,:,round(para.foi(1,1)):round(para.foi(1,2))),sa.L_coarse);
    elseif strcmp(para.grid,'cortex')
      A1 = pconn_beamformer(para.cs(:,:,round(para.foi(1,1)):round(para.foi(1,2))),sa.L_coarse);
    elseif strcmp(para.grid,'aal')
      A1 = pconn_beamformer(para.cs(:,:,round(para.foi(1,1)):round(para.foi(1,2))),sa.L_aal);      
    elseif strcmp(para.grid,'medium')
      A1 = pconn_beamformer(para.cs(:,:,round(para.foi(1,1)):round(para.foi(1,2))),sa.L_medium);      
    end
  end
elseif ndims(para.cs)==2
  % LCMV WITH CODE FROM GUIDO
  if strcmp(para.filt,'lcmv')
    if strcmp(para.grid,'xcoarse')
      load(para.sa);
      [~,A1] = mkfilt_lcmv(sa.L_xcoarse,para.cs);
    elseif strcmp(para.grid,'coarse')
      load(para.sa);
      [~,A1] = mkfilt_lcmv(sa.L_coarse,para.cs);
    elseif strcmp(para.grid,'cortex')
      load(para.sa);
      [~,A1] = mkfilt_lcmv(sa.L_coarse,para.cs);
      A1 = getdipdir(para.cs,A1);
    elseif strcmp(para.grid,'aal')
      load(para.sa);
      [~,A1] = mkfilt_lcmv(sa.L_aal,para.cs);
      A1 = getdipdir(para.cs,A1);    
    elseif strcmp(para.grid,'aal')
      load(para.sa);
      [~,A1] = mkfilt_lcmv(sa.L_medium,para.cs);
      A1 = getdipdir(para.cs,A1);    
    end

  % ELORETA WITH CODE FROM GUIDO
  elseif strcmp(para.filt,'eloreta')
    if isstr(para.sa)
      load(para.sa);
    else
      sa = para.sa;
    end
    if strcmp(para.grid,'xcoarse')
      % v2
      A1 = mkfilt_eloreta_v2(sa.L_xcoarse);
      A1 = getdipdir(para.cs,A1);
    elseif strcmp(para.grid,'coarse')
      % v1
      A1 = mkfilt_eloreta_v2(sa.L_coarse);
      A1 = getdipdir(para.cs,A1);
    elseif strcmp(para.grid,'cortex')
      % v3
      % watch out! L_coarse is correct even for cortex grid
      A1 = mkfilt_eloreta_v2(sa.L_coarse);
      A1 = getdipdir(para.cs,A1);
    elseif strcmp(para.grid,'aal')
%       load(para.sa);
      [~,A1] = mkfilt_lcmv(sa.L_aal,para.cs);
      A1 = getdipdir(para.cs,A1);  
    elseif strcmp(para.grid,'medium')
%       load(para.sa);
      [~,A1] = mkfilt_lcmv(sa.L_medium,para.cs);
      A1 = getdipdir(para.cs,A1);  
    elseif strcmp(para.grid,'L_aal_4mm')
      [~,A1] = mkfilt_lcmv(sa.L_aal_4mm,para.cs);
      A1 = getdipdir(para.cs,A1);  
    end

    % LCMV WITH CODE FROM JOERG
  elseif strcmp(para.filt,'jh_lcmv')
    load(para.sa);
    if strcmp(para.grid,'xcoarse')
      A1 = pconn_beamformer(para.cs,sa.L_xcoarse);
    elseif strcmp(para.grid,'coarse')
      A1 = pconn_beamformer(para.cs,sa.L_coarse);
    elseif strcmp(para.grid,'cortex')
      A1 = pconn_beamformer(para.cs,sa.L_coarse);
    elseif strcmp(para.grid,'aal')
      A1 = pconn_beamformer(para.cs,sa.L_aal);
    elseif strcmp(para.grid,'medium')
      A1 = pconn_beamformer(para.cs,sa.L_medium);
    end
  end
end
  
