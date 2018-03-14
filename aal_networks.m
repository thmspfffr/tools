function [new_matrix,all_ROI,networks_mask,net_def,net_mat]=aal_networks(matrix)

%% DMN - default mode network
net_def.dmn{1,1} = [31,32,...  % Cingulum Ant
                    21,22,...  % Olfactory
                    23,24,...  % Frontal Sup Medial
                    25,26,...  % Frontal Med Orb
                    15,16,...  % Frontal Inf Orb
                    89,90,...  % Temporal Inf
                    51,52,...  % Occiptal Mid   
                    65,66,...  % Angular
                    67,68,...  % Precuneus
                    35,36];    % Cingulum Post 
                net_def.dmn{1,2} = 1;

net_mat.dmn = nan(length(net_def.dmn{1,1})+1,length(net_def.dmn{1,1})+1);

for j = 1:length(net_def.dmn{1,1})-1
    a = net_def.dmn{1,1}(1,j);
    for k = j+1:length(net_def.dmn{1,1})
        a2 = net_def.dmn{1,1}(1,k);
        if a < a2
            net_mat.dmn(j,k) = matrix(a,a2);
        elseif a > a2
            net_mat.dmn(j,k) = matrix(a2,a);
        end
    end
end

%% FPC - frontoparietal executive / control network / WORKING MEMORY
net_def.fpc{1,1} = [47,48,...  % Lingual
                    61,62,...  % Parietal Inf
                    63,64,...  % SupraMarginal
                    65,66,...  % Angular *DMN
                    13,14,...  % Frontal Inf Tri
                    15,16,...  % Frontal Inf Orb * DMN
                    5,6,...    % Frontal Sup Orb
                    9,10,...   % Frontal Mid Orb
                    11,12,...  % Frontal Inf Oper
                    7,8];      % Frontal Mid
                net_def.fpc{1,2} = 2;
      
net_mat.fpc = nan(length(net_def.fpc{1,1})+1,length(net_def.fpc{1,1})+1);

for j = 1:length(net_def.fpc{1,1})-1
    a = net_def.fpc{1,1}(j);
    for k = j+1:length(net_def.fpc{1,1})
        a2 = net_def.fpc{1,1}(k);
        if a < a2
            net_mat.fpc(j,k) = matrix(a,a2);
        elseif a > a2
            net_mat.fpc(j,k) = matrix(a2,a);
        end
    end
end

%% SALIENCE NETWORK
net_def.salience{1,1} = [27,28,... % Rectus
                        31,32,...  % Cingulum Ant
                        33,34,...  % Cingulum Mid
                        59,60,...  % Parietal Sup
                        67,68,...  % Precuneus 
                        45,46];    % Cuneus
                     net_def.salience{1,2} = 3;
                    
net_mat.salience = nan(length(net_def.salience{1,1} )+1,length(net_def.salience{1,1} )+1);

for j = 1:length(net_def.salience{1,1} )-1
    a = net_def.salience{1,1} (j);
    for k = j+1:length(net_def.salience{1,1} )
        a2 = net_def.salience{1,1} (k);     
        if a < a2
            net_mat.salience(j,k) = matrix(a,a2);
        elseif a > a2
            net_mat.salience(j,k) = matrix(a2,a);
        end
    end
end


%% VAN - ventral attention network
net_def.van{1,1} = [83,84,...  % Temporal Pole Sup
                    87,88,...  % Temporal Pole Mid
                    85,86,...  % Temporal Mid
                    89,90,...  % Temporal Inf
                    43,44,...  % Calcarine
                    3,4,...    % Frontal Sup
                    23,24,...  % Frontal Sup Medial
                    55,56];    % Fusiform
                net_def.van{1,2} = 4;

net_mat.van = nan(length(net_def.van{1,1})+1,length(net_def.van{1,1})+1);

for j = 1:length(net_def.van{1,1})-1
    a = net_def.van{1,1}(j);
    for k = j+1:length(net_def.van{1,1})
        a2 = net_def.van{1,1}(k);     
        if a < a2
            net_mat.van(j,k) = matrix(a,a2);
        elseif a > a2
            net_mat.van(j,k) = matrix(a2,a);
        end
    end
end

%% DAN - dorsal attention network
net_def.dan{1,1} = [59,60,...  % Parietal Sup
                    61,62,...  % Parietal Inf
                    1,2,...    % Precentral 
                    53,54];    % Occipital Inf
                net_def.dan{1,2} = 5;

net_mat.dan = nan(length(net_def.dan{1,1})+1,length(net_def.dan{1,1})+1);

for j = 1:length(net_def.dan{1,1})-1
    a = net_def.dan{1,1}(j);
    for k = j+1:length(net_def.dan{1,1})
        a2 = net_def.dan{1,1}(k);     
        if a < a2
            net_mat.dan(j,k) = matrix(a,a2);
        elseif a > a2
            net_mat.dan(j,k) = matrix(a2,a);
        end
    end
end

%% SMN - somatosensory network
net_def.smn{1,1} = [19,20,... % SAM
                    1,2,...   % Precentral 
                    57,58,... % Postcentral
                    69,70];   % Paracentral Lobule 
                net_def.smn{1,2} = 6;

net_mat.smn = nan(length(net_def.smn{1,1})+1,length(net_def.smn{1,1})+1);

for j = 1:length(net_def.smn{1,1})-1
    a = net_def.smn{1,1}(j);
    for k = j+1:length(net_def.smn{1,1})
        a2 = net_def.smn{1,1}(k);     
        if a < a2
            net_mat.smn(j,k) = matrix(a,a2);
        elseif a > a2
            net_mat.smn(j,k) = matrix(a2,a);
        end
    end
end

%% VIS - visual network
net_def.vis{1,1} = [55,56,... % Fusiform
                    53,54,... % Occipital Inf
                    51,52,... % Occipital Mid
                    49,50,... % Occipital Sup                    
                    89,90];   % Temporal Inf (geh?rt zu Fusiform) 
                net_def.vis{1,2} = 7;

net_mat.vis = nan(length(net_def.vis{1,1})+1,length(net_def.vis{1,1})+1);

for j = 1:length(net_def.vis{1,1})-1
    a = net_def.vis{1,1}(j);
    for k = j+1:length(net_def.vis{1,1})
        a2 = net_def.vis{1,1}(k);     
        if a < a2
            net_mat.vis (j,k) = matrix(a,a2);
        elseif a > a2
            net_mat.vis (j,k) = matrix(a2,a);
        end
    end
end
%% Auditory
net_def.aud{1,1} = [17,18,... % Rolandic Oper
                    57,58,... % Fusiform
                    79,80,... % Heschl
                    81,82,... % Temporal Sup
                    29,30];   % Insula
                net_def.aud{1,2} = 8;


net_mat.aud = nan(length(net_def.aud{1,1})+1,length(net_def.aud{1,1})+1);

for j = 1:length(net_def.aud{1,1})-1
    a = net_def.aud{1,1}(j);
    for k = j+1:length(net_def.aud{1,1})
        a2 = net_def.aud{1,1}(k);     
        if a < a2
            net_mat.aud (j,k) = matrix(a,a2);
        elseif a > a2
            net_mat.aud (j,k) = matrix(a2,a);
        end
    end
end

%% Basalganglia

net_def.BG{1,1} = [75,76,71,72,73,74]; % 'Pallidum','Caudate','Putamen'
net_def.BG{1,2} = 9;

net_mat.BG = nan(length(net_def.BG{1,1})+1,length(net_def.BG{1,1})+1);

for j = 1:length(net_def.BG{1,1})-1
    a = net_def.BG{1,1}(j);
    for k = j+1:length(net_def.BG{1,1})
        a2 = net_def.BG{1,1}(k);     
        if a < a2
            net_mat.BG(j,k) = matrix(a,a2);
        elseif a > a2
            net_mat.BG(j,k) = matrix(a2,a);
        end
    end
end

%% SORTING

b = [net_def.dmn{1,1},net_def.fpc{1,1},net_def.salience{1,1},net_def.van{1,1},net_def.dan{1,1},net_def.smn{1,1},net_def.vis{1,1},net_def.aud{1,1},net_def.BG{1,1}];


[A,m,n]=unique(b);
b_new = zeros(100,1);

for i = 1:length(A)
    b_new(m(i)) = A(i);
end

b = b_new;
b = b';
b(b==0) = [];
   
bs = sort(b);

other = [37,38,... % Hippocampus
    39,40,...      % Parahippocampal
    41,42,...      % Amygdala
    77,78];        % Thalamus

net_def.hip{1,1} = 83:86;
net_def.thal{1,1} = 87:88;
net_def.amyg{1,1} = 89:90;

net_def.hip{1,2} = 10;
net_def.thal{1,2} = 11;
net_def.hamygp{1,2} = 12;

all_ROI = [b other];

new_matrix = nan(90,90);

for j = 1:length(all_ROI)-1
    a = all_ROI(j);
    for k = j+1:length(all_ROI)
        a2 = all_ROI(k);
        if a < a2
            new_matrix(j,k) = matrix(a,a2);
        elseif a > a2
            new_matrix(j,k) = matrix(a2,a);
        end
    end
end
 
%% make network mask
networks_mask = nan(90,90);

% ----------------------------- dmn
for j = 1:length(net_def.dmn{1,1})
    for k = 1:length(net_def.dmn{1,1})
        networks_mask(j,k) = 1;
    end
end

% -----------------------------  fpc % 20 / 16
for j = 21:36
    for k = 21:36
        networks_mask(j,k) = 2;
    end
end

a = find(all_ROI==65 | all_ROI==66 | all_ROI==15 | all_ROI==16);
for j = 1:length(a)
    ja = a(j);
    for k = 21:36
    networks_mask(ja,k) = 2;
    networks_mask(k,ja) = 2;
    end
end

% -----------------------------  salience 12 / 8
for j = 37:44
    for k = 37:44
        networks_mask(j,k) = 3;
    end
end

a = find(all_ROI==31 | all_ROI==32 | all_ROI==67 | all_ROI==68);
for j = 1:length(a)
    ja = a(j);
    for k = 37:44
    networks_mask(ja,k) = 3;
    networks_mask(k,ja) = 3;
    end
end

% -----------------------------  van 16 / 12
for j = 45:56
    for k = 45:56
        networks_mask(j,k) = 4;
    end
end

a = find(all_ROI==89 | all_ROI==90 | all_ROI==23 | all_ROI==24);
for j = 1:length(a)
    ja = a(j);
    for k = 45:56
    networks_mask(ja,k) = 4;
    networks_mask(k,ja) = 4;
    end
end

% -----------------------------  dan 8 / 4
for j = 57:60
    for k = 57:60
        networks_mask(j,k) = 5;
    end
end

a = find(all_ROI==59 | all_ROI==60 | all_ROI==61 | all_ROI==62);
for j = 1:length(a)
    ja = a(j);
    for k = 57:60
    networks_mask(ja,k) = 5;
    networks_mask(k,ja) = 5;
    end
end

% ----------------------------- smn 8 / 6
for j = 61:66
    for k = 61:66
        networks_mask(j,k) = 6;
    end
end

a = find(all_ROI==1 | all_ROI==2);
for j = 1:length(a)
    ja = a(j);
    for k = 61:66
    networks_mask(ja,k) = 6;
    networks_mask(k,ja) = 6;
    end
end

% ----------------------------- vis 10 / 2
for j = 67:68
    for k = 67:68
        networks_mask(j,k) = 7;
    end
end

a = find(all_ROI==55 | all_ROI==56 | all_ROI==53 | all_ROI==54 | all_ROI==51 | all_ROI==52 | all_ROI==89 | all_ROI==90);
for j = 1:length(a)
    ja = a(j);
    for k = 67:68
    networks_mask(ja,k) = 7;
    networks_mask(k,ja) = 7;
    end
end

% ----------------------------- aud 10 / 8
for j = 69:76
    for k = 69:76
        networks_mask(j,k) = 8;
    end
end

a = find(all_ROI==57 | all_ROI==58);
for j = 1:length(a)
    ja = a(j);
    for k = 69:76
    networks_mask(ja,k) = 8;
    networks_mask(k,ja) = 8;
    end
end

% ----------------------------- BG 6
for j = 77:82
    for k = 77:82
        networks_mask(j,k) = 9;
    end
end

% ----------------------------- HIP
for j = 83:86
    for k = 83:86
        networks_mask(j,k) = 10;
    end
end

% ----------------------------- THAL
for j = 87:88
    for k = 87:88
        networks_mask(j,k) = 11;
    end
end

% ----------------------------- AMYG
for j = 89:90
    for k = 89:90
        networks_mask(j,k) = 12;
    end
end