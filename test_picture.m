% Plot the figure of encoding picture ``peppers.png''
% The test of qiskit is based on the unitary synthesis: https://docs.quantum.ibm.com/api/qiskit/synthesis.

clc; clear; close all;
addpath( "bitble-qclab" ) ;
addpath( "test_fable" ) ;

offset = 0 ;
logging = 1 ;
compr_type = 'cutoff' ;%'percentage'; 
compr_val = 1e-8 ;
circuit_sim = false ;

record = struct();

img = imread('peppers.png');
img1 = img(:,:,1) ;
img2 = img(:,:,2) ;
img3 = img(:,:,3) ;

%% loading data for qiskit
pepper_mat = struct() ;
for i = 1:3
    if i == 1
        A = img1;
    elseif i == 2
        A = img2;
    else
        A = img3;
    end
    n = ceil(log(length(A))/log(2)) ;
    A(pow2(n),pow2(n)) = 0 ;
    [BlockEncodingMatrix, normalization_factor] = BlockEncoding(double(A)) ;
    pepper_mat(i).M = BlockEncodingMatrix;
end

save("pepper_mat.mat","pepper_mat") ;


% BlockEncoding()
function [BlockEncodingMatrix, normalized_factor] = BlockEncoding(A)
% Generate the (1.0001*norm(A),1,0)-block encoding 'BlockEncodingA' of 'A'
% input: A \in \mathbb{C}^{n\times n}
% output: normlized_factor = 1.0001*norm(A); A = A ./ normlized_factor;
%         BlockEncodingMatrix = [A, sqrt(eye(n)-A*A'); sqrt(eye(n)-A'*A), A']
    
    assert( size(A,1) == size(A,2) );
    n = size(A,1);
    normalized_factor = 1.0001*norm(A);
    A = A ./ normalized_factor;
    [uA12, sA12] = eig( eye(n)-A*A' );
    UA12 = uA12 * diag( sqrt(diag(sA12)) ) * uA12';
    [uA21, sA21] = svd( eye(n)-A'*A );
    UA21 = uA21 * diag( sqrt(diag(sA21)) ) * uA21';
    BlockEncodingMatrix = [ A,  UA12; -UA21, A' ];
end % end of BlockEncoding


%% test bitble and fable
for i = 1:3
    disp("Fro-norm     Spectra-norm   ratio") ;
    if i == 1
        A = img1 ;
    elseif i == 2 
        A = img2 ;
    else
        A = img3 ;
    end
    A = double(A);
    n = floor(log(max(size(A)))/log(2)) ;
    if min(size(A)) ~= pow2(n)
        A(pow2(n),pow2(n)) = 0;
    end
        
    fprintf(" %e       %e       %e \n",norm(A,'fro'),norm(A,2),norm(A,'fro')/norm(A,2)) ;
    record(i).F_norm = norm(A,'fro');
    record(i).Spectral_norm = norm(A,2);
    record(i).norm_rotia = norm(A,'fro')/norm(A,2) ;
    record(i).max_element = max(max(abs(A))) ;

    %% Simulate the FABLE quantum circuit 
    
    fprintf("\n\nFABLE Block Encoding \n");
    fprintf("------------------------------------------------------------ \n");
    fprintf("parameter computing... \n");
    t1 = clock ;
    [circuit, OA, alpha, info] = fable( A, compr_type, compr_val, logging, circuit_sim) ;
    t2 = clock ;
    time3 = etime(t2,t1) ;
    record(i).fable_nRY = info.circ.nRY ;
    if ~isreal(A)
    record(i).fable_nRZ = info.circ.nRZ ;
    end
    subnormalization_factor = max(max(abs(A))) .* pow2(n) ;
    record(i).fable_nCNOT = info.circ.nCNOT ;
    record(i).fable_subnormalized_factor = subnormalization_factor ;
    record(i).fable_time = time3 ;
    if ~real(A)
        record(i).fable_nG_merit = subnormalization_factor .* (info.circ.nRY + info.circ.nRZ) ;
    else
        record(i).fable_nG_merit = subnormalization_factor .* info.circ.nRY ;
    end
    record(i).fable_nCNOT_merit = subnormalization_factor .* info.circ.nCNOT ;
    
    if logging
        info.circ
    end 
    nRG = info.circ.nRY ;
    if ~isreal(A)
        nRG = nRG + info.circ.nRZ ;
    end
    fprintf( "number of RY + RZ gates * subnormalization factor = %e \n", subnormalization_factor .* nRG ) ;
    fprintf( "number of CNOT gates * subnormalization factor = %e \n", subnormalization_factor .* info.circ.nCNOT ) ;


     %% Simulate the BITBLE quantum circuit 
    fprintf("\n\nBITBLE Block Encoding \n");
    fprintf("------------------------------------------------------------ \n");
    fprintf("parameter computing... \n");
    t1 = clock ;
    [~, subnormalization_factor, info] = bitble( A, compr_type, compr_val, logging, offset, circuit_sim ) ;
    t2 = clock ;
    time1 = etime(t2,t1) ;
    record(i).dimension = n ;
    record(i).bitble_nRY = info.circ.nRY ;
    if ~isreal(A)
    record(i).bitble_nRZ = info.circ.nRZ ;
    end
    record(i).bitble_nCNOT = info.circ.nCNOT ;
    record(i).bitble_subnormalized_factor = subnormalization_factor ;
    record(i).bitble_time = time1 ;
    record(i).bitble_nG_merit = subnormalization_factor .* (info.circ.nRY + info.circ.nRZ) ;
    record(i).bitble_nCNOT_merit = subnormalization_factor .* info.circ.nCNOT ;
    fprintf( "number of RY + RZ gates * subnormalization factor = %e \n", subnormalization_factor .* (info.circ.nRY + info.circ.nRZ) ) ;
    fprintf( "number of CNOT gates * subnormalization factor = %e \n", subnormalization_factor .* info.circ.nCNOT ) ;
    
    %% Simulate the BITBLE2 quantum circuit 
    
    fprintf("\n\nBITBLE2 Block Encoding \n");
    fprintf("------------------------------------------------------------ \n");
    fprintf("parameter computing... \n");
    t1 = clock ;
    [~, subnormalization_factor, info] = bitble2( A, compr_type, compr_val, logging, offset, circuit_sim ) ;
    t2 = clock ;
    time2 = etime(t2,t1) ;
    record(i).bitble2_nRY = info.circ.nRY ;
    if ~isreal(A)
    record(i).bitble2_nRZ = info.circ.nRZ ;
    end
    record(i).bitble2_nCNOT = info.circ.nCNOT ;
    record(i).bitble2_subnormalized_factor = subnormalization_factor ;
    record(i).bitble2_time = time2 ;
    record(i).bitble2_nG_merit = subnormalization_factor .* (info.circ.nRY + info.circ.nRZ) ;
    record(i).bitble2_nCNOT_merit = subnormalization_factor .* info.circ.nCNOT ;
    fprintf( "number of RY + RZ gates * subnormalization factor = %e \n", subnormalization_factor .* (info.circ.nRY + info.circ.nRZ) ) ;
    fprintf( "number of CNOT gates * subnormalization factor = %e \n", subnormalization_factor .* info.circ.nCNOT ) ;

    %% Simulate the BITBLE3 quantum circuit 
    fprintf("\n\nBITBLE3 Block Encoding \n");
    fprintf("------------------------------------------------------------ \n");
    fprintf("parameter computing... \n");
    t1 = clock ;
    [~, subnormalization_factor, info] = bitble3( A, 0.5, compr_type, compr_val, logging, offset, circuit_sim ) ;
    t2 = clock ;
    time1 = etime(t2,t1) ;
    record(i).dimension = n ;
    record(i).bitble3_nRY = info.circ.nRY ;
    if ~isreal(A)
    record(i).bitble3_nRZ = info.circ.nRZ ;
    end
    record(i).bitble3_nCNOT = info.circ.nCNOT ;
    record(i).bitble3_subnormalized_factor = subnormalization_factor ;
    record(i).bitble3_time = time1 ;
    record(i).bitble3_nG_merit = subnormalization_factor .* (info.circ.nRY + info.circ.nRZ) ;
    record(i).bitble3_nCNOT_merit = subnormalization_factor .* info.circ.nCNOT ;
    
    fprintf( "number of RY + RZ gates * subnormalization factor = %e \n", subnormalization_factor .* (info.circ.nRY + info.circ.nRZ) ) ;
    fprintf( "number of CNOT gates * subnormalization factor = %e \n", subnormalization_factor .* info.circ.nCNOT ) ;
    
end
record_peppers = record;
save("record_peppers.mat","record_peppers") ;

%% Plot bar
clc;
close all;
load("record_peppers.mat") ;

multiple = 1.3 ;
y = [record_peppers(3).bitble_nCNOT_merit./(multiple*record_peppers(3).fable_nCNOT_merit), record_peppers(3).bitble2_nCNOT_merit/(multiple*record_peppers(3).fable_nCNOT_merit), 1/multiple;
    record_peppers(3).bitble_nG_merit./(multiple*record_peppers(3).fable_nG_merit), record_peppers(3).bitble2_nG_merit./(multiple*record_peppers(3).fable_nG_merit), 1/multiple;
    record_peppers(3).bitble_time./(multiple*record_peppers(3).fable_time), record_peppers(3).bitble2_time./(multiple*record_peppers(3).fable_time), 1/multiple ];
str = ["CNOT's size metric"; "RY's size metric"; "Time(s)"];
text_labels = {num2str(record_peppers(3).bitble_nCNOT_merit,"%.3g"), num2str(record_peppers(3).bitble2_nCNOT_merit,"%.3g"), num2str(record_peppers(3).fable_nCNOT_merit,"%.2g");
    num2str(record_peppers(3).bitble_nG_merit,"%.3g"), num2str(record_peppers(3).bitble2_nG_merit,"%.3g"), num2str(record_peppers(3).fable_nG_merit,"%.2g");
    num2str(record_peppers(3).bitble_time,"%.2f"), num2str(record_peppers(3).bitble2_time,"%.2f"), num2str(record_peppers(3).fable_time,"%.2f")};


bar_fig = bar(y, "GroupWidth", 0.92) ;
set(gca, 'XTickLabel',str, 'XTick',1:numel(str)) ;
bar_fig(1).FaceColor = 'r';
bar_fig(2).FaceColor = 'b';
bar_fig(3).FaceColor = [.8 .8 .8];
legend("BITBLE^1","BITBLE^2","BITBLE^3","FABLE",'Location','northwest');
set(gca,'ytick',[],'yticklabel',[]);

% Label BITBLE
xtips = bar_fig(1).XEndPoints;
subxtips = xtips(2) ;
ytips = bar_fig(1).YEndPoints;
subytips = ytips(2) ;
labels = text_labels(2,1);
t = text(subxtips+0.05,subytips + 0.09,labels,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom','FontWeight','bold') ;
set(t,'Rotation',90);
subxtips = xtips(1) ;
subytips = ytips(1) ;
labels = text_labels(1,1);
t = text(subxtips+0.05,subytips + 0.09,labels,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom') ;
set(t,'Rotation',90);

% Label BITBLE2
xtips = bar_fig(2).XEndPoints;
ytips = bar_fig(2).YEndPoints;
subxtips = xtips(1) ;
subytips = ytips(1) ;
labels = text_labels(1,2);
t = text(subxtips+0.05,subytips + 0.09,labels,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom','FontWeight','bold') ;
set(t,'Rotation',90);
subxtips = xtips(2) ;
subytips = ytips(2) ;
labels = text_labels(2,2);
t = text(subxtips+0.05,subytips + 0.09,labels,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom','FontWeight','bold') ;
set(t,'Rotation',90);

% Label FABLE
xtips = bar_fig(3).XEndPoints;
ytips = bar_fig(3).YEndPoints;
subxtips = xtips(1:2) ;
subytips = ytips(1:2) ;
labels = text_labels(1:2,3);
t = text(subxtips+0.05,subytips + 0.09,labels,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom') ;
set(t,'Rotation',90);


for i=1:length(bar_fig)
    xtips = bar_fig(i).XEndPoints;
    xtips = xtips(3) ;
    ytips = bar_fig(i).YEndPoints;
    ytips = ytips(3) ;
    labels = text_labels(3,i);
    if i == 1
        t = text(xtips,ytips,labels,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom','FontWeight','bold') ;
    else
        t = text(xtips,ytips,labels,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom') ;
    end
end



ylim([0,1]);



%% PLOT the 'MERIT', CNOT's gate count $\times$ normalization factor, of "Peppers.png"

clc;clear;close all;

MarkerSize1 = 8 ;
MarkerSize2 = 6 ;
LineWidth = 0.5 ;

color_shallow_blue = [0.8, 0.8, 1.0];
color_gray = [0.7 0.7 0.7] ;
color_gray2 = [0.9 0.9 0.9] ;

load("record_peppers.mat") ;
len = length(record_peppers) ;

fable_nCNOT_merit = zeros(1,len) ;
fable_nG_merit = zeros(1,len) ;
fable_time = zeros(1,len) ;
bitble_nCNOT_merit = zeros(1,len) ;
bitble_nG_merit = zeros(1,len) ;
bitble_time = zeros(1,len) ;
bitble2_nCNOT_merit = zeros(1,len) ;
bitble2_nG_merit = zeros(1,len) ;
bitble2_time = zeros(1,len) ;
bitble3_nCNOT_merit = zeros(1,len) ;
bitble3_nG_merit = zeros(1,len) ;
bitble3_time = zeros(1,len) ;
qiskit_time = [670.7749, 893.5286, 782.1517] ;
qiskit_nCNOT_merit = [435371, 435371, 435371] .* record_peppers(3).Spectral_norm .* 1.001;
qiskit_nG_merit = [369834, 369835, 369834] .* record_peppers(3).Spectral_norm .* 1.001;

for i = 1:len
    fable_nCNOT_merit(i) = record_peppers(i).fable_nCNOT_merit ;
    fable_nG_merit(i) = record_peppers(i).fable_nG_merit ;
    fable_time(i) = record_peppers(i).fable_time ;
    bitble_nCNOT_merit(i) = record_peppers(i).bitble_nCNOT_merit ;
    bitble_nG_merit(i) = record_peppers(i).bitble_nG_merit ;
    bitble_time(i) = record_peppers(i).bitble_time ;
    bitble2_nCNOT_merit(i) = record_peppers(i).bitble2_nCNOT_merit ;
    bitble2_nG_merit(i) = record_peppers(i).bitble2_nG_merit ;
    bitble2_time(i) = record_peppers(i).bitble2_time ;
    bitble3_nCNOT_merit(i) = record_peppers(i).bitble3_nCNOT_merit ;
    bitble3_nG_merit(i) = record_peppers(i).bitble3_nG_merit ;
    bitble3_time(i) = record_peppers(i).bitble3_time ;
end

figure;
yyaxis right ;
y = zeros(len,5);
labels = cell(len,5); 
for i = 1:len
    y(i,:) = [record_peppers(i).bitble_time, record_peppers(i).bitble2_time, record_peppers(i).bitble3_time, record_peppers(i).fable_time, qiskit_time(i) ];
    labels(i,:) = {num2str(record_peppers(i).bitble_time,"%.1f"), num2str(record_peppers(i).bitble2_time,"%.1f"), num2str(record_peppers(i).bitble3_time,"%.1f"), num2str(record_peppers(i).fable_time,"%.1f"), num2str(qiskit_time(i),"%.f")};
end


bar_fig = bar(y, "GroupWidth", 0.92) ;
set(gca,'YScale','log')
bar_fig(1).FaceColor = 'r' ;
bar_fig(2).FaceColor = 'b' ;
bar_fig(3).FaceColor = color_shallow_blue ;
bar_fig(4).FaceColor = color_gray ;
bar_fig(5).FaceColor = color_gray2 ;

for i=1:length(bar_fig)
    xtips = bar_fig(i).XEndPoints;
    xtips = xtips(1:len) ;
    ytips = bar_fig(i).YEndPoints;
    ytips = ytips(1:len) ;
    sublabels = labels(:,i);
    if i == 1
        t = text(xtips+0.05,ytips + 0.3,sublabels,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom','FontWeight','bold') ;
    elseif i < 5
        t = text(xtips+0.05,ytips + 0.8,sublabels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom') ;
    else
        t = text(xtips+0.05,ytips + 200,sublabels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom') ;
    end
    set(t,'Rotation',90);
end

ylabel('time(s)','FontWeight','normal');
ylim([0,3.*max(qiskit_time)]) ;
axis on;


hold on;
yyaxis left ;

semilogy(1:len,bitble_nCNOT_merit,'-o','LineWidth',LineWidth,...
    'Color','r','MarkerSize',MarkerSize1,'MarkerFaceColor','r'); 
hold on;
semilogy(1:len,bitble2_nCNOT_merit,'-o','LineWidth',LineWidth,...
    'Color','b','MarkerSize',MarkerSize2,'MarkerFaceColor','b'); 
hold on;
semilogy(1:len,bitble3_nCNOT_merit,'-o','LineWidth',LineWidth,...
    'Color',color_shallow_blue,'MarkerSize',MarkerSize2,'MarkerFaceColor',color_shallow_blue); 
hold on;
semilogy(1:len,fable_nCNOT_merit,'-o','LineWidth',LineWidth,...
    'Color',color_gray,'MarkerSize',MarkerSize1,'MarkerFaceColor',color_gray);
hold on;
semilogy(1:len,qiskit_nG_merit,'-o','LineWidth',LineWidth,...
    'Color','k','MarkerSize',MarkerSize1,'MarkerFaceColor',color_gray2);

x_labels = {"Red channel","Green channel","Blue channel"};
ylabel("CNOT's size metric",'Interpreter','latex');
ylim([1e8,5*1e10]) ;
set(gca,'YScale','linear')


h2 = legend(["BITBLE$^1$","BITBLE$^2$","BITBLE$^3$","FABLE","Qiskit"],'Interpreter','latex','Location','north','Orientation','horizontal');
set(h2,'Box','off');
set(h2,'FontSize',8);
ax = gca;
set(gca, 'XTick', 1:len);
set(gca, 'XTickLabel', x_labels);

ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';


%% Plot RGB channels

% clc;clear;close all;
image = imread('peppers.png');

% RGB
red_channel = image(:, :, 1);
green_channel = image(:, :, 2);
blue_channel = image(:, :, 3);
red_image = cat(3, red_channel, zeros(size(red_channel)), zeros(size(red_channel)));
green_image = cat(3, zeros(size(green_channel)), green_channel, zeros(size(green_channel)));
blue_image = cat(3, zeros(size(blue_channel)), zeros(size(blue_channel)), blue_channel);

figure;
subplot(2, 2, 1);
imshow(image);
xlabel('Original Image');

% Red
subplot(2, 2, 2);
imshow(red_image);
xlabel('Red Channel');
% Green
subplot(2, 2, 3);
imshow(green_image);
xlabel('Green Channel');
% Blue
subplot(2, 2, 4);
imshow(blue_image);
xlabel('Blue Channel');

