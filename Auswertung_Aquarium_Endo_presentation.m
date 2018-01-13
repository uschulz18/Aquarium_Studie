%Aquaiium Studie
% Die Tiefe, die Winkel und das abzufahrende Volumen sind zu bestimmen aus
% aufgenommen 3D daten
%--------------------------------------------------------------------------

function auswertung_aquarium()

%% Einlesen aller Möglichen Files die zu einer Variante(Gelenk) gehört

[filename,path] = uigetfile('C:\Users\schulzu\Desktop\Aquarium\Files\*.trc');

% % suchen nach allen .trc die dem selben Werkzeug angehören 
% trenn = regexp(filename,'_');
% link = filename(1:trenn(1)-1);
% mat = dir([path link '*.trc']);
% unnamed raussuchen und löschen 
% files = struct2cell(mat);
% files = files(1,:)';
files = filename;
% un = regexp(files, 'Unnamed');
% 
% for i = length(un):-1:1
%     if isempty(un{i}) == 0 
%         files(i,:) = [];
%     end
% end
% 
% clear i mat un filename;

%% die Files auswerten
for i = 1:length(files)    
   
    % einlesen der Daten 
    %----------------------------------------------------------------------
    %dlmread liest alle Zahlenwerte ein
    data = dlmread([path files],'\t',6,0);
    fid= fopen([path files]);
    
    %mit textscan werden die Überschriften eingelesen
    strings = textscan(fid,'%s');%, 22,'delimiter','\t');
    
    if str2num(strings{1}{16}) == 6
        marker  = [cellstr(strings{1}{23}) cellstr(strings{1}{24}) cellstr(strings{1}{25})...
           cellstr(strings{1}{26}) cellstr(strings{1}{27}) cellstr(strings{1}{28})];
        ueberschrift = [cellstr(strings{1}{21}) cellstr(strings{1}{22}) cellstr(strings{1}{29})...
                cellstr(strings{1}{30}) cellstr(strings{1}{31}) cellstr(strings{1}{32}) ...
                cellstr(strings{1}{33}) cellstr(strings{1}{34}) cellstr(strings{1}{35}) ...
                cellstr( strings{1}{36}) cellstr(strings{1}{37}) cellstr(strings{1}{38}) ...
                cellstr(strings{1}{39}) cellstr(strings{1}{40}) cellstr(strings{1}{41})...
                cellstr(strings{1}{42}) cellstr(strings{1}{43}) cellstr(strings{1}{44})...
                cellstr(strings{1}{45}) cellstr(strings{1}{46})];

    elseif str2num(strings{1}{16}) == 10
        marker  = [cellstr(strings{1}{23}) cellstr(strings{1}{24}) cellstr(strings{1}{25})...
           cellstr(strings{1}{26}) cellstr(strings{1}{27}) cellstr(strings{1}{28}) ...
           cellstr(strings{1}{29}) cellstr(strings{1}{30}) cellstr(strings{1}{31}) cellstr(strings{1}{32})];
       
       ueberschrift = [cellstr(strings{1}{21}) cellstr(strings{1}{22}) cellstr(strings{1}{29})...
                cellstr(strings{1}{30}) cellstr(strings{1}{31}) cellstr(strings{1}{32}) ...
                cellstr(strings{1}{33}) cellstr(strings{1}{34}) cellstr(strings{1}{35}) ...
                cellstr( strings{1}{36}) cellstr(strings{1}{37}) cellstr(strings{1}{38}) ...
                cellstr(strings{1}{39}) cellstr(strings{1}{40}) cellstr(strings{1}{41})...
                cellstr(strings{1}{42}) cellstr(strings{1}{43}) cellstr(strings{1}{44})...
                cellstr(strings{1}{45}) cellstr(strings{1}{46}) cellstr(strings{1}{47})...
                cellstr(strings{1}{48}) cellstr(strings{1}{49}) cellstr(strings{1}{50}) ...
                cellstr(strings{1}{51}) cellstr(strings{1}{52}) cellstr(strings{1}{53}) ...
                cellstr(strings{1}{54}) cellstr(strings{1}{55}) cellstr(strings{1}{56}) ...
                cellstr(strings{1}{57}) cellstr(strings{1}{58}) cellstr(strings{1}{59}) ...
                cellstr(strings{1}{60}) cellstr(strings{1}{61}) ];
    end
    
    % Den Titel zerlegen um zu erkennen welche Methode angewandt werden
    % muss
    trenn = regexp(files,'_');
    if length(trenn) < 3 
        tool = files(1:trenn(1)-1);
        anzahlwerkzeug = files(trenn(1)+1:trenn(2)-1);
        option = files(trenn(2)+1:length(files)-5);
        messung =  files(end-4);
        name_messung = [tool anzahlwerkzeug option messung];
    elseif length(trenn) == 3 
        tool = files(1:trenn(1)-1);
        anzahlwerkzeug = files(trenn(1)+1:trenn(2)-1);
        option = files(trenn(3)+1:length(files)-5);
        extra = files(trenn(2)+1:trenn(3)-1);
        messung =  files(end-4);
        name_messung = [tool anzahlwerkzeug option extra messung];
    end
      
    
    %fclose([path files{i}]);
    
    for j=1:length(marker)
        daten.Rohdaten.(name_messung).(marker{j}).time = data(:,2);
        daten.Rohdaten.(name_messung).(marker{j}).coordinates = data(:,j*3:j*3+2);
        daten.Rohdaten.(name_messung).(marker{j}).ueberschrift = ueberschrift(j*3-2:j*3);
    end

    clear j strings data fid
    
    % aussuchen welche Parameter ausgewertet wird
    %----------------------------------------------------------------------    
    % zur einfachen handhabung für die einzelnen Marker Namen zu
    % ordnen
    if length(marker) == 6
        daten_stab_innen = daten.Rohdaten.(name_messung).(marker{4}).coordinates;
        daten_stab_aussen = daten.Rohdaten.(name_messung).(marker{5}).coordinates;
        daten_boxunten = daten.Rohdaten.(name_messung).(marker{3}).coordinates;
        daten_boxoben_vorne = daten.Rohdaten.(name_messung).(marker{1}).coordinates;
        daten_boxoben_hinten = daten.Rohdaten.(name_messung).(marker{2}).coordinates;
        
    elseif length(marker) == 10
        daten_stab_innen = daten.Rohdaten.(name_messung).(marker{5}).coordinates;
        daten_stab_aussen = daten.Rohdaten.(name_messung).(marker{6}).coordinates;
        daten_boxunten = daten.Rohdaten.(name_messung).(marker{3}).coordinates;
        daten_boxoben_vorne = daten.Rohdaten.(name_messung).(marker{2}).coordinates;
        daten_boxoben_hinten = daten.Rohdaten.(name_messung).(marker{1}).coordinates;
    end
    
    name = name_messung;
    path = path; 
    
    %switch wählt aus welche Funktion aufgerufen wird
    switch option
        case 'Tiefe'
            daten.Tiefe.ergebnis.name(i,1) = files;          
            daten.Tiefe.ergebnis.wert(i,1:3) = tiefe(daten_stab_innen, daten_boxunten);
        
        case 'Winkel'   
            z = str2num(files{i}(end-4));
            win = winkel(daten_stab_innen, daten_stab_aussen, daten_boxunten, files{i});
            if isfield(daten,'Winkel') == 0
                daten.Winkel.ergebniss.name(1,1) = files;
                daten.Winkel.ergebniss.wert(1,1) = win;
            elseif z == 1
                daten.Winkel.ergebniss.name(size(daten.Winkel.ergebniss.name,1)+1,1) = files;
                daten.Winkel.ergebniss.wert(size(daten.Winkel.ergebniss.wert,1)+1,1) = win;
            elseif z ~= 1    
                daten.Winkel.ergebniss.name(size(daten.Winkel.ergebniss.name,1),z) = files;
                daten.Winkel.ergebniss.wert(size(daten.Winkel.ergebniss.wert,1),z) = win;
            end
            clear win z
            
        case 'Volumen'
            daten.Volumen.ergebniss.name{i} = files;
            
            volumen = volumen(daten_stab_innen,daten_stab_aussen,daten_boxunten,daten_boxoben_vorne,daten_boxoben_hinten,name,path);
            daten.Volumen.ergebniss.netz.(name_messung) = volumen.gesamt;
            daten.Volumen.ergebniss.vol(i,1) = volumen.TotalVol;      
            daten.Volumen.ergebniss.berech_vol(i,1)  = volumen.groberrechnet;
    end
    
end


save([path '\Ergebnise\' tool '_ergebnisse.mat'],'daten');
disp('das mat-File wurde gespeichert');
end


function[tiefe] =  tiefe(daten_stab_innen, daten_boxunten)
% Tiefe bestimmen
%--------------------------------------------------------------------------
%abstand zwischen oben und unteren Boxenmarker
    tiefe = mean(daten_stab_innen) - mean(daten_boxunten)+ 20 ; % 20mm sind der abstand des Markers zur Oberfläche des Aquariums

end

function[winkel] = winkel(daten_stab_innen,daten_stab_aussen,daten_boxunten, filename)
% Gerade berechnen die die Punkte innen und außen aufspannen 
% Geradengleichung P = Pstart + g*(Pend-Pstart)
    innen_mean = mean(daten_stab_innen);
    aussen_mean = mean(daten_stab_aussen);
    boxunten_mean = mean(daten_boxunten);

    faktor = (0:0.001:1);
    gx = aussen_mean(1) + (faktor*(innen_mean(1)) - aussen_mean(1));
    gy = aussen_mean(2) + (faktor*(innen_mean(2)) - aussen_mean(2));
    gz = aussen_mean(3) + (faktor*(innen_mean(3)) - aussen_mean(3));
    gerade = [gx' gy' gz'];

    %finden den Punkt der am in Z-richtung am nächten an den Wert des unteren
    %Markes (Z3) der Box kommt
    [I Zeile] = min(abs(gerade(:,3)- boxunten_mean(3) + 20));
    schnittpunkt = gerade(Zeile,:);

    %Längen der gegenkathede und Ankathede berechnen 
    %wenn die Messungen 1 und 3 analysiert werden muss ank auf der x-achs
    %eangeben werden wenn die Messungen 2 und 4 dran sind dann auf der y-achse
     
    u = str2num(filename(end-4));
    switch u 
       case {1,3}    
            lange_ank = abs(innen_mean(1) - schnittpunkt(1));
            lange_geg = abs(innen_mean(3) - schnittpunkt(3));
       case {2,4}
           lange_ank = abs(innen_mean(2) - schnittpunkt(2));
           lange_geg = abs(innen_mean(3) - schnittpunkt(3));
    end


    % winkel berechnen über den Arctan; mit atand wird der winkel in grad
    % ausgegeben
    winkel = atand(lange_geg/lange_ank);
    
%     % 3DPlot
%     %--------------------------------------------------------------------------
%     plot3(innen_mean_x, innen_mean_y, innen_mean_z,'ob');
%     hold on 
%     plot3(ausse_mean_x, ausse_mean_y, ausse_mean_z,'*m');
%     plot3(gx,gy,gz,'k');
%     grid on;
%     plot3(boxunten_mean_x,boxunten_mean_y,boxunten_mean_z + 15,'+r');
%     plot3(schnittpunkt(1),schnittpunkt(2),schnittpunkt(3),'*g');
%     plot3(innen_mean_x, schnittpunkt(2),schnittpunkt(3),'*g');
% 
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
% 
%     %fläche
%     fx = [mean(daten_boxunten(:,1)), mean(daten_boxunten(:,1))+1000, mean(daten_boxunten(:,1))+1000, mean(daten_boxunten(:,1))];   
%     fy = [mean(daten_boxunten(:,2)), mean(daten_boxunten(:,2)), mean(daten_boxunten(:,2))+1000, mean(daten_boxunten(:,2))+1000];    
%     fz = [mean(daten_boxunten(:,3))+15 ,mean(daten_boxunten(:,3))+15,mean(daten_boxunten(:,3))+15, mean(daten_boxunten(:,3))+15 ];
% 
%     patch(fx,fy,fz,[.7 .7 .7]);           
%     hold off   
end

function volumen = volumen(daten_stab_innen,daten_stab_aussen,daten_boxunten,daten_boxoben_vorne,daten_boxoben_hinten,name,path)

%Messdaten vom inneren des Stabes ploten
fig = figure;
plot3(daten_stab_innen(:,1),daten_stab_innen(:,2),daten_stab_innen(:,3),'b.');
xlabel('x');
ylabel('y');
zlabel('z');
hold on
plot3(daten_stab_innen(1,1),daten_stab_innen(1,2),daten_stab_innen(1,3),'ro');
plot3(daten_stab_innen(end,1),daten_stab_innen(end,2),daten_stab_innen(end,3),'go');

%Auswahl von Start und Endpunkt des Kreise sowie dem obersten Punkt zur
%Angebe des Drehpunkts
datacursormode on;
dcm_obj = datacursormode(fig);
disp('Start des Kreises :');
pause
c_info = getCursorInfo(dcm_obj);
panfang=c_info.DataIndex; 
disp(panfang);

disp('Ende des Kreises :');
pause
c_info = getCursorInfo(dcm_obj);
pende=c_info.DataIndex; 
disp(pende);

disp('Anfangspunkt :');
pause
c_info = getCursorInfo(dcm_obj);
pcenter=c_info.DataIndex; 
disp(pcenter);

datacursormode off;

clear dcm_obj c_info;


% plot3(daten_stab_innen(panfang,1),daten_stab_innen(panfang,2),daten_stab_innen(panfang,3),'m*');
% plot3(daten_stab_innen(pende,1),daten_stab_innen(pende,2),daten_stab_innen(pende,3),'m*');
hold off

k1 = min(panfang,pende);
k2 = max(panfang,pende);
c = pcenter;

close all

%% Sphere einpassen
%__________________________________________________________________________
% Schnittebne festelgen 
% mit dem ausgewählten Punkt wird eine Gerade errechnet und das Center ist
% an der Schnittstelle vom unteren Marker 
Schnittebene = mean(daten_boxunten(:,3)) +20;

%[c,r] = sphereFit(daten_stab_innen(k1:k2,:));
oben = daten_stab_innen(c,:);
unten = daten_stab_aussen(c,:);

faktor = (0:0.001:1);
cx = oben(1) + (faktor*(unten(1) - oben(1)));
cy = oben(2) + (faktor*(unten(2) - oben(2)));
cz = oben(3) + (faktor*(unten(3) - oben(3)));

[~, z] = min(abs(cz - Schnittebene));
center(1) = cx(1,z);
center(2) = cy(1,z);
center(3) = cz(1,z);

clear oben unten faktor cx cy cz wert z 

hold on
plot3(center(1),center(2),center(3),'co');
xlabel('x');
ylabel('y');
zlabel('z');
% erstellen der Punktewolke von k1 bis zum letzten Punkt vor k2
%--------------------------------------------------------------------------
v = 1;
for p = k1:5:k2
    phi_punkt(p) = atan2d(daten_stab_innen(p,2)-center(2),(daten_stab_innen(p,1)-center(1)));
    lokales_r = sqrt(((daten_stab_innen(p,1)-center(1))^2)+((daten_stab_innen(p,2)-center(2))^2)+((daten_stab_innen(p,3)-center(3))^2)) ;
    
    %theta_punkt(p) = acosd(sqrt((daten_stab_innen(p,1)-center(1))^2+(daten_stab_innen(p,2)-center(2))^2)/r);
    theta_punkt_lokalesr(p) = acosd(sqrt((daten_stab_innen(p,1)-center(1))^2+(daten_stab_innen(p,2)-center(2))^2)/lokales_r);
    
    %theta = theta_punkt(p);
    theta_lr = theta_punkt_lokalesr(p);
  
    x_punkte_lr = center(1) + (lokales_r.*cosd(theta_lr)*cosd(phi_punkt(p)));
    y_punkte_lr = center(2) + (lokales_r.*cosd(theta_lr)*sind(phi_punkt(p)));
    z_punkte_lr = center(3) + lokales_r.*sind(theta_lr);
    
    
    hold on 
%     plot3(daten_stab_innen(p,1),daten_stab_innen(p,2),daten_stab_innen(p,3),'m*')
%     plot3(x_punkte,y_punkte,z_punkte,'k*');
    
    %bis auf 90 grad die Punktewolke in 2 Gradschritten hinauf rechenen
    while theta_lr <= 90   

        
        x_punkte_lr(end+1) = center(1) + (lokales_r.*cosd(theta_lr)*cosd(phi_punkt(p)));
        y_punkte_lr(end+1) = center(2) + (lokales_r.*cosd(theta_lr)*sind(phi_punkt(p)));
        z_punkte_lr(end+1) = center(3) + lokales_r.*sind(theta_lr);
        
        theta_lr = round(theta_lr)+2;
        plot3(x_punkte_lr,y_punkte_lr ,z_punkte_lr,'m*');
    end

    % die errechneten Punkte werden Spaltenweise abgespeichert
    volumen.netz_kugel.x(1:length(x_punkte_lr),v) = x_punkte_lr';
    volumen.netz_kugel.y(1:length(y_punkte_lr),v) = y_punkte_lr';
    volumen.netz_kugel.z(1:length(z_punkte_lr),v) = z_punkte_lr';
    
    v = v+1;
    clear x_punkte_lr y_punkte_lr z_punkte_lr
end

% erstellend er Punktewolke von k2, da bei der aufsummierung nicht garantiert 
% ist das der Punkt k2 getroffen wird 
phi_punkt(k2) = atan2d(daten_stab_innen(k2,2)-center(2),(daten_stab_innen(k2,1)-center(1)));
lokales_r = sqrt(((daten_stab_innen(k2,1)-center(1))^2)+((daten_stab_innen(k2,2)-center(2))^2)+((daten_stab_innen(k2,3)-center(3))^2)) ;
theta_punkt_lokalesr(k2) = acosd(sqrt((daten_stab_innen(k2,1)-center(1))^2+(daten_stab_innen(k2,2)-center(2))^2)/lokales_r); 
theta = theta_punkt_lokalesr(k2);

x_punkte_lr(1) = center(1) + (lokales_r.*cosd(theta)*cosd(phi_punkt(k2)));
y_punkte_lr(1) = center(2) + (lokales_r.*cosd(theta)*sind(phi_punkt(k2)));
z_punkte_lr(1) = center(3) + lokales_r.*sind(theta);
%plot3(x_punkte,y_punkte,z_punkte,'k*')

while theta <= 90   
    x_punkte_lr(end+1)= center(1) + (lokales_r.*cosd(theta)*cosd(phi_punkt(p)));
    y_punkte_lr(end+1) = center(2) + (lokales_r.*cosd(theta)*sind(phi_punkt(p)));
    z_punkte_lr(end+1) = center(3) + lokales_r.*sind(theta);
    
    theta = round(theta)+5;
    plot3(x_punkte_lr,y_punkte_lr ,z_punkte_lr,'r*');
    
end 

volumen.netz_kugel.x(1:length(x_punkte_lr),end+1) = x_punkte_lr';
volumen.netz_kugel.y(1:length(y_punkte_lr),end+1) = y_punkte_lr';
volumen.netz_kugel.z(1:length(z_punkte_lr),end+1) = z_punkte_lr';
    
clear x_punkte y_punkte z_punkte


% Schließen der Kugelfläche - jener Bereich der mit dem 
%--------------------------------------------------------------------------
% der Winkel phi von Punkt k1 und k2 wird erechnet und sowei auch Theta und
% dann wird die jeweilige Differenz durch 40 geteilt, was meine
% Schrittweite vorgibt

differenz_winkelphi = phi_punkt(k1)- phi_punkt(k2);
lokales_rk2 = sqrt(((daten_stab_innen(k2,1)-center(1))^2)+((daten_stab_innen(k2,2)-center(2))^2)+((daten_stab_innen(k2,3)-center(3))^2)) ;
lokales_rk1 = sqrt(((daten_stab_innen(k1,1)-center(1))^2)+((daten_stab_innen(k1,2)-center(2))^2)+((daten_stab_innen(k1,3)-center(3))^2)) ;
% der mittelwert von R(k1) und R(k2), da mit einem Radius gerechent werden
% muss aber ich durch die Interpolation 
mean_r = mean([lokales_rk1,lokales_rk2]);

theta_punkt_k1 = acosd(sqrt((daten_stab_innen(k1,1)-center(1))^2+(daten_stab_innen(k1,2)-center(2))^2)/lokales_rk1);
theta_punkt_k2 = acosd(sqrt((daten_stab_innen(k2,1)-center(1))^2+(daten_stab_innen(k2,2)-center(2))^2)/lokales_rk2);
differenz_winkeltheta = theta_punkt_k1 - theta_punkt_k2;
deltaPhi = differenz_winkelphi/40;
deltaTheta = differenz_winkeltheta/40;


v = size(volumen.netz_kugel.x,2)+1;
for o =1:40
    theta =  theta_punkt_k2 + (o*deltaTheta);
    theta_punkt(k2+(o*20)) = theta;
    phi = phi_punkt(k2) + (o*deltaPhi);
    phi_punkt(k2+(o*20)) = phi;
    
    x_punkte_lr(1)= center(1) + (mean_r.*cosd(theta)*cosd(phi));
    y_punkte_lr(1) = center(2) + (mean_r.*cosd(theta)*sind(phi));
    z_punkte_lr(1) = center(3) + mean_r.*sind(theta);
    
    interpolierte_anfangspunkte(o,:) =  [x_punkte_lr(1), y_punkte_lr(1),z_punkte_lr(1)];
    lokales_r = sqrt(((x_punkte_lr(1)-center(1))^2)+((y_punkte_lr(1)-center(2))^2)+((z_punkte_lr(1)-center(3))^2)) ;
    
    plot3(x_punkte_lr(1),y_punkte_lr(1),z_punkte_lr(1),'b*');
    
    while theta <= 90   
        x_punkte_lr(end+1)= center(1) + (lokales_r.*cosd(theta)*cosd(phi));
        y_punkte_lr(end+1) = center(2) + (lokales_r.*cosd(theta)*sind(phi));
        z_punkte_lr(end+1) = center(3) + lokales_r.*sind(theta);
        theta = round(theta)+2;
      plot3(x_punkte_lr(end),y_punkte_lr(end),z_punkte_lr(end),'b*')
    end 
    
    volumen.netz_kugel.x(1:length(x_punkte_lr),v) = x_punkte_lr';
    volumen.netz_kugel.y(1:length(y_punkte_lr),v) = y_punkte_lr';
    volumen.netz_kugel.z(1:length(z_punkte_lr),v) = z_punkte_lr';
    
    v = v+1;
    
    clear x_punkte_l y_punkte_l z_punkte_l
end

hold on 
%plot3(daten_stab_innen(k1:k2,1),daten_stab_innen(k1:k2,2),daten_stab_innen(k1:k2,3),'k^');

% da es sein kann das Punkte durch die rotationsentstehung wo sind, wo aber
% kein volumen abgefahren werden kann, müssen diese punkte mittels
% inpolygon() ermittelt und aussotiert werden. 
%--------------------------------------------------------------------------

%die Markerdaten mit den interpolierten zu einer geschlossenen Kurve führen
grundkurve(:,1) = vertcat(daten_stab_innen(k1:k2,1), interpolierte_anfangspunkte(:,1));
grundkurve(:,2) = vertcat(daten_stab_innen(k1:k2,2), interpolierte_anfangspunkte(:,2));

plot(grundkurve(:,1),grundkurve(:,2),'k');
%die Daten des Meshes zu einerm x,y,z struktur zusammen führen
gesamt = [];
for m = 1:size(volumen.netz_kugel.x,2)
    gesamt(length(gesamt)+1:length(gesamt)+length(volumen.netz_kugel.x(:,m)),1:3) = [volumen.netz_kugel.x(1:end,m),volumen.netz_kugel.y(1:end,m),volumen.netz_kugel.z(1:end,m) ];
end
gesamt = gesamt(gesamt(:,1)~=0,:);

hold off
plot3(gesamt(:,1),gesamt(:,2),gesamt(:,3),'m+');
hold on 

%hold on 
% plot(gesamt(:,1),gesamt(:,2),'bo');
[in,~] = inpolygon(gesamt(:,1),gesamt(:,2),grundkurve(:,1),grundkurve(:,2));

% plot(gesamt(in,1),gesamt(in,2),'bo');
gesamt= gesamt(in,:);

clear o l differenz_winkelphi theta_punkt_k1 theta_punkt_k2 differenz_winkeltheta deltaPhi deltaTheta theta v x_punkte y_punkte z_punkte phi_punkt theta_punkt in on lokales_rk2 lokales_rk1

%% Kegel erstellen 
%--------------------------------------------------------------------------
% Gerade berechnen die die Punkte innen und außen aufspanne 
% Geradengleichung P = Pstart + g*(Pend-Pstart)
%--------------------------------------------------------------------------

%finden den Punkt der am in Z-richtung am nächten an den Wert des unteren
%Markes (Z3) der Box kommt
faktor = (0:0.01:1);
v = 1;
for m = k1:5:k2       
    gx = daten_stab_innen(m,1) + (faktor*(daten_stab_aussen(m,1) - daten_stab_innen(m,1)));
    gy = daten_stab_innen(m,2) + (faktor*(daten_stab_aussen(m,2) - daten_stab_innen(m,2)));
    gz = daten_stab_innen(m,3) + (faktor*(daten_stab_aussen(m,3) - daten_stab_innen(m,3)));

%Zelle finden aus jeder Spalte in der der abstand von der Schnittebene am
%kleinsten ist   
    [~, z] = min(abs(gz - Schnittebene));
    gx = gx(1:z);
    gy = gy(1:z);
    gz = gz(1:z);
    
    plot3(gx,gy,gz,'g.');
    hold on 
    volumen.netz_kegel.x(1:length(gx),v) = gx';
    volumen.netz_kegel.y(1:length(gy),v) = gy';
    volumen.netz_kegel.z(1:length(gz),v) = gz';
      
    untererkreis_x(v) = gx(end);
    untererkreis_y(v) = gy(end);
    untererkreis_z(v) = gz(end);
    
    v = v+1;
    
end

gx = daten_stab_innen(k2,1) + (faktor*(daten_stab_aussen(k2,1) - daten_stab_innen(k2,1)));
gy = daten_stab_innen(k2,2) + (faktor*(daten_stab_aussen(k2,2) - daten_stab_innen(k2,2)));
gz = daten_stab_innen(k2,3) + (faktor*(daten_stab_aussen(k2,3) - daten_stab_innen(k2,3)));

%Zelle finden aus jeder Spalte in der der abstand von der Schnittebene am 
%kleinsten ist   
[~, z] = min(abs(gz - Schnittebene));
gx = gx(1:z);
gy = gy(1:z);
gz = gz(1:z);
plot3(gx,gy,gz,'g.');

volumen.netz_kegel.x(1:length(gx),end+1) = gx';
volumen.netz_kegel.y(1:length(gy),end+1) = gy';
volumen.netz_kegel.z(1:length(gz),end+1) = gz';

untererkreis_x(end+1) = gx(end);
untererkreis_y(end+1) = gy(end);
untererkreis_z(end+1) = gz(end);

unterkreis = [untererkreis_x', untererkreis_y', untererkreis_z'];
% plot3(untererkreis_x,untererkreis_y,untererkreis_z,'y.');

% schließen des Kreise am unteren Ende
%--------------------------------------------------------------------------
% plot3(untererkreis_x(end,1),untererkreis_y(end,1),untererkreis_z(end,1),'yx');
% plot3(untererkreis_x(end,end),untererkreis_y(end,end),untererkreis_z(end,end),'yx');

runten_k1 = sqrt(((untererkreis_x(end,1)-center(1))^2)+((untererkreis_y(end,1)-center(2))^2)+((untererkreis_z(end,1)-center(3))^2));
runten_k2 = sqrt(((untererkreis_x(end,end)-center(1))^2)+((untererkreis_y(end,end)-center(2))^2)+((untererkreis_z(end,end)-center(3))^2));
mean_r = mean([runten_k1,runten_k2]);

phi_unten_k1 = atan2d(untererkreis_y(end,1)-center(2),untererkreis_x(end,1)-center(1));
phi_unten_k2 = atan2d(untererkreis_y(end,end)-center(2),untererkreis_x(end,end)-center(1));
theta_unten_k1 = acosd(sqrt((untererkreis_x(end,1)-center(1))^2+(untererkreis_y(end,1)-center(2))^2)/runten_k1);  
theta_unten_k2 = acosd(sqrt((untererkreis_x(end,end)-center(1))^2+(untererkreis_y(end,end)-center(2))^2)/runten_k2);  


% wenn ein Komplexer Winkel entsteht weil die Daten nichts andres zulassen
% muss zwischen zwei VOrgehensweisen unterschieden werden. 
% 
% wenn ein komplexer Winkel besteht (if- fall) dann wird aus den unten
% ermittelntn Punkt ein Mean-Wert berechnet und auf diesen Hin die fehlende
% Segment berechnet
%
% wenn kein komplexer Winkel besteht (else -fall) wird so vorgegangen das
% das Fehlnde Kreissemgment unten in gelich viele Elemetne unterteilt wird
% wie oben und diese Punkte dann verbunden werden
if isreal(theta_unten_k1) ==0 || isreal(theta_unten_k2) ==0
   untererpunkt = mean(unterkreis);
   v = size(volumen.netz_kegel.x,2);
    for n = 1:40   
        gx = interpolierte_anfangspunkte(n,1) + (faktor*(untererpunkt(1) - interpolierte_anfangspunkte(n,1)));
        gy = interpolierte_anfangspunkte(n,2) + (faktor*(untererpunkt(2) - interpolierte_anfangspunkte(n,2)));
        gz = interpolierte_anfangspunkte(n,3) + (faktor*(untererpunkt(3) - interpolierte_anfangspunkte(n,3)));

        volumen.netz_kegel.x(1:length(gx),v) = gx';
        volumen.netz_kegel.y(1:length(gy),v) = gy';
        volumen.netz_kegel.z(1:length(gz),v) = gz';

%          plot3(gx,gy,gz,'g.');
         v=v+1;
    end
else 
  %Schließen der nicht aufgenommen Daten
    differenz_winkelphi_unten = phi_unten_k1 - phi_unten_k2;
    differenz_winkeltheta_unten = theta_unten_k1 - theta_unten_k2;

    deltaPhi = differenz_winkelphi_unten/40;
    deltaTheta = differenz_winkeltheta_unten/40;
    
    %die Anfagnspunkte
    for n = 1:40
        theta =  theta_unten_k2 + (n*deltaTheta);
        phi = phi_unten_k2 + (n*deltaPhi);

        x_punkt(n)= center(1) - (mean_r.*cosd(theta)*cosd(phi));
        y_punkt(n) = center(2) - (mean_r.*cosd(theta)*sind(phi));
        z_punkt(n) = center(3) -+ mean_r.*sind(theta);

    end

    plot3(x_punkt,y_punkt,z_punkt,'g.');
    v = size(volumen.netz_kegel.x,2);
    for n = 1:40   
        gx = interpolierte_anfangspunkte(n,1) + (faktor*(x_punkt(n) - interpolierte_anfangspunkte(n,1)));
        gy = interpolierte_anfangspunkte(n,2) + (faktor*(y_punkt(n) - interpolierte_anfangspunkte(n,2)));
        gz = interpolierte_anfangspunkte(n,3) + (faktor*(z_punkt(n) - interpolierte_anfangspunkte(n,3)));

        volumen.netz_kegel.x(1:length(gx),v) = gx';
        volumen.netz_kegel.y(1:length(gy),v) = gy';
        volumen.netz_kegel.z(1:length(gz),v) = gz';

       
       plot3(gx,gy,gz,'g.');
        v=v+1;
    end  
end

clear v gx gy gz theta phi deltaPhi deltaTheta differenz_winkeltheta_unten differenz_winkelphi_unten theta_unten_k2 theta_unten_k1 phi_unten_k2 phi_unten_k1 c1 r1 unterkreis untererkreis_x untererkreis_y untererkreis_z Schnittebene faktor; 
clear x_punkt y_punkt z_punkt wert z untererpunkt

%abfragen ob um volumen.netz_kegel NANs sind und rauslöschen
[zeile,spalte] = find(isnan(volumen.netz_kegel.x));
volumen.netz_kegel.x(:,spalte) = [];
volumen.netz_kegel.y(:,spalte) = [];
volumen.netz_kegel.z(:,spalte) = [];

clear zeile spalte
%% Gesamtes Mesh zusammensetzten
for m = 1:size(volumen.netz_kegel.x,2)
    gesamt(length(gesamt)+1:length(gesamt)+length(volumen.netz_kegel.x(:,m)),1:3) = [volumen.netz_kegel.x(1:end,m),volumen.netz_kegel.y(1:end,m),volumen.netz_kegel.z(1:end,m) ];
end

%alle Nuller -koordinaten raus löschen
volumen.gesamt = gesamt(gesamt(:,1)~=0,:);
%clear gesamt m

clear sph;
% ein Grid mit polygons erstellen und anschließend das volumen berechnen
sph = alphaShape(volumen.gesamt);
volumen.TotalVol = volume(sph);

clear sph;

%% Theoretisches Volumen berechnen 

%kleinsten Wert des Inneren Stabs finden 
[~, stabunten] = min(daten_stab_innen(:,3));

% Kegel bis zum kleinsten Wert der Messdaten 
% V = (1/3)*pi*r^2*h
kegelhohe = daten_stab_innen(stabunten,3) - daten_boxunten(1,3) - 20;
kegelradius = sqrt((daten_stab_innen(stabunten,1)-center(1))^2 + (daten_stab_innen(stabunten,2)-center(2))^2 ); 
volumenkegel = (pi*kegelradius^2*kegelhohe)/3;

% Kugelsegment bis zum kleinsten Wert der Messdaten 
% V = (h^2*pi/3)*(3r-h)
%niedrigsten Wert Z finden von den Kugeldaten  
kugelsegmentradius = sqrt((daten_stab_innen(stabunten,1)-center(1))^2+(daten_stab_innen(stabunten,2)-center(2))^2+(daten_stab_innen(stabunten,3)-center(3))^2);

kugelhohe = center(3) + kugelsegmentradius - daten_stab_innen(stabunten,3);
volumenkugel = (((kugelhohe^2)*pi)/3)*((3*kugelsegmentradius)-kugelhohe);

volumen.groberrechnet = volumenkegel + volumenkugel;
clear w stabunten kegelhohe kegelradius volumenkegel kugelhohe volumenkugel kugelsegmentradius


%% Volumen graphisch darstellen 
box_1 = [daten_boxunten(1,1) daten_boxunten(1,2) daten_boxunten(1,3)+20];
box_2 = [daten_boxunten(1,1) daten_boxunten(1,2) daten_boxoben_hinten(1,3)-20];
box_3 = [daten_boxunten(1,1) daten_boxoben_vorne(1,2) daten_boxoben_vorne(1,3)-20];
box_4 = [daten_boxunten(1,1) daten_boxoben_vorne(1,2) daten_boxunten(1,3)+20];


box_5 = [daten_boxoben_hinten(1,1) daten_boxunten(1,2) daten_boxunten(1,3)+20];
box_7 = [daten_boxoben_vorne(1,1) daten_boxoben_vorne(1,2) daten_boxoben_vorne(1,3)-20];
box_6 = [daten_boxoben_hinten(1,1) daten_boxoben_hinten(1,2) daten_boxoben_hinten(1,3)-20];
box_8 = [daten_boxoben_vorne(1,1) daten_boxoben_vorne(1,2) daten_boxunten(1,3)+20];


fig1 = figure
plot3(volumen.gesamt(:,1),volumen.gesamt(:,2),volumen.gesamt(:,3),'c.');
hold on 
plot3(daten_stab_innen(:,1),daten_stab_innen(:,2),daten_stab_innen(:,3),'m.');
grid on;
xlabel('x');
ylabel('y');
zlabel('z');
title(name)
%aquarium volumen einzeichnen 
patch([box_1(1) box_2(1) box_3(1) box_4(1)],[box_1(2) box_2(2) box_3(2) box_4(2)],[box_1(3) box_2(3) box_3(3) box_4(3)],[0.3,0.3,0.3],'faceAlpha',0.3);
patch([box_4(1) box_8(1) box_7(1) box_3(1)],[box_4(2) box_8(2) box_7(2) box_3(2)],[box_4(3) box_8(3) box_7(3) box_3(3)],[0.3,0.3,0.3],'faceAlpha',0.3);
patch([box_5(1) box_6(1) box_7(1) box_8(1)],[box_5(2) box_6(2) box_7(2) box_8(2)],[box_5(3) box_6(3) box_7(3) box_8(3)],[0.3,0.3,0.3],'faceAlpha',0.3);
patch([box_1(1) box_5(1) box_6(1) box_2(1)],[box_1(2) box_5(2) box_6(2) box_2(2)],[box_1(3) box_5(3) box_6(3) box_2(3)],[0.3,0.3,0.3],'faceAlpha',0.3);
patch([box_2(1) box_6(1) box_7(1) box_3(1)],[box_2(2) box_6(2) box_7(2) box_3(2)],[box_2(3) box_6(3) box_7(3) box_3(3)],[0.3,0.3,0.3],'faceAlpha',0.3);
patch([box_1(1) box_5(1) box_8(1) box_4(1)],[box_1(2) box_5(2) box_8(2) box_4(2)],[box_1(3) box_5(3) box_8(3) box_4(3)],[0.3,0.3,0.3],'faceAlpha',0.3);
hold off  

if isdir([path '\Ergebnise']) == 0
    mkdir([path '\Ergebnise'])
end

save 
namegraphik1 = [path '\Ergebnise\' name '_pre.png'];
namegraphik2 = [path '\Ergebnise\' name '_pre.fig'];
%     %für png speichern
X = 42.0;                  %# A3 paper size

Y = 29.7;                  %# A3 paper size
set(gcf,'PaperUnits','centimeters','PaperType','A3')
print(fig1, '-dpng', '-r400',  namegraphik1)
savefig(namegraphik2);

clear fig1 namegraphik1 X Y

close all

end



