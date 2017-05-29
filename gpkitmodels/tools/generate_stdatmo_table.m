% Generates a lookup table of standard atmospheric data. SI units.

h = (0:100:90000)';
[rho,a,T,P,kvisc] = stdatmo(h,0,'SI');

fileID = fopen('stdatmo_table.txt','w');
fprintf(fileID,'Standard-Atmospheric Data Table (SI Units)');
fprintf(fileID,'\n');
fprintf(fileID,'h (m)\trho (kg/m^3)\ta (m/s)\tT (K)\tP (Pa)\tkvisc (m^2/s)');
fprintf(fileID,'\n');
for i = 1:1:length(h)
    fprintf(fileID,'%0.1f\t%0.10f\t%0.2f\t%0.2f\t%6.5f\t%0.9f',...
        h(i),rho(i),a(1),T(i),P(i),kvisc(i));
    fprintf(fileID,'\n');
end

fclose(fileID);


