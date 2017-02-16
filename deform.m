% data preparation: the raw data of the penetration length should be
% without the first line of the asigned names of each column. also, the
% first line the xlsx file should be the undeformed cell's circular shape
% parameters. the output files are text files and xls files.
clear;
cirpos=17;
presspos=1;
penpos=12;
kslide=0.08667;
outwidth=32;
unitf=6.7e-11;

modu0=0;
% this script is to calculate the elastic modulus of the cells.

global timeunit;
%please revise these two parameters mannually!!!!!!
timeunit=0.1;
type='3t3';
%dcell=14.4*1e-6
widthcol=12;
global pressure;
pressure=0.2; %the pressure of the experiment
%anglefactor=0.043333;


global pen10;
global pen15;
global pen20;
global pen25;
global pen30;
global xlen10;
global xlen15;
global xlen20;
global xlen25;
global xlen30;
global ylen10;
global ylen15;
global ylen20;
global ylen25;
global ylen30;
global force;
global deformwidth;
global filecount;
global circ;
stdpath='G:\deform\140204-CellDiaPosNForce-verB.xlsx';
pen10=xlsread(stdpath,'total F by Flow','a3:f8');
pen15=xlsread(stdpath,'total F by Flow','a12:f17');
pen20=xlsread(stdpath,'total F by Flow','a21:f29');
pen25=xlsread(stdpath,'total F by Flow','a33:f40');
pen30=xlsread(stdpath,'total F by Flow','a44:f51');
[xlen10 ylen10]=size(pen10);
[xlen15 ylen15]=size(pen15);
[xlen20 ylen20]=size(pen20);
[xlen25 ylen25]=size(pen25);
[xlen30 ylen30]=size(pen30);
% 
% pen10(:,4)
deformfiles = dir('*.xlsx');
filecount = length(deformfiles);
force=cell(filecount);
deformwidth=cell(filecount);
i = 0;
%ubfile=fopen('unit_area_bond.txt','w');

%cbfile=fopen('unit_cell_bond.txt','w');
ftotal=fopen('total.txt','w');
for i=1:filecount
    a=xlsread(deformfiles(i).name);
    [xlen,ylen]=size(a)
%    cir=a(1,cirpos);
 f=fopen(strcat(num2str(i),'.txt'),'w'); 
   for j=1:(xlen/2)
    cirtemp=a(2*j-1,penpos);
    cirpen=a(2*j,penpos);
    cir1=32-0.08667*cirpen;
    cir=(1.5*cirtemp*cirtemp*cir1-0.5*cir1*cir1*cir1).^(1/3);
    %circ(i)=cir;
    test(i)=outwidth-kslide*a(2,penpos)-cir;
      
    f0=exforce(cir,cirpen);
    modu0=pressure*f0/(8*0.04333*sqrt(cir*(cir-cir1)*(cir-cir1)*(cir-cir1))/9);
    bility=cirtemp/cir1/f0/pressure;
    fprintf(f,'%f\t%f\t%f\t%f\n',cir,cirpen,bility,modu0);
    fprintf(ftotal,'%f\t%f\t%f\t%f\n',cir,cirpen,bility,modu0);
%     pressep=a(2*j-1,1)/1.5;
%     penf=300-sqrt(2/3*(cir*cir*cir/cir1-cir1*cir1))/2;
%     f1=pressep*exforce(cir,penf);
%     cir2=30-0.08667*penf;
%     f0=modu0*8*0.04333*sqrt(cir*(cir-cir2)*(cir-cir2)*(cir-cir2))/9;
%     df=f1-f0;
%     contact=3.1416*(cir*cir*cir/6/cir2-cir2*cir2/6)*2;
%     bondc=df*1e-9/unitf;
%     unitbond=bondc/contact;
%     cellbond=unitbond*3.1416*cir*cir;
%     fprintf(ubfile,'%f\t%f\t%f\n',cir,modu0,unitbond);
%     fprintf(cbfile,'%f\t%f\t%f\n',cir,modu0,cellbond);
   end
%     for j=1:(xlen-1)
%         pen=a(j+1,penpos);
%         force{i}(j)=exforce(cir,pen);
%        deformwidth{i}(j,1)=outwidth-kslide*pen;
%         fprintf(f,'%f\n',force{i,1}(j));
%     end
    %xlswrite(strcat(num2str(i),'\','force',num2str(i),'.xls'),force{i});
    %xlswrite(strcat(num2str(i),'\','deformwidth',num2str(i),'.xlsx'),deformwidth{i,1});
    fclose(f);

    
    
end
fclose(ftotal);
% ubfile=fopen('unit_area_bond.txt','w');
% cbfile=fopen('unit_cell_bond.txt','w');
% for i=1:filecount
%     f=force{i};
%     cir=circ(i);
%     a=xlsread(deformfiles(i).name);
%     pen=a(2,penpos);
%     width0=30-kslide*pen;
%     f0=8*0.04333/2*modu0*sqrt(cir*(cir-width0)*(cir-width0)*(cir-width0))/9;
%     addf=f*pressure-f0;
%     contact=3.1416*(cir*cir*cir/6/width0-width0*width0/6);
%     bondc=addf*1e-9/unitf;
%     unitbond=bondc/contact;
%     cellbond=unitbond*3.1416*cir*cir;
%     fprintf(ubfile,'%f\n',unitbond);
%     fprintf(cbfile,'%f\n',cellbond);
% end

% fclose(ubfile);
% fclose(cbfile);


    