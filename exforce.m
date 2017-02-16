function y=exforce(cir,pen)

stdpath='G:\deform\140204-CellDiaPosNForce-verB.xlsx';
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

forcecol=4;
rowstart=[3 12 21 33 44];
rowend=[8 17 29 40 51];

pen10(:,forcecol)
if cir<=15
   
    %calculate pen10
    indict10=1;
    for i=1:xlen10
       if pen<=pen10(i,1) 
           break;
       end
       indict10=indict10+1;
    end
    
    if indict10==1
        indict10=2
    end
    
    if indict10==7
        indict10=6
    end
    force10=ex(pen,pen10(indict10-1,1),pen10(indict10,1),pen10(indict10-1,forcecol),pen10(indict10,forcecol));
    
    %calculate pen15
    indict15=1;
    for i=1:xlen15
       if pen<=pen15(i,1) 
           break;
       end
       indict15=indict15+1;
    end
    if indict15==1
        indict15=2
    end
    if indict15>6
        indict15=6
    end
    force15=ex(pen,pen15(indict15-1,1),pen15(indict15,1),pen15(indict15-1,forcecol),pen15(indict15,forcecol));
    
    %calculate final
    y=ex(cir,10,15,force10,force15);
    
end

if cir>15&&cir<=20
    
    %calculate pen15
    indict15=1;
    for i=1:xlen15
       if pen<=pen15(i,1) 
           break;
       end
       indict15=indict15+1;
    end
    if indict15==1
        indict15=2
    end
    if indict15==7
        indict15=6
    end
    force15=ex(pen,pen15(indict15-1,1),pen15(indict15,1),pen15(indict15-1,forcecol),pen15(indict15,forcecol));
    
     %calculate pen20
    indict20=1;
    for i=1:xlen20
       if pen<=pen20(i,1) 
           break;
       end
       indict20=indict20+1;
    end
    if indict20==1
        indict20=2
    end
    if indict20>9
        indict20=9
    end
    force20=ex(pen,pen20(indict20-1,1),pen20(indict20,1),pen20(indict20-1,forcecol),pen20(indict20,forcecol));
    
    %calculate final
    y=ex(cir,15,20,force15,force20);
    
end

if cir>20 && cir<=25
    
     %calculate pen20
    indict20=1;
    for i=1:xlen20
       if pen<=pen20(i,1) 
           break;
       end
       indict20=indict20+1;
    end
    if indict20==1
        indict20=2
    end
    if indict20==10
        indict20=9
    end
    force20=ex(pen,pen20(indict20-1,1),pen20(indict20,1),pen20(indict20-1,forcecol),pen20(indict20,forcecol));
    
     %calculate pen25
    indict25=1;
    for i=1:xlen25
       if pen<=pen25(i,1) 
           break;
       end
       indict25=indict25+1;
    end
    if indict25==1
        indict25=2
    end
    if indict25>8
        indict25=8
    end
    force25=ex(pen,pen25(indict25-1,1),pen25(indict25,1),pen25(indict25-1,forcecol),pen25(indict25,forcecol));
    
    %calculate final
    y=ex(cir,20,25,force20,force25);
    
end


if cir>25
    %calculate pen25
    indict25=1;
    for i=1:xlen25
       if pen<=pen25(i,1) 
           break;
       end
       indict25=indict25+1;
    end
    if indict25==1
        indict25=2;
    end
    if indict25==9
        indict25=8;
    end
    force25=ex(pen,pen25(indict25-1,1),pen25(indict25,1),pen25(indict25-1,forcecol),pen25(indict25,forcecol));
    
     %calculate pen30
    indict30=1;
    for i=1:xlen30
       if pen<=pen30(i,1) 
           break;
       end
       indict30=indict30+1;
    end
    if indict30==1
        indict30=2
    end
    if indict30==9
        indict30=8
    end
    force30=ex(pen,pen30(indict30-1,1),pen30(indict30,1),pen30(indict30-1,forcecol),pen30(indict30,forcecol));
    
    %calculate final
    y=ex(cir,25,30,force25,force30);
end


end

    