function PlotBox_b(Box)
Box_vertex=[
            0       0       0      
            Box(1)        0       0        
            Box(1)         Box(2)       0        
            0        Box(2)       0       
            0       0       Box(3)   
            Box(1)        0       Box(3)   
            Box(1)         Box(2)       Box(3)   
            0        Box(2)       Box(3)   
        ];
plot3( Box_vertex(1:4,1),Box_vertex(1:4,2),Box_vertex(1:4,3),'k','LineWidth',2 );hold on;
plot3( Box_vertex(5:8,1),Box_vertex(5:8,2),Box_vertex(5:8,3),'k','LineWidth',2 );hold on;
for i=1:4
    plot3( Box_vertex([i,i+4],1),Box_vertex([i,i+4],2),Box_vertex([i,i+4],3),'k','LineWidth',2 );hold on;
end
for i=[1,5]
    plot3( Box_vertex([i,i+3],1),Box_vertex([i,i+3],2),Box_vertex([i,i+3],3),'k','LineWidth',2 );hold on;
end
%geometry of simulation box
