function PlotBox(Box)
Box_vertex=[
            -Box(1)       -Box(2)       -Box(3)      
            Box(1)        -Box(2)       -Box(3)        
            Box(1)         Box(2)       -Box(3)        
            -Box(1)        Box(2)       -Box(3)       
            -Box(1)       -Box(2)       Box(3)   
            Box(1)        -Box(2)       Box(3)   
            Box(1)         Box(2)       Box(3)   
            -Box(1)        Box(2)       Box(3)   
        ];
plot3( Box_vertex(1:4,1),Box_vertex(1:4,2),Box_vertex(1:4,3) );hold on;
plot3( Box_vertex(5:8,1),Box_vertex(5:8,2),Box_vertex(5:8,3) );hold on;
for i=1:4
    plot3( Box_vertex([i,i+4],1),Box_vertex([i,i+4],2),Box_vertex([i,i+4],3) );hold on;
end
for i=[1,5]
    plot3( Box_vertex([i,i+3],1),Box_vertex([i,i+3],2),Box_vertex([i,i+3],3) );hold on;
end
%geometry of simulation box
