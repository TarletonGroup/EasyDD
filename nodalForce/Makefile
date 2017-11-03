#all: sq dc qdc
#sq: SegQuadTrianForces_final_2.2_quad_debugg2.c
#	gcc -o sq.out SegQuadTrianForces_final_2.2_quad_debugg2.c -lquadmath
#dc: nodal_surface_force_quad_triangle.c
#	gcc -o dc.out nodal_surface_force_quad_triangle.c -lm
#qdc: nodal_surface_force_quad_triangle.c
#	gcc -Wall -o qdc.out nodal_surface_force_quad_triangle_qp.c -lquadmath
#clean:
#	rm -f *.out
#-Wunused-but-set-variable
#all: qdc
#qdc: nodal_surface_force_quad_triangle.c
#	gcc -Wall -Wno-format -o qdc.out nodal_surface_force_quad_triangle_qp.c -lquadmath
debugflags := -Ddebug2 #-Ddebug2 #-Ddebug
NVCCFLAGS	:= -lineinfo -arch=sm_50 --ptxas-options=-v --use_fast_math
gccflags :=  -lm -Wall
optflags := -O3
#-arch=sm_35

all:	c_lin_rect lin_rect
c_lin_rect:	cuda_nodal_surface_force_linear_rectangle.cu Makefile
	nvcc cuda_nodal_surface_force_linear_rectangle.cu -o c_lin_rect.out $(NVCCFLAGS) $(debugflags) $(optflags)
lin_rect: nodal_surface_force_linear_rectangle.c Makefile
	gcc nodal_surface_force_linear_rectangle.c -o lin_rect.out $(gccflags) $(debugflags) $(optflags)
clean:
	rm -f *.o *.out
