#include "cgnslib.h"
#include <stdio.h>
#include <stdlib.h>
#include "string.h"

int main(int argc, char *argv[])
{
  float time;
  
  int*ijk;
  int ijk_0;
  int i,j,k,int_iMeshPoints,int_jMeshPoints,int_kMeshPoints;

  float *pressureSave;
  float *pressureSurfaceSave;
	int noPoints;
	noPoints=10;
	ijk = (int*) calloc(noPoints,sizeof(int));

	int_iMeshPoints=2048;
	int_jMeshPoints=256;
	int_kMeshPoints=64;
	
	j=50;
	k=32;	
	
	i=654;
	ijk[0]=k*int_jMeshPoints*int_iMeshPoints+j*int_iMeshPoints+i;
	i=687;
	ijk[1]=k*int_jMeshPoints*int_iMeshPoints+j*int_iMeshPoints+i;
	i=722;
	ijk[2]=k*int_jMeshPoints*int_iMeshPoints+j*int_iMeshPoints+i;
	i=766;
	ijk[3]=k*int_jMeshPoints*int_iMeshPoints+j*int_iMeshPoints+i;
	i=879;
	ijk[4]=k*int_jMeshPoints*int_iMeshPoints+j*int_iMeshPoints+i;
	i=1188;
	ijk[5]=k*int_jMeshPoints*int_iMeshPoints+j*int_iMeshPoints+i;
	i=1498;
	ijk[6]=k*int_jMeshPoints*int_iMeshPoints+j*int_iMeshPoints+i;
	i=1734;
	ijk[7]=k*int_jMeshPoints*int_iMeshPoints+j*int_iMeshPoints+i;
	i=1827;
	ijk[8]=k*int_jMeshPoints*int_iMeshPoints+j*int_iMeshPoints+i;
	i=1913;
	ijk[9]=k*int_jMeshPoints*int_iMeshPoints+j*int_iMeshPoints+i;
  
	FILE * file0;
	FILE * file1;
	char output_file[200];
	char output_file_2[200];
	int t;

	sprintf(output_file,"extracted_PressureHistory_BAC_VG.dat");
	sprintf(output_file_2,"surfacePressure_BAC_VG.dat");


	file0=fopen(output_file,"w");
	fprintf(file0,"TITLE = \"PressureHistory\"\n");
	fprintf(file0,"VARIABLES = \"Time (S)\" \"Pressure\"\n");
	
	file1=fopen(output_file_2,"w");
	fprintf(file1,"TITLE = \"PressureHistory\"\n");
	fprintf(file1,"VARIABLES = \"CoordinateX\" \"Pressure\"\n");
	
	
	int buffer;

	float* x;
	float* y;
	float* z;

	// 	float* u;
	// 	float* v;
	// 	float* w;
	float* p;
	// 	float* rho;

  
	char *text;
	char name[33];
	int ndescriptors,D;
	char BaseIterName[30];
	char zonename[30];
	char chr_FilmPath[200];
	int celldim,int_meshDimensions;
	int index_file,index_base,index_zone;
	int int_actualIteration;
	double dbl_time_dim;
	int index_flow,max_index_flow;
	int actual_FilmFile,numberFilmFiles;
	
	numberFilmFiles=21;

	int_meshDimensions=3;
	
	for(actual_FilmFile=1;actual_FilmFile<=numberFilmFiles;actual_FilmFile++)
	{
	      sprintf(chr_FilmPath,"/home/mag/RESULTS/VG/BAC_VG/16384CPU/%dk-%dk/Film_3d_bac_2048x256x128_VG_Re1E6_periodic_ct_16384CPU_sdrc.cgns",(125+(actual_FilmFile-1)*30),(125+(actual_FilmFile)*30));
	      printf("File: %s\n",chr_FilmPath);
	      cg_open(chr_FilmPath,CG_MODE_READ,&index_file);
	      index_base=1;
	      cg_cell_dim(index_file, index_base, &celldim);

	      cgsize_t isize[3][celldim],irmin[celldim],irmax[celldim];
	      
	      index_zone=1;

	      cg_zone_read(index_file,index_base,index_zone,zonename,isize[0]);

	      /* lower range index */
	      irmin[0]=1;
	      irmin[1]=1;
	      /* upper range index of vertices */
	      irmax[0]=isize[0][0];
	      irmax[1]=isize[0][1];

	      if (int_meshDimensions==3)
	      {
		      irmin[2]=1;
		      irmax[2]=isize[0][2];
	      }
	      
	      if (int_meshDimensions==2)
	      {
		      buffer=irmax[0]*irmax[1];
	      }
	      else
	      {
		      buffer=irmax[0]*irmax[1]*irmax[2];
	      }
	      

// 	      cg_goto(index_file,index_base,"end");
// 	      cg_ndescriptors(&ndescriptors);
// 	      for(D=1;D<=ndescriptors;D++)
// 	      {
// 	      cg_descriptor_read(D, name, &text);
// 	      if (strcmp(name,"Iterations")==0)
// 	      {
// 		int_actualIteration = atoi( text );
// 	      }
// 	      if (strcmp(name,"ActualTime")==0)
// 	      {
// 		dbl_time_dim = strtod( text,NULL );
// 	      }
// 	      }
	      
	      cg_biter_read(index_file, index_base, BaseIterName, &max_index_flow);
	      printf("Zone:%s bei t=%g und Iteration=%d mit %d timesteps wird geladen.\n",zonename,dbl_time_dim,int_actualIteration,max_index_flow);
	      
	      p = (float*) calloc(buffer,sizeof(float));
		x = (float*) calloc(buffer,sizeof(float));
		y = (float*) calloc(buffer,sizeof(float));
		z = (float*) calloc(buffer,sizeof(float));
	// 	u = (float*) calloc(buffer,sizeof(float));
	// 	v = (float*) calloc(buffer,sizeof(float));
	// 	w = (float*) calloc(buffer,sizeof(float));
	// 	rho = (float*) calloc(buffer,sizeof(float));	
		
	      if(actual_FilmFile==1)
	      {		
		pressureSave = (float*) calloc(noPoints*numberFilmFiles*max_index_flow,sizeof(float));
		pressureSurfaceSave = (float*) calloc((1912-654),sizeof(float));
	      }
	      
	      /* read grid coordinates */
      // 	cg_field_read(index_file,index_base,index_zone,index_flow,"VelocityX",
      // 	RealSingle,irmin,irmax,u);
      // 	cg_field_read(index_file,index_base,index_zone,index_flow,"VelocityY",
      // 	RealSingle,irmin,irmax,v);
      // 	cg_field_read(index_file,index_base,index_zone,index_flow,"VelocityZ",
      // 	RealSingle,irmin,irmax,w);
	      
	      cg_coord_read(index_file,index_base,index_zone,"CoordinateX",RealSingle,irmin,irmax,x);
	      cg_coord_read(index_file,index_base,index_zone,"CoordinateY",RealSingle,irmin,irmax,y);
	      cg_coord_read(index_file,index_base,index_zone,"CoordinateZ",RealSingle,irmin,irmax,z);

	      for (index_flow=1;index_flow<=max_index_flow;index_flow++)
	      {
		cg_field_read(index_file,index_base,index_zone,index_flow,"Pressure",RealSingle,irmin,irmax,p);
		for(i=0;i<noPoints;i++)
		{
 		  pressureSave[i*max_index_flow*numberFilmFiles+(actual_FilmFile-1)*max_index_flow+(index_flow-1)]=p[ijk[i]];
		}
		
 		fprintf(file1,"ZONE T=\"%d\", F=POINT, I=%d, DT=(SINGLE)\n",((actual_FilmFile-1)*max_index_flow+index_flow),(1911-654+1));
		for(i=654;i<=1911;i++)
		{
		  ijk_0=k*int_jMeshPoints*int_iMeshPoints+j*int_iMeshPoints+i;
 		  fprintf(file1,"%e %e\n",x[ijk_0],p[ijk_0]);
		  pressureSurfaceSave[i-654]+=p[ijk_0];
		}
		
		
	      }

	      
    
	      cg_close(index_file);

      // 	free(u);
      // 	free(v);
      // 	free(w);
      // 	free(rho);
	}
	
	for(i=0;i<noPoints;i++)
	{
	  fprintf(file0,"ZONE T=\"Point%d(x=%f, y=%f, z=%f)\", F=POINT, I=%d, DT=(SINGLE)\n",i,x[ijk[i]],y[ijk[i]],z[ijk[i]],(max_index_flow*numberFilmFiles));
	for(actual_FilmFile=1;actual_FilmFile<=numberFilmFiles;actual_FilmFile++)
	{
	  for (index_flow=1;index_flow<=max_index_flow;index_flow++)
	 {
	  ijk_0=k*int_jMeshPoints*int_iMeshPoints+j*int_iMeshPoints+i;
	  time=((actual_FilmFile-1)*30000.+index_flow*3000.)*0.000008*0.08/263.8634;
	   fprintf(file0,"%f %f\n",time,pressureSave[i*max_index_flow*numberFilmFiles+(actual_FilmFile-1)*max_index_flow+(index_flow-1)]);
	 }
	}
	}
	
	
	fprintf(file1,"ZONE T=\"%d\", F=POINT, I=%d, DT=(SINGLE)\n",((actual_FilmFile-1)*max_index_flow+index_flow),(1911-654+1));
	for(i=654;i<=1911;i++)
	{
	  ijk_0=k*int_jMeshPoints*int_iMeshPoints+j*int_iMeshPoints+i;
	  fprintf(file1,"%e %e\n",x[ijk_0],(pressureSurfaceSave[i-654]/(max_index_flow*numberFilmFiles)));
	}
	
  fclose(file0);
   fclose(file1);
}

