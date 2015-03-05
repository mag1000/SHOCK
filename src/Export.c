#include "Import.h"
#include "SHOCK.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cgnslib.h"
#include "mpi.h"
#include "string.h"
//make and remove directories
#include <unistd.h>
#include <sys/types.h>  /* Linux/UNIX */
#include <sys/stat.h>   /* Linux/UNIX */

void SnapshotExportEachCPU(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep)
{
	char * actual_file;
	actual_file=malloc(50*sizeof(char));
	sprintf(actual_file,"SnapShotGC_%d.dat",pnt_config->MPI_rank);

	int i,j,k,ijk;
	FILE *file0;
	file0=fopen(actual_file,"w");
	fprintf(file0,"variables = x, y, z, rho, u, v, w, p, e, gradRho \n");

	if(MESHDIMENSIONS==2)
	{
		fprintf(file0,"zone t='1', i= %d, j= %d, f=point \n",pnt_config->int_iMeshPointsGhostCells,pnt_config->int_jMeshPointsGhostCells);
	}
	if(MESHDIMENSIONS==3)
	{
		fprintf(file0,"zone t='1', i= %d, j= %d, k= %d, f=point \n",pnt_config->int_iMeshPointsGhostCells,pnt_config->int_jMeshPointsGhostCells,pnt_config->int_kMeshPointsGhostCells);
	}

	for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				fprintf(file0," %le %le %le %le %le %le %le %le %le %le\n",
						pnt_mesh->x[ijk],
						pnt_mesh->y[ijk],
						pnt_mesh->z[ijk],
						pnt_U_lastStep->rho[ijk],
						pnt_U_lastStep->u[ijk],
						pnt_U_lastStep->v[ijk],
						pnt_U_lastStep->w[ijk],
						pnt_U_lastStep->p[ijk],
						pnt_U_lastStep->e[ijk],
						pnt_U_lastStep->gradRho[ijk]);
			}
		}
	}

	fclose(file0);
	free(actual_file);
}



void MeshMetricExport(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
//	double *buffer = calloc(pnt_config->int_iStartGhosts * pnt_config->int_jStartGhosts, sizeof(double));


	char actual_file[200];

	FILE * file0;
	int i,j,k,ijk;


	char foldername[200];
	sprintf(foldername,"%sMetric/",pnt_config->chr_folder);
	if (pnt_config->MPI_rank==0)
	{
//		rmdir(foldername);
		mkdir(foldername,0755);
	}
	MPI_Barrier(pnt_config->MPI_comm);


	sprintf(actual_file,"%s/Metric/MeshMetric_Zone%d.dat",pnt_config->chr_folder,(pnt_config->MPI_rank+1));

	file0=fopen(actual_file,"w");


	fprintf(file0,"variables = x, y, z, xi_x, xi_y, xi_z, eta_x, eta_y, eta_z, zeta_x, zeta_y, zeta_z, jacobian, u, v, w, p, rho\n");
	fprintf(file0,"zone t='1', k= %d, j= %d, i= %d, f=point \n",pnt_config->int_kMeshPointsGhostCells,pnt_config->int_jMeshPointsGhostCells,pnt_config->int_iMeshPointsGhostCells);


	for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				fprintf(file0," %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",
						pnt_mesh->x[ijk],
						pnt_mesh->y[ijk],
						pnt_mesh->z[ijk],
						pnt_mesh->xi_x[ijk],
						pnt_mesh->xi_y[ijk],
						pnt_mesh->xi_z[ijk],
						pnt_mesh->eta_x[ijk],
						pnt_mesh->eta_y[ijk],
						pnt_mesh->eta_z[ijk],
						pnt_mesh->zeta_x[ijk],
						pnt_mesh->zeta_y[ijk],
						pnt_mesh->zeta_z[ijk],
						pnt_mesh->jacobian[ijk],
						pnt_U->u[ijk],
						pnt_U->v[ijk],
						pnt_U->w[ijk],
						pnt_U->p[ijk],
						pnt_U->rho[ijk]);

			}
		}
	}
	fclose(file0);


//	sprintf(actual_file,"%s/Metric/MeshExtrapolate_Zone%d.dat",pnt_config->chr_folder,(pnt_config->MPI_rank+1));
//
//	file0=fopen(actual_file,"w");
//	fprintf(file0,"variables = x_extrapolate, y_extrapolate, z_extrapolate\n");
//	fprintf(file0,"zone t='1', k= %d, j= %d, i= %d, f=point \n",pnt_config->int_kMeshPointsGhostCells_extrapolate,pnt_config->int_jMeshPointsGhostCells_extrapolate,pnt_config->int_iMeshPointsGhostCells_extrapolate);
//
//
//	for (k=pnt_config->int_kStartGhosts_extrapolate; k <= pnt_config->int_kEndGhosts_extrapolate; k++)
//	{
//		for (j=pnt_config->int_jStartGhosts_extrapolate; j <= pnt_config->int_jEndGhosts_extrapolate; j++)
//		{
//			for (i=pnt_config->int_iStartGhosts_extrapolate; i <= pnt_config->int_iEndGhosts_extrapolate; i++)
//			{
//				ijk=i*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+j*pnt_config->int_kMeshPointsGhostCells_extrapolate+k;
//
//				fprintf(file0," %le %le %le\n",
//						pnt_mesh->x_extrapolate[ijk],
//						pnt_mesh->y_extrapolate[ijk],
//						pnt_mesh->z_extrapolate[ijk]);
//
//			}
//		}
//	}
//
//	fclose(file0);

}


void CGNS_PressureHistoryValuesExportParallel(
		struct strct_configuration * pnt_config,
		int p_index)
{
	int index_file,icelldim,iphysdim,index_base;
	int index_zone,index_coord,index_field,index_flow;


	char basename[33];
	char zonename[33];


	int t,nsteps;
	int buffer=1;

	cgsize_t idata[2],nuse;


	//    timedepended CGNS-DATA
	double time[pnt_config->int_TotalIterations];
	char sn[pnt_config->int_TotalIterations][33];
	char solname[pnt_config->int_TotalIterations*32+1];  /* need an extra byte for the terminating 0 */

	strcpy(solname,"");

	double* x;
	double* y;
	double* z;
	double* p;

//Speicherallokierung für die dynamischen-1D Arrays, worin die 3D-Lösungen gespeichert werden
	x = (double*) calloc(buffer,sizeof(double));
	y = (double*) calloc(buffer,sizeof(double));
	z = (double*) calloc(buffer,sizeof(double));
	p = (double*) calloc(buffer,sizeof(double));

	for(t=0;t<pnt_config->int_TotalIterations;t++)
	{
		time[t]=pnt_config->start_Time+(double)(t+1)*pnt_config->int_IterationsBetweenSamples*
				(pnt_config->dbl_L0_dim*pnt_config->dbl_numericalTau/pnt_config->dbl_u0_dim);

		sprintf(sn[t],"FlowSolution%d",t);
		sprintf(solname,"%s%-32s",solname,sn[t]);
	}

		x[0] = pnt_config->PressureHistory_x_P_real[p_index];
		y[0] = pnt_config->PressureHistory_y_P_real[p_index];
		z[0] = pnt_config->PressureHistory_z_P_real[p_index];

		icelldim=pnt_config->int_meshDimensions;
	   cgsize_t isize[3][icelldim];

	/* WRITE X, Y, Z GRID POINTS TO CGNS FILE */
	/* open CGNS file for write */

	/* create base (user can give any name) */
	   sprintf(basename,"Base");

	   iphysdim=3;
	   if(p_index==0)
	   {
		   cg_open(pnt_config->chr_PressureHistoryPath,CG_MODE_WRITE,&index_file);
		   cg_base_write(index_file,basename,icelldim,iphysdim,&index_base);
	   }
	   else
	   {
		   cg_open(pnt_config->chr_PressureHistoryPath,CG_MODE_MODIFY,&index_file);
		   index_base=1;
	   }


	   cg_simulation_type_write(index_file,index_base,TimeAccurate);

	/* vertex size */
	   isize[0][0]=1;
	   isize[0][1]=1;

	/* cell size */
	   isize[1][0]=isize[0][0]-1;
	   isize[1][1]=isize[0][1]-1;

	/* boundary vertex size (always zero for structured grids) */
	   isize[2][0]=0;
	   isize[2][1]=0;

	   if (pnt_config->int_meshDimensions==3)
	   {
		   isize[0][2]=1;
		   isize[1][2]=isize[0][2]-1;
		   isize[2][2]=0;
	   }

	   sprintf(zonename,"Point%d",p_index);
	   /* create zone */
	   cg_zone_write(index_file,index_base,zonename,*isize,Structured,&index_zone);
	/* write grid coordinates (user must use SIDS-standard names here) */
	   cg_coord_write(index_file,index_base,index_zone,RealSingle,"CoordinateX",
		   x,&index_coord);
	   cg_coord_write(index_file,index_base,index_zone,RealSingle,"CoordinateY",
		   y,&index_coord);
	   cg_coord_write(index_file,index_base,index_zone,RealSingle,"CoordinateZ",
		   z,&index_coord);

	for(t=0;t<pnt_config->int_TotalIterations;t++)
	{
		p[0] = pnt_config->PressureHistory_pressure[p_index][t];

		cg_sol_write(index_file,index_base,index_zone,sn[t],Vertex,&index_flow);

	   /* write flow solution (user must use SIDS-standard names here) */
	   cg_field_write(index_file,index_base,index_zone,index_flow,
			   RealSingle,"Pressure",p,&index_field);
	}


	   /* create BaseIterativeData */
		   nsteps=pnt_config->int_TotalIterations;
		   cg_biter_write(index_file,index_base,"TimeIterValues",nsteps);
	   /* go to BaseIterativeData level and write time values */
		   cg_goto(index_file,index_base,"BaseIterativeData_t",1,"end");
		   nuse=pnt_config->int_TotalIterations;
		   cg_array_write("TimeValues",RealSingle,1,&nuse,&time);
	   /* create ZoneIterativeData */
		   cg_ziter_write(index_file,index_base,index_zone,"ZoneIterativeData");
	   /* go to ZoneIterativeData level and give info telling which */
	   /* flow solution corresponds with which time (solname(1) corresponds */
	   /* with time(1), solname(2) with time(2), and solname(3) with time(3)) */
		   cg_goto(index_file,index_base,"Zone_t",index_zone,"ZoneIterativeData_t",1,"end");
		   idata[0]=32;
		   idata[1]=pnt_config->int_TotalIterations;
		   cg_array_write("FlowSolutionPointers",Character,2,idata,solname);



	free(x);
	free(y);
	free(z);
	free(p);

	cg_close(index_file);



	}

void ASCii_PressureHistoryValuesExportParallel(
		struct strct_configuration * pnt_config,
		int p_index)
{
	char actual_file[200];

	FILE * file0;
	int t;

	sprintf(actual_file,"%sPressureHistory_Points_Iteration%d.dat",pnt_config->chr_folder,(pnt_config->int_actualIteration-1));

	if(p_index==0)
	{
		file0=fopen(actual_file,"w");
		fprintf(file0,"TITLE = \"PressureHistory von %s\"\n",
				pnt_config->chr_MeshFile);
		fprintf(file0,"VARIABLES = \"Time (S)\" \"Pressure\"\n");
	}
	else
	{
		file0=fopen(actual_file,"a");
	}


	fprintf(file0,"ZONE T=\"Point%d(x=%f, y=%f, z=%f)\", F=POINT, I=%d, DT=(SINGLE)\n",
			p_index,
			pnt_config->PressureHistory_x_P_real[p_index],
			pnt_config->PressureHistory_y_P_real[p_index],
			pnt_config->PressureHistory_z_P_real[p_index],
			pnt_config->int_TotalIterations);

	for(t=0;t<pnt_config->int_TotalIterations;t++)
	{
		fprintf(file0,"%le %le\n",
				pnt_config->PressureHistory_time[t],
				pnt_config->PressureHistory_pressure[p_index][t]);
	}


	fclose(file0);
	}


void ASCii_VelocityHistoryValuesExportParallel(
		struct strct_configuration * pnt_config,
		int p_index)
{
	char actual_file[200];

	FILE * file0;
	int t;

	sprintf(actual_file,"%sVelocityHistory_Points_Iteration%d.dat",pnt_config->chr_folder,(pnt_config->int_actualIteration-1));

	if(p_index==0)
	{
		file0=fopen(actual_file,"w");
		fprintf(file0,"TITLE = \"VelocityHistory von %s\"\n",
				pnt_config->chr_MeshFile);
		fprintf(file0,"VARIABLES = \"Time (S)\" \"VelocityX\" \"VelocityY\" \"VelocityZ\"\n");
	}
	else
	{
		file0=fopen(actual_file,"a");
	}


	fprintf(file0,"ZONE T=\"Point%d(x=%f, y=%f, z=%f)\", F=POINT, I=%d, DT=(SINGLE)\n",
			p_index,
			pnt_config->VelocityHistory_x_P_real[p_index],
			pnt_config->VelocityHistory_y_P_real[p_index],
			pnt_config->VelocityHistory_z_P_real[p_index],
			pnt_config->int_TotalIterations);

	for(t=0;t<pnt_config->int_TotalIterations;t++)
	{
		fprintf(file0,"%le %le %le %le\n",
				pnt_config->VelocityHistory_time[t],
				pnt_config->VelocityHistory_VelocityX[p_index][t],
				pnt_config->VelocityHistory_VelocityY[p_index][t],
				pnt_config->VelocityHistory_VelocityZ[p_index][t]);
	}


	fclose(file0);
	}

void NANExport(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
//	double *buffer = calloc(pnt_config->int_iStartGhosts * pnt_config->int_jStartGhosts, sizeof(double));


	char actual_file[200];

	FILE * file0;
	int i,j,k,ijk;


	sprintf(actual_file,"%s/NAN-Export_Zone%d.dat",pnt_config->chr_folder,(pnt_config->MPI_rank+1));

	file0=fopen(actual_file,"w");


	fprintf(file0,"variables = x, y, z, xi_x, xi_y, xi_z, eta_x, eta_y, eta_z, zeta_x, zeta_y, zeta_z, jacobian, u, v, w, p, rho\n");
	fprintf(file0,"zone t='1', k= %d, j= %d, i= %d, f=point \n",pnt_config->int_kMeshPointsGhostCells,pnt_config->int_jMeshPointsGhostCells,pnt_config->int_iMeshPointsGhostCells);


	for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				fprintf(file0," %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",
						pnt_mesh->x[ijk],
						pnt_mesh->y[ijk],
						pnt_mesh->z[ijk],
						pnt_mesh->xi_x[ijk],
						pnt_mesh->xi_y[ijk],
						pnt_mesh->xi_z[ijk],
						pnt_mesh->eta_x[ijk],
						pnt_mesh->eta_y[ijk],
						pnt_mesh->eta_z[ijk],
						pnt_mesh->zeta_x[ijk],
						pnt_mesh->zeta_y[ijk],
						pnt_mesh->zeta_z[ijk],
						pnt_mesh->jacobian[ijk],
						pnt_U->u[ijk],
						pnt_U->v[ijk],
						pnt_U->w[ijk],
						pnt_U->p[ijk],
						pnt_U->rho[ijk]);

			}
		}
	}
	fclose(file0);

}
