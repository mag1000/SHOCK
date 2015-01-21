#include "Import.h"
#include "saiwenos.h"
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

void SnapshotExport(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep)
{
	char * actual_file;

	actual_file=malloc(50*sizeof(char));
	sprintf(actual_file,"Snapshot.dat");

	int i,j,k,ijk;
	FILE *file0;
	file0=fopen(actual_file,"w");
	fprintf(file0,"variables = x, y, z, rho, u, v, w, p, e, gradRho \n");

	if(MESHDIMENSIONS==2)
	{
		fprintf(file0,"zone t='1', i= %d, j= %d, f=point \n",pnt_config->int_iMeshPoints,pnt_config->int_jMeshPoints);
	}
	if(MESHDIMENSIONS==3)
	{
		fprintf(file0,"zone t='1', i= %d, j= %d, k= %d, f=point \n",pnt_config->int_iMeshPoints,pnt_config->int_jMeshPoints,pnt_config->int_kMeshPoints);
	}

	for (k=pnt_config->int_kStartReal_original; k <= pnt_config->int_kEndReal_original; k++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
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


void CGNS_SnapshotExportParallel(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep)
{
	int buffer=pnt_config->int_iMeshPoints*pnt_config->int_jMeshPoints*pnt_config->int_kMeshPoints;
	int i,j,k,ijk,ijk2,i2,j2,k2;

	int index_file,icelldim,iphysdim,index_base;
	int index_zone_out,index_coord,index_field,index_flow;

	char basename[33],solname[33];

	char foldername[200];
	sprintf(foldername,"%sSnapshot/",pnt_config->chr_folder);
	if ((pnt_config->MPI_rank==0)&&(pnt_config->int_initializeType!=1))
	{
//		rmdir(foldername);
		mkdir(foldername,0755);
	}
	MPI_Barrier(pnt_config->MPI_comm);

	float* x;
	float* y;
	float* z;
	float* u;
	float* v;
	float* w;
	float* p;
	float* rho;
	float* gradrho;
	float* mach;
	float* Lambda2;

//Speicherallokierung für die dynamischen-1D Arrays, worin die 3D-Lösungen gespeichert werden
	x = (float*) calloc(buffer,sizeof(float));
	y = (float*) calloc(buffer,sizeof(float));
	z = (float*) calloc(buffer,sizeof(float));
	u = (float*) calloc(buffer,sizeof(float));
	v = (float*) calloc(buffer,sizeof(float));
	w = (float*) calloc(buffer,sizeof(float));
	p = (float*) calloc(buffer,sizeof(float));
	rho = (float*) calloc(buffer,sizeof(float));
	gradrho = (float*) calloc(buffer,sizeof(float));
	mach = (float*) calloc(buffer,sizeof(float));
	Lambda2 = (float*) calloc(buffer,sizeof(float));

	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal_original; k <= pnt_config->int_kEndReal_original; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//				ijk=k*pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells+j*pnt_config->int_iMeshPointsGhostCells+i;

				i2=i-pnt_config->int_iStartReal;
				j2=j-pnt_config->int_jStartReal;
				k2=k-pnt_config->int_kStartReal_original;

//				ijk2=i2*int_jPoints*int_kPoints+j2*int_kPoints+k2;
				ijk2=k2*pnt_config->int_jMeshPoints*pnt_config->int_iMeshPoints+j2*pnt_config->int_iMeshPoints+i2;


				x[ijk2] = pnt_mesh->x[ijk];
				y[ijk2] = pnt_mesh->y[ijk];
				z[ijk2] = pnt_mesh->z[ijk];
				u[ijk2] = pnt_U_lastStep->u[ijk];
				v[ijk2] = pnt_U_lastStep->v[ijk];
				w[ijk2] = pnt_U_lastStep->w[ijk];
				p[ijk2] = pnt_U_lastStep->p[ijk];
				rho[ijk2] = pnt_U_lastStep->rho[ijk];
				gradrho[ijk2] = pnt_U_lastStep->gradRho[ijk];
				mach[ijk2] = pnt_U_lastStep->MachNumber[ijk];
				Lambda2[ijk2] = pnt_U_lastStep->Lambda2[ijk];
			}
		}
	}


	icelldim=pnt_config->int_meshDimensions;
	   cgsize_t isize[3][icelldim];

	/* WRITE X, Y, Z GRID POINTS TO CGNS FILE */
	/* open CGNS file for write */

	/* create base (user can give any name) */
	   sprintf(basename,"Base");

	   iphysdim=3;

	   if(pnt_config->int_initializeType==1)
	   {
		   cg_open(pnt_config->chr_SnapshotPath,CG_MODE_MODIFY,&index_file);
	   }
	   else
	   {
		   cg_open(pnt_config->chr_SnapshotPath,CG_MODE_WRITE,&index_file);
	   }

	   if(pnt_config->int_initializeType!=1)
	   {
		   cg_base_write(index_file,basename,icelldim,iphysdim,&index_base);
	   }
	   else
	   {
		   index_base=1;
	   }

	/* vertex size */
	   isize[0][0]=pnt_config->int_iMeshPoints;
	   isize[0][1]=pnt_config->int_jMeshPoints;

	/* cell size */
	   isize[1][0]=isize[0][0]-1;
	   isize[1][1]=isize[0][1]-1;

	/* boundary vertex size (always zero for structured grids) */
	   isize[2][0]=0;
	   isize[2][1]=0;

	   if (pnt_config->int_meshDimensions==3)
	   {
		   isize[0][2]=pnt_config->int_kMeshPoints;
		   isize[1][2]=isize[0][2]-1;
		   isize[2][2]=0;
	   }

	   if(pnt_config->int_initializeType!=1)
	   {
		/* create zone */
		   cg_zone_write(index_file,index_base,pnt_config->Zonename,*isize,Structured,&index_zone_out);
		/* write grid coordinates (user must use SIDS-standard names here) */
		   cg_coord_write(index_file,index_base,index_zone_out,RealSingle,"CoordinateX",
			   x,&index_coord);
		   cg_coord_write(index_file,index_base,index_zone_out,RealSingle,"CoordinateY",
			   y,&index_coord);
		   cg_coord_write(index_file,index_base,index_zone_out,RealSingle,"CoordinateZ",
			   z,&index_coord);
	   }
	   else
	   {
		   index_zone_out=1;
	   }


   /* define flow solution node name (user can give any name) */
	   sprintf(solname,"FlowSolution");
   /* create flow solution node (NOTE USE OF CellCenter HERE) */
	   cg_sol_write(index_file,index_base,index_zone_out,solname,Vertex,&index_flow);
   /* go to position within tree at FlowSolution_t node */
	   cg_goto(index_file,index_base,"Zone_t",index_zone_out,"FlowSolution_t",index_flow,"end");

   /* write flow solution (user must use SIDS-standard names here) */
	   cg_field_write(index_file,index_base,index_zone_out,index_flow,
		   RealSingle,"VelocityX",u,&index_field);
	   cg_field_write(index_file,index_base,index_zone_out,index_flow,
		   RealSingle,"VelocityY",v,&index_field);
	   cg_field_write(index_file,index_base,index_zone_out,index_flow,
		   RealSingle,"VelocityZ",w,&index_field);
	   cg_field_write(index_file,index_base,index_zone_out,index_flow,
		   RealSingle,"Density",rho,&index_field);
	   cg_field_write(index_file,index_base,index_zone_out,index_flow,
		   RealSingle,"Pressure",p,&index_field);
	   cg_field_write(index_file,index_base,index_zone_out,index_flow,
		   RealSingle,"DensityGradient",gradrho,&index_field);
	   cg_field_write(index_file,index_base,index_zone_out,index_flow,
		   RealSingle,"MachNumber",mach,&index_field);
	   cg_field_write(index_file,index_base,index_zone_out,index_flow,
		   RealSingle,"Lambda2",Lambda2,&index_field);

		char iteration[33];
		char name[33];
		sprintf(name,"Iterations");
		sprintf(iteration,"%d",(pnt_config->int_actualIteration-1));
		if (cg_goto(index_file,index_base,"end")) cg_error_exit();
		if (cg_descriptor_write(name, iteration)) cg_error_exit();

		char actualTime[33];
		sprintf(name,"ActualTime");
		sprintf(actualTime,"%g",(pnt_config->dbl_time_dim));
		if (cg_goto(index_file,index_base,"end")) cg_error_exit();
		if (cg_descriptor_write(name, actualTime)) cg_error_exit();

	   if(pnt_config->int_initializeType!=1)
	   {
		   //	   SDRC-TAG-COPY
			int index_file_in;
			int index_base_in;
			int index_zone_in;
			int ndescriptors;
			int D;
			char *text;
			cg_open(pnt_config->chr_MeshPath,CG_MODE_READ,&index_file_in);
			index_base_in=1;
			if((pnt_config->flag_SplitMeshFile==2)||((pnt_config->int_initializeType==1)))
			{
				index_zone_in=1;
			}
			else
			{
				index_zone_in=pnt_config->MPI_rank+1;
			}
			cg_goto(index_file_in,index_base_in,"Zone_t",index_zone_in,"end");
			cg_ndescriptors(&ndescriptors);
			for(D=1;D<=ndescriptors;D++)
			{
			  cg_descriptor_read(D, name, &text);
			  if (strcmp(name,"SDDC_TAG")==0)
			  {
				  cg_goto(index_file,index_base,"Zone_t",index_zone_out,"end");
				  cg_descriptor_write(name, text);
			  }
			}
			cg_free(text);

			//		READ & WRITE INTERFACES
			int NumberInterfaces,interface_in,interface_out;
			int *TransformMatrixOfInterface;
			cgsize_t *RangeOfInterface;
			cgsize_t *DonorRangeOfInterface;
			char donorname[33];
			char connectname[33];
			float RotationCenter[pnt_config->int_meshDimensions];
			float RotationAngle[pnt_config->int_meshDimensions];
			float Translation[pnt_config->int_meshDimensions];

			TransformMatrixOfInterface= (int *)calloc(pnt_config->int_meshDimensions, sizeof(int));
			RangeOfInterface= (cgsize_t *)calloc(2*pnt_config->int_meshDimensions, sizeof(cgsize_t));
			DonorRangeOfInterface= (cgsize_t *)calloc(2*pnt_config->int_meshDimensions, sizeof(cgsize_t));
			if (cg_n1to1(index_file_in, index_base_in,index_zone_in, &NumberInterfaces)) cg_error_exit();
			for(interface_in=1;interface_in<=NumberInterfaces;interface_in++)
			{
				if (cg_1to1_read(index_file_in, index_base_in,index_zone_in,interface_in, connectname, donorname,
						RangeOfInterface, DonorRangeOfInterface, TransformMatrixOfInterface)) cg_error_exit();
				if (cg_1to1_write(index_file,index_base,index_zone_out,connectname,donorname,
						RangeOfInterface,DonorRangeOfInterface,TransformMatrixOfInterface,&interface_out)) cg_error_exit();

				if(cg_1to1_periodic_read(index_file_in, index_base_in,index_zone_in,interface_in,
						RotationCenter, RotationAngle, Translation)!=CG_NODE_NOT_FOUND)
				{
					if(cg_1to1_periodic_write(index_file, index_base,index_zone_out,interface_out,
							RotationCenter, RotationAngle, Translation)) cg_error_exit();
				}
			}


	   //		READ & WRITE BC

			int nbocos,ib;
			int normalindex[3];
			int ndataset;
			int normallist;
			char boconame[33];
			char famname[33];

			BCType_t ibocotype;
			PointSetType_t iptset;
			DataType_t normaldatatype;
			cgsize_t ipnts[2*pnt_config->int_meshDimensions];
			cgsize_t npts,normallistflag;


			int nfamilies;
			char dummy[33];
			int nFamBC;
			int nGeo;
			int BC;
			int Fam;
			char unspecified[30];


			sprintf(unspecified,"Unspecified");


			cg_nfamilies(index_file_in,index_base_in,&nfamilies);

			BCType_t BCType;
			for(i=1;i<=nfamilies;i++)
			{
				cg_family_read(index_file_in,index_base_in, i, dummy,&nFamBC, &nGeo);
				if (strcmp(dummy,unspecified)!=0)
				{
					cg_family_write(index_file,index_base, dummy,&Fam);

					cg_fambc_read(index_file_in,index_base_in, i, 1,dummy, &BCType);
					cg_fambc_write(index_file,index_base, Fam, dummy, BCType,&BC);
				}
			}

		/* find out number of BCs that exist under this zone */
			cg_nbocos(index_file_in,index_base_in,index_zone_in,&nbocos);
		/* do loop over the total number of BCs */
			for (ib=1; ib <= nbocos; ib++)
			{
			/* get BC info */
				  cg_boco_info(index_file_in,index_base_in,index_zone_in,ib,boconame,&ibocotype,
							   &iptset,&npts,normalindex,&normallistflag,&normaldatatype,&ndataset);
				  cg_goto(index_file_in,index_base_in,"Zone_t",index_zone_in,"ZoneBC_t",1,"BC_t",ib,"end");
				  cg_famname_read(famname);
				  cg_boco_read(index_file_in,index_base_in,index_zone_in,ib,ipnts,&normallist);

	   //			  cg_fambc_write(index_file,index_base, 1, "FamilyName", BCType,&BC);
				  cg_boco_write(index_file,index_base, index_zone_out, boconame, ibocotype,
						  iptset,npts, ipnts, &BC);
				  cg_goto(index_file,index_base,"Zone_t",index_zone_out,"ZoneBC_t",1,"BC_t",BC,"end");
				  cg_famname_write(famname);
	   //			  cg_fambc_write(index_file_in, index_base_in, Fam, FamBCName,
	   //			        BCType, BC);
			}

			cg_close(index_file_in);
	   }

	   cg_close(index_file);

    free(x);
	free(y);
	free(z);
	free(u);
	free(v);
	free(w);
	free(p);
	free(rho);
	free(gradrho);
	free(mach);
	free(Lambda2);
}


void CGNS_FilmExportParallel(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_Film * pnt_Film)
{
	int index_file,icelldim,iphysdim,index_base;
	int index_zone,index_coord,index_field,index_flow;


	char basename[33];

	char foldername[200];
	sprintf(foldername,"%sFilm/",pnt_config->chr_folder);
	if (pnt_config->MPI_rank==0)
	{
//		rmdir(foldername);
		mkdir(foldername,0755);
	}
	MPI_Barrier(pnt_config->MPI_comm);


	int buffer=pnt_config->int_iMeshPoints*pnt_config->int_jMeshPoints*pnt_config->int_kMeshPoints;
	int i,j,k,ijkFilm,ijk,ijk2,i2,j2,k2,t,nsteps;

	cgsize_t idata[2],nuse;


	//    timedepended CGNS-DATA
	float time[pnt_config->int_Samples];
	char sn[pnt_config->int_Samples][33];
	char solname[pnt_config->int_Samples*32+1];  /* need an extra byte for the terminating 0 */

	strcpy(solname,"");

	float* x;
	float* y;
	float* z;
	float* u;
	float* v;
	float* w;
	float* p;
	float* rho;
	float* mach;
	float* Lambda2;
	float* gradrho;

//Speicherallokierung für die dynamischen-1D Arrays, worin die 3D-Lösungen gespeichert werden
	x = (float*) calloc(buffer,sizeof(float));
	y = (float*) calloc(buffer,sizeof(float));
	z = (float*) calloc(buffer,sizeof(float));
	u = (float*) calloc(buffer,sizeof(float));
	v = (float*) calloc(buffer,sizeof(float));
	w = (float*) calloc(buffer,sizeof(float));
	p = (float*) calloc(buffer,sizeof(float));
	rho = (float*) calloc(buffer,sizeof(float));
	gradrho = (float*) calloc(buffer,sizeof(float));
	Lambda2 = (float*) calloc(buffer,sizeof(float));
	mach = (float*) calloc(buffer,sizeof(float));

	for(t=0;t<pnt_config->int_Samples;t++)
	{
		time[t]=pnt_config->start_Time+(float)(t+1)*pnt_config->int_IterationsBetweenSamples*
				(pnt_config->dbl_L0_dim*pnt_config->dbl_numericalTau/pnt_config->dbl_u0_dim);
		sprintf(sn[t],"FlowSolution%d",t);
		sprintf(solname,"%s%-32s",solname,sn[t]);
	}

	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal_original; k <= pnt_config->int_kEndReal_original; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//				ijk=k*pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells+j*pnt_config->int_iMeshPointsGhostCells+i;

				i2=i-pnt_config->int_iStartReal;
				j2=j-pnt_config->int_jStartReal;
				k2=k-pnt_config->int_kStartReal_original;

//				ijk2=i2*int_jPoints*int_kPoints+j2*int_kPoints+k2;
				ijk2=k2*pnt_config->int_jMeshPoints*pnt_config->int_iMeshPoints+j2*pnt_config->int_iMeshPoints+i2;


				x[ijk2] = pnt_mesh->x[ijk];
				y[ijk2] = pnt_mesh->y[ijk];
				z[ijk2] = pnt_mesh->z[ijk];
			}
		}
	}


	icelldim=pnt_config->int_meshDimensions;
	   cgsize_t isize[3][icelldim];

	/* WRITE X, Y, Z GRID POINTS TO CGNS FILE */
	/* open CGNS file for write */

	/* create base (user can give any name) */
	   sprintf(basename,"Base");

	   iphysdim=3;
	   cg_open(pnt_config->chr_FilmPath,CG_MODE_WRITE,&index_file);

	   cg_base_write(index_file,basename,icelldim,iphysdim,&index_base);

	   cg_simulation_type_write(index_file,index_base,TimeAccurate);

		char iteration[33];
		char name[33];
		sprintf(name,"Iterations");
		sprintf(iteration,"%d",(pnt_config->int_actualIteration-1));
		if (cg_goto(index_file,index_base,"end")) cg_error_exit();
		if (cg_descriptor_write(name, iteration)) cg_error_exit();

	/* vertex size */

	   isize[0][0]=pnt_config->int_iMeshPoints;
	   isize[0][1]=pnt_config->int_jMeshPoints;

	/* cell size */
	   isize[1][0]=isize[0][0]-1;
	   isize[1][1]=isize[0][1]-1;

	/* boundary vertex size (always zero for structured grids) */
	   isize[2][0]=0;
	   isize[2][1]=0;

	   if (pnt_config->int_meshDimensions==3)
	   {
		   isize[0][2]=pnt_config->int_kMeshPoints;
		   isize[1][2]=isize[0][2]-1;
		   isize[2][2]=0;
	   }
	/* create zone */
	   cg_zone_write(index_file,index_base,pnt_config->Zonename,*isize,Structured,&index_zone);
	/* write grid coordinates (user must use SIDS-standard names here) */
	   cg_coord_write(index_file,index_base,index_zone,RealSingle,"CoordinateX",
	       x,&index_coord);
	   cg_coord_write(index_file,index_base,index_zone,RealSingle,"CoordinateY",
	       y,&index_coord);
	   cg_coord_write(index_file,index_base,index_zone,RealSingle,"CoordinateZ",
	       z,&index_coord);

	for(t=0;t<pnt_config->int_Samples;t++)
	{


		for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
		{
			for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
			{
				for (k=pnt_config->int_kStartReal_original; k <= pnt_config->int_kEndReal_original; k++)
				{
					i2=i-pnt_config->int_iStartReal;
					j2=j-pnt_config->int_jStartReal;
					k2=k-pnt_config->int_kStartReal_original;
					ijk2=k2*pnt_config->int_jMeshPoints*pnt_config->int_iMeshPoints+j2*pnt_config->int_iMeshPoints+i2;

					ijkFilm=(i-pnt_config->int_iStartReal)*pnt_config->int_jMeshPoints*pnt_config->int_kMeshPoints+(j-pnt_config->int_jStartReal)*pnt_config->int_kMeshPoints+(k-pnt_config->int_kStartReal_original)+
							t*pnt_config->int_iMeshPoints*pnt_config->int_jMeshPoints*pnt_config->int_kMeshPoints;


					u[ijk2] = pnt_Film->u[ijkFilm];
					v[ijk2] = pnt_Film->v[ijkFilm];
					w[ijk2] = pnt_Film->w[ijkFilm];
					p[ijk2] = pnt_Film->p[ijkFilm];
					rho[ijk2] = pnt_Film->rho[ijkFilm];
					gradrho[ijk2] = pnt_Film->gradRho[ijkFilm];
					mach[ijk2] = pnt_Film->MachNumber[ijkFilm];
					Lambda2[ijk2] = pnt_Film->Lambda2[ijkFilm];
				}
			}
		}

		cg_sol_write(index_file,index_base,index_zone,sn[t],Vertex,&index_flow);

	   /* write flow solution (user must use SIDS-standard names here) */
		   cg_field_write(index_file,index_base,index_zone,index_flow,
			   RealSingle,"VelocityX",u,&index_field);
		   cg_field_write(index_file,index_base,index_zone,index_flow,
			   RealSingle,"VelocityY",v,&index_field);
		   cg_field_write(index_file,index_base,index_zone,index_flow,
			   RealSingle,"VelocityZ",w,&index_field);
		   cg_field_write(index_file,index_base,index_zone,index_flow,
			   RealSingle,"Density",rho,&index_field);
		   cg_field_write(index_file,index_base,index_zone,index_flow,
			   RealSingle,"Pressure",p,&index_field);
		   cg_field_write(index_file,index_base,index_zone,index_flow,
			   RealSingle,"DensityGradient",gradrho,&index_field);
		   cg_field_write(index_file,index_base,index_zone,index_flow,
			   RealSingle,"MachNumber",mach,&index_field);
		   cg_field_write(index_file,index_base,index_zone,index_flow,
			   RealSingle,"Lambda2",Lambda2,&index_field);
	}

		/* create BaseIterativeData */
		   nsteps=pnt_config->int_Samples;
		   cg_biter_write(index_file,index_base,"TimeIterValues",nsteps);
	   /* go to BaseIterativeData level and write time values */
		   cg_goto(index_file,index_base,"BaseIterativeData_t",1,"end");
		   nuse=pnt_config->int_Samples;
		   cg_array_write("TimeValues",RealSingle,1,&nuse,&time);
	   /* create ZoneIterativeData */
		   cg_ziter_write(index_file,index_base,index_zone,"ZoneIterativeData");
	   /* go to ZoneIterativeData level and give info telling which */
	   /* flow solution corresponds with which time (solname(1) corresponds */
	   /* with time(1), solname(2) with time(2), and solname(3) with time(3)) */
		   cg_goto(index_file,index_base,"Zone_t",index_zone,"ZoneIterativeData_t",1,"end");
		   idata[0]=32;
		   idata[1]=pnt_config->int_Samples;
		   cg_array_write("FlowSolutionPointers",Character,2,idata,solname);




//	   SDRC-TAG-COPY
	   int index_file_in;
	   int index_base_in;
	   int index_zone_in;
	   int ndescriptors;
	   int D;
	   char *text;
	   cg_open(pnt_config->chr_MeshPath,CG_MODE_READ,&index_file_in);
	   index_base_in=1;
		if((pnt_config->flag_SplitMeshFile==2)||((pnt_config->int_initializeType==1)))
		{
			index_zone_in=1;
		}
		else
		{
			index_zone_in=pnt_config->MPI_rank+1;
		}

	  cg_goto(index_file_in,index_base_in,"Zone_t",index_zone_in,"end");
	  cg_ndescriptors(&ndescriptors);
	  for(D=1;D<=ndescriptors;D++)
	  {
		  cg_descriptor_read(D, name, &text);
		  if (strcmp(name,"SDDC_TAG")==0)
		  {
			  cg_goto(index_file,index_base,"Zone_t",index_zone,"end");
			  cg_descriptor_write(name, text);
		  }
	  }

		//		READ & WRITE INTERFACES
		int NumberInterfaces,interface_in,interface_out;
		int *TransformMatrixOfInterface;
		cgsize_t *RangeOfInterface;
		cgsize_t *DonorRangeOfInterface;
		char donorname[33];
		char connectname[33];
		float RotationCenter[pnt_config->int_meshDimensions];
		float RotationAngle[pnt_config->int_meshDimensions];
		float Translation[pnt_config->int_meshDimensions];

		TransformMatrixOfInterface= (int *)calloc(pnt_config->int_meshDimensions, sizeof(int));
		RangeOfInterface= (cgsize_t *)calloc(2*pnt_config->int_meshDimensions, sizeof(cgsize_t));
		DonorRangeOfInterface= (cgsize_t *)calloc(2*pnt_config->int_meshDimensions, sizeof(cgsize_t));
		if (cg_n1to1(index_file_in, index_base_in,index_zone_in, &NumberInterfaces)) cg_error_exit();
		for(interface_in=1;interface_in<=NumberInterfaces;interface_in++)
		{
			if (cg_1to1_read(index_file_in, index_base_in,index_zone_in,interface_in, connectname, donorname,
					RangeOfInterface, DonorRangeOfInterface, TransformMatrixOfInterface)) cg_error_exit();
			if (cg_1to1_write(index_file,index_base,index_zone,connectname,donorname,
					RangeOfInterface,DonorRangeOfInterface,TransformMatrixOfInterface,&interface_out)) cg_error_exit();

			if(cg_1to1_periodic_read(index_file_in, index_base_in,index_zone_in,interface_in,
					RotationCenter, RotationAngle, Translation)!=CG_NODE_NOT_FOUND)
			{
				if(cg_1to1_periodic_write(index_file, index_base,index_zone,interface_out,
						RotationCenter, RotationAngle, Translation)) cg_error_exit();
			}
		}
	  cg_close(index_file_in);
	  cg_close(index_file);



	free(x);
	free(y);
	free(z);
	free(u);
	free(v);
	free(w);
	free(p);
	free(rho);
	free(gradrho);
	free(Lambda2);
	free(mach);


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
	float time[pnt_config->int_TotalIterations];
	char sn[pnt_config->int_TotalIterations][33];
	char solname[pnt_config->int_TotalIterations*32+1];  /* need an extra byte for the terminating 0 */

	strcpy(solname,"");

	float* x;
	float* y;
	float* z;
	float* p;

//Speicherallokierung für die dynamischen-1D Arrays, worin die 3D-Lösungen gespeichert werden
	x = (float*) calloc(buffer,sizeof(float));
	y = (float*) calloc(buffer,sizeof(float));
	z = (float*) calloc(buffer,sizeof(float));
	p = (float*) calloc(buffer,sizeof(float));

	for(t=0;t<pnt_config->int_TotalIterations;t++)
	{
		time[t]=pnt_config->start_Time+(float)(t+1)*pnt_config->int_IterationsBetweenSamples*
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
