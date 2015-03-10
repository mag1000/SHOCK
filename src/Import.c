#include "Import.h"
#include "SHOCK.h"
#include "Functions.h"
#include "BC.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cgnslib.h"
#include <string.h>
#include "mpi.h"
#include <unistd.h>
#include "iniparser.h"




void ConfigImport(
		struct strct_configuration * pnt_config)
{
	int i,p;
	char pressureHistory_actualName[50];
	char velocityHistory_actualName[50];
	//Konfigurationsdatei wird importiert. Hier befinden sich sämtliche Informationen, die der Solver benötigt.
	//Die Datei "config.sws", die importiert wird, enthält Beschreibungen zu den einzelnen Größen.
	//Es handelt sich um Einstellungen bezüglich des Netzes, der MPI-Verteilung, der Physik und der Numerik.

    dictionary  *   ini ;

    ini = iniparser_load(pnt_config->chr_configPath);
    if (ini==NULL) {
        fprintf(stderr, "cannot parse config file: %s\n", pnt_config->chr_configPath);
        MPI_Finalize();
    }

//	General
	pnt_config->int_initializeType = iniparser_getint(ini, "general:initializationtype", 1);
	sprintf(pnt_config->chr_MeshPath, "%s", iniparser_getstring(ini, "general:meshpath", NULL));
	pnt_config->int_TotalIterations= iniparser_getint(ini, "general:iterations", -1);
	pnt_config->int_EndIteration=pnt_config->int_TotalIterations;
	pnt_config->int_Samples = iniparser_getint(ini, "general:samples", -1);
	pnt_config->int_StartSampling = iniparser_getint(ini, "general:StartSampling", 0);
	pnt_config->int_NumberBackups = iniparser_getint(ini, "general:Backup-Exports", 0);
	pnt_config->flag_TauAccelerator = iniparser_getboolean(ini, "general:TauAccelerator", 0);
	pnt_config->flag_swapDivisionFile = iniparser_getint(ini, "general:SwapDivFile", 0);	

//	Numerics
	pnt_config->int_SpaceOrder=SPACEORDER;
	pnt_config->int_TimeOrder = iniparser_getint(ini, "numerics:TimeOrder", 4);
	pnt_config->numericalTau = iniparser_getdouble(ini, "numerics:tau", -1);
	pnt_config->numericalTauStart=pnt_config->numericalTau;

//	Tau
	pnt_config->TauAccelerator_factor = iniparser_getdouble(ini, "tau:tauaccelerator_factor", 1.1);
	pnt_config->TauDecelerator_factor = iniparser_getdouble(ini, "tau:taudecelerator_factor", 0.7);
	pnt_config->int_NumberResets = iniparser_getint(ini, "tau:NumberResets", 10);
	pnt_config->int_distanceNAN = iniparser_getdouble(ini, "tau:DistanceNAN", 5000);
	pnt_config->int_distanceForward = iniparser_getdouble(ini, "tau:DistanceForward", 1000);
	pnt_config->int_distanceBackward = iniparser_getdouble(ini, "tau:DistanceBackward", 200);
	pnt_config->int_distanceForwardStart=pnt_config->int_distanceForward;

//	BoundaryConditions
	pnt_config->flag_BC_option_inflow_normal_sub=iniparser_getboolean(ini, "boundaryconditions:Inflow-Normal-Subsonic", -1);
	pnt_config->flag_BC_option_inflow_riemann_sub=iniparser_getboolean(ini, "boundaryconditions:Inflow-Riemann-Subsonic", -1);
	pnt_config->flag_BC_option_inflow_isentrop_sub=iniparser_getboolean(ini, "boundaryconditions:Inflow-Isentrop-Subsonic", -1);
	pnt_config->flag_BC_option_inflow_normal_super=iniparser_getboolean(ini, "boundaryconditions:Inflow-Normal-Supersonic", -1);

	if((pnt_config->flag_BC_option_inflow_normal_sub+
			pnt_config->flag_BC_option_inflow_riemann_sub+
			pnt_config->flag_BC_option_inflow_isentrop_sub+
			pnt_config->flag_BC_option_inflow_normal_super)!=1)
	{
		if(pnt_config->MPI_rank==0){printf("SHOCK: ERROR: Fehlerhafte Eingabe fuer Inflow-Randbedingung\n");}
		MPI_Barrier(pnt_config->MPI_comm);
		MPI_Abort(pnt_config->MPI_comm,13371);
	}

	pnt_config->flag_BC_option_outflow_normal_sub=iniparser_getboolean(ini, "boundaryconditions:Outflow-Normal-Subsonic", -1);
	pnt_config->flag_BC_option_outflow_riemann_sub=iniparser_getboolean(ini, "boundaryconditions:Outflow-Riemann-Subsonic", -1);
	pnt_config->flag_BC_option_outflow_rudy_sub=iniparser_getboolean(ini, "boundaryconditions:Outflow-Rudy-Subsonic", -1);

	if((pnt_config->flag_BC_option_outflow_normal_sub+
			pnt_config->flag_BC_option_outflow_riemann_sub+
			pnt_config->flag_BC_option_outflow_rudy_sub)!=1)
	{
		if(pnt_config->MPI_rank==0){printf("SHOCK: ERROR: Fehlerhafte Eingabe fuer Outflow-Randbedingung\n");}
		MPI_Barrier(pnt_config->MPI_comm);
		MPI_Abort(pnt_config->MPI_comm,13372);
	}

	pnt_config->rho_inflow = iniparser_getdouble(ini, "boundaryconditions:rho_inflow", 1);
	pnt_config->p_inflow = iniparser_getdouble(ini, "boundaryconditions:p_inflow", 1);
	pnt_config->AoA = iniparser_getdouble(ini, "boundaryconditions:AoA", 0);

	pnt_config->p_out = iniparser_getdouble(ini, "boundaryconditions:p_outflow", 1);

	pnt_config->T_wall = iniparser_getdouble(ini, "boundaryconditions:T_wall", 1);

//	Fluid Properties
	pnt_config->machNumber = iniparser_getdouble(ini, "fluidproperties:mach", -1);
	pnt_config->reynoldsNumber = iniparser_getdouble(ini, "fluidproperties:reynolds", -1);
	pnt_config->prandtlNumber = iniparser_getdouble(ini, "fluidproperties:Prandtl", -1);
	pnt_config->gammaNumber = iniparser_getdouble(ini, "fluidproperties:Gamma", -1);
	pnt_config->gasConstantNumber = iniparser_getdouble(ini, "fluidproperties:R", -1);
	pnt_config->T0_dim = iniparser_getdouble(ini, "fluidproperties:T", -1);
	pnt_config->L0_dim = iniparser_getdouble(ini, "fluidproperties:L", -1);

//	Export
	pnt_config->flag_exportMetric = iniparser_getboolean(ini, "export:metric", -1);
	pnt_config->flag_ReducedExport = iniparser_getboolean(ini, "export:reduced", 0);

//	Options
	pnt_config->int_specialInitializeType = iniparser_getint(ini, "options:SpecialInitializeType", 0);
	pnt_config->flag_constantZValues = iniparser_getboolean(ini, "options:3Dto2D", 0);
	pnt_config->flag_IBC = iniparser_getboolean(ini, "options:IBC", 0);
	pnt_config->flag_LaminarBoundary = iniparser_getboolean(ini, "options:LaminarBoundary", 0);
	pnt_config->flag_PressureHistory = iniparser_getboolean(ini, "options:PressureHistory", 0);
	pnt_config->flag_VelocityHistory = iniparser_getboolean(ini, "options:VelocityHistory", 0);
	pnt_config->flag_PressureWaves = iniparser_getboolean(ini, "options:PressureWaves", 0);
	pnt_config->flag_Inviscid = iniparser_getboolean(ini, "options:Inviscid", 0);
	pnt_config->flag_Vortex = iniparser_getboolean(ini, "options:Vortex", 0);
	pnt_config->flag_rotation_symmetric = iniparser_getboolean(ini, "options:2D-Rotation-Symmetric", 0);
	pnt_config->flag_ManufacturedSolution = iniparser_getboolean(ini, "options:ManufacturedSolution", 0);


//	LaminarBoundary
	pnt_config->LaminarBoundary_xStart = iniparser_getdouble(ini, "LaminarBoundary:x-Startposition", -1);

//	PressureWaves
	pnt_config->pw_UseBC = iniparser_getboolean(ini, "PressureWaves:UseBC", 0);
	pnt_config->pw_UseFlowAverage = iniparser_getboolean(ini, "PressureWaves:UseFlowAverage", 0);
	pnt_config->pw_numberSources = iniparser_getint(ini, "PressureWaves:NumberSources", 1);
	pnt_config->pw_x0 = iniparser_getdouble(ini, "PressureWaves:x0-Location", 0);
	pnt_config->pw_y0 = iniparser_getdouble(ini, "PressureWaves:y0-Location", 0);
	pnt_config->pw_z0 = iniparser_getdouble(ini, "PressureWaves:z0-Location", 0);
	pnt_config->pw_x1 = iniparser_getdouble(ini, "PressureWaves:x1-Location", 0);
	pnt_config->pw_y1 = iniparser_getdouble(ini, "PressureWaves:y1-Location", 0);
	pnt_config->pw_z1 = iniparser_getdouble(ini, "PressureWaves:z1-Location", 0);
	pnt_config->pw_r0 = iniparser_getdouble(ini, "PressureWaves:Radius", 1);
	pnt_config->pw_amplitude = iniparser_getdouble(ini, "PressureWaves:Amplitude", 1);
	pnt_config->pw_frequency = iniparser_getdouble(ini, "PressureWaves:Frequency", 1);

//	Vortex
	pnt_config->Vortex_x_wirb_zentr = iniparser_getdouble(ini, "Vortex:x-Location", -1);
	pnt_config->Vortex_y_wirb_zentr = iniparser_getdouble(ini, "Vortex:y-Location", -1);
	pnt_config->Vortex_r_wirb_max = iniparser_getdouble(ini, "Vortex:Radius", -1);
	pnt_config->Vortex_faktor_quer = iniparser_getdouble(ini, "Vortex:f", -1);
	pnt_config->Vortex_beta = iniparser_getdouble(ini, "Vortex:beta", -1);

//	PressureHistory
	if(pnt_config->flag_PressureHistory==1)
	{
		pnt_config->PressureHistory_No = iniparser_getint(ini, "PressureHistory:NumberLocations", 4);
		pnt_config->PressureHistory_x_P = (FLT *)calloc(pnt_config->PressureHistory_No, sizeof(FLT ));
		pnt_config->PressureHistory_y_P = (FLT *)calloc(pnt_config->PressureHistory_No, sizeof(FLT ));
		pnt_config->PressureHistory_z_P = (FLT *)calloc(pnt_config->PressureHistory_No, sizeof(FLT ));
		pnt_config->PressureHistory_x_P_real = (FLT *)calloc(pnt_config->PressureHistory_No, sizeof(FLT ));
		pnt_config->PressureHistory_y_P_real = (FLT *)calloc(pnt_config->PressureHistory_No, sizeof(FLT ));
		pnt_config->PressureHistory_z_P_real = (FLT *)calloc(pnt_config->PressureHistory_No, sizeof(FLT ));
		pnt_config->PressureHistory_time=(FLT *)calloc(pnt_config->int_TotalIterations, sizeof(FLT ));
		pnt_config->PressureHistory_pressure=(FLT **)calloc(pnt_config->PressureHistory_No, sizeof(FLT *));
		for(p=0;p<pnt_config->PressureHistory_No;p++)
		{
			pnt_config->PressureHistory_pressure[p]=(FLT *)calloc(pnt_config->int_TotalIterations, sizeof(FLT));
		}
		pnt_config->ijk_PressureHistory_P = (int *)calloc(pnt_config->PressureHistory_No, sizeof(int ));
		pnt_config->flag_PressureHistory_P = (int *)calloc(pnt_config->PressureHistory_No, sizeof(int ));
		for(i=0;i<pnt_config->PressureHistory_No;i++)
		{
			sprintf(pressureHistory_actualName,"PressureHistory:P%d_x-Location",i);
			pnt_config->PressureHistory_x_P[i] = iniparser_getdouble(ini, pressureHistory_actualName, 0);
			sprintf(pressureHistory_actualName,"PressureHistory:P%d_y-Location",i);
			pnt_config->PressureHistory_y_P[i] = iniparser_getdouble(ini, pressureHistory_actualName, 0);
			sprintf(pressureHistory_actualName,"PressureHistory:P%d_z-Location",i);
			pnt_config->PressureHistory_z_P[i] = iniparser_getdouble(ini, pressureHistory_actualName, 0);
		}
	}

//	VelocityHistory
	if(pnt_config->flag_VelocityHistory==1)
	{
		pnt_config->VelocityHistory_No = iniparser_getint(ini, "VelocityHistory:NumberLocations", 4);
		pnt_config->VelocityHistory_x_P = (FLT *)calloc(pnt_config->VelocityHistory_No, sizeof(FLT ));
		pnt_config->VelocityHistory_y_P = (FLT *)calloc(pnt_config->VelocityHistory_No, sizeof(FLT ));
		pnt_config->VelocityHistory_z_P = (FLT *)calloc(pnt_config->VelocityHistory_No, sizeof(FLT ));
		pnt_config->VelocityHistory_x_P_real = (FLT *)calloc(pnt_config->VelocityHistory_No, sizeof(FLT ));
		pnt_config->VelocityHistory_y_P_real = (FLT *)calloc(pnt_config->VelocityHistory_No, sizeof(FLT ));
		pnt_config->VelocityHistory_z_P_real = (FLT *)calloc(pnt_config->VelocityHistory_No, sizeof(FLT ));
		pnt_config->VelocityHistory_time=(FLT *)calloc(pnt_config->int_TotalIterations, sizeof(FLT ));
		pnt_config->VelocityHistory_VelocityX=(FLT **)calloc(pnt_config->VelocityHistory_No, sizeof(FLT *));
		pnt_config->VelocityHistory_VelocityY=(FLT **)calloc(pnt_config->VelocityHistory_No, sizeof(FLT *));
		pnt_config->VelocityHistory_VelocityZ=(FLT **)calloc(pnt_config->VelocityHistory_No, sizeof(FLT *));
		for(p=0;p<pnt_config->VelocityHistory_No;p++)
		{
			pnt_config->VelocityHistory_VelocityX[p]=(FLT *)calloc(pnt_config->int_TotalIterations, sizeof(FLT));
			pnt_config->VelocityHistory_VelocityY[p]=(FLT *)calloc(pnt_config->int_TotalIterations, sizeof(FLT));
			pnt_config->VelocityHistory_VelocityZ[p]=(FLT *)calloc(pnt_config->int_TotalIterations, sizeof(FLT));
		}
		pnt_config->ijk_VelocityHistory_P = (int *)calloc(pnt_config->VelocityHistory_No, sizeof(int ));
		pnt_config->flag_VelocityHistory_P = (int *)calloc(pnt_config->VelocityHistory_No, sizeof(int ));
		for(i=0;i<pnt_config->VelocityHistory_No;i++)
		{
			sprintf(velocityHistory_actualName,"VelocityHistory:V%d_x-Location",i);
			pnt_config->VelocityHistory_x_P[i] = iniparser_getdouble(ini, velocityHistory_actualName, 0);
			sprintf(velocityHistory_actualName,"VelocityHistory:V%d_y-Location",i);
			pnt_config->VelocityHistory_y_P[i] = iniparser_getdouble(ini, velocityHistory_actualName, 0);
			sprintf(velocityHistory_actualName,"VelocityHistory:V%d_z-Location",i);
			pnt_config->VelocityHistory_z_P[i] = iniparser_getdouble(ini, velocityHistory_actualName, 0);
		}
	}

//	IBC
	pnt_config->flag_IBC_Moving = iniparser_getboolean(ini, "IBC:MovingBC", 0);
	pnt_config->IBC_Type = iniparser_getint(ini, "IBC:Type", 0);
	pnt_config->IBC_StartpositionX = iniparser_getdouble(ini, "IBC:StartpositionX", 0);
	pnt_config->IBC_StartpositionY = iniparser_getdouble(ini, "IBC:StartpositionY", 0);
	pnt_config->IBC_StartpositionZ = iniparser_getdouble(ini, "IBC:StartpositionZ", 0);
	pnt_config->IBC_SizeX = iniparser_getdouble(ini, "IBC:SizeX", 0);
	pnt_config->IBC_SizeY = iniparser_getdouble(ini, "IBC:SizeY", 0);
	pnt_config->IBC_SizeZ = iniparser_getdouble(ini, "IBC:SizeZ", 0);
	pnt_config->IBC_MovingSpeed = iniparser_getdouble(ini, "IBC:Speed", 0);
	pnt_config->IBC_MovingType = iniparser_getint(ini, "IBC:Movement", 0);
	pnt_config->IBC_SpeedFactor = iniparser_getdouble(ini, "IBC:SpeedFactor", 1.0);
	pnt_config->IBC_MovingStepsize = iniparser_getint(ini, "IBC:Stepsize", 100);

	//InitializeValues
	pnt_config->InitializeValues_u0 = iniparser_getdouble(ini, "InitializeValues:u0", 1);
	pnt_config->InitializeValues_p0 = iniparser_getdouble(ini, "InitializeValues:p0", 1);
	pnt_config->InitializeValues_rho0 = iniparser_getdouble(ini, "InitializeValues:rho0", 1);
	pnt_config->InitializeValues_p1 = iniparser_getdouble(ini, "InitializeValues:p1", 1);
	pnt_config->InitializeValues_rho1 = iniparser_getdouble(ini, "InitializeValues:rho1", 1);
	pnt_config->InitializeValues_u1 = iniparser_getdouble(ini, "InitializeValues:u1", 1);
	pnt_config->InitializeValues_xBorder = iniparser_getdouble(ini, "InitializeValues:xBorder", 0.);	
	pnt_config->u_inflow = cos(pnt_config->AoA/360*2*MY_PI)*pnt_config->InitializeValues_u0;
	pnt_config->v_inflow = sin(pnt_config->AoA/360*2*MY_PI)*pnt_config->InitializeValues_u0;
	pnt_config->w_inflow = 0.0;

	//ManufacturedSolution
	pnt_config->ManufacturedSolution_case = iniparser_getint(ini, "ManufacturedSolution:case", 1);


    iniparser_freedict(ini);

	if((pnt_config->flag_rotation_symmetric==1)&&(pnt_config->int_meshDimensions==3))
	{
		if(pnt_config->MPI_rank==0){printf("SHOCK: ERROR: Rotationssymmetrisch und 3D schliesst sich aus\n");}
		MPI_Barrier(pnt_config->MPI_comm);
		MPI_Abort(pnt_config->MPI_comm,13373);
	}

	if(pnt_config->int_Samples>pnt_config->int_TotalIterations)
	{
		pnt_config->int_Samples=pnt_config->int_TotalIterations;
	}

	if(pnt_config->int_StartSampling>pnt_config->int_TotalIterations)
	{
		pnt_config->int_StartSampling=0;
	}

	if((pnt_config->int_TotalIterations==0)||(pnt_config->int_Samples==0))
	{
		pnt_config->int_Samples=1;
	}

//	Definiere für jeden Prozess die Filenamen
	DefineFilesPath(
			pnt_config);

}






void MeshImport_CGNS(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh)
{
	/*
	  dimension statements (note that tri-dimensional arrays
	  x,y,z must be dimensioned exactly as [N][17][21] (N>=9)
	  for this particular case or else they will be read from
	  the CGNS file incorrectly!  Other options are to use 1-D
	  arrays, use dynamic memory, or pass index values to a
	  subroutine and dimension exactly there):
	*/

		int index_file,index_base,index_zone;
		int celldim=pnt_config->int_meshDimensions;
		cgsize_t isize[3][celldim],irmin[celldim],irmax[celldim];
		char zonename[33];

	//Speicherallokierung für die dynamischen-1D Arrays, worin die 3D-Lösungen gespeichert werden
	    int buffer=pnt_config->int_iMeshPoints*pnt_config->int_jMeshPoints*pnt_config->int_kMeshPoints;
		int i,j,k,ijk,ijk2,nzone;
		int i2,j2,k2;

		FLT* x;
		FLT* y;
		FLT* z;
		x = (FLT*) calloc(buffer,sizeof(FLT));
		y = (FLT*) calloc(buffer,sizeof(FLT));
		z = (FLT*) calloc(buffer,sizeof(FLT));

	/* READ X, Y, Z GRID POINTS FROM CGNS FILE */
	/* open CGNS file for read-only */
	    cg_open(pnt_config->chr_MeshPath,CG_MODE_READ,&index_file);
	/* we know there is only one base (real working code would check!) */
	    index_base=1;
	    //Sofern es nur eine Mesh-File gibt, bestimmt der rank die zone-ID.
	    //Bei einer gesplitteten Mesh-File ist diese stets 1
	    if(pnt_config->int_initializeType==1)
	    {
	    	index_zone=1;
	    }
	    else
	    {
	    	index_zone=pnt_config->MPI_rank+1;

		    cg_nzones(index_file,index_base,&nzone);
		    if (nzone != pnt_config->MPI_size)
		    {
		      printf("SHOCK: Error. Es werden %d Zonen erwartet. %d sind vorhanden.\n",pnt_config->MPI_size,nzone);
		    }
	    }

	/* get zone size (and name - although not needed here) */
	    cg_zone_read(index_file,index_base,index_zone,zonename,isize[0]);


//	    printf("My rank is %d and my zonename is %s\n",pnt_config->MPI_rank, zonename);

	/* lower range index */
	    irmin[0]=1;
	    irmin[1]=1;
	/* upper range index of vertices */
	    irmax[0]=isize[0][0];
	    irmax[1]=isize[0][1];

	    if (pnt_config->int_meshDimensions==3)
	    {
		    irmin[2]=1;
	    	irmax[2]=isize[0][2];
	    }

	/* read grid coordinates */
	    cg_coord_read(index_file,index_base,index_zone,"CoordinateX",
	                  RealSingle,irmin,irmax,x);
	    cg_coord_read(index_file,index_base,index_zone,"CoordinateY",
	                  RealSingle,irmin,irmax,y);
	    cg_coord_read(index_file,index_base,index_zone,"CoordinateZ",
	                  RealSingle,irmin,irmax,z);

		for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
		{
			for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
			{
				for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
				{
					ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
	//				ijk=k*pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells+j*pnt_config->int_iMeshPointsGhostCells+i;

					i2=i-pnt_config->int_iStartReal;
					j2=j-pnt_config->int_jStartReal;
					k2=k-pnt_config->int_kStartReal;

					ijk2=k2*pnt_config->int_jMeshPoints*pnt_config->int_iMeshPoints+j2*pnt_config->int_iMeshPoints+i2;

					pnt_mesh->x[ijk]=x[ijk2];
					pnt_mesh->y[ijk]=y[ijk2];
					if(MESHDIMENSIONS==3)
					{
						pnt_mesh->z[ijk]=z[ijk2];
					}
					else
					{
						pnt_mesh->z[ijk]=0.0;
					}

				}
			}
		}
	    free(x);
	    free(y);
	    free(z);


	/* close CGNS file */
	    cg_close(index_file);
}


void BCImport_CGNS(
		struct strct_configuration * pnt_config)
{
	int index_file,index_base,index_zone,nbocos,ib;
	int normalindex[3],ndataset;
	int normallist;
	char boconame[33];
	char famname[33];
	int dim=pnt_config->int_meshDimensions;
	BCType_t ibocotype;
	PointSetType_t iptset;
	DataType_t normaldatatype;
	cgsize_t ipnts[2*dim];
	cgsize_t npts,normallistflag;

/* READ BOUNDARY CONDITIONS FROM EXISTING CGNS FILE */
/* open CGNS file for read-only */
	cg_open(pnt_config->chr_MeshPath,CG_MODE_READ,&index_file);
/* we know there is only one base (real working code would check!) */
	index_base=1;

    if(pnt_config->int_initializeType==1)
    {
    	index_zone=1;
    }
    else
    {
    	index_zone=pnt_config->MPI_rank+1;
    }


	int nfamilies,i;
	char dummy[30];
	cg_nfamilies(index_file,index_base,&nfamilies);
	char FamilyName[nfamilies][30];
	char BCName[nfamilies][30];
	int nFamBC;
	int nGeo;
	char BCFarfield[30];
	char BCInflow[30];
	char BCOutflow[30];
	char BCOutflowSubsonic[30];
	char BCWallInviscid[30];
	char BCWallViscous[30];
	char BCInflowSupersonic[30];
	char BCInflowSubsonic[30];

	int intern_i;
	intern_i=0;

	strcpy(BCFarfield,"BCFarfield");
	strcpy(BCInflow,"BCInflow");
	strcpy(BCOutflow,"BCOutflow");
	strcpy(BCOutflowSubsonic,"BCOutflowSubsonic");
	strcpy(BCWallInviscid,"BCWallInviscid");
	strcpy(BCWallViscous,"BCWallViscous");
	strcpy(BCInflowSupersonic,"BCInflowSupersonic");
	strcpy(BCInflowSubsonic,"BCInflowSubsonic");


	BCType_t BCType;

	if (pnt_config->MPI_rank==0){printf("SHOCK: Folgende Randbedingungen wurden gesetzt: ");}

	for(i=1;i<=nfamilies;i++)
	{
		cg_goto(index_file,index_base,"Family_t",i,"end");
		cg_famname_read(famname);
		cg_fambc_read(index_file,index_base, i, 1,dummy, &BCType);
		cg_family_read(index_file,index_base, i, dummy,&nFamBC, &nGeo);

		if((strcmp(dummy,"Unspecified")!=0)&&(strcmp(dummy,"UnspecifiedBC")!=0))
		{
			strcpy(FamilyName[intern_i],dummy);
			strcpy(BCName[intern_i],BCTypeName[BCType]);

			if (pnt_config->MPI_rank==0)
			{
				printf("%s ",BCName[intern_i]);
			}
			intern_i++;
		}
	}
	if (pnt_config->MPI_rank==0)
	{
		printf("\n");
	}

/* find out number of BCs that exist under this zone */
	cg_nbocos(index_file,index_base,index_zone,&nbocos);
/* do loop over the total number of BCs */
	for (ib=1; ib <= nbocos; ib++)
	{
	/* get BC info */
		  cg_boco_info(index_file,index_base,index_zone,ib,boconame,&ibocotype,
					   &iptset,&npts,normalindex,&normallistflag,&normaldatatype,&ndataset);
		  cg_goto(index_file,index_base,"Zone_t",index_zone,"ZoneBC_t",1,"BC_t",ib,"end");
		  cg_famname_read(famname);
//		  printf("famname: %s\n",famname);

		  i=-1;
		  do
		  {
			  i++;
		  }while(strcmp(famname,FamilyName[i])!=0);

		  cg_boco_read(index_file,index_base,index_zone,ib,ipnts,&normallist);

		  if (
				  (strcmp(BCName[i],BCFarfield)==0)||
				  (strcmp(BCName[i],BCInflow)==0)||
				  (strcmp(BCName[i],BCOutflow)==0)||
				  (strcmp(BCName[i],BCOutflowSubsonic)==0)||
				  (strcmp(BCName[i],BCWallInviscid)==0)||
				  (strcmp(BCName[i],BCWallViscous)==0)||
				  (strcmp(BCName[i],BCInflowSupersonic)==0)||
				  (strcmp(BCName[i],BCInflowSubsonic)==0)
				  )
		  {
			if((ipnts[0]==1)&&(ipnts[0+dim]==1)){strcpy(pnt_config->BC_Left,BCName[i]);}
			if((ipnts[0]!=1)&&(ipnts[0+dim]!=1)){strcpy(pnt_config->BC_Right,BCName[i]);}
			if((ipnts[1]==1)&&(ipnts[1+dim]==1)){strcpy(pnt_config->BC_Bottom,BCName[i]);}
			if((ipnts[1]!=1)&&(ipnts[1+dim]!=1)){strcpy(pnt_config->BC_Top,BCName[i]);}
			if(dim==3)
			{
				if((ipnts[2]==1)&&(ipnts[2+dim]==1)){strcpy(pnt_config->BC_Behind,BCName[i]);}
				if((ipnts[2]!=1)&&(ipnts[2+dim]!=1)){strcpy(pnt_config->BC_InFront,BCName[i]);}
			}
		  }

	}

//	if(pnt_config->MPI_rank==0)
//	{
//		printf("zone: %s | left: %s | right: %s | top: %s | bottom: %s | behind: %s | infront: %s | \n",
//				pnt_config->Zonename,
//				pnt_config->BC_Left,
//				pnt_config->BC_Right,
//				pnt_config->BC_Top,
//				pnt_config->BC_Bottom,
//				pnt_config->BC_Behind,
//				pnt_config->BC_InFront);
//	}


/* close CGNS file */
	cg_close(index_file);
}


