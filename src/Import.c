#include "Import.h"
#include "saiwenos.h"
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
	pnt_config->int_IOType = iniparser_getint(ini, "general:IOtype", 0);
	sprintf(pnt_config->chr_MeshPath, "%s", iniparser_getstring(ini, "general:meshpath", NULL));
	pnt_config->int_TotalIterations= iniparser_getint(ini, "general:iterations", -1);
	pnt_config->int_EndIteration=pnt_config->int_TotalIterations;
	pnt_config->int_Samples = iniparser_getint(ini, "general:samples", -1);
	pnt_config->int_StartSampling = iniparser_getint(ini, "general:StartSampling", 0);
	pnt_config->int_NumberBackups = iniparser_getint(ini, "general:Backup-Exports", 0);
	pnt_config->flag_TauAccelerator = iniparser_getboolean(ini, "general:TauAccelerator", 0);

//	Numerics
	pnt_config->int_SpaceOrder=SPACEORDER;
	pnt_config->int_TimeOrder = iniparser_getint(ini, "numerics:TimeOrder", 4);
	pnt_config->dbl_numericalTau = iniparser_getdouble(ini, "numerics:tau", -1);
	pnt_config->dbl_numericalTauStart=pnt_config->dbl_numericalTau;

//	Tau
	pnt_config->dbl_TauAccelerator_factor = iniparser_getdouble(ini, "tau:tauaccelerator_factor", 1.1);
	pnt_config->dbl_TauDecelerator_factor = iniparser_getdouble(ini, "tau:taudecelerator_factor", 0.7);
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
		if(pnt_config->MPI_rank==0){printf("saiWENOs: ERROR: Fehlerhafte Eingabe fuer Inflow-Randbedingung\n");}
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
		if(pnt_config->MPI_rank==0){printf("saiWENOs: ERROR: Fehlerhafte Eingabe fuer Outflow-Randbedingung\n");}
		MPI_Barrier(pnt_config->MPI_comm);
		MPI_Abort(pnt_config->MPI_comm,13372);
	}

	pnt_config->dbl_rho_inflow = iniparser_getdouble(ini, "boundaryconditions:rho_inflow", 1);
	pnt_config->dbl_p_inflow = iniparser_getdouble(ini, "boundaryconditions:p_inflow", 1);
	pnt_config->dbl_AoA = iniparser_getdouble(ini, "boundaryconditions:AoA", 0);

	pnt_config->dbl_p_out = iniparser_getdouble(ini, "boundaryconditions:p_outflow", 1);

	pnt_config->dbl_T_wall = iniparser_getdouble(ini, "boundaryconditions:T_wall", 1);

//	Fluid Properties
	pnt_config->dbl_machNumber = iniparser_getdouble(ini, "fluidproperties:mach", -1);
	pnt_config->dbl_reynoldsNumber = iniparser_getdouble(ini, "fluidproperties:reynolds", -1);
	pnt_config->dbl_prandtlNumber = iniparser_getdouble(ini, "fluidproperties:Prandtl", -1);
	pnt_config->dbl_gammaNumber = iniparser_getdouble(ini, "fluidproperties:Gamma", -1);
	pnt_config->dbl_gasConstantNumber = iniparser_getdouble(ini, "fluidproperties:R", -1);
	pnt_config->dbl_T0_dim = iniparser_getdouble(ini, "fluidproperties:T", -1);
	pnt_config->dbl_L0_dim = iniparser_getdouble(ini, "fluidproperties:L", -1);

//	Export
	pnt_config->flag_exportMetric = iniparser_getboolean(ini, "export:metric", -1);
	pnt_config->flag_exportSnapshot = iniparser_getboolean(ini, "export:snapshot", -1);
	pnt_config->flag_exportFilm = iniparser_getboolean(ini, "export:film", -1);
	pnt_config->flag_ReducedExport = iniparser_getboolean(ini, "export:reduced", 0);

//	Options
	pnt_config->flag_constantZValues = iniparser_getboolean(ini, "options:3Dto2D", 0);
	pnt_config->flag_IBC = iniparser_getboolean(ini, "options:IBC", 0);
	pnt_config->flag_LaminarBoundary = iniparser_getboolean(ini, "options:LaminarBoundary", 0);
	pnt_config->flag_PressureHistory = iniparser_getboolean(ini, "options:PressureHistory", 0);
	pnt_config->flag_VelocityHistory = iniparser_getboolean(ini, "options:VelocityHistory", 0);
	pnt_config->flag_PressureWaves = iniparser_getboolean(ini, "options:PressureWaves", 0);
	pnt_config->flag_SplitMeshFile = iniparser_getint(ini, "options:MeshSplitting", 0);
	pnt_config->flag_Inviscid = iniparser_getboolean(ini, "options:Inviscid", 0);
	pnt_config->flag_Vortex = iniparser_getboolean(ini, "options:Vortex", 0);
	pnt_config->flag_rotation_symmetric = iniparser_getboolean(ini, "options:2D-Rotation-Symmetric", 0);
	pnt_config->flag_BC_Changer = iniparser_getboolean(ini, "options:BC-Changer", 0);


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
		pnt_config->PressureHistory_x_P = (double *)calloc(pnt_config->PressureHistory_No, sizeof(double ));
		pnt_config->PressureHistory_y_P = (double *)calloc(pnt_config->PressureHistory_No, sizeof(double ));
		pnt_config->PressureHistory_z_P = (double *)calloc(pnt_config->PressureHistory_No, sizeof(double ));
		pnt_config->PressureHistory_x_P_real = (double *)calloc(pnt_config->PressureHistory_No, sizeof(double ));
		pnt_config->PressureHistory_y_P_real = (double *)calloc(pnt_config->PressureHistory_No, sizeof(double ));
		pnt_config->PressureHistory_z_P_real = (double *)calloc(pnt_config->PressureHistory_No, sizeof(double ));
		pnt_config->PressureHistory_time=(double *)calloc(pnt_config->int_TotalIterations, sizeof(double ));
		pnt_config->PressureHistory_pressure=(double **)calloc(pnt_config->PressureHistory_No, sizeof(double *));
		for(p=0;p<pnt_config->PressureHistory_No;p++)
		{
			pnt_config->PressureHistory_pressure[p]=(double *)calloc(pnt_config->int_TotalIterations, sizeof(double));
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
		pnt_config->VelocityHistory_x_P = (double *)calloc(pnt_config->VelocityHistory_No, sizeof(double ));
		pnt_config->VelocityHistory_y_P = (double *)calloc(pnt_config->VelocityHistory_No, sizeof(double ));
		pnt_config->VelocityHistory_z_P = (double *)calloc(pnt_config->VelocityHistory_No, sizeof(double ));
		pnt_config->VelocityHistory_x_P_real = (double *)calloc(pnt_config->VelocityHistory_No, sizeof(double ));
		pnt_config->VelocityHistory_y_P_real = (double *)calloc(pnt_config->VelocityHistory_No, sizeof(double ));
		pnt_config->VelocityHistory_z_P_real = (double *)calloc(pnt_config->VelocityHistory_No, sizeof(double ));
		pnt_config->VelocityHistory_time=(double *)calloc(pnt_config->int_TotalIterations, sizeof(double ));
		pnt_config->VelocityHistory_VelocityX=(double **)calloc(pnt_config->VelocityHistory_No, sizeof(double *));
		pnt_config->VelocityHistory_VelocityY=(double **)calloc(pnt_config->VelocityHistory_No, sizeof(double *));
		pnt_config->VelocityHistory_VelocityZ=(double **)calloc(pnt_config->VelocityHistory_No, sizeof(double *));
		for(p=0;p<pnt_config->VelocityHistory_No;p++)
		{
			pnt_config->VelocityHistory_VelocityX[p]=(double *)calloc(pnt_config->int_TotalIterations, sizeof(double));
			pnt_config->VelocityHistory_VelocityY[p]=(double *)calloc(pnt_config->int_TotalIterations, sizeof(double));
			pnt_config->VelocityHistory_VelocityZ[p]=(double *)calloc(pnt_config->int_TotalIterations, sizeof(double));
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
	pnt_config->flag_IBC_ApplyBC = iniparser_getboolean(ini, "IBC:ApplyBC", 0);

	//InitializeValues
	pnt_config->InitializeValues_u0 = iniparser_getdouble(ini, "InitializeValues:u0", 1);
	pnt_config->InitializeValues_p0 = iniparser_getdouble(ini, "InitializeValues:p0", 1);
	pnt_config->InitializeValues_rho0 = iniparser_getdouble(ini, "InitializeValues:rho0", 1);
	pnt_config->InitializeValues_p1 = iniparser_getdouble(ini, "InitializeValues:p1", 1);
	pnt_config->InitializeValues_rho1 = iniparser_getdouble(ini, "InitializeValues:rho1", 1);
	pnt_config->InitializeValues_u1 = iniparser_getdouble(ini, "InitializeValues:u1", 1);
	pnt_config->dbl_u_inflow = cos(pnt_config->dbl_AoA/360*2*M_PI)*pnt_config->InitializeValues_u0;
	pnt_config->dbl_v_inflow = sin(pnt_config->dbl_AoA/360*2*M_PI)*pnt_config->InitializeValues_u0;
	pnt_config->dbl_w_inflow = 0.0;


    iniparser_freedict(ini);

	if((pnt_config->flag_rotation_symmetric==1)&&(pnt_config->int_meshDimensions==3))
	{
		if(pnt_config->MPI_rank==0){printf("saiWENOs: ERROR: Rotationssymmetrisch und 3D schliesst sich aus\n");}
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

		float* x;
		float* y;
		float* z;
		x = (float*) calloc(buffer,sizeof(float));
		y = (float*) calloc(buffer,sizeof(float));
		z = (float*) calloc(buffer,sizeof(float));

	/* READ X, Y, Z GRID POINTS FROM CGNS FILE */
	/* open CGNS file for read-only */
	    cg_open(pnt_config->chr_MeshPath,CG_MODE_READ,&index_file);
	/* we know there is only one base (real working code would check!) */
	    index_base=1;
	    //Sofern es nur eine Mesh-File gibt, bestimmt der rank die zone-ID.
	    //Bei einer gesplitteten Mesh-File ist diese stets 1
	    if((pnt_config->flag_SplitMeshFile==2)||(pnt_config->int_initializeType==1))
	    {
	    	index_zone=1;
	    }
	    else
	    {
	    	index_zone=pnt_config->MPI_rank+1;

		    cg_nzones(index_file,index_base,&nzone);
		    if (nzone != pnt_config->MPI_size)
		    {
		      printf("saiWENOS: Error. Es werden %d Zonen erwartet. %d sind vorhanden.\n",pnt_config->MPI_size,nzone);
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


void ResultImport_CGNS(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep)
{
	/*
	  dimension statements (note that tri-dimensional arrays
	  x,y,z must be dimensioned exactly as [N][17][21] (N>=9)
	  for this particular case or else they will be read from
	  the CGNS file incorrectly!  Other options are to use 1-D
	  arrays, use dynamic memory, or pass index values to a
	  subroutine and dimension exactly there):
	*/
	    int celldim=pnt_config->int_meshDimensions;
		cgsize_t isize[3][celldim],irmin[celldim],irmax[celldim];
	    int index_file,index_base,index_zone;
	    char zonename[33];
	    char actual_path[500];

	//Speicherallokierung für die dynamischen-1D Arrays, worin die 3D-Lösungen gespeichert werden
	    int buffer=pnt_config->int_iMeshPoints*pnt_config->int_jMeshPoints*pnt_config->int_kMeshPoints;
		int i,j,k,ijk,ijk2;
		int i2,j2,k2;
		int index_flow;

		float* u;
		float* v;
		float* w;
		float* p;
		float* rho;
		u = (float*) calloc(buffer,sizeof(float));
		v = (float*) calloc(buffer,sizeof(float));
		w = (float*) calloc(buffer,sizeof(float));
		p = (float*) calloc(buffer,sizeof(float));
		rho = (float*) calloc(buffer,sizeof(float));


	/* we know there is only one base (real working code would check!) */
	    index_base=1;

	    //Sofern es nur eine Mesh-File gibt, bestimmt der rank die zone-ID.
	    //Bei einer gesplitteten Mesh-File ist diese stets 1
    	cg_open(pnt_config->chr_SnapshotPath,CG_MODE_READ,&index_file);

    	strcpy(actual_path,pnt_config->chr_SnapshotPath);
    	index_zone=1;

//    	int nzone;
//	    if(pnt_config->flag_SplitMeshFile==2)
//	    {
//	    	cg_open(pnt_config->chr_SnapshotPath,CG_MODE_READ,&index_file);
//	    	strcpy(actual_path,pnt_config->chr_SnapshotPath);
//	    	index_zone=1;
//	    }
//	    else
//	    {
//	    	cg_open(pnt_config->chr_MeshPath,CG_MODE_READ,&index_file);
//	    	strcpy(actual_path,pnt_config->chr_MeshPath);
//	    	index_zone=pnt_config->MPI_rank+1;
//		    cg_nzones(index_file,index_base,&nzone);
//		    if (nzone != pnt_config->MPI_size)
//		    {
//		      printf("saiWENOS: Error. Es werden %d Zonen erwartet. %d sind vorhanden.\n",pnt_config->MPI_size,nzone);
//		    }
//	    }
	    cg_zone_read(index_file,index_base,index_zone,zonename,isize[0]);

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

		char *text;
		char name[33];
		int ndescriptors,D;
		cg_goto(index_file,index_base,"end");
		cg_ndescriptors(&ndescriptors);
		for(D=1;D<=ndescriptors;D++)
		{
		  cg_descriptor_read(D, name, &text);
		  if (strcmp(name,"Iterations")==0)
		  {
			  pnt_config->int_StartIteration = atoi( text );
			  pnt_config->int_actualIteration=pnt_config->int_StartIteration;
			  pnt_config->int_EndIteration=pnt_config->int_StartIteration+pnt_config->int_TotalIterations;
		  }
		  if (strcmp(name,"ActualTime")==0)
		  {
			  pnt_config->start_Time = strtod( text,NULL );
			  pnt_config->dbl_time_dim=pnt_config->start_Time;
			  pnt_config->dbl_time_dim_lastAction=pnt_config->start_Time;

		  }
		}

	    SimulationType_t SimulationType;
	    cg_simulation_type_read(index_file,index_base,&SimulationType);
	    char BaseIterName[30];
	    if(SimulationType==TimeAccurate)
		{
		    cg_biter_read(index_file, index_base, BaseIterName, &index_flow);
			cg_goto(index_file,index_base,"BaseIterativeData_t",1,"end");
			double time[index_flow];
			cg_array_read_as(1,RealDouble,&time);
			pnt_config->start_Time=time[index_flow-1];
		    if(pnt_config->MPI_rank==0)
		    {
		    	printf("saiWENOs: Es gibt %d timesteps in Datei %s. Der letzte (t=%le) wird importiert.\n",
		    			index_flow,actual_path,pnt_config->start_Time);
		    }
		}
	    else
	    {
	    	index_flow=1;
		    if(pnt_config->MPI_rank==0)
		    {
		    	printf("saiWENOs: Es gibt nur einen timestep in Datei %s. Dieser wird importiert.\n",
		    			actual_path);
		    }
	    }



	/* read grid coordinates */
	    cg_field_read(index_file,index_base,index_zone,index_flow,"VelocityX",
	    		RealSingle,irmin,irmax,u);
	    cg_field_read(index_file,index_base,index_zone,index_flow,"VelocityY",
	    		RealSingle,irmin,irmax,v);
	    cg_field_read(index_file,index_base,index_zone,index_flow,"VelocityZ",
	    		RealSingle,irmin,irmax,w);
	    cg_field_read(index_file,index_base,index_zone,index_flow,"Pressure",
	    		RealSingle,irmin,irmax,p);
	    cg_field_read(index_file,index_base,index_zone,index_flow,"Density",
	    		RealSingle,irmin,irmax,rho);

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

					pnt_U_lastStep->u[ijk]=u[ijk2];
					pnt_U_lastStep->v[ijk]=v[ijk2];
					pnt_U_lastStep->w[ijk]=w[ijk2];
					pnt_U_lastStep->rho[ijk]=rho[ijk2];
					pnt_U_lastStep->p[ijk]=p[ijk2];

					pnt_U_lastStep->e[ijk]=(0.5*((pnt_U_lastStep->u[ijk]*pnt_U_lastStep->u[ijk])+(pnt_U_lastStep->v[ijk]*pnt_U_lastStep->v[ijk])+(pnt_U_lastStep->w[ijk]*pnt_U_lastStep->w[ijk]))+
							pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk]/(pnt_config->dbl_gammaNumber-1.0)*pnt_config->dbl_Upsilon);

				}
			}
		}


	/* close CGNS file */
	    cg_close(index_file);

		free(u);
		free(v);
		free(w);
		free(rho);
		free(p);
}


void MeshConfig_CGNS(
		struct strct_configuration * pnt_config)
{
	/*
	  dimension statements (note that tri-dimensional arrays
	  x,y,z must be dimensioned exactly as [N][17][21] (N>=9)
	  for this particular case or else they will be read from
	  the CGNS file incorrectly!  Other options are to use 1-D
	  arrays, use dynamic memory, or pass index values to a
	  subroutine and dimension exactly there):
	*/
	int i;
	int index_file,index_base,index_zone;
	int cellDim;


//	if (pnt_config->int_initializeType==1)
//	{
//		cg_open(pnt_config->chr_SnapshotPath,CG_MODE_READ,&index_file);
//	}
//	else
//	{
//		cg_open(pnt_config->chr_MeshPath,CG_MODE_READ,&index_file);
//	}

	cg_open(pnt_config->chr_MeshPath,CG_MODE_READ,&index_file);

    index_base=1;

    if((pnt_config->flag_SplitMeshFile==2)||(pnt_config->int_initializeType==1))
    {
    	index_zone=1;
    }
	else
	{
		index_zone=pnt_config->MPI_rank+1;
	}

	cg_cell_dim(index_file, index_base, &cellDim);
	pnt_config->int_meshDimensions=cellDim;
	if((MESHDIMENSIONS!=pnt_config->int_meshDimensions)&&(pnt_config->MPI_rank==0))
	{
		printf("ERROR: Es wird versucht ein %dD-Loeser auf einem %dD-Gitter anzuwenden!\n",MESHDIMENSIONS,pnt_config->int_meshDimensions);
		MPI_Abort(pnt_config->MPI_comm,13373);
	}


	cgsize_t isize[3][cellDim];
	cg_zone_read(index_file,index_base,index_zone,pnt_config->Zonename,isize[0]);
	cg_close(index_file);

	pnt_config->int_meshDimensions=MESHDIMENSIONS;

	pnt_config->int_iMeshPoints=isize[0][0];
	pnt_config->int_jMeshPoints=isize[0][1];
	pnt_config->int_iMeshPointsGhostCells=pnt_config->int_iMeshPoints+(pnt_config->int_SpaceOrder+1);
	pnt_config->int_jMeshPointsGhostCells=pnt_config->int_jMeshPoints+(pnt_config->int_SpaceOrder+1);

	if(cellDim==2)
	{
		pnt_config->int_kMeshPoints=1;
		pnt_config->int_kMeshPointsGhostCells=1;
	}
	else
	{
		pnt_config->int_kMeshPoints=isize[0][2];
		pnt_config->int_kMeshPointsGhostCells=isize[0][2]+(pnt_config->int_SpaceOrder+1);
	}

	pnt_config->int_iStartGhosts=0;
	pnt_config->int_jStartGhosts=0;
	pnt_config->int_kStartGhosts=0;
	pnt_config->int_iEndGhosts=pnt_config->int_iStartGhosts+pnt_config->int_iMeshPoints+pnt_config->int_SpaceOrder;
	pnt_config->int_jEndGhosts=pnt_config->int_jStartGhosts+pnt_config->int_jMeshPoints+pnt_config->int_SpaceOrder;
	pnt_config->int_kEndGhosts=pnt_config->int_kStartGhosts+pnt_config->int_kMeshPoints+pnt_config->int_SpaceOrder;

	pnt_config->int_iStartReal=(pnt_config->int_SpaceOrder-1)/2+1;
	pnt_config->int_jStartReal=(pnt_config->int_SpaceOrder-1)/2+1;
	pnt_config->int_kStartReal=(pnt_config->int_SpaceOrder-1)/2+1;
	pnt_config->int_iEndReal=pnt_config->int_iStartReal+pnt_config->int_iMeshPoints-1;
	pnt_config->int_jEndReal=pnt_config->int_jStartReal+pnt_config->int_jMeshPoints-1;
	pnt_config->int_kEndReal=pnt_config->int_kStartReal+pnt_config->int_kMeshPoints-1;

	pnt_config->int_iMid=(int)(0.5*(pnt_config->int_iEndGhosts+pnt_config->int_iStartGhosts));
	pnt_config->int_jMid=(int)(0.5*(pnt_config->int_jEndGhosts+pnt_config->int_jStartGhosts));
	pnt_config->int_kMid=(int)(0.5*(pnt_config->int_kEndGhosts+pnt_config->int_kStartGhosts));

	if(MESHDIMENSIONS==2)
	{
		pnt_config->int_kStartReal=0;
		pnt_config->int_kEndReal=0;

		pnt_config->int_kStartGhosts=0;
		pnt_config->int_kEndGhosts=0;

		pnt_config->int_kMid=0;
	}
	pnt_config->int_kStartReal_original=pnt_config->int_kStartReal;
	pnt_config->int_kEndReal_original=pnt_config->int_kEndReal;

//	Extrapolate
	pnt_config->int_iStartReal_extrapolate=pnt_config->int_iStartReal+1;
	pnt_config->int_jStartReal_extrapolate=pnt_config->int_jStartReal+1;
	pnt_config->int_kStartReal_extrapolate=pnt_config->int_kStartReal+1;
	pnt_config->int_iEndReal_extrapolate=pnt_config->int_iEndReal+1;
	pnt_config->int_jEndReal_extrapolate=pnt_config->int_jEndReal+1;
	pnt_config->int_kEndReal_extrapolate=pnt_config->int_kEndReal+1;

	pnt_config->int_iStartGhosts_extrapolate=pnt_config->int_iStartGhosts;
	pnt_config->int_jStartGhosts_extrapolate=pnt_config->int_jStartGhosts;
	pnt_config->int_kStartGhosts_extrapolate=pnt_config->int_kStartGhosts;
	pnt_config->int_iEndGhosts_extrapolate=pnt_config->int_iEndGhosts+2;
	pnt_config->int_jEndGhosts_extrapolate=pnt_config->int_jEndGhosts+2;
	pnt_config->int_kEndGhosts_extrapolate=pnt_config->int_kEndGhosts+2;

	pnt_config->int_iMeshPointsGhostCells_extrapolate=pnt_config->int_iMeshPointsGhostCells+2;
	pnt_config->int_jMeshPointsGhostCells_extrapolate=pnt_config->int_jMeshPointsGhostCells+2;
	pnt_config->int_kMeshPointsGhostCells_extrapolate=pnt_config->int_kMeshPointsGhostCells+2;

	if(MESHDIMENSIONS==2)
	{
		pnt_config->int_kStartReal_extrapolate=0;
		pnt_config->int_kEndReal_extrapolate=0;

		pnt_config->int_kStartGhosts_extrapolate=0;
		pnt_config->int_kEndGhosts_extrapolate=0;

		pnt_config->int_kMeshPointsGhostCells_extrapolate=1;

	}


	//		Transfer Interface-Informations
	getInterfaceInformations(
			pnt_config);

//	pnt_config->bufferSendFlowLeft = (double *)calloc(5*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPoints*pnt_config->int_kMeshPoints,sizeof(double));
//	pnt_config->bufferSendFlowRight= (double *)calloc(5*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPoints*pnt_config->int_kMeshPoints,sizeof(double));
//	pnt_config->bufferSendFlowBottom = (double *)calloc(5*pnt_config->int_iMeshPoints * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPoints,sizeof(double));
//	pnt_config->bufferSendFlowTop= (double *)calloc(5*pnt_config->int_iMeshPoints * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPoints,sizeof(double));
//	pnt_config->bufferRecieveFlowLeft = (double *)calloc(5*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPoints*pnt_config->int_kMeshPoints,sizeof(double));
//	pnt_config->bufferRecieveFlowRight= (double *)calloc(5*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPoints*pnt_config->int_kMeshPoints,sizeof(double));
//	pnt_config->bufferRecieveFlowBottom = (double *)calloc(5*pnt_config->int_iMeshPoints * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPoints,sizeof(double));
//	pnt_config->bufferRecieveFlowTop= (double *)calloc(5*pnt_config->int_iMeshPoints * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPoints,sizeof(double));

	pnt_config->bufferSendMeshLeft = (double *)calloc(13*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells,sizeof(double));
	pnt_config->bufferSendMeshRight = (double *)calloc(13*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells,sizeof(double));
	pnt_config->bufferSendMeshBottom = (double *)calloc(13*pnt_config->int_iMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPointsGhostCells,sizeof(double));
	pnt_config->bufferSendMeshTop = (double *)calloc(13*pnt_config->int_iMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPointsGhostCells,sizeof(double));
	pnt_config->bufferRecieveMeshLeft = (double *)calloc(13*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells,sizeof(double));
	pnt_config->bufferRecieveMeshRight = (double *)calloc(13*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells,sizeof(double));
	pnt_config->bufferRecieveMeshBottom = (double *)calloc(13*pnt_config->int_iMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPointsGhostCells,sizeof(double));
	pnt_config->bufferRecieveMeshTop = (double *)calloc(13*pnt_config->int_iMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPointsGhostCells,sizeof(double));

	if(MESHDIMENSIONS==2)
	{
		pnt_config->bufferSendFlowWithGhostsLeft = (double *)calloc(4*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells,sizeof(double));
		pnt_config->bufferSendFlowWithGhostsRight= (double *)calloc(4*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells,sizeof(double));
		pnt_config->bufferSendFlowWithGhostsBottom = (double *)calloc(4*pnt_config->int_iMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPointsGhostCells,sizeof(double));
		pnt_config->bufferSendFlowWithGhostsTop= (double *)calloc(4*pnt_config->int_iMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPointsGhostCells,sizeof(double));

		pnt_config->bufferRecieveFlowWithGhostsLeft = (double *)calloc(4*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells,sizeof(double));
		pnt_config->bufferRecieveFlowWithGhostsRight= (double *)calloc(4*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells,sizeof(double));
		pnt_config->bufferRecieveFlowWithGhostsBottom = (double *)calloc(4*pnt_config->int_iMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPointsGhostCells,sizeof(double));
		pnt_config->bufferRecieveFlowWithGhostsTop= (double *)calloc(4*pnt_config->int_iMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPointsGhostCells,sizeof(double));
	}
	if(MESHDIMENSIONS==3)
	{
		pnt_config->bufferSendFlowWithGhostsLeft = (double *)calloc(5*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells,sizeof(double));
		pnt_config->bufferSendFlowWithGhostsRight= (double *)calloc(5*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells,sizeof(double));
		pnt_config->bufferSendFlowWithGhostsBottom = (double *)calloc(5*pnt_config->int_iMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPointsGhostCells,sizeof(double));
		pnt_config->bufferSendFlowWithGhostsTop= (double *)calloc(5*pnt_config->int_iMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPointsGhostCells,sizeof(double));

		pnt_config->bufferRecieveFlowWithGhostsLeft = (double *)calloc(5*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells,sizeof(double));
		pnt_config->bufferRecieveFlowWithGhostsRight= (double *)calloc(5*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells,sizeof(double));
		pnt_config->bufferRecieveFlowWithGhostsBottom = (double *)calloc(5*pnt_config->int_iMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPointsGhostCells,sizeof(double));
		pnt_config->bufferRecieveFlowWithGhostsTop= (double *)calloc(5*pnt_config->int_iMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPointsGhostCells,sizeof(double));

		pnt_config->bufferSendFlowWithGhostsBehind = (double *)calloc(5*pnt_config->int_iMeshPointsGhostCells * pnt_config->int_jMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2),sizeof(double));
		pnt_config->bufferSendFlowWithGhostsInFront= (double *)calloc(5*pnt_config->int_iMeshPointsGhostCells * pnt_config->int_jMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2),sizeof(double));
		pnt_config->bufferRecieveFlowWithGhostsBehind = (double *)calloc(5*pnt_config->int_iMeshPointsGhostCells * pnt_config->int_jMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2),sizeof(double));
		pnt_config->bufferRecieveFlowWithGhostsInFront= (double *)calloc(5*pnt_config->int_iMeshPointsGhostCells * pnt_config->int_jMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2),sizeof(double));

		pnt_config->bufferSendMeshBehind = (double *)calloc(13*pnt_config->int_iMeshPointsGhostCells * pnt_config->int_jMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2),sizeof(double));
		pnt_config->bufferSendMeshInFront = (double *)calloc(13*pnt_config->int_iMeshPointsGhostCells * pnt_config->int_jMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2),sizeof(double));
		pnt_config->bufferRecieveMeshBehind = (double *)calloc(13*pnt_config->int_iMeshPointsGhostCells * pnt_config->int_jMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2),sizeof(double));
		pnt_config->bufferRecieveMeshInFront= (double *)calloc(13*pnt_config->int_iMeshPointsGhostCells * pnt_config->int_jMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2),sizeof(double));
	}


	/* get neighbour */
	pnt_config->MPI_rankNeighbours= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_tag= (long *)calloc(pnt_config->NumberInterfaces, sizeof(long));
	getNeighbour(
			pnt_config);


	pnt_config->MPI_intTransformation_IMax= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intTransformation_JMax= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intTransformation_KMax= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intTransformation_IMax_Mesh= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intTransformation_JMax_Mesh= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intTransformation_KMax_Mesh= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intTransformation_flag_I0_I= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intTransformation_flag_I0_J= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intTransformation_flag_I0_K= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intTransformation_flag_J0_I= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intTransformation_flag_J0_J= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intTransformation_flag_J0_K= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intTransformation_flag_K0_I= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intTransformation_flag_K0_J= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intTransformation_flag_K0_K= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intTransformation_Offset_I= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intTransformation_Offset_J= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intTransformation_Offset_K= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intTransformation_Offset_I_Ghosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intTransformation_Offset_J_Ghosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intTransformation_Offset_K_Ghosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));

	pnt_config->MPI_dblTransformation_xi_x= (double *)calloc(pnt_config->NumberInterfaces, sizeof(double));
	pnt_config->MPI_dblTransformation_xi_y= (double *)calloc(pnt_config->NumberInterfaces, sizeof(double));
	pnt_config->MPI_dblTransformation_xi_z= (double *)calloc(pnt_config->NumberInterfaces, sizeof(double));
	pnt_config->MPI_dblTransformation_eta_x= (double *)calloc(pnt_config->NumberInterfaces, sizeof(double));
	pnt_config->MPI_dblTransformation_eta_y= (double *)calloc(pnt_config->NumberInterfaces, sizeof(double));
	pnt_config->MPI_dblTransformation_eta_z= (double *)calloc(pnt_config->NumberInterfaces, sizeof(double));
	pnt_config->MPI_dblTransformation_zeta_x= (double *)calloc(pnt_config->NumberInterfaces, sizeof(double));
	pnt_config->MPI_dblTransformation_zeta_y= (double *)calloc(pnt_config->NumberInterfaces, sizeof(double));
	pnt_config->MPI_dblTransformation_zeta_z= (double *)calloc(pnt_config->NumberInterfaces, sizeof(double));

	//		pnt_config->MPI_charNeighbours= (char **)calloc(pnt_config->NumberInterfaces, sizeof(char*));
	for(i=0;i<pnt_config->NumberInterfaces;i++)
	{
	//			pnt_config->RankNeighbours[i]= (int *)calloc(pnt_config->int_meshDimensions+1, sizeof(int));
	//			pnt_config->MPI_charNeighbours[i]= (char *)calloc(30, sizeof(char));

		pnt_config->MPI_dblTransformation_xi_x[i]=1.0;
		pnt_config->MPI_dblTransformation_xi_y[i]=1.0;
		pnt_config->MPI_dblTransformation_xi_z[i]=1.0;
		pnt_config->MPI_dblTransformation_eta_x[i]=1.0;
		pnt_config->MPI_dblTransformation_eta_y[i]=1.0;
		pnt_config->MPI_dblTransformation_eta_z[i]=1.0;
		pnt_config->MPI_dblTransformation_zeta_x[i]=1.0;
		pnt_config->MPI_dblTransformation_zeta_y[i]=1.0;
		pnt_config->MPI_dblTransformation_zeta_z[i]=1.0;
	}

	pnt_config->MPI_intTransferSizeMesh=(int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intTransferSizeFlow_WithGhosts=(int *)calloc(pnt_config->NumberInterfaces, sizeof(int));

	pnt_config->MPI_SendBufferMesh= (double **)calloc(pnt_config->NumberInterfaces, sizeof(double *));
	pnt_config->MPI_RecieveBufferMesh= (double **)calloc(pnt_config->NumberInterfaces, sizeof(double *));
	pnt_config->MPI_SendBufferFlowWithGhosts= (double **)calloc(pnt_config->NumberInterfaces, sizeof(double *));
	pnt_config->MPI_RecieveBufferFlowWithGhosts= (double **)calloc(pnt_config->NumberInterfaces, sizeof(double *));

	pnt_config->MPI_intIStartSend= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intIEndSend= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intJStartSend= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intJEndSend= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intKStartSend= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intKEndSend= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));

	pnt_config->MPI_intIStartRecieve= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intIEndRecieve= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intJStartRecieve= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intJEndRecieve= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intKStartRecieve= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intKEndRecieve= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));

	pnt_config->MPI_intIStartSend_WithGhosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intIEndSend_WithGhosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intJStartSend_WithGhosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intJEndSend_WithGhosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intKStartSend_WithGhosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intKEndSend_WithGhosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));

	pnt_config->MPI_intIStartRecieve_WithGhosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intIEndRecieve_WithGhosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intJStartRecieve_WithGhosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intJEndRecieve_WithGhosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intKStartRecieve_WithGhosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
	pnt_config->MPI_intKEndRecieve_WithGhosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));




	//		Teile den Interfaces, an denen der einzelne Rank beteiligt ist, die Informationen zu,
	//		ob dieses Interface aus seiner Sicht links, rechts... ist.
	for(i=0;i<pnt_config->NumberInterfaces;i++)
	{
			if(i==pnt_config->InterfaceNeighbourLeft)
			{
//				pnt_config->MPI_SendBufferMesh[i]=pnt_config->bufferSendMeshI;
//				pnt_config->MPI_RecieveBufferMesh[i]=pnt_config->bufferRecieveMeshI;
//				pnt_config->MPI_SendBufferFlow[i]=pnt_config->bufferSendFlowI;
//				pnt_config->MPI_RecieveBufferFlow[i]=pnt_config->bufferRecieveFlowI;
//				pnt_config->MPI_SendBufferViscid[i]=pnt_config->bufferSendViscidI;
//				pnt_config->MPI_RecieveBufferViscid[i]=pnt_config->bufferRecieveViscidI;

				pnt_config->MPI_SendBufferMesh[i]=pnt_config->bufferSendMeshLeft;
				pnt_config->MPI_RecieveBufferMesh[i]=pnt_config->bufferRecieveMeshLeft;
				pnt_config->MPI_SendBufferFlowWithGhosts[i]=pnt_config->bufferSendFlowWithGhostsLeft;
				pnt_config->MPI_RecieveBufferFlowWithGhosts[i]=pnt_config->bufferRecieveFlowWithGhostsLeft;

			   pnt_config->MPI_intIStartSend[i]=pnt_config->int_iStartReal;
			   pnt_config->MPI_intIEndSend[i]=pnt_config->int_iStartReal+((pnt_config->int_SpaceOrder+1)/2-1);
			   pnt_config->MPI_intJStartSend[i]=pnt_config->int_jStartReal;
			   pnt_config->MPI_intJEndSend[i]=pnt_config->int_jEndReal;
			   pnt_config->MPI_intKStartSend[i]=pnt_config->int_kStartReal;
			   pnt_config->MPI_intKEndSend[i]=pnt_config->int_kEndReal;

			   pnt_config->MPI_intIStartRecieve[i]=pnt_config->int_iStartGhosts;
			   pnt_config->MPI_intIEndRecieve[i]=pnt_config->int_iStartReal-1;
			   pnt_config->MPI_intJStartRecieve[i]=pnt_config->int_jStartReal;
			   pnt_config->MPI_intJEndRecieve[i]=pnt_config->int_jEndReal;
			   pnt_config->MPI_intKStartRecieve[i]=pnt_config->int_kStartReal;
			   pnt_config->MPI_intKEndRecieve[i]=pnt_config->int_kEndReal;


			   pnt_config->MPI_intIStartSend_WithGhosts[i]=pnt_config->int_iStartReal;
			   pnt_config->MPI_intIEndSend_WithGhosts[i]=pnt_config->int_iStartReal+((pnt_config->int_SpaceOrder+1)/2-1);
			   pnt_config->MPI_intJStartSend_WithGhosts[i]=pnt_config->int_jStartGhosts;
			   pnt_config->MPI_intJEndSend_WithGhosts[i]=pnt_config->int_jEndGhosts;
			   pnt_config->MPI_intKStartSend_WithGhosts[i]=pnt_config->int_kStartGhosts;
			   pnt_config->MPI_intKEndSend_WithGhosts[i]=pnt_config->int_kEndGhosts;

			   pnt_config->MPI_intIStartRecieve_WithGhosts[i]=pnt_config->int_iStartGhosts;
			   pnt_config->MPI_intIEndRecieve_WithGhosts[i]=pnt_config->int_iStartReal-1;
			   pnt_config->MPI_intJStartRecieve_WithGhosts[i]=pnt_config->int_jStartGhosts;
			   pnt_config->MPI_intJEndRecieve_WithGhosts[i]=pnt_config->int_jEndGhosts;
			   pnt_config->MPI_intKStartRecieve_WithGhosts[i]=pnt_config->int_kStartGhosts;
			   pnt_config->MPI_intKEndRecieve_WithGhosts[i]=pnt_config->int_kEndGhosts;
			   if(MESHDIMENSIONS==2)
			   {
				   pnt_config->MPI_intKStartSend_WithGhosts[i]=0;
				   pnt_config->MPI_intKEndSend_WithGhosts[i]=0;
				   pnt_config->MPI_intKStartRecieve_WithGhosts[i]=0;
				   pnt_config->MPI_intKEndRecieve_WithGhosts[i]=0;

				   pnt_config->MPI_intTransferSizeFlow_WithGhosts[i]=
						   (pnt_config->int_SpaceOrder+1)/2*4*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells;
			   }
			   else
			   {
				   pnt_config->MPI_intTransferSizeFlow_WithGhosts[i]=
					   (pnt_config->int_SpaceOrder+1)/2*5*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells;
			   }
			   pnt_config->MPI_intTransferSizeMesh[i]=
					   (pnt_config->int_SpaceOrder+1)/2*13*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells;



			}
				if(i==pnt_config->InterfaceNeighbourRight)
			{
//				pnt_config->MPI_SendBufferMesh[i]=pnt_config->bufferSendMeshI;
//				pnt_config->MPI_RecieveBufferMesh[i]=pnt_config->bufferRecieveMeshI;
//				pnt_config->MPI_SendBufferFlow[i]=pnt_config->bufferSendFlowI;
//				pnt_config->MPI_RecieveBufferFlow[i]=pnt_config->bufferRecieveFlowI;
//				pnt_config->MPI_SendBufferViscid[i]=pnt_config->bufferSendViscidI;
//				pnt_config->MPI_RecieveBufferViscid[i]=pnt_config->bufferRecieveViscidI;

				pnt_config->MPI_SendBufferMesh[i]=pnt_config->bufferSendMeshRight;
				pnt_config->MPI_RecieveBufferMesh[i]=pnt_config->bufferRecieveMeshRight;
				pnt_config->MPI_SendBufferFlowWithGhosts[i]=pnt_config->bufferSendFlowWithGhostsRight;
				pnt_config->MPI_RecieveBufferFlowWithGhosts[i]=pnt_config->bufferRecieveFlowWithGhostsRight;

			   pnt_config->MPI_intIStartSend[i]=pnt_config->int_iEndReal-((pnt_config->int_SpaceOrder+1)/2-1);
			   pnt_config->MPI_intIEndSend[i]=pnt_config->int_iEndReal;
			   pnt_config->MPI_intJStartSend[i]=pnt_config->int_jStartReal;
			   pnt_config->MPI_intJEndSend[i]=pnt_config->int_jEndReal;
			   pnt_config->MPI_intKStartSend[i]=pnt_config->int_kStartReal;
			   pnt_config->MPI_intKEndSend[i]=pnt_config->int_kEndReal;

			   pnt_config->MPI_intIStartRecieve[i]=pnt_config->int_iEndReal+1;
			   pnt_config->MPI_intIEndRecieve[i]=pnt_config->int_iEndGhosts;
			   pnt_config->MPI_intJStartRecieve[i]=pnt_config->int_jStartReal;
			   pnt_config->MPI_intJEndRecieve[i]=pnt_config->int_jEndReal;
			   pnt_config->MPI_intKStartRecieve[i]=pnt_config->int_kStartReal;
			   pnt_config->MPI_intKEndRecieve[i]=pnt_config->int_kEndReal;


			   pnt_config->MPI_intIStartSend_WithGhosts[i]=pnt_config->int_iEndReal-((pnt_config->int_SpaceOrder+1)/2-1);
			   pnt_config->MPI_intIEndSend_WithGhosts[i]=pnt_config->int_iEndReal;
			   pnt_config->MPI_intJStartSend_WithGhosts[i]=pnt_config->int_jStartGhosts;
			   pnt_config->MPI_intJEndSend_WithGhosts[i]=pnt_config->int_jEndGhosts;
			   pnt_config->MPI_intKStartSend_WithGhosts[i]=pnt_config->int_kStartGhosts;
			   pnt_config->MPI_intKEndSend_WithGhosts[i]=pnt_config->int_kEndGhosts;

			   pnt_config->MPI_intIStartRecieve_WithGhosts[i]=pnt_config->int_iEndReal+1;
			   pnt_config->MPI_intIEndRecieve_WithGhosts[i]=pnt_config->int_iEndGhosts;
			   pnt_config->MPI_intJStartRecieve_WithGhosts[i]=pnt_config->int_jStartGhosts;
			   pnt_config->MPI_intJEndRecieve_WithGhosts[i]=pnt_config->int_jEndGhosts;
			   pnt_config->MPI_intKStartRecieve_WithGhosts[i]=pnt_config->int_kStartGhosts;
			   pnt_config->MPI_intKEndRecieve_WithGhosts[i]=pnt_config->int_kEndGhosts;
			   if(MESHDIMENSIONS==2)
			   {
				   pnt_config->MPI_intKStartSend_WithGhosts[i]=0;
				   pnt_config->MPI_intKEndSend_WithGhosts[i]=0;
				   pnt_config->MPI_intKStartRecieve_WithGhosts[i]=0;
				   pnt_config->MPI_intKEndRecieve_WithGhosts[i]=0;
				   pnt_config->MPI_intTransferSizeFlow_WithGhosts[i]=
				   					   (pnt_config->int_SpaceOrder+1)/2*4*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells;
			   }
			   else
			   {
				   pnt_config->MPI_intTransferSizeFlow_WithGhosts[i]=
				   					   (pnt_config->int_SpaceOrder+1)/2*5*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells;
			   }

			   pnt_config->MPI_intTransferSizeMesh[i]=
					   (pnt_config->int_SpaceOrder+1)/2*13*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells;

			}
				if(i==pnt_config->InterfaceNeighbourBottom)
			{
//				pnt_config->MPI_SendBufferMesh[i]=pnt_config->bufferSendMeshJ;
//				pnt_config->MPI_RecieveBufferMesh[i]=pnt_config->bufferRecieveMeshJ;
//				pnt_config->MPI_SendBufferFlow[i]=pnt_config->bufferSendFlowJ;
//				pnt_config->MPI_RecieveBufferFlow[i]=pnt_config->bufferRecieveFlowJ;
//				pnt_config->MPI_SendBufferViscid[i]=pnt_config->bufferSendViscidJ;
//				pnt_config->MPI_RecieveBufferViscid[i]=pnt_config->bufferRecieveViscidJ;

				pnt_config->MPI_SendBufferMesh[i]=pnt_config->bufferSendMeshBottom;
				pnt_config->MPI_RecieveBufferMesh[i]=pnt_config->bufferRecieveMeshBottom;
				pnt_config->MPI_SendBufferFlowWithGhosts[i]=pnt_config->bufferSendFlowWithGhostsBottom;
				pnt_config->MPI_RecieveBufferFlowWithGhosts[i]=pnt_config->bufferRecieveFlowWithGhostsBottom;

			   pnt_config->MPI_intJStartSend[i]=pnt_config->int_jStartReal;
			   pnt_config->MPI_intJEndSend[i]=pnt_config->int_jStartReal+((pnt_config->int_SpaceOrder+1)/2-1);
			   pnt_config->MPI_intIStartSend[i]=pnt_config->int_iStartReal;
			   pnt_config->MPI_intIEndSend[i]=pnt_config->int_iEndReal;
			   pnt_config->MPI_intKStartSend[i]=pnt_config->int_kStartReal;
			   pnt_config->MPI_intKEndSend[i]=pnt_config->int_kEndReal;

			   pnt_config->MPI_intJStartRecieve[i]=pnt_config->int_jStartGhosts;
			   pnt_config->MPI_intJEndRecieve[i]=pnt_config->int_jStartReal-1;
			   pnt_config->MPI_intIStartRecieve[i]=pnt_config->int_iStartReal;
			   pnt_config->MPI_intIEndRecieve[i]=pnt_config->int_iEndReal;
			   pnt_config->MPI_intKStartRecieve[i]=pnt_config->int_kStartReal;
			   pnt_config->MPI_intKEndRecieve[i]=pnt_config->int_kEndReal;


			   pnt_config->MPI_intIStartSend_WithGhosts[i]=pnt_config->int_iStartGhosts;
			   pnt_config->MPI_intIEndSend_WithGhosts[i]=pnt_config->int_iEndGhosts;
			   pnt_config->MPI_intJStartSend_WithGhosts[i]=pnt_config->int_jStartReal;
			   pnt_config->MPI_intJEndSend_WithGhosts[i]=pnt_config->int_jStartReal+((pnt_config->int_SpaceOrder+1)/2-1);
			   pnt_config->MPI_intKStartSend_WithGhosts[i]=pnt_config->int_kStartGhosts;
			   pnt_config->MPI_intKEndSend_WithGhosts[i]=pnt_config->int_kEndGhosts;

			   pnt_config->MPI_intIStartRecieve_WithGhosts[i]=pnt_config->int_iStartGhosts;
			   pnt_config->MPI_intIEndRecieve_WithGhosts[i]=pnt_config->int_iEndGhosts;
			   pnt_config->MPI_intJStartRecieve_WithGhosts[i]=pnt_config->int_jStartGhosts;
			   pnt_config->MPI_intJEndRecieve_WithGhosts[i]=pnt_config->int_jStartReal-1;
			   pnt_config->MPI_intKStartRecieve_WithGhosts[i]=pnt_config->int_kStartGhosts;
			   pnt_config->MPI_intKEndRecieve_WithGhosts[i]=pnt_config->int_kEndGhosts;
			   if(MESHDIMENSIONS==2)
			   {
				   pnt_config->MPI_intKStartSend_WithGhosts[i]=0;
				   pnt_config->MPI_intKEndSend_WithGhosts[i]=0;
				   pnt_config->MPI_intKStartRecieve_WithGhosts[i]=0;
				   pnt_config->MPI_intKEndRecieve_WithGhosts[i]=0;

				   pnt_config->MPI_intTransferSizeFlow_WithGhosts[i]=
						   (pnt_config->int_SpaceOrder+1)/2*4*pnt_config->int_iMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells;
			   }
			   else
			   {
				   pnt_config->MPI_intTransferSizeFlow_WithGhosts[i]=
						   (pnt_config->int_SpaceOrder+1)/2*5*pnt_config->int_iMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells;
			   }

			   pnt_config->MPI_intTransferSizeMesh[i]=
					   (pnt_config->int_SpaceOrder+1)/2*13*pnt_config->int_iMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells;
			}
				if(i==pnt_config->InterfaceNeighbourTop)
			{
//				pnt_config->MPI_SendBufferMesh[i]=pnt_config->bufferSendMeshJ;
//				pnt_config->MPI_RecieveBufferMesh[i]=pnt_config->bufferRecieveMeshJ;
//				pnt_config->MPI_SendBufferFlow[i]=pnt_config->bufferSendFlowJ;
//				pnt_config->MPI_RecieveBufferFlow[i]=pnt_config->bufferRecieveFlowJ;
//				pnt_config->MPI_SendBufferViscid[i]=pnt_config->bufferSendViscidJ;
//				pnt_config->MPI_RecieveBufferViscid[i]=pnt_config->bufferRecieveViscidJ;

				pnt_config->MPI_SendBufferMesh[i]=pnt_config->bufferSendMeshTop;
				pnt_config->MPI_RecieveBufferMesh[i]=pnt_config->bufferRecieveMeshTop;
				pnt_config->MPI_SendBufferFlowWithGhosts[i]=pnt_config->bufferSendFlowWithGhostsTop;
				pnt_config->MPI_RecieveBufferFlowWithGhosts[i]=pnt_config->bufferRecieveFlowWithGhostsTop;

			   pnt_config->MPI_intJStartSend[i]=pnt_config->int_jEndReal-((pnt_config->int_SpaceOrder+1)/2-1);
			   pnt_config->MPI_intJEndSend[i]=pnt_config->int_jEndReal;
			   pnt_config->MPI_intIStartSend[i]=pnt_config->int_iStartReal;
			   pnt_config->MPI_intIEndSend[i]=pnt_config->int_iEndReal;
			   pnt_config->MPI_intKStartSend[i]=pnt_config->int_kStartReal;
			   pnt_config->MPI_intKEndSend[i]=pnt_config->int_kEndReal;

			   pnt_config->MPI_intJStartRecieve[i]=pnt_config->int_jEndReal+1;
			   pnt_config->MPI_intJEndRecieve[i]=pnt_config->int_jEndGhosts;
			   pnt_config->MPI_intIStartRecieve[i]=pnt_config->int_iStartReal;
			   pnt_config->MPI_intIEndRecieve[i]=pnt_config->int_iEndReal;
			   pnt_config->MPI_intKStartRecieve[i]=pnt_config->int_kStartReal;
			   pnt_config->MPI_intKEndRecieve[i]=pnt_config->int_kEndReal;


			   pnt_config->MPI_intIStartSend_WithGhosts[i]=pnt_config->int_iStartGhosts;
			   pnt_config->MPI_intIEndSend_WithGhosts[i]=pnt_config->int_iEndGhosts;
			   pnt_config->MPI_intJStartSend_WithGhosts[i]=pnt_config->int_jEndReal-((pnt_config->int_SpaceOrder+1)/2-1);
			   pnt_config->MPI_intJEndSend_WithGhosts[i]=pnt_config->int_jEndReal;
			   pnt_config->MPI_intKStartSend_WithGhosts[i]=pnt_config->int_kStartGhosts;
			   pnt_config->MPI_intKEndSend_WithGhosts[i]=pnt_config->int_kEndGhosts;

			   pnt_config->MPI_intIStartRecieve_WithGhosts[i]=pnt_config->int_iStartGhosts;
			   pnt_config->MPI_intIEndRecieve_WithGhosts[i]=pnt_config->int_iEndGhosts;
			   pnt_config->MPI_intJStartRecieve_WithGhosts[i]=pnt_config->int_jEndReal+1;
			   pnt_config->MPI_intJEndRecieve_WithGhosts[i]=pnt_config->int_jEndGhosts;
			   pnt_config->MPI_intKStartRecieve_WithGhosts[i]=pnt_config->int_kStartGhosts;
			   pnt_config->MPI_intKEndRecieve_WithGhosts[i]=pnt_config->int_kEndGhosts;
			   if(MESHDIMENSIONS==2)
			   {
				   pnt_config->MPI_intKStartSend_WithGhosts[i]=0;
				   pnt_config->MPI_intKEndSend_WithGhosts[i]=0;
				   pnt_config->MPI_intKStartRecieve_WithGhosts[i]=0;
				   pnt_config->MPI_intKEndRecieve_WithGhosts[i]=0;

				   pnt_config->MPI_intTransferSizeFlow_WithGhosts[i]=
						   (pnt_config->int_SpaceOrder+1)/2*4*pnt_config->int_iMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells;
			   }
			   else
			   {
				   pnt_config->MPI_intTransferSizeFlow_WithGhosts[i]=
						   (pnt_config->int_SpaceOrder+1)/2*5*pnt_config->int_iMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells;
			   }

			   pnt_config->MPI_intTransferSizeMesh[i]=
					   (pnt_config->int_SpaceOrder+1)/2*13*pnt_config->int_iMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells;

			}
				if(i==pnt_config->InterfaceNeighbourBehind)
			{
//				pnt_config->MPI_SendBufferMesh[i]=pnt_config->bufferSendMeshK;
//				pnt_config->MPI_RecieveBufferMesh[i]=pnt_config->bufferRecieveMeshK;
//				pnt_config->MPI_SendBufferFlow[i]=pnt_config->bufferSendFlowK;
//				pnt_config->MPI_RecieveBufferFlow[i]=pnt_config->bufferRecieveFlowK;
//				pnt_config->MPI_SendBufferViscid[i]=pnt_config->bufferSendViscidK;
//				pnt_config->MPI_RecieveBufferViscid[i]=pnt_config->bufferRecieveViscidK;

				pnt_config->MPI_SendBufferMesh[i]=pnt_config->bufferSendMeshBehind;
				pnt_config->MPI_RecieveBufferMesh[i]=pnt_config->bufferRecieveMeshBehind;
				pnt_config->MPI_SendBufferFlowWithGhosts[i]=pnt_config->bufferSendFlowWithGhostsBehind;
				pnt_config->MPI_RecieveBufferFlowWithGhosts[i]=pnt_config->bufferRecieveFlowWithGhostsBehind;

			   pnt_config->MPI_intKStartSend[i]=pnt_config->int_kStartReal;
			   pnt_config->MPI_intKEndSend[i]=pnt_config->int_kStartReal+((pnt_config->int_SpaceOrder+1)/2-1);
			   pnt_config->MPI_intIStartSend[i]=pnt_config->int_iStartReal;
			   pnt_config->MPI_intIEndSend[i]=pnt_config->int_iEndReal;
			   pnt_config->MPI_intJStartSend[i]=pnt_config->int_jStartReal;
			   pnt_config->MPI_intJEndSend[i]=pnt_config->int_jEndReal;

			   pnt_config->MPI_intKStartRecieve[i]=pnt_config->int_kStartGhosts;
			   pnt_config->MPI_intKEndRecieve[i]=pnt_config->int_kStartReal-1;
			   pnt_config->MPI_intIStartRecieve[i]=pnt_config->int_iStartReal;
			   pnt_config->MPI_intIEndRecieve[i]=pnt_config->int_iEndReal;
			   pnt_config->MPI_intJStartRecieve[i]=pnt_config->int_jStartReal;
			   pnt_config->MPI_intJEndRecieve[i]=pnt_config->int_jEndReal;


			   pnt_config->MPI_intIStartSend_WithGhosts[i]=pnt_config->int_iStartGhosts;
			   pnt_config->MPI_intIEndSend_WithGhosts[i]=pnt_config->int_iEndGhosts;
			   pnt_config->MPI_intJStartSend_WithGhosts[i]=pnt_config->int_jStartGhosts;
			   pnt_config->MPI_intJEndSend_WithGhosts[i]=pnt_config->int_jEndGhosts;
			   pnt_config->MPI_intKStartSend_WithGhosts[i]=pnt_config->int_kStartReal;
			   pnt_config->MPI_intKEndSend_WithGhosts[i]=pnt_config->int_kStartReal+((pnt_config->int_SpaceOrder+1)/2-1);

			   pnt_config->MPI_intIStartRecieve_WithGhosts[i]=pnt_config->int_iStartGhosts;
			   pnt_config->MPI_intIEndRecieve_WithGhosts[i]=pnt_config->int_iEndGhosts;
			   pnt_config->MPI_intJStartRecieve_WithGhosts[i]=pnt_config->int_jStartGhosts;
			   pnt_config->MPI_intJEndRecieve_WithGhosts[i]=pnt_config->int_jEndGhosts;
			   pnt_config->MPI_intKStartRecieve_WithGhosts[i]=pnt_config->int_kStartGhosts;
			   pnt_config->MPI_intKEndRecieve_WithGhosts[i]=pnt_config->int_kStartReal-1;

			   pnt_config->MPI_intTransferSizeMesh[i]=
					   (pnt_config->int_SpaceOrder+1)/2*13*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_iMeshPointsGhostCells;
			   pnt_config->MPI_intTransferSizeFlow_WithGhosts[i]=
					   (pnt_config->int_SpaceOrder+1)/2*5*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_iMeshPointsGhostCells;

			}
				if(i==pnt_config->InterfaceNeighbourInFront)
			{
//				pnt_config->MPI_SendBufferMesh[i]=pnt_config->bufferSendMeshK;
//				pnt_config->MPI_RecieveBufferMesh[i]=pnt_config->bufferRecieveMeshK;
//				pnt_config->MPI_SendBufferFlow[i]=pnt_config->bufferSendFlowK;
//				pnt_config->MPI_RecieveBufferFlow[i]=pnt_config->bufferRecieveFlowK;
//				pnt_config->MPI_SendBufferViscid[i]=pnt_config->bufferSendViscidK;
//				pnt_config->MPI_RecieveBufferViscid[i]=pnt_config->bufferRecieveViscidK;

				pnt_config->MPI_SendBufferMesh[i]=pnt_config->bufferSendMeshInFront;
				pnt_config->MPI_RecieveBufferMesh[i]=pnt_config->bufferRecieveMeshInFront;
				pnt_config->MPI_SendBufferFlowWithGhosts[i]=pnt_config->bufferSendFlowWithGhostsInFront;
				pnt_config->MPI_RecieveBufferFlowWithGhosts[i]=pnt_config->bufferRecieveFlowWithGhostsInFront;

			   pnt_config->MPI_intKStartSend[i]=pnt_config->int_kEndReal-((pnt_config->int_SpaceOrder+1)/2-1);
			   pnt_config->MPI_intKEndSend[i]=pnt_config->int_kEndReal;
			   pnt_config->MPI_intIStartSend[i]=pnt_config->int_iStartReal;
			   pnt_config->MPI_intIEndSend[i]=pnt_config->int_iEndReal;
			   pnt_config->MPI_intJStartSend[i]=pnt_config->int_jStartReal;
			   pnt_config->MPI_intJEndSend[i]=pnt_config->int_jEndReal;

			   pnt_config->MPI_intKStartRecieve[i]=pnt_config->int_kEndReal+1;
			   pnt_config->MPI_intKEndRecieve[i]=pnt_config->int_kEndGhosts;
			   pnt_config->MPI_intIStartRecieve[i]=pnt_config->int_iStartReal;
			   pnt_config->MPI_intIEndRecieve[i]=pnt_config->int_iEndReal;
			   pnt_config->MPI_intJStartRecieve[i]=pnt_config->int_jStartReal;
			   pnt_config->MPI_intJEndRecieve[i]=pnt_config->int_jEndReal;


			   pnt_config->MPI_intIStartSend_WithGhosts[i]=pnt_config->int_iStartGhosts;
			   pnt_config->MPI_intIEndSend_WithGhosts[i]=pnt_config->int_iEndGhosts;
			   pnt_config->MPI_intJStartSend_WithGhosts[i]=pnt_config->int_jStartGhosts;
			   pnt_config->MPI_intJEndSend_WithGhosts[i]=pnt_config->int_jEndGhosts;
			   pnt_config->MPI_intKStartSend_WithGhosts[i]=pnt_config->int_kEndReal-((pnt_config->int_SpaceOrder+1)/2-1);
			   pnt_config->MPI_intKEndSend_WithGhosts[i]=pnt_config->int_kEndReal;

			   pnt_config->MPI_intIStartRecieve_WithGhosts[i]=pnt_config->int_iStartGhosts;
			   pnt_config->MPI_intIEndRecieve_WithGhosts[i]=pnt_config->int_iEndGhosts;
			   pnt_config->MPI_intJStartRecieve_WithGhosts[i]=pnt_config->int_jStartGhosts;
			   pnt_config->MPI_intJEndRecieve_WithGhosts[i]=pnt_config->int_jEndGhosts;
			   pnt_config->MPI_intKStartRecieve_WithGhosts[i]=pnt_config->int_kEndReal+1;
			   pnt_config->MPI_intKEndRecieve_WithGhosts[i]=pnt_config->int_kEndGhosts;

			   pnt_config->MPI_intTransferSizeMesh[i]=
					   (pnt_config->int_SpaceOrder+1)/2*13*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_iMeshPointsGhostCells;
			   pnt_config->MPI_intTransferSizeFlow_WithGhosts[i]=
					   (pnt_config->int_SpaceOrder+1)/2*5*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_iMeshPointsGhostCells;

			}
				check_TransformationMatrix(
					i,
					pnt_config);
	}


//	printf("saiWENOS: Rank: %d (%s)",pnt_config->MPI_rank,pnt_config->Zonename);
//	printf(": Left: %d",pnt_config->MPI_rankNeighbours[pnt_config->InterfaceNeighbourLeft]);
//	printf(": Right: %d ",pnt_config->MPI_rankNeighbours[pnt_config->InterfaceNeighbourRight]);
//	printf(": Bottom: %d",pnt_config->MPI_rankNeighbours[pnt_config->InterfaceNeighbourBottom]);
//	printf(": Top: %d",pnt_config->MPI_rankNeighbours[pnt_config->InterfaceNeighbourTop]);
//	printf(": Behind: %d ",pnt_config->MPI_rankNeighbours[pnt_config->InterfaceNeighbourBehind]);
//	printf(": InFront: %d\n",pnt_config->MPI_rankNeighbours[pnt_config->InterfaceNeighbourInFront]);
//
//	for (i=0;i<pnt_config->NumberInterfaces;i++)
//	{
//		printf("myrank: %d | rank: %d | interface:%d | transform: %d %d %d | donorname: %s\n",
//				pnt_config->MPI_rank,
//				 pnt_config->MPI_rankNeighbours[i],i,
//				pnt_config->TransformMatrixOfInterface[i][0],
//				pnt_config->TransformMatrixOfInterface[i][1],
//				pnt_config->TransformMatrixOfInterface[i][2],
//				pnt_config->Donorname[i]);
//
//	}


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

    if((pnt_config->flag_SplitMeshFile==2)||(pnt_config->int_initializeType==1))
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

	if (pnt_config->MPI_rank==0){printf("saiWENOS: Folgende Randbedingungen wurden gesetzt: ");}

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


