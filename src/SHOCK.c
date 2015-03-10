#define _ISOC99_SOURCE  1

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include "string.h"
#include <math.h>
#include "Load.h"
#include "Save.h"

#include "SHOCK.h"
#include "Functions.h"
#include "Import.h"
#include "Export.h"
#include "WENO.h"
#include "ZD.h"
#include "BC.h"
#include "iniparser.h"
#include <hdf5.h>
#include "ManufacturedSolution.h"

int int_interationsStart;
int int_interationsEnd;

double t0,t1,t2;

int int_myCPUID;
int actualID,actualP;
int int_MaxNumberCPUs;

int i,j,k,ijk;

FLT leftEigenvector[5][5];
FLT rightEigenvector[5][5];


int int_helpValue1;
FLT helpValue2;

int int_iterationCounterStdOut;
int int_iterationCounterBackupOut;
int int_iterationCounterSamples;

int int_ResetCounter;

int main(int argc, char *argv[])
{
	MPI_Init(
			&argc,
			&argv);


	//Import of Configuration
	struct strct_configuration configuration;
	//Erzeugung der Strukturen für Gitter(strct_mesh), Erhaltungsgrößen(strct_U) und Flüssen(strct_Flux)
	struct strct_mesh mesh;
	//Erhaltungsgrößen des letzten Schrittes
	struct strct_U U_lastStep;
	//Erhaltungsgrößen innerhalb der RK-Schritte
	struct strct_U U_RK;
	//Flussvektor aller Flüsse
	struct strct_Flux Q;

	struct strct_Flux Flux;
	struct strct_Flux Flux_PlusHalf;
	//Flussvektor aller Flüsse über mehrere RK-Schritte
	struct strct_Flux Q_sum;

	//	Filmstruktur (U(u,v,w,rho,p) von mehreren Zeitschritten)
	struct strct_Film Film;

	//Backupspeicher für Erhaltungsgrößen
	struct strct_U U_backup1;
	struct strct_U U_backup2;


	configuration.MPI_comm=MPI_COMM_WORLD;

	MPI_Barrier(configuration.MPI_comm);
	t0 = MPI_Wtime();

	//Zuweisung Koordinate <-> MyID//
	MPI_Comm_size(
			configuration.MPI_comm,
			&configuration.MPI_size);
	MPI_Comm_rank(
			configuration.MPI_comm,
			&configuration.MPI_rank);


	if(configuration.MPI_rank==0){printf("SHOCK: ####################################\n");}
	if(configuration.MPI_rank==0){printf("SHOCK:                Start\n");}
	if(configuration.MPI_rank==0){printf("SHOCK: (git-ID: %s)\n",GITID);}
	if(configuration.MPI_rank==0){printf("SHOCK: (Precision: %s (sizeof: float=%lu, double=%lu, long double=%lu, __float128=%lu)\n",FLT_name,sizeof(float),sizeof(double),sizeof(long double),sizeof(__float128));}
	if(configuration.MPI_rank==0){printf("SHOCK: (MPI_Version: %d.%d.%d)\n",OMPI_MAJOR_VERSION,OMPI_MINOR_VERSION,OMPI_RELEASE_VERSION);}
	if(configuration.MPI_rank==0){printf("SHOCK: (CGNS_Version: %d)\n",CGNS_VERSION);}
	if(configuration.MPI_rank==0){printf("SHOCK: (HDF5_Version: %d.%d.%d)\n",H5_VERS_MAJOR,H5_VERS_MINOR,H5_VERS_RELEASE);}
	if(configuration.MPI_rank==0){printf("SHOCK: ####################################\n");}
	if(argv[1]!=NULL)
	{
		strcpy(configuration.chr_configPath,argv[1]);
	}

	if(configuration.MPI_rank==0){printf("SHOCK: Einstellungen werden geladen aus: %s \n",configuration.chr_configPath);}

	if(check_ConfigFile(&configuration)){
		ConfigImport(
				&configuration);
	}
	else
	{
		goto exit;
	}

	if(check_CGNSFile(&configuration)==0)
	{
		goto exit;
	}

	int_iterationCounterSamples=0;
	int_iterationCounterStdOut=0;
	int_iterationCounterBackupOut=0;
	configuration.int_actualSample=0;

	int_ResetCounter=0;

	//Hier werden einige Variablen gesetzt, die später gebraucht werden
	DefineParameters(
			&configuration);
	if(configuration.MPI_rank==0){printf("SHOCK: Definition der Variablen fertig!\n");}
	double ts;

	ts = MPI_Wtime( );

	if(configuration.MPI_rank==0){printf("SHOCK: Parallel Loading of CGNS-file...!\n");}

	vta x = { NULL },y = { NULL },z = { NULL },u = { NULL },v = { NULL },w = { NULL },rho = { NULL },p = { NULL };

	loadFile( &configuration,&x,&y,&z,&u,&v,&w,&rho,&p );

	postprocessLoad(
		&x,&y,&z,&u,&v,&w,&rho,&p,
		&configuration,
		&mesh,
		&U_lastStep,
		&U_RK,
		&Flux,
		&Flux_PlusHalf,
		&Q,
		&Q_sum,
		&Film,
		&U_backup1,
		&U_backup2);

	freeVTA(&x,&y,&z,&u,&v,&w,&rho,&p);

	MPI_Barrier( MPI_COMM_WORLD );

	if(configuration.MPI_rank==0){printf("SHOCK: Loading %s fertig (%f min.)!\n",configuration.chr_MeshFile,(MPI_Wtime( )-ts)/60. );}

//	printf("rank: %d left %d right %d top %d bottom %d behind %d inFront %d \n",
//			configuration.MPI_rank,
//			configuration.InterfaceNeighbourLeft,
//			configuration.InterfaceNeighbourRight,
//			configuration.InterfaceNeighbourTop,
//			configuration.InterfaceNeighbourBottom,
//			configuration.InterfaceNeighbourBehind,
//			configuration.InterfaceNeighbourInFront);
	ExtrapolateGhostCells(
			&configuration,
			&mesh);
	if(configuration.MPI_rank==0){printf("SHOCK: GhostCells-Extrapolation fertig!\n");}


	//Berechnung der Metrik
//	Transfer der Gitterpunkte(Metrik ist überall noch 0)
	TransferMeshParameter(
			&configuration,
			&mesh);
	TransferMeshParameter(
			&configuration,
			&mesh);
	if (MESHDIMENSIONS==3)
	{
		TransferMeshParameter(
				&configuration,
				&mesh);
	}

//	Berechnung der Metrik
	CreateMetric(
			&configuration,
			&mesh);

//	Transfer damit die auch in den GhostCells verfügbar ist
	TransferMeshParameter(
			&configuration,
			&mesh);
	TransferMeshParameter(
			&configuration,
			&mesh);
	if (MESHDIMENSIONS==3)
	{
		TransferMeshParameter(
				&configuration,
				&mesh);
	}

	//	Entgültige Überprüfung der Metrik
	check_Metric(
			&configuration,
			&mesh);

//	An Wänden leglicher Art wird die Metrik gespiegelt. Bei entsprechender Einstellung der Randbedingung
//	fuer theta und u fuehrt dies zu antimetrischen fluessen fuer masse und energie
//	und symmetrischen fluessen fuer die impulse. die impulsfluesse sind symmetrisch,
//	da dort stets zwei geschwindigkeiten multipliziert werden (minuszeichen verschwindet)
//	MitHilfe der BC_Corrector-Variablen wird dieses Problem beseitigt
	mirrorGhostCellsMetric(
			&configuration,
			&mesh);

//	Berechnung der Hilfs-Metrik fuer viskose Fluesse
	CreateViscidMetric(
			&configuration,
			&mesh);

	if(configuration.MPI_rank==0){printf("SHOCK: Metrik-Erstellung fertig!\n");}

	if(configuration.flag_IBC==1)
	{
		IBC_prepare(
				&configuration,
				&mesh);
		if(configuration.MPI_rank==0){printf(">>>>> Options: Immerged BC fertig\n");}
	}

	SetAllBoundaryConditions(
			&configuration,
			&mesh,
			&U_lastStep,
			&U_lastStep);


	//	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	//	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	//	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> OPTIONS<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	//	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	//	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	if(configuration.MPI_rank==0){printf("SHOCK: Moegliche Optionen werden gesetzt...\n");}
	setOptions(
			&configuration,
			&mesh,
			&U_lastStep,
			&U_RK,
			&Q,
			&Q_sum,
			&Flux,
			&Flux_PlusHalf,
			&Film);
	MPI_Barrier(configuration.MPI_comm);
	if(configuration.MPI_rank==0){printf("SHOCK: Optionen sind gesetzt!\n");}



	if(configuration.MPI_rank==0){printf("\n");}
	if(configuration.MPI_rank==0){printf("SHOCK: ####################################\n");}
	if(configuration.MPI_rank==0){printf("SHOCK: Mesh: '%s'\n",configuration.chr_MeshPath);}
	if(configuration.MPI_rank==0){printf("SHOCK: %d Prozesse werden %d Iterationen rechnen!\n",configuration.MPI_size,configuration.int_TotalIterations);}
	if(configuration.MPI_rank==0){printf("SHOCK: Angle of Attack: %g [Degree] | Ma: %g | Re: %g | Pr: %g \n",
			(double)configuration.AoA,(double)configuration.machNumber,(double)configuration.reynoldsNumber*abs(configuration.flag_Inviscid-1),(double)configuration.prandtlNumber);}
	if(configuration.MPI_rank==0){printf("SHOCK: SpaceOrder: %d | TimeOrder: %d | numerisches delta t (Tau): %g | Epsilon: %g\n",
				configuration.int_SpaceOrder,configuration.int_TimeOrder,(double)configuration.numericalTau,(double)configuration.wenoEpsilon);}
	if(configuration.MPI_rank==0){printf("SHOCK: ####################################\n");}
	if(configuration.MPI_rank==0){printf("\n");}

	t1 = MPI_Wtime();
	if(configuration.MPI_rank==0){print_memusage_c();}
	if(configuration.MPI_rank==0){printf("SHOCK: Die Initialisierung fertig (%f min.)\n",(t1-t0)/60.);}

	//########################################
	//########################################
	//###########SIMULATION STARTEN###########
	//########################################
	//########################################
	//########################################

	//	Schreibe Daten in den Tmp_U-Array U_RK
	WriteValuesFromU1ToU2(
			&configuration,
			&U_lastStep,
			&U_RK);

	if(configuration.MPI_rank==0){printf("SHOCK: %dD Simulation wird gestartet...\n",configuration.int_meshDimensions);}
	startSimulation(
			&configuration,
			&mesh,
			&U_lastStep,
			&U_RK,
			&Q,
			&Q_sum,
			&Flux,
			&Flux_PlusHalf,
			&Film,
			&U_backup1,
			&U_backup2);
	MPI_Barrier(configuration.MPI_comm);
	if(configuration.MPI_rank==0){printf("SHOCK: %dD Simulation fertig (%g min.)!\n",
			configuration.int_meshDimensions,
			((MPI_Wtime())-t1)/60.);}
	if(configuration.MPI_rank==0){printf("SHOCK: Gesamtkommunikationsdauer: %g min.\n",
			configuration.comm_time/60.);}
	if(configuration.MPI_rank==0){printf("SHOCK: Letztes Tau: %g!\n",
			(double)configuration.numericalTau);}
	if(configuration.MPI_rank==0){printf("SHOCK: Letztes DistanceForwad: %d!\n",
			configuration.int_distanceForward);}

	if(configuration.flag_ManufacturedSolution==1)
	{
		ErrorManufacturedSolution(
				&configuration,
				&mesh,
				&U_lastStep,
				1);
	}

	if(configuration.int_TotalIterations==0)
	{
		if(configuration.MPI_rank==0){printf("SHOCK: Initialisierungswerte für Export speichern.\n");}

		CalcValuesForPost(
			&configuration,
			&mesh,
			&U_lastStep);

		WriteValuesFromUToFilm(
			&configuration,
			&U_lastStep,
			&Film,
			&mesh);
	}


	if(configuration.flag_constantZValues==1)
	{
		WriteConstantZValues(
				&configuration,
				&U_lastStep);
	}


//	Berechne Größen, die für den Export relevant sind
	//		Strömungsparameter (rho,u,v,w) werden übertragen bzw. Randbedinungen gesetzt
	TransferFlowParameterWithGhosts(
			&configuration,
			&mesh,
			&U_lastStep);
	TransferFlowParameterWithGhosts(
			&configuration,
			&mesh,
			&U_lastStep);
	TransferFlowParameterWithGhosts(
			&configuration,
			&mesh,
			&U_lastStep);
	CalcValuesForPost(
			&configuration,
			&mesh,
			&U_lastStep);

	if(configuration.flag_constantZValues==1)
	{
		WriteConstantZValues(
				&configuration,
				&U_lastStep);
	}

	//	Gitter wird herausgeschrieben
	if(configuration.flag_exportMetric==1)
	{
		MeshMetricExport(
				&configuration,
				&mesh,
				&U_lastStep);
		if(configuration.MPI_rank==0){printf("SHOCK: Metric-Export fertig! \n");}
	}
	//	Ergebnisse werden herausgeschrieben
	ts = MPI_Wtime( );
	saveFile( &configuration,&Film );

	if(configuration.MPI_rank==0){printf("SHOCK: CGNS-Export (%d Samples) fertig (%f min.)!\n",configuration.int_Samples,((MPI_Wtime( ))-ts)/60.);}


	if(configuration.flag_PressureHistory==1)
	{
		ts = MPI_Wtime( );
		for(actualP=0;actualP<configuration.PressureHistory_No;actualP++)
		{
//			for(actualID=0;actualID<configuration.MPI_size;actualID++)
//			{
//				if((configuration.MPI_rank==actualID)&&(configuration.flag_PressureHistory_P[actualP]==1))
//				{
//					CGNS_PressureHistoryValuesExportParallel(
//							&configuration,
//							actualP);
//				}
//				MPI_Barrier(configuration.MPI_comm);
//			}
			if(configuration.flag_PressureHistory_P[actualP]==1)
			{
				ASCii_PressureHistoryValuesExportParallel(
							&configuration,
							actualP);
			}
			MPI_Barrier(configuration.MPI_comm);
		}
		if(configuration.MPI_rank==0){printf("SHOCK: PressureHistory-Export fertig (%f min.)!\n",((MPI_Wtime( ))-ts)/60.);}
	}

	if(configuration.flag_VelocityHistory==1)
	{
		ts = MPI_Wtime( );
		for(actualP=0;actualP<configuration.VelocityHistory_No;actualP++)
		{
//			for(actualID=0;actualID<configuration.MPI_size;actualID++)
//			{
//				if((configuration.MPI_rank==actualID)&&(configuration.flag_VelocityHistory_P[actualP]==1))
//				{
//					CGNS_VelocityHistoryValuesExportParallel(
//							&configuration,
//							actualP);
//				}
//				MPI_Barrier(configuration.MPI_comm);
//			}
			if(configuration.flag_VelocityHistory_P[actualP]==1)
			{
				ASCii_VelocityHistoryValuesExportParallel(
							&configuration,
							actualP);
			}
			MPI_Barrier(configuration.MPI_comm);
		}
		if(configuration.MPI_rank==0){printf("SHOCK: VelocityHistory-Export fertig (%f min.)!\n",((MPI_Wtime( ))-ts)/60.);}
	}



	ts = MPI_Wtime( );
	//Allokierten Speicher freigeben
	FreeMemory(
			&configuration,
			&mesh,
			&U_lastStep,
			&U_RK,
			&Flux,
			&Flux_PlusHalf,
			&Q,
			&Q_sum,
			&Film);

if(configuration.MPI_rank==0){printf("SHOCK: Speicher freigegeben (%f min.)!\n",((MPI_Wtime( ))-ts)/60.);}

//	if(configuration.MPI_rank==0){printf("Glaettungsindikatoren: Durchschnitt: %le   Minimum: %le  Maximum: %le,\n",
//			(configuration.is_avrg/configuration.is_avrg_counter),configuration.is_minimum,configuration.is_maximum);}

	exit:
	if(configuration.MPI_rank==0){printf("SHOCK: Exit!\n");}
	MPI_Finalize();
	return 0;
}














//########################################
//########################################
//###########  2D  SIMULATION  ###########
//########################################
//########################################
//########################################

void startSimulation(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep,
		struct strct_U * pnt_U_RK,
		struct strct_Flux * pnt_Q,
		struct strct_Flux * pnt_Q_sum,
		struct strct_Flux * pnt_Flux,
		struct strct_Flux * pnt_Flux_PlusHalf,
		struct strct_Film * pnt_Film,
		struct strct_U * pnt_U_backup1,
		struct strct_U * pnt_U_backup2)
{
// Nur noetig fur system_command, der aber auf juqueen nicht funktioniert
//  char chr_SnapshotPath_normal[500];
//  char chr_SnapshotPath_new[500];
//  char command[500];

////		printf("SHOCK: Rank: %d (%s)",pnt_config->MPI_rank,pnt_config->Zonename);
////		printf(": Left: %d",pnt_config->MPI_rankNeighbours[pnt_config->InterfaceNeighbourLeft]);
////		printf(": Right: %d ",pnt_config->MPI_rankNeighbours[pnt_config->InterfaceNeighbourRight]);
////		printf(": Bottom: %d",pnt_config->MPI_rankNeighbours[pnt_config->InterfaceNeighbourBottom]);
////		printf(": Top: %d",pnt_config->MPI_rankNeighbours[pnt_config->InterfaceNeighbourTop]);
////		printf(": Behind: %d ",pnt_config->MPI_rankNeighbours[pnt_config->InterfaceNeighbourBehind]);
////		printf(": InFront: %d\n",pnt_config->MPI_rankNeighbours[pnt_config->InterfaceNeighbourInFront]);
//
//		for (i=0;i<pnt_config->NumberInterfaces;i++)
//		{
//			printf("myrank: %d | rank: %d | interface:%d | transform: %d %d | range: %d %d | donorrange: %d %d \n",
//					pnt_config->MPI_rank,
//					 pnt_config->MPI_rankNeighbours[i],i,
//					pnt_config->TransformMatrixOfInterface[i][0],
//					pnt_config->TransformMatrixOfInterface[i][1],
//					(int)pnt_config->RangeOfInterface[i][0],
//					(int)pnt_config->RangeOfInterface[i][1],
//					(int)pnt_config->DonorRangeOfInterface[i][0],
//					(int)pnt_config->DonorRangeOfInterface[i][1]);
//
//		}
//		if(pnt_config->MPI_rank==0)
//		{
//			printf("zone: %s | left: %s | right: %s | top: %s | bottom: %s | behind: %s | infront: %s | \n",
//					pnt_config->Zonename,
//					pnt_config->BC_Left,
//					pnt_config->BC_Right,
//					pnt_config->BC_Top,
//					pnt_config->BC_Bottom,
//					pnt_config->BC_Behind,tec
//					pnt_config->BC_InFront);
//		}


	WriteValuesFromU1ToU2(
		pnt_config,
		pnt_U_lastStep,
		pnt_U_backup1);
	pnt_config->time_dim_backup1=pnt_config->time_dim;
	pnt_config->int_actualIteration_backup1=1;
	pnt_config->int_actualIteration_backup2=999999;
	pnt_config->int_actualSample_backup1=pnt_config->int_actualSample;
	pnt_config->int_iterationCounterSamples_backup1=int_iterationCounterSamples;
	pnt_config->int_iterationCounterBackupOut_backup1=int_iterationCounterBackupOut;
	pnt_config->int_iterationCounterStdOut_backup1=int_iterationCounterStdOut;

	for(pnt_config->int_actualIteration=pnt_config->int_StartIteration+1;
			pnt_config->int_actualIteration<=pnt_config->int_EndIteration;
			pnt_config->int_actualIteration++)
	{
		pnt_config->time_dim+=
				pnt_config->numericalTau
				*pnt_config->L0_dim
				/pnt_config->u0_dim;
		CalcRungeKutta(
				pnt_config,
				pnt_mesh,
				pnt_U_lastStep,
				pnt_U_RK,
				pnt_Q,
				pnt_Q_sum,
				pnt_Flux,
				pnt_Flux_PlusHalf);


//		Nach jedem Iterationsschritt wird das Gebiet nach NaN abgesucht
		if(checkNAN(pnt_config,pnt_U_RK,pnt_mesh)==1)
		{
			if(int_ResetCounter>=pnt_config->int_NumberResets)
			{
				pnt_config->int_actualIteration=pnt_config->int_EndIteration;
				if(pnt_config->MPI_rank==0){printf("SHOCK: Maximal Number of Resets reached. Aborting Simulation!\n");}
				if(pnt_config->flag_NAN==1){NANExport( pnt_config,pnt_mesh,pnt_U_RK);}
//				MPI_Abort(pnt_config->MPI_comm,1337);
				break;
			}
			int_ResetCounter++;

			pnt_config->numericalTau*=pnt_config->TauDecelerator_factor;
			pnt_config->int_distanceForward=pnt_config->int_distanceForwardStart;

			if(pnt_config->int_actualIteration_backup1<pnt_config->int_actualIteration_backup2)
			{
				if(pnt_config->MPI_rank==0){printf("SHOCK: ----> Backup von Iterationsschritt %d wird geladen (Neues Tau = %g).\n",
						pnt_config->int_actualIteration_backup1,
						(double)pnt_config->numericalTau);}
				WriteValuesFromU1ToU2(
						pnt_config,
						pnt_U_backup1,
						pnt_U_RK);
				pnt_config->time_dim_lastAction=pnt_config->time_dim;
				pnt_config->time_dim=pnt_config->time_dim_backup1;
				pnt_config->int_actualIteration=pnt_config->int_actualIteration_backup1;
				pnt_config->int_actualSample=pnt_config->int_actualSample_backup1;
				int_iterationCounterSamples=pnt_config->int_iterationCounterSamples_backup1;
				int_iterationCounterBackupOut=pnt_config->int_iterationCounterBackupOut_backup1;
				int_iterationCounterStdOut=pnt_config->int_iterationCounterStdOut_backup1;
			}
			else
			{
				if(pnt_config->MPI_rank==0){printf("SHOCK: ----> Backup von Iterationsschritt %d wird geladen (Neues Tau = %g).\n",
						pnt_config->int_actualIteration_backup2,
						(double)pnt_config->numericalTau);}
				WriteValuesFromU1ToU2(
						pnt_config,
						pnt_U_backup2,
						pnt_U_RK);
				pnt_config->time_dim_lastAction=pnt_config->time_dim;
				pnt_config->time_dim=pnt_config->time_dim_backup2;
				pnt_config->int_actualIteration=pnt_config->int_actualIteration_backup2;
				pnt_config->int_actualSample=pnt_config->int_actualSample_backup2;
				int_iterationCounterSamples=pnt_config->int_iterationCounterSamples_backup2;
				int_iterationCounterBackupOut=pnt_config->int_iterationCounterBackupOut_backup2;
				int_iterationCounterStdOut=pnt_config->int_iterationCounterStdOut_backup2;
			}
		}
//		Tau-Resetter
		else if(
			((pnt_config->time_dim-pnt_config->time_dim_lastAction)>
		pnt_config->int_distanceNAN*pnt_config->numericalTau*pnt_config->L0_dim/pnt_config->u0_dim)
			&&(pnt_config->flag_TauAccelerator!=0)
			&&(pnt_config->int_actualIteration+pnt_config->int_distanceBackward<pnt_config->int_EndIteration)
			&&(pnt_config->numericalTau<pnt_config->numericalTauStart)
			)
		{
			pnt_config->int_distanceForward=pnt_config->int_distanceForwardStart;
			pnt_config->numericalTau=pnt_config->numericalTauStart;
			pnt_config->time_dim_lastAction=pnt_config->time_dim;
			if(pnt_config->MPI_rank==0){printf("SHOCK: ----> Tau-Resetting bei Iterationsschritt %d (Neues Tau = %g)\n",
					pnt_config->int_actualIteration,
					(double)pnt_config->numericalTau);}
		}
//		Tau-Beschleuniger
		else if(
			((pnt_config->time_dim-pnt_config->time_dim_lastAction)>
		pnt_config->int_distanceForward*pnt_config->numericalTau*pnt_config->L0_dim/pnt_config->u0_dim)
			&&(pnt_config->flag_TauAccelerator!=0)
			&&(pnt_config->int_actualIteration+pnt_config->int_distanceBackward<pnt_config->int_EndIteration)
			&&(pnt_config->numericalTau>=pnt_config->numericalTauStart)
			)
		{
			pnt_config->int_distanceForward*=4.;
			pnt_config->numericalTau*=pnt_config->TauAccelerator_factor;
			pnt_config->time_dim_lastAction=pnt_config->time_dim;
			if(pnt_config->MPI_rank==0){printf("SHOCK: Tau-Erhoehung bei Iterationsschritt %d (Tau = %g)\n",
					pnt_config->int_actualIteration,
					(double)pnt_config->numericalTau);}
		}

		//	In U_RK sind nach dem vierten RK-Schritt die aktuellen Ergebnisse des Zeitschrittes n+1
		//	Diese werden nun übertragen
		WriteValuesFromU1ToU2(
				pnt_config,
				pnt_U_RK,
				pnt_U_lastStep);


		if((pnt_config->int_actualIteration)%pnt_config->int_distanceBackward==0)
		{
			WriteValuesFromU1ToU2(
				pnt_config,
				pnt_U_lastStep,
				pnt_U_backup1);
			pnt_config->time_dim_backup1=pnt_config->time_dim;
			pnt_config->int_actualIteration_backup1=pnt_config->int_actualIteration;
			pnt_config->int_actualSample_backup1=pnt_config->int_actualSample;
			pnt_config->int_iterationCounterSamples_backup1=int_iterationCounterSamples;
			pnt_config->int_iterationCounterBackupOut_backup1=int_iterationCounterBackupOut;
			pnt_config->int_iterationCounterStdOut_backup1=int_iterationCounterStdOut;
		}else if((pnt_config->int_actualIteration+pnt_config->int_distanceBackward/2)%pnt_config->int_distanceBackward==0)
		{
			WriteValuesFromU1ToU2(
				pnt_config,
				pnt_U_lastStep,
				pnt_U_backup2);
			pnt_config->time_dim_backup2=pnt_config->time_dim;
			pnt_config->int_actualIteration_backup2=pnt_config->int_actualIteration;
			pnt_config->int_actualSample_backup2=pnt_config->int_actualSample;
			pnt_config->int_iterationCounterSamples_backup2=int_iterationCounterSamples;
			pnt_config->int_iterationCounterBackupOut_backup2=int_iterationCounterBackupOut;
			pnt_config->int_iterationCounterStdOut_backup2=int_iterationCounterStdOut;
		}




		int_iterationCounterSamples++;
		if((int_iterationCounterSamples%pnt_config->int_IterationsBetweenSamples==0)&&(pnt_config->int_actualSample<pnt_config->int_Samples)&&(pnt_config->int_actualIteration>=(pnt_config->int_StartSampling)))
		{
			if((pnt_config->MPI_rank==0)&&(pnt_config->int_actualSample==0)){printf("SHOCK: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");}
			if((pnt_config->MPI_rank==0)&&(pnt_config->int_actualSample==0)){printf("SHOCK:           Start Sampling.\n");}
			if((pnt_config->MPI_rank==0)&&(pnt_config->int_actualSample==0)){printf("SHOCK: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");}

			CalcValuesForPost(
					 pnt_config,
					 pnt_mesh,
					 pnt_U_lastStep);

			WriteValuesFromUToFilm(
				 pnt_config,
				 pnt_U_lastStep,
				 pnt_Film,
				 pnt_mesh);

			pnt_config->int_actualSample++;
		}


		if(pnt_config->flag_PressureHistory==1)
		{
			WriteValuesForPressureHistory(
				pnt_config,
				pnt_U_lastStep);
		}

		if(pnt_config->flag_VelocityHistory==1)
		{
			WriteValuesForVelocityHistory(
				pnt_config,
				pnt_U_lastStep);
		}

		int_iterationCounterStdOut++;
		if(int_iterationCounterStdOut==(pnt_config->int_TotalIterations/100)||(pnt_config->int_EndIteration<100))
		{
			int_iterationCounterStdOut=0;
			if(pnt_config->MPI_rank==0){printf("SHOCK: Iteration: %d von %d (Zeit: %g [sec])\n",pnt_config->int_actualIteration,pnt_config->int_EndIteration,(double)pnt_config->time_dim);}
			if(pnt_config->MPI_rank==0){printf("\t rho[x=%g, y=%g, z=%g] = %g \n",
					(double)pnt_mesh->x[pnt_config->int_iMid*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+pnt_config->int_jMid*pnt_config->int_kMeshPointsGhostCells+pnt_config->int_kMid],
					(double)pnt_mesh->y[pnt_config->int_iMid*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+pnt_config->int_jMid*pnt_config->int_kMeshPointsGhostCells+pnt_config->int_kMid],
					(double)pnt_mesh->z[pnt_config->int_iMid*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+pnt_config->int_jMid*pnt_config->int_kMeshPointsGhostCells+pnt_config->int_kMid],
					(double)pnt_U_lastStep->rho[pnt_config->int_iMid*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+pnt_config->int_jMid*pnt_config->int_kMeshPointsGhostCells+pnt_config->int_kMid]);}
			if((pnt_config->flag_IBC_Moving==1)&&(pnt_config->flag_IBC==1))
			{
				if(pnt_config->MPI_rank==0){printf("SHOCK: MovingWall at x-Position %g!\n",
						(double)pnt_config->IBC_MovingActualPosition);}
			}
			if(pnt_config->flag_ManufacturedSolution==1)
			{
				if(pnt_config->MPI_rank==0){printf("SHOCK: Recent L2-Delta(rho): %.8Le)\n",
						pnt_config->ManufacturedSolution_L2_Delta);}
			}
			if(pnt_config->MPI_rank==0){printf("SHOCK: --------\n");}
		}
		if(pnt_config->flag_ManufacturedSolution==1)
		{
			ErrorManufacturedSolution(
					pnt_config,
					pnt_mesh,
					pnt_U_lastStep,
					0);
		}
	}
}





void setOptions(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep,
		struct strct_U * pnt_U_RK,
		struct strct_Flux * pnt_Q,
		struct strct_Flux * pnt_Q_sum,
		struct strct_Flux * pnt_Flux,
		struct strct_Flux * pnt_Flux_PlusHalf,
		struct strct_Film * pnt_Film)
{
	//	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Inviscid<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	if(pnt_config->flag_Inviscid==1){
		if(pnt_config->MPI_rank==0){printf(">>>>> Options: Inviscid (Euler)\n");}}


	//	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>PressureWaves<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	if(pnt_config->flag_PressureWaves==1){
		preparePressureWaves(
				pnt_config,
				pnt_mesh,
				pnt_U_lastStep);
		if(pnt_config->MPI_rank==0){printf(">>>>> Options: Druckwellen werden bei x=%f,y=%f,z=%f mit einem Radius von %f induziert\n",
				(double)pnt_config->pw_x0,
				(double)pnt_config->pw_y0,
				(double)pnt_config->pw_z0,
				(double)pnt_config->pw_r0);}}

	//	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>PressureHistory<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	if(pnt_config->flag_PressureHistory==1)	{
		if(pnt_config->MPI_rank==0){printf(">>>>> Options: PressureHistory (%d locations)\n",pnt_config->PressureHistory_No);}
		MPI_Barrier(pnt_config->MPI_comm);
		preparePressureHistoryValues(
				pnt_config,
				pnt_mesh);
		MPI_Barrier(pnt_config->MPI_comm);	}

	//	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>VelocityHistory<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	if(pnt_config->flag_VelocityHistory==1)	{
		if(pnt_config->MPI_rank==0){printf(">>>>> Options: VelocityHistory (%d locations)\n",pnt_config->VelocityHistory_No);}
		MPI_Barrier(pnt_config->MPI_comm);
		prepareVelocityHistoryValues(
				pnt_config,
				pnt_mesh);
		MPI_Barrier(pnt_config->MPI_comm);	}

	//	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>SpecialInitialisation<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	if((pnt_config->int_specialInitializeType>0)&&(pnt_config->int_initializeType==0)){
		InitializeSpecialConditions(
							pnt_config,
							pnt_mesh,
							pnt_U_lastStep);
		if(pnt_config->MPI_rank==0){printf(">>>>> Options: SpecialInitialisation: %d aktiviert.\n",pnt_config->int_specialInitializeType);}}
		
	//	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>LaminarBoundary<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	if(pnt_config->flag_LaminarBoundary==1){
		InitializeLaminarBoundary(
							pnt_config,
							pnt_mesh,
							pnt_U_lastStep);
		if(pnt_config->MPI_rank==0){printf(">>>>> Options: Laminare Grenzschicht wurd ab x=%f generiert.\n",(double)pnt_config->LaminarBoundary_xStart);}}

	//	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Vortex<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	if(pnt_config->flag_Vortex==1){
		InitializeVortex(
				pnt_config,
				pnt_mesh,
				pnt_U_lastStep);
		if(pnt_config->MPI_rank==0){printf(">>>>> Options: Vortex wurde generiert (x=%lf, y=%lf, r=%lf, f=%lf, beta=%lf).\n",
				(double)pnt_config->Vortex_x_wirb_zentr,(double)pnt_config->Vortex_y_wirb_zentr,(double)pnt_config->Vortex_r_wirb_max,(double)pnt_config->Vortex_faktor_quer,(double)pnt_config->Vortex_beta);}
		}

	//	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>movingWall<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	if((pnt_config->flag_IBC_Moving==1)&&(pnt_config->flag_IBC==1)){
		if(pnt_config->MPI_rank==0){printf(">>>>> Options: Moving Wall wurde eingestellt.\n");}
		if(pnt_config->MPI_rank==0){printf(">>>>> Options: Start: %g %g %g, Speed: %g, Movement: %d, SpeedFactor: %g.\n",(double)pnt_config->IBC_StartpositionX,(double)pnt_config->IBC_StartpositionY,(double)pnt_config->IBC_StartpositionZ,(double)pnt_config->IBC_MovingSpeed,pnt_config->IBC_MovingType, (double)pnt_config->IBC_SpeedFactor);}
		if(pnt_config->MPI_rank==0){printf(">>>>> Options: Aktuelle Position %g\n",(double)pnt_config->IBC_MovingActualPosition);}
	}

	//	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>rotation_symmetric<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	if(pnt_config->flag_rotation_symmetric==1){
		if(pnt_config->MPI_rank==0){printf(">>>>> Options: 2D-Rotationssymmetrisch\n");}}

	//	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>ManufacturedSolution<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	if(pnt_config->flag_ManufacturedSolution==1){
		ConfigureManufacturedSolution(
				pnt_config);
		if(pnt_config->MPI_rank==0){printf(">>>>> Options: ManufacturedSolution konfiguriert.\n");}

		if(pnt_config->MPI_rank==0){printf(">>>>> Options: ManufacturedSolution-Case: %d (0: sub-, 1: supersonic).\n",
				pnt_config->ManufacturedSolution_case);}

		InitializeManufacturedSolution(
				pnt_config,
				pnt_mesh,
				pnt_U_lastStep);
		if(pnt_config->MPI_rank==0){printf(">>>>> Options: ManufacturedSolution initialisiert.\n");}

		SetAllBoundaryConditions(
				pnt_config,
				pnt_mesh,
				pnt_U_lastStep,
				pnt_U_lastStep);
		if(pnt_config->MPI_rank==0){printf(">>>>> Options: ManufacturedSolution BoundaryConditions aktualisiert.\n");}
	}

}


