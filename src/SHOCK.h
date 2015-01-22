#include "mpi.h"
#include "cgnslib.h"

#define NO_NEIGHBOUR -1

#ifndef SPACEORDER
#define SPACEORDER 5
#endif

#ifndef GITID
#define GITID "00000"
#endif

#ifndef MESHDIMENSIONS
#define MESHDIMENSIONS 2
#endif

#ifndef SHOCK_H
#define SHOCK_H

#ifndef OMPI_MAJOR_VERSION
#define OMPI_MAJOR_VERSION 0
#define OMPI_MINOR_VERSION 0
#define OMPI_RELEASE_VERSION 0
#endif


#define HDF5_HAVE_MULTI_DATASETS_ 0
#if (H5_VERS_MAJOR > 1) || (H5_VERS_MAJOR == 1 && H5_VERS_MINOR > 8)
#	define HDF5_HAVE_MULTI_DATASETS_ 1
#else
#	if (H5_VERS_MAJOR == 1 && H5_VERS_MINOR > 7 && H5_VERS_RELEASE > 13)
#		define HDF5_HAVE_MULTI_DATASETS_ 1
#	endif
#endif

struct strct_configuration
{
	//	Rudy-BC
	double AlphaNonRef;

	double dbl_is_minimum;
	double dbl_is_maximum;
	double dbl_is_avrg;
	double dbl_is_avrg_counter;


	int flag_ReducedExport;
	int flag_Inviscid;
	int flag_NAN;
	int flag_constantZValues;
	int flag_LaminarBoundary;
	int flag_PressureHistory;
	int *flag_PressureHistory_P;
	int flag_VelocityHistory;
	int *flag_VelocityHistory_P;


	int PressureHistory_No;
	double *PressureHistory_x_P,*PressureHistory_y_P,*PressureHistory_z_P;
	double *PressureHistory_x_P_real,*PressureHistory_y_P_real,*PressureHistory_z_P_real;
	double *PressureHistory_time;
	double **PressureHistory_pressure;
	int *ijk_PressureHistory_P;

	int VelocityHistory_No;
	double *VelocityHistory_x_P,*VelocityHistory_y_P,*VelocityHistory_z_P;
	double *VelocityHistory_x_P_real,*VelocityHistory_y_P_real,*VelocityHistory_z_P_real;
	double *VelocityHistory_time;
	double **VelocityHistory_VelocityX;
	double **VelocityHistory_VelocityY;
	double **VelocityHistory_VelocityZ;
	int *ijk_VelocityHistory_P;

	double LaminarBoundary_xStart;

	double global_lambdaMax;

	int int_CGNS_Sample_File;
	int int_CGNS_Sample_Base;
	int int_CGNS_Sample_Zone;
	int int_CGNS_Sample_FlowSolution;

	char chr_configPath[500];
	char chr_DivisionFile[200];
	char chr_DivisionPath[500];
	char chr_MeshFile[200];
	char chr_MeshPath[500];
	char chr_MeshPathOriginal[500];
	char chr_PressureHistoryFile[200];
	char chr_PressureHistoryPath[500];
	char chr_VelocityHistoryFile[200];
	char chr_VelocityHistoryPath[500];
	char chr_folder[500];

	int int_NoCPUsI;
	int int_NoCPUsJ;
	int int_NoCPUsK;

	int int_initializeType;
	int int_EndIteration;
	int int_StartIteration;
	int int_StartSampling;
	int int_Samples;
	int int_actualIteration;
	int int_actualIteration_backup1;
	int int_actualIteration_backup2;
	int int_iterationCounterSamples_backup1;
	int int_iterationCounterSamples_backup2;
	int int_iterationCounterBackupOut_backup1;
	int int_iterationCounterBackupOut_backup2;
	int int_iterationCounterStdOut_backup1;
	int int_iterationCounterStdOut_backup2;
	int int_TotalIterations;
	int int_actualSample;
	int int_actualSample_backup1;
	int int_actualSample_backup2;
	int int_IterationsBetweenSamples;

	int int_meshDimensions;
	int int_conservationEquations;

	int flag_extrapolate3Dimension;
	int flag_exportMetric;


	int int_NumberBackups;

	int int_iMeshPoints;
	int int_jMeshPoints;
	int int_kMeshPoints;

	int int_iMeshPointsGhostCells;
	int int_jMeshPointsGhostCells;
	int int_kMeshPointsGhostCells;

	int int_SpaceOrder;
	int int_TimeOrder;

	int int_iMid;
	int int_iStartReal;
	int int_iEndReal;
	int int_iStartGhosts;
	int int_iEndGhosts;

	int int_jMid;
	int int_jStartReal;
	int int_jEndReal;
	int int_jStartGhosts;
	int int_jEndGhosts;

	int int_kMid;
	int int_kStartReal;
	int int_kEndReal;
	int int_kStartGhosts;
	int int_kEndGhosts;

	int ijkMid;

	int int_kStartReal_original;
	int int_kEndReal_original;
	int int_kStartGhosts_original;
	int int_kEndGhosts_original;

	double comm_time;
	double dbl_p_out;
	double dbl_T_wall;

	int flag_BC_option_inflow_normal_sub;
	int flag_BC_option_inflow_riemann_sub;
	int flag_BC_option_inflow_isentrop_sub;
	int flag_BC_option_inflow_normal_super;

	int flag_BC_option_outflow_normal_sub;
	int flag_BC_option_outflow_riemann_sub;
	int flag_BC_option_outflow_rudy_sub;

	double dbl_p_inflow;
	double dbl_rho_inflow;
	double dbl_u_inflow;
	double dbl_v_inflow;
	double dbl_w_inflow;

	double dbl_deltaXi;
	double dbl_deltaEta;
	double dbl_deltaZeta;

	double dbl_AoA;
	double dbl_numericalTau;
	double dbl_numericalTauStart;
	double dbl_machNumber;
	double dbl_reynoldsNumber;
	double dbl_prandtlNumber;
	double dbl_gammaNumber;
	double dbl_gasConstantNumber;

	double dbl_wenoP;
	double dbl_wenoEpsilon;
	double dbl_wenoOptimalerKoeffizient_W9[5];
	double dbl_wenoOptimalerKoeffizient_W5[3];

	//fuer neue ZD-Berechnung
	double *dbl_ZD_Interpolation_Koeffizient;
	double *dbl_ZD_ZweiteAbleitungZwischenPunkt_Koeffizient;
	double *dbl_ZD_Ableitung_Koeffizient;
	double *dbl_ZD_ZweiteAbleitung_Koeffizient;

	double dbl_RK_U_n_Faktor[4];
	double dbl_RK_U_ABC_Faktor[4];
	double dbl_RK_Q_Faktor[4];
	double dbl_RK_Q_Summe_Flag[4];
	double dbl_RK_Q_Summe_Faktor[4];

	double dbl_Upsilon;
	double dbl_Psi;
	double *dbl_Gamma; //ist jetzt eine lokale Größe, um den Wärmefluss (lokales Lambda) unabhängig vom Reibungsterm (nü) zu bestimmen

	double dbl_SutherlandConstant;
	double dbl_T0_dim;

	int * MPI_intArray_NoCPUs;
	int MPI_rank;
	int MPI_size;
	int MPI_intMyDatatypeSize;
	int *MPI_intTransferSizeMesh;
	int *MPI_intTransferSizeFlow_WithGhosts;
	int InterfaceNeighbourLeft;
	int InterfaceNeighbourRight;
	int InterfaceNeighbourTop;
	int InterfaceNeighbourBottom;
	int InterfaceNeighbourInFront;
	int InterfaceNeighbourBehind;


	char Zonename[33];
	cgsize_t zonesize[ MESHDIMENSIONS ];
	cgsize_t offset[ MESHDIMENSIONS ];
	char **ZonenameAll;
	char **Donorname;
	char **Interfacename;
	cgsize_t **RangeOfInterface;
	cgsize_t **DonorRangeOfInterface;
	int **TransformMatrixOfInterface;
	float **RotationCenter;
	float **RotationAngle;
	float **Translation;

	int NumberInterfaces;
	long *MPI_tag;
	int *MPI_rankNeighbours;
	int *MPI_intIStartSend;
	int *MPI_intIEndSend;
	int *MPI_intJStartSend;
	int *MPI_intJEndSend;
	int *MPI_intKStartSend;
	int *MPI_intKEndSend;
	int *MPI_intIStartRecieve;
	int *MPI_intIEndRecieve;
	int *MPI_intJStartRecieve;
	int *MPI_intJEndRecieve;
	int *MPI_intKStartRecieve;
	int *MPI_intKEndRecieve;

	int *MPI_intIStartSend_WithGhosts;
	int *MPI_intIEndSend_WithGhosts;
	int *MPI_intJStartSend_WithGhosts;
	int *MPI_intJEndSend_WithGhosts;
	int *MPI_intKStartSend_WithGhosts;
	int *MPI_intKEndSend_WithGhosts;
	int *MPI_intIStartRecieve_WithGhosts;
	int *MPI_intIEndRecieve_WithGhosts;
	int *MPI_intJStartRecieve_WithGhosts;
	int *MPI_intJEndRecieve_WithGhosts;
	int *MPI_intKStartRecieve_WithGhosts;
	int *MPI_intKEndRecieve_WithGhosts;

	double **MPI_SendBufferMesh;
	double **MPI_RecieveBufferMesh;
//	double **MPI_SendBufferFlow;
//	double **MPI_RecieveBufferFlow;
	double **MPI_SendBufferFlowWithGhosts;
	double **MPI_RecieveBufferFlowWithGhosts;


	int *MPI_intTransformation_IMax;
	int *MPI_intTransformation_JMax;
	int *MPI_intTransformation_KMax;
	int *MPI_intTransformation_IMax_Mesh;
	int *MPI_intTransformation_JMax_Mesh;
	int *MPI_intTransformation_KMax_Mesh;
	int *MPI_intTransformation_flag_I0_I;
	int *MPI_intTransformation_flag_I0_J;
	int *MPI_intTransformation_flag_I0_K;
	int *MPI_intTransformation_flag_J0_I;
	int *MPI_intTransformation_flag_J0_J;
	int *MPI_intTransformation_flag_J0_K;
	int *MPI_intTransformation_flag_K0_I;
	int *MPI_intTransformation_flag_K0_J;
	int *MPI_intTransformation_flag_K0_K;
	int *MPI_intTransformation_Offset_I;
	int *MPI_intTransformation_Offset_J;
	int *MPI_intTransformation_Offset_K;
	int *MPI_intTransformation_Offset_I_Ghosts;
	int *MPI_intTransformation_Offset_J_Ghosts;
	int *MPI_intTransformation_Offset_K_Ghosts;


	double *MPI_dblTransformation_xi_x;
	double *MPI_dblTransformation_xi_y;
	double *MPI_dblTransformation_xi_z;
	double *MPI_dblTransformation_eta_x;
	double *MPI_dblTransformation_eta_y;
	double *MPI_dblTransformation_eta_z;
	double *MPI_dblTransformation_zeta_x;
	double *MPI_dblTransformation_zeta_y;
	double *MPI_dblTransformation_zeta_z;

    MPI_Status MPI_status;
    MPI_Comm MPI_comm;
    MPI_Info MPI_info;

//    double * bufferSendFlowLeft;
//	double * bufferSendFlowRight;
//	double * bufferSendFlowBottom;
//	double * bufferSendFlowTop;
//	double * bufferSendFlowBehind;
//	double * bufferSendFlowInFront;
//	double * bufferRecieveFlowLeft;
//	double * bufferRecieveFlowRight;
//	double * bufferRecieveFlowBottom;
//	double * bufferRecieveFlowTop;
//	double * bufferRecieveFlowBehind;
//	double * bufferRecieveFlowInFront;

    double * bufferSendFlowWithGhostsLeft;
	double * bufferSendFlowWithGhostsRight;
	double * bufferSendFlowWithGhostsBottom;
	double * bufferSendFlowWithGhostsTop;
	double * bufferSendFlowWithGhostsBehind;
	double * bufferSendFlowWithGhostsInFront;
	double * bufferRecieveFlowWithGhostsLeft;
	double * bufferRecieveFlowWithGhostsRight;
	double * bufferRecieveFlowWithGhostsBottom;
	double * bufferRecieveFlowWithGhostsTop;
	double * bufferRecieveFlowWithGhostsBehind;
	double * bufferRecieveFlowWithGhostsInFront;

	double * bufferSendMeshLeft;
	double * bufferSendMeshRight;
	double * bufferSendMeshBottom;
	double * bufferSendMeshTop;
	double * bufferSendMeshBehind;
	double * bufferSendMeshInFront;
	double * bufferRecieveMeshLeft;
	double * bufferRecieveMeshRight;
	double * bufferRecieveMeshBottom;
	double * bufferRecieveMeshTop;
	double * bufferRecieveMeshBehind;
	double * bufferRecieveMeshInFront;

    double dbl_L0_dim;
    double dbl_u0_dim;
    double dbl_c0_dim;
	double dbl_time_dim;
	double dbl_time_dim_lastAction;
	double dbl_time_dim_backup1;
	double dbl_time_dim_backup2;

//    double dbl_c0;


//	IBC
	int flag_IBC;
	int IBC_Type;
    double IBC_yKolben;
	double IBC_alphaKolben;

//	IBC:Moving
	int flag_IBC_Moving;
	double IBC_StartpositionX;
	double IBC_StartpositionY;
	double IBC_StartpositionZ;
	double IBC_SizeX;
	double IBC_SizeY;
	double IBC_SizeZ;
	double IBC_MovingSpeed;
	int IBC_MovingType;
	double IBC_SpeedFactor;
	int IBC_MovingStepsize;
	double IBC_MovingActualPosition;
	double IBC_MovingLastPosition;

//	IBC:Apply Boundary Conditions at Walls of IBC
	int flag_IBC_ApplyBC;


    //    VortexGenerator
	int flag_Vortex;
	double Vortex_x_wirb_zentr;
	double Vortex_beta;
	double Vortex_y_wirb_zentr;
	double Vortex_faktor_quer;
	double Vortex_r_wirb_max;

	double start_Time;


	//	Pressure Waves
	int flag_PressureWaves;
	int pw_UseBC;
	int pw_UseFlowAverage;
	int pw_numberSources;
	double pw_amplitude;
	double pw_frequency;
	double pw_x0;
	double pw_y0;
	double pw_z0;
	double pw_x1;
	double pw_y1;
	double pw_z1;
	double pw_r0;

//	BoundaryConditions
	char BC_Left[30];
	char BC_Right[30];
	char BC_Top[30];
	char BC_Bottom[30];
	char BC_Behind[30];
	char BC_InFront[30];

	char BCFarfield[30];
	char BCInflow[30];
	char BCOutflow[30];
	char BCOutflowSubsonic[30];
	char BCWallInviscid[30];
	char BCWallViscous[30];
	char BCInflowSupersonic[30];
	char BCInflowSubsonic[30];
	char BCWallViscousIsothermal[30];


//		InitializeValues
	double InitializeValues_u0;
	double InitializeValues_p0;
	double InitializeValues_rho0;
	double InitializeValues_u1;
	double InitializeValues_p1;
	double InitializeValues_rho1;


	//	2D-Rotation-Symmetric
	int flag_rotation_symmetric;

	//Tau
	int flag_TauAccelerator;
	int flag_reinitialization;
	double dbl_TauAccelerator_factor;
	double dbl_TauDecelerator_factor;
	int int_distanceNAN;
	int int_distanceForward;
	int int_distanceForwardStart;
	int int_distanceBackward;
	int int_NumberResets;


//	Extrapolate
	int int_iMeshPointsGhostCells_extrapolate;
	int int_jMeshPointsGhostCells_extrapolate;
	int int_kMeshPointsGhostCells_extrapolate;

	int int_iStartReal_extrapolate;
	int int_iEndReal_extrapolate;
	int int_iStartGhosts_extrapolate;
	int int_iEndGhosts_extrapolate;

	int int_jStartReal_extrapolate;
	int int_jEndReal_extrapolate;
	int int_jStartGhosts_extrapolate;
	int int_jEndGhosts_extrapolate;

	int int_kStartReal_extrapolate;
	int int_kEndReal_extrapolate;
	int int_kStartGhosts_extrapolate;
	int int_kEndGhosts_extrapolate;
};


struct strct_mesh
{
	double *BC_Corrector;

	double *x;
	double *y;
	double *z;

	double *x_extrapolate;
	double *y_extrapolate;
	double *z_extrapolate;

	double *xi_x;
	double *xi_y;
	double *xi_z;

	double *eta_x;
	double *eta_y;
	double *eta_z;

	double *zeta_x;
	double *zeta_y;
	double *zeta_z;

	double *jacobian;

	//    ImmergedBC
	int* flag_IBC;
	int* flag_IBC_last;

	//	Pressure-Waves
	int* flag_PressureWaves;
	double* startPressure_PressureWaves;
	double* startDensity_PressureWaves;

//	Metric fuer viskose Fluesse
	double *xiFluss_Faktor;
	double *etaFluss_Faktor;
	double *zetaFluss_Faktor;
};

struct strct_U
{
	double *rho;
	double *u;
	double *v;
	double *w;
	double *p;
	double *e;
	double *theta1;
	double *theta2;
	double *theta3;
	double *c;
	double *gradRho;
	double *Lambda2;
	double *MachNumber;
	double *T;
	double *mue;

	double *u_xi;
	double *u_eta;
	double *u_zeta;
	double *v_xi;
	double *v_eta;
	double *v_zeta;
	double *w_xi;
	double *w_eta;
	double *w_zeta;
	double *T_xi;
	double *T_eta;
	double *T_zeta;
};

struct strct_Film
{
	double *rho;
	double *u;
	double *v;
	double *w;
	double *p;
	double *gradRho;
	double *Lambda2;
	double *MachNumber;
	double *dbl_time_dim;
};

struct strct_Flux
{
	double *Mass;
	double *xiMomentum;
	double *etaMomentum;
	double *zetaMomentum;
	double *Energy;
};

//Derzeit nicht unter Verwendung
struct strct_ZD
{
	double *u_xi;
	double *u_eta;
	double *u_zeta;
	double *v_xi;
	double *v_eta;
	double *v_zeta;
	double *w_xi;
	double *w_eta;
	double *w_zeta;
	double *T_xi;
	double *T_eta;
	double *T_zeta;
	double *mue_xi;
	double *mue_eta;
	double *mue_zeta;

	double *tau_xx;
	double *tau_yy;
	double *tau_zz;
	double *tau_xy;
	double *tau_xz;
	double *tau_yz;
	double *q_x;
	double *q_y;
	double *q_z;
};

//Variables
extern double dbl_leftEigenvector[5][5];
extern double dbl_rightEigenvector[5][5];

extern int int_interations;
extern int int_interationsStart;
extern int int_interationsEnd;

extern int int_myCPUID;
extern int int_MaxNumberCPUs;

extern struct strct_configuration configuration;
extern struct strct_mesh *mesh;

extern double dbl_helpValue1;
extern double dbl_helpValue2;

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
		struct strct_U * pnt_U_backup2);

void setOptions(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep,
		struct strct_U * pnt_U_RK,
		struct strct_Flux * pnt_Q,
		struct strct_Flux * pnt_Q_sum,
		struct strct_Flux * pnt_Flux,
		struct strct_Flux * pnt_Flux_PlusHalf,
		struct strct_Film * pnt_Film);

//Variables of MeshImport

int flag;
#endif /*SHOCK_H */

