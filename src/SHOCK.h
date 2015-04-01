#include "mpi.h"
#include "cgnslib.h"
#define NO_NEIGHBOUR -1

#ifndef PRECISION
#define PRECISION 3
#endif
#if PRECISION == 1
	#define FLT_name "float"
	#define FLT float
	#define MY_FLT_MIN FLT_MIN
	#define MPI_FLT MPI_FLOAT
	#define MY_PI 3.14159265358979323846264338327950288419716939937510
	#define CONV_ERROR	1.0E-10L
#elif PRECISION == 2
	#define FLT_name "double"
	#define FLT double
	#define MY_FLT_MIN DBL_MIN
	#define MPI_FLT MPI_DOUBLE
	#define MY_PI 3.14159265358979323846264338327950288419716939937510
	#define CONV_ERROR	1.0E-18L
#elif PRECISION == 3
	#define FLT_name "long double"
	#define FLT long double
	#define MY_FLT_MIN DBL_MIN
	#define MPI_FLT MPI_LONG_DOUBLE
	#define MY_PI 3.14159265358979323846264338327950288419716939937510L
	#define CONV_ERROR	1.0E-18L
#elif PRECISION == 4
	#define FLT_name "quad"
	#define FLT __float128
	#define MY_FLT_MIN DBL_MIN
	#define MPI_FLT MPI_LONG_DOUBLE
	#define MY_PI 3.14159265358979323846264338327950288419716939937510Q
	#define CONV_ERROR	1.0E-20L
#endif

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
	FLT AlphaNonRef;

	FLT is_minimum;
	FLT is_maximum;
	FLT is_avrg;
	FLT is_avrg_counter;


	int flag_ManufacturedSolution;
	int flag_ReducedExport;
	int flag_Inviscid;
	int flag_NAN;
	int flag_constantZValues;
	int flag_LaminarBoundary;
	int int_specialInitializeType;
	int flag_PressureHistory;
	int *flag_PressureHistory_P;
	int flag_VelocityHistory;
	int *flag_VelocityHistory_P;


	int PressureHistory_No;
	FLT *PressureHistory_x_P,*PressureHistory_y_P,*PressureHistory_z_P;
	FLT *PressureHistory_x_P_real,*PressureHistory_y_P_real,*PressureHistory_z_P_real;
	FLT *PressureHistory_time;
	FLT **PressureHistory_pressure;
	int *ijk_PressureHistory_P;

	int VelocityHistory_No;
	FLT *VelocityHistory_x_P,*VelocityHistory_y_P,*VelocityHistory_z_P;
	FLT *VelocityHistory_x_P_real,*VelocityHistory_y_P_real,*VelocityHistory_z_P_real;
	FLT *VelocityHistory_time;
	FLT **VelocityHistory_VelocityX;
	FLT **VelocityHistory_VelocityY;
	FLT **VelocityHistory_VelocityZ;
	int *ijk_VelocityHistory_P;

	FLT LaminarBoundary_xStart;

	FLT global_lambdaMax;

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

	int flag_swapDivisionFile;

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
	FLT p_out;
	FLT T_wall;

	int flag_BC_option_inflow_normal_sub;
	int flag_BC_option_inflow_riemann_sub;
	int flag_BC_option_inflow_isentrop_sub;
	int flag_BC_option_inflow_normal_super;

	int flag_BC_option_outflow_normal_sub;
	int flag_BC_option_outflow_riemann_sub;
	int flag_BC_option_outflow_rudy_sub;

	FLT p_inflow;
	FLT rho_inflow;
	FLT u_inflow;
	FLT v_inflow;
	FLT w_inflow;

	FLT deltaXi;
	FLT deltaEta;
	FLT deltaZeta;

	FLT AoA;
	FLT numericalTau;
	FLT numericalTauStart;
	FLT machNumber;
	FLT reynoldsNumber;
	FLT prandtlNumber;
	FLT gammaNumber;
	FLT gasConstantNumber;

	FLT wenoP;
	FLT wenoEpsilon;
	FLT wenoOptimalerKoeffizient_W9[5];
	FLT wenoOptimalerKoeffizient_W5[3];

	//fuer neue ZD-Berechnung
	FLT *ZD_Interpolation_Koeffizient;
	FLT *ZD_AbleitungZwischenPunkt_Koeffizient;
	FLT *ZD_Ableitung_Koeffizient;
	FLT *ZD_ZweiteAbleitung_Koeffizient;

	FLT RK_U_n_Faktor[4];
	FLT RK_U_ABC_Faktor[4];
	FLT RK_Q_Faktor[4];
	FLT RK_Q_Summe_Flag[4];
	FLT RK_Q_Summe_Faktor[4];

	FLT Upsilon;
	FLT Psi;
	FLT *Gamma; //ist jetzt eine lokale Größe, um den Wärmefluss (lokales Lambda) unabhängig vom Reibungsterm (nü) zu bestimmen

	FLT SutherlandConstant;
	FLT T0_dim;

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

	FLT **MPI_SendBufferMesh;
	FLT **MPI_RecieveBufferMesh;
//	FLT **MPI_SendBufferFlow;
//	FLT **MPI_RecieveBufferFlow;
	FLT **MPI_SendBufferFlowWithGhosts;
	FLT **MPI_RecieveBufferFlowWithGhosts;


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


	FLT *MPI_dblTransformation_xi_x;
	FLT *MPI_dblTransformation_xi_y;
	FLT *MPI_dblTransformation_xi_z;
	FLT *MPI_dblTransformation_eta_x;
	FLT *MPI_dblTransformation_eta_y;
	FLT *MPI_dblTransformation_eta_z;
	FLT *MPI_dblTransformation_zeta_x;
	FLT *MPI_dblTransformation_zeta_y;
	FLT *MPI_dblTransformation_zeta_z;

    MPI_Status MPI_status;
    MPI_Comm MPI_comm;
    MPI_Info MPI_info;

//    FLT * bufferSendFlowLeft;
//	FLT * bufferSendFlowRight;
//	FLT * bufferSendFlowBottom;
//	FLT * bufferSendFlowTop;
//	FLT * bufferSendFlowBehind;
//	FLT * bufferSendFlowInFront;
//	FLT * bufferRecieveFlowLeft;
//	FLT * bufferRecieveFlowRight;
//	FLT * bufferRecieveFlowBottom;
//	FLT * bufferRecieveFlowTop;
//	FLT * bufferRecieveFlowBehind;
//	FLT * bufferRecieveFlowInFront;

    FLT * bufferSendFlowWithGhostsLeft;
	FLT * bufferSendFlowWithGhostsRight;
	FLT * bufferSendFlowWithGhostsBottom;
	FLT * bufferSendFlowWithGhostsTop;
	FLT * bufferSendFlowWithGhostsBehind;
	FLT * bufferSendFlowWithGhostsInFront;
	FLT * bufferRecieveFlowWithGhostsLeft;
	FLT * bufferRecieveFlowWithGhostsRight;
	FLT * bufferRecieveFlowWithGhostsBottom;
	FLT * bufferRecieveFlowWithGhostsTop;
	FLT * bufferRecieveFlowWithGhostsBehind;
	FLT * bufferRecieveFlowWithGhostsInFront;

	FLT * bufferSendMeshLeft;
	FLT * bufferSendMeshRight;
	FLT * bufferSendMeshBottom;
	FLT * bufferSendMeshTop;
	FLT * bufferSendMeshBehind;
	FLT * bufferSendMeshInFront;
	FLT * bufferRecieveMeshLeft;
	FLT * bufferRecieveMeshRight;
	FLT * bufferRecieveMeshBottom;
	FLT * bufferRecieveMeshTop;
	FLT * bufferRecieveMeshBehind;
	FLT * bufferRecieveMeshInFront;

    FLT L0_dim;
    FLT u0_dim;
    FLT c0_dim;
	FLT time_dim;
	FLT time_dim_lastAction;
	FLT time_dim_backup1;
	FLT time_dim_backup2;

	//ManufacturedSolution
	int ManufacturedSolution_case;
	char BCManufacturedSolution[30];
	long double ManufacturedSolution_L2_last;
	long double ManufacturedSolution_L2_last_pressure;
	int ManufacturedSolution_L2_Converged[15];
	long double ManufacturedSolution_L2_Delta;
	char ManufacturedSolution_L2_Delta_name[32];
	long double all_L2_norm_rho;
	long double all_L2_norm_pressure;
	long double all_Linf_norm_rho;
	long double all_Linf_norm_pressure;
	long double ManufacturedSolution_param_rho[7];
	long double ManufacturedSolution_param_u[7];
	long double ManufacturedSolution_param_v[7];
	long double ManufacturedSolution_param_w[7];
	long double ManufacturedSolution_param_p[7];
	int ManufacturedSolution_L2_counter;


//    FLT c0;


//	IBC
	int flag_IBC;
	int IBC_Type;
    FLT IBC_yKolben;
	FLT IBC_alphaKolben;

//	IBC:Moving
	int flag_IBC_Moving;
	FLT IBC_StartpositionX;
	FLT IBC_StartpositionY;
	FLT IBC_StartpositionZ;
	FLT IBC_SizeX;
	FLT IBC_SizeY;
	FLT IBC_SizeZ;
	FLT IBC_MovingSpeed;
	int IBC_MovingType;
	FLT IBC_SpeedFactor;
	int IBC_MovingStepsize;
	FLT IBC_MovingActualPosition;
	FLT IBC_MovingLastPosition;



    //    VortexGenerator
	int flag_Vortex;
	FLT Vortex_x_wirb_zentr;
	FLT Vortex_beta;
	FLT Vortex_y_wirb_zentr;
	FLT Vortex_faktor_quer;
	FLT Vortex_r_wirb_max;

	FLT start_Time;


	//	Pressure Waves
	int flag_PressureWaves;
	int pw_UseBC;
	int pw_UseFlowAverage;
	int pw_numberSources;
	FLT pw_amplitude;
	FLT pw_frequency;
	FLT pw_x0;
	FLT pw_y0;
	FLT pw_z0;
	FLT pw_x1;
	FLT pw_y1;
	FLT pw_z1;
	FLT pw_r0;

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
	FLT InitializeValues_u0;
	FLT InitializeValues_p0;
	FLT InitializeValues_rho0;
	FLT InitializeValues_u1;
	FLT InitializeValues_p1;
	FLT InitializeValues_rho1;
	FLT InitializeValues_xBorder;


	//	2D-Rotation-Symmetric
	int flag_rotation_symmetric;

	//Tau
	int flag_TauAccelerator;
	int flag_reinitialization;
	FLT TauAccelerator_factor;
	FLT TauDecelerator_factor;
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
	FLT *BC_Corrector_xiMomentum;
	FLT *BC_Corrector_etaMomentum;
	FLT *BC_Corrector_zetaMomentum;

	FLT *x;
	FLT *y;
	FLT *z;

	FLT *x_extrapolate;
	FLT *y_extrapolate;
	FLT *z_extrapolate;

	FLT *xi_x;
	FLT *xi_y;
	FLT *xi_z;

	FLT *eta_x;
	FLT *eta_y;
	FLT *eta_z;

	FLT *zeta_x;
	FLT *zeta_y;
	FLT *zeta_z;

	FLT *jacobian;

	//    ImmergedBC
	int* flag_IBC;
	int* flag_IBC_last;

	//	Pressure-Waves
	int* flag_PressureWaves;
	FLT* startPressure_PressureWaves;
	FLT* startDensity_PressureWaves;

//	Metric fuer viskose Fluesse
	FLT *xiFluss_Faktor;
	FLT *etaFluss_Faktor;
	FLT *zetaFluss_Faktor;
};

struct strct_U
{
	FLT *rho;
	FLT *u;
	FLT *v;
	FLT *w;
	FLT *p;
	FLT *e;
	FLT *theta1;
	FLT *theta2;
	FLT *theta3;
	FLT *c;
	FLT *gradRho;
	FLT *Lambda2;
	FLT *MachNumber;
	FLT *T;
	FLT *mue;

	FLT *u_xi;
	FLT *u_eta;
	FLT *u_zeta;
	FLT *v_xi;
	FLT *v_eta;
	FLT *v_zeta;
	FLT *w_xi;
	FLT *w_eta;
	FLT *w_zeta;
	FLT *T_xi;
	FLT *T_eta;
	FLT *T_zeta;
};

struct strct_Film
{
	FLT *rho;
	FLT *u;
	FLT *v;
	FLT *w;
	FLT *p;
	FLT *gradRho;
	FLT *Lambda2;
	FLT *MachNumber;
	FLT *time_dim;
};

struct strct_Flux
{
	FLT *Mass;
	FLT *xiMomentum;
	FLT *etaMomentum;
	FLT *zetaMomentum;
	FLT *Energy;
};

//Variables
extern FLT leftEigenvector[5][5];
extern FLT rightEigenvector[5][5];

extern int int_interations;
extern int int_interationsStart;
extern int int_interationsEnd;

extern int int_myCPUID;
extern int int_MaxNumberCPUs;

extern struct strct_configuration configuration;
extern struct strct_mesh *mesh;

extern FLT helpValue1;
extern FLT helpValue2;

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

