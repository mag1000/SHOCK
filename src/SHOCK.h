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
	float AlphaNonRef;

	float flt_is_minimum;
	float flt_is_maximum;
	float flt_is_avrg;
	float flt_is_avrg_counter;


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
	float *PressureHistory_x_P,*PressureHistory_y_P,*PressureHistory_z_P;
	float *PressureHistory_x_P_real,*PressureHistory_y_P_real,*PressureHistory_z_P_real;
	float *PressureHistory_time;
	float **PressureHistory_pressure;
	int *ijk_PressureHistory_P;

	int VelocityHistory_No;
	float *VelocityHistory_x_P,*VelocityHistory_y_P,*VelocityHistory_z_P;
	float *VelocityHistory_x_P_real,*VelocityHistory_y_P_real,*VelocityHistory_z_P_real;
	float *VelocityHistory_time;
	float **VelocityHistory_VelocityX;
	float **VelocityHistory_VelocityY;
	float **VelocityHistory_VelocityZ;
	int *ijk_VelocityHistory_P;

	float LaminarBoundary_xStart;

	float global_lambdaMax;

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
	float flt_p_out;
	float flt_T_wall;

	int flag_BC_option_inflow_normal_sub;
	int flag_BC_option_inflow_riemann_sub;
	int flag_BC_option_inflow_isentrop_sub;
	int flag_BC_option_inflow_normal_super;

	int flag_BC_option_outflow_normal_sub;
	int flag_BC_option_outflow_riemann_sub;
	int flag_BC_option_outflow_rudy_sub;

	float flt_p_inflow;
	float flt_rho_inflow;
	float flt_u_inflow;
	float flt_v_inflow;
	float flt_w_inflow;

	float flt_deltaXi;
	float flt_deltaEta;
	float flt_deltaZeta;

	float flt_AoA;
	float flt_numericalTau;
	float flt_numericalTauStart;
	float flt_machNumber;
	float flt_reynoldsNumber;
	float flt_prandtlNumber;
	float flt_gammaNumber;
	float flt_gasConstantNumber;

	float flt_wenoP;
	float flt_wenoEpsilon;
	float flt_wenoOptimalerKoeffizient_W9[5];
	float flt_wenoOptimalerKoeffizient_W5[3];

	//fuer neue ZD-Berechnung
	float *flt_ZD_Interpolation_Koeffizient;
	float *flt_ZD_ZweiteAbleitungZwischenPunkt_Koeffizient;
	float *flt_ZD_Ableitung_Koeffizient;
	float *flt_ZD_ZweiteAbleitung_Koeffizient;

	float flt_RK_U_n_Faktor[4];
	float flt_RK_U_ABC_Faktor[4];
	float flt_RK_Q_Faktor[4];
	float flt_RK_Q_Summe_Flag[4];
	float flt_RK_Q_Summe_Faktor[4];

	float flt_Upsilon;
	float flt_Psi;
	float *flt_Gamma; //ist jetzt eine lokale Größe, um den Wärmefluss (lokales Lambda) unabhängig vom Reibungsterm (nü) zu bestimmen

	float flt_SutherlandConstant;
	float flt_T0_dim;

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

	float **MPI_SendBufferMesh;
	float **MPI_RecieveBufferMesh;
//	float **MPI_SendBufferFlow;
//	float **MPI_RecieveBufferFlow;
	float **MPI_SendBufferFlowWithGhosts;
	float **MPI_RecieveBufferFlowWithGhosts;


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


	float *MPI_fltTransformation_xi_x;
	float *MPI_fltTransformation_xi_y;
	float *MPI_fltTransformation_xi_z;
	float *MPI_fltTransformation_eta_x;
	float *MPI_fltTransformation_eta_y;
	float *MPI_fltTransformation_eta_z;
	float *MPI_fltTransformation_zeta_x;
	float *MPI_fltTransformation_zeta_y;
	float *MPI_fltTransformation_zeta_z;

    MPI_Status MPI_status;
    MPI_Comm MPI_comm;
    MPI_Info MPI_info;

//    float * bufferSendFlowLeft;
//	float * bufferSendFlowRight;
//	float * bufferSendFlowBottom;
//	float * bufferSendFlowTop;
//	float * bufferSendFlowBehind;
//	float * bufferSendFlowInFront;
//	float * bufferRecieveFlowLeft;
//	float * bufferRecieveFlowRight;
//	float * bufferRecieveFlowBottom;
//	float * bufferRecieveFlowTop;
//	float * bufferRecieveFlowBehind;
//	float * bufferRecieveFlowInFront;

    float * bufferSendFlowWithGhostsLeft;
	float * bufferSendFlowWithGhostsRight;
	float * bufferSendFlowWithGhostsBottom;
	float * bufferSendFlowWithGhostsTop;
	float * bufferSendFlowWithGhostsBehind;
	float * bufferSendFlowWithGhostsInFront;
	float * bufferRecieveFlowWithGhostsLeft;
	float * bufferRecieveFlowWithGhostsRight;
	float * bufferRecieveFlowWithGhostsBottom;
	float * bufferRecieveFlowWithGhostsTop;
	float * bufferRecieveFlowWithGhostsBehind;
	float * bufferRecieveFlowWithGhostsInFront;

	float * bufferSendMeshLeft;
	float * bufferSendMeshRight;
	float * bufferSendMeshBottom;
	float * bufferSendMeshTop;
	float * bufferSendMeshBehind;
	float * bufferSendMeshInFront;
	float * bufferRecieveMeshLeft;
	float * bufferRecieveMeshRight;
	float * bufferRecieveMeshBottom;
	float * bufferRecieveMeshTop;
	float * bufferRecieveMeshBehind;
	float * bufferRecieveMeshInFront;

    float flt_L0_dim;
    float flt_u0_dim;
    float flt_c0_dim;
	float flt_time_dim;
	float flt_time_dim_lastAction;
	float flt_time_dim_backup1;
	float flt_time_dim_backup2;

//    float flt_c0;


//	IBC
	int flag_IBC;
	int IBC_Type;
    float IBC_yKolben;
	float IBC_alphaKolben;

//	IBC:Moving
	int flag_IBC_Moving;
	float IBC_StartpositionX;
	float IBC_StartpositionY;
	float IBC_StartpositionZ;
	float IBC_SizeX;
	float IBC_SizeY;
	float IBC_SizeZ;
	float IBC_MovingSpeed;
	int IBC_MovingType;
	float IBC_SpeedFactor;
	int IBC_MovingStepsize;
	float IBC_MovingActualPosition;
	float IBC_MovingLastPosition;

//	IBC:Apply Boundary Conditions at Walls of IBC
	int flag_IBC_ApplyBC;


    //    VortexGenerator
	int flag_Vortex;
	float Vortex_x_wirb_zentr;
	float Vortex_beta;
	float Vortex_y_wirb_zentr;
	float Vortex_faktor_quer;
	float Vortex_r_wirb_max;

	float start_Time;


	//	Pressure Waves
	int flag_PressureWaves;
	int pw_UseBC;
	int pw_UseFlowAverage;
	int pw_numberSources;
	float pw_amplitude;
	float pw_frequency;
	float pw_x0;
	float pw_y0;
	float pw_z0;
	float pw_x1;
	float pw_y1;
	float pw_z1;
	float pw_r0;

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
	float InitializeValues_u0;
	float InitializeValues_p0;
	float InitializeValues_rho0;
	float InitializeValues_u1;
	float InitializeValues_p1;
	float InitializeValues_rho1;
	float InitializeValues_xBorder;


	//	2D-Rotation-Symmetric
	int flag_rotation_symmetric;

	//Tau
	int flag_TauAccelerator;
	int flag_reinitialization;
	float flt_TauAccelerator_factor;
	float flt_TauDecelerator_factor;
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
	float *BC_Corrector;

	float *x;
	float *y;
	float *z;

	float *x_extrapolate;
	float *y_extrapolate;
	float *z_extrapolate;

	float *xi_x;
	float *xi_y;
	float *xi_z;

	float *eta_x;
	float *eta_y;
	float *eta_z;

	float *zeta_x;
	float *zeta_y;
	float *zeta_z;

	float *jacobian;

	//    ImmergedBC
	int* flag_IBC;
	int* flag_IBC_last;

	//	Pressure-Waves
	int* flag_PressureWaves;
	float* startPressure_PressureWaves;
	float* startDensity_PressureWaves;

//	Metric fuer viskose Fluesse
	float *xiFluss_Faktor;
	float *etaFluss_Faktor;
	float *zetaFluss_Faktor;
};

struct strct_U
{
	float *rho;
	float *u;
	float *v;
	float *w;
	float *p;
	float *e;
	float *theta1;
	float *theta2;
	float *theta3;
	float *c;
	float *gradRho;
	float *Lambda2;
	float *MachNumber;
	float *T;
	float *mue;

	float *u_xi;
	float *u_eta;
	float *u_zeta;
	float *v_xi;
	float *v_eta;
	float *v_zeta;
	float *w_xi;
	float *w_eta;
	float *w_zeta;
	float *T_xi;
	float *T_eta;
	float *T_zeta;
};

struct strct_Film
{
	float *rho;
	float *u;
	float *v;
	float *w;
	float *p;
	float *gradRho;
	float *Lambda2;
	float *MachNumber;
	float *flt_time_dim;
};

struct strct_Flux
{
	float *Mass;
	float *xiMomentum;
	float *etaMomentum;
	float *zetaMomentum;
	float *Energy;
};

//Derzeit nicht unter Verwendung
struct strct_ZD
{
	float *u_xi;
	float *u_eta;
	float *u_zeta;
	float *v_xi;
	float *v_eta;
	float *v_zeta;
	float *w_xi;
	float *w_eta;
	float *w_zeta;
	float *T_xi;
	float *T_eta;
	float *T_zeta;
	float *mue_xi;
	float *mue_eta;
	float *mue_zeta;

	float *tau_xx;
	float *tau_yy;
	float *tau_zz;
	float *tau_xy;
	float *tau_xz;
	float *tau_yz;
	float *q_x;
	float *q_y;
	float *q_z;
};

//Variables
extern float flt_leftEigenvector[5][5];
extern float flt_rightEigenvector[5][5];

extern int int_interations;
extern int int_interationsStart;
extern int int_interationsEnd;

extern int int_myCPUID;
extern int int_MaxNumberCPUs;

extern struct strct_configuration configuration;
extern struct strct_mesh *mesh;

extern float flt_helpValue1;
extern float flt_helpValue2;

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

