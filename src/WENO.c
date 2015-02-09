#include "WENO.h"
#include "SHOCK.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"

#define CONSVEQU 5


float Phi_Function_W9(
        struct strct_configuration * pnt_config,
        float * pnt_flux,
        float * pnt_deltaFlux)
{
    float flt_is0, flt_is1, flt_is2, flt_is3, flt_is4;
    float flt_omegaTmp0, flt_omegaTmp1, flt_omegaTmp2, flt_omegaTmp3, flt_omegaTmp4;
    float flt_omega0, flt_omega2, flt_omega3, flt_omega4;
    float flt_q1, flt_q2, flt_q3, flt_q4;
    float flt_result;

//  Berechnung der Glättungsindikatoren basierend auf FLUESSEN
//    Eine abweichende Berechnung auf Basis von DeltaFlüssen exisitert ebenfalls in der Literatur und wir von
//    Klioutchnikov bei WENO5 verwendet
    flt_is0 = pnt_flux[0] *(22658.*pnt_flux[0] - 208501.*pnt_flux[1] +  364863.*pnt_flux[2] -  288007.*pnt_flux[3] +  86329.*pnt_flux[4])
                     +pnt_flux[1] *(482963.*pnt_flux[1] - 1704396.*pnt_flux[2] + 1358458.*pnt_flux[3] - 411487.*pnt_flux[4])
                                  +pnt_flux[2] *(1521393.*pnt_flux[2] - 2462076.*pnt_flux[3] + 758823.*pnt_flux[4])
                                                +pnt_flux[3] *(1020563.*pnt_flux[3] - 649501.*pnt_flux[4])
+pnt_flux[4] *(107918.*pnt_flux[4]);

    flt_is1 = pnt_flux[1] *(6908.*pnt_flux[1] -  60871.*pnt_flux[2] +  99213.*pnt_flux[3] -  70237.*pnt_flux[4] +  18079.*pnt_flux[5])
                    +pnt_flux[2] *(138563.*pnt_flux[2] - 464976.*pnt_flux[3] + 337018.*pnt_flux[4] -  88297.*pnt_flux[5])
                                 +pnt_flux[3] *(406293.*pnt_flux[3] - 611976.*pnt_flux[4] + 165153.*pnt_flux[5])
                                              +pnt_flux[4] *(242723.*pnt_flux[4] - 140251.*pnt_flux[5])
+pnt_flux[5] *( 22658.*pnt_flux[5]);

    flt_is2 = pnt_flux[2] *(6908.*pnt_flux[2] -  51001.*pnt_flux[3] +  67923.*pnt_flux[4] -  38947.*pnt_flux[5] +  8209.*pnt_flux[6])
                    +pnt_flux[3] *(104963.*pnt_flux[3] - 299076.*pnt_flux[4] + 179098.*pnt_flux[5] - 38947.*pnt_flux[6])
                                 +pnt_flux[4] *(231153.*pnt_flux[4] - 299076.*pnt_flux[5] + 67923.*pnt_flux[6])
                                              +pnt_flux[5] *(104963.*pnt_flux[5] - 51001.*pnt_flux[6])
+pnt_flux[6] *( 6908.*pnt_flux[6]);

    flt_is3 = pnt_flux[3] *(22658.*pnt_flux[3] - 140251.*pnt_flux[4] + 165153.*pnt_flux[5] -  88297.*pnt_flux[6] + 18079.*pnt_flux[7])
                     +pnt_flux[4] *(242723.*pnt_flux[4] - 611976.*pnt_flux[5] + 337018.*pnt_flux[6] - 70237.*pnt_flux[7])
                                  +pnt_flux[5] *(406293.*pnt_flux[5] - 464976.*pnt_flux[6] + 99213.*pnt_flux[7])
                                               +pnt_flux[6] *(138563.*pnt_flux[6] - 60871.*pnt_flux[7])
+pnt_flux[7] *( 6908.*pnt_flux[7]);

    flt_is4 = pnt_flux[4] *(107918.*pnt_flux[4] - 649501.*pnt_flux[5] +  758823.*pnt_flux[6] -  411487.*pnt_flux[7] + 86329.*pnt_flux[8])
                      +pnt_flux[5] *(1020563.*pnt_flux[5] - 2462076.*pnt_flux[6] + 1358458.*pnt_flux[7] - 288007.*pnt_flux[8])
                                    +pnt_flux[6] *(1521393.*pnt_flux[6] - 1704396.*pnt_flux[7] + 364863.*pnt_flux[8])
                                                  +pnt_flux[7] *( 482963.*pnt_flux[7] - 208501.*pnt_flux[8])
+pnt_flux[8] *( 22658.*pnt_flux[8]);

//    if (flt_is0>pnt_config->flt_is_maximum)
//    	pnt_config->flt_is_maximum=flt_is0;
//    if (flt_is1>pnt_config->flt_is_maximum)
//    	pnt_config->flt_is_maximum=flt_is1;
//    if (flt_is2>pnt_config->flt_is_maximum)
//    	pnt_config->flt_is_maximum=flt_is2;
//    if (flt_is3>pnt_config->flt_is_maximum)
//    	pnt_config->flt_is_maximum=flt_is3;
//    if (flt_is4>pnt_config->flt_is_maximum)
//    	pnt_config->flt_is_maximum=flt_is4;
//
//    if (flt_is0<pnt_config->flt_is_minimum)
//    	pnt_config->flt_is_minimum=flt_is0;
//    if (flt_is1<pnt_config->flt_is_minimum)
//    	pnt_config->flt_is_minimum=flt_is1;
//    if (flt_is2<pnt_config->flt_is_minimum)
//    	pnt_config->flt_is_minimum=flt_is2;
//    if (flt_is3<pnt_config->flt_is_minimum)
//    	pnt_config->flt_is_minimum=flt_is3;
//    if (flt_is4<pnt_config->flt_is_minimum)
//    	pnt_config->flt_is_minimum=flt_is4;
//
//    pnt_config->flt_is_avrg=pnt_config->flt_is_avrg+flt_is0+flt_is1+flt_is2+flt_is3+flt_is4;
//    pnt_config->flt_is_avrg_counter+=5;


//    Berechnung der Gewichtungsfaktoren flt_omegaTmp und darauf aufbauend die resultierenden finalen Gewichtungsfaktoren flt_omega
//    flt_omegaTmp0 = pnt_config->flt_wenoOptimalerKoeffizient[0]/pow((flt_is0 + pnt_config->flt_wenoEpsilon),pnt_config->flt_wenoP);
//    flt_omegaTmp1 = pnt_config->flt_wenoOptimalerKoeffizient[1]/pow((flt_is1 + pnt_config->flt_wenoEpsilon),pnt_config->flt_wenoP);
//    flt_omegaTmp2 = pnt_config->flt_wenoOptimalerKoeffizient[2]/pow((flt_is2 + pnt_config->flt_wenoEpsilon),pnt_config->flt_wenoP);
//    flt_omegaTmp3 = pnt_config->flt_wenoOptimalerKoeffizient[3]/pow((flt_is3 + pnt_config->flt_wenoEpsilon),pnt_config->flt_wenoP);
//    flt_omegaTmp4 = pnt_config->flt_wenoOptimalerKoeffizient[4]/pow((flt_is4 + pnt_config->flt_wenoEpsilon),pnt_config->flt_wenoP);

    flt_omegaTmp0 = pnt_config->flt_wenoOptimalerKoeffizient_W9[0]/((flt_is0 + pnt_config->flt_wenoEpsilon)*(flt_is0 + pnt_config->flt_wenoEpsilon));
    flt_omegaTmp1 = pnt_config->flt_wenoOptimalerKoeffizient_W9[1]/((flt_is1 + pnt_config->flt_wenoEpsilon)*(flt_is1 + pnt_config->flt_wenoEpsilon));
    flt_omegaTmp2 = pnt_config->flt_wenoOptimalerKoeffizient_W9[2]/((flt_is2 + pnt_config->flt_wenoEpsilon)*(flt_is2 + pnt_config->flt_wenoEpsilon));
    flt_omegaTmp3 = pnt_config->flt_wenoOptimalerKoeffizient_W9[3]/((flt_is3 + pnt_config->flt_wenoEpsilon)*(flt_is3 + pnt_config->flt_wenoEpsilon));
    flt_omegaTmp4 = pnt_config->flt_wenoOptimalerKoeffizient_W9[4]/((flt_is4 + pnt_config->flt_wenoEpsilon)*(flt_is4 + pnt_config->flt_wenoEpsilon));

    flt_omega0 = flt_omegaTmp0/(flt_omegaTmp0 + flt_omegaTmp1 + flt_omegaTmp2 + flt_omegaTmp3 + flt_omegaTmp4);
    flt_omega2 = flt_omegaTmp2/(flt_omegaTmp0 + flt_omegaTmp1 + flt_omegaTmp2 + flt_omegaTmp3 + flt_omegaTmp4);
    flt_omega3 = flt_omegaTmp3/(flt_omegaTmp0 + flt_omegaTmp1 + flt_omegaTmp2 + flt_omegaTmp3 + flt_omegaTmp4);
    flt_omega4 = flt_omegaTmp4/(flt_omegaTmp0 + flt_omegaTmp1 + flt_omegaTmp2 + flt_omegaTmp3 + flt_omegaTmp4);


//    Berechnung der kombinierten Flüsse (dient nur der Übersicht)
    flt_q1 = pnt_deltaFlux[0] - 4.* pnt_deltaFlux[1] + 6.* pnt_deltaFlux[2] - 4.* pnt_deltaFlux[3] + pnt_deltaFlux[4];
    flt_q2 = pnt_deltaFlux[1] - 4.* pnt_deltaFlux[2] + 6.* pnt_deltaFlux[3] - 4.* pnt_deltaFlux[4] + pnt_deltaFlux[5];
    flt_q3 = pnt_deltaFlux[2] - 4.* pnt_deltaFlux[3] + 6.* pnt_deltaFlux[4] - 4.* pnt_deltaFlux[5] + pnt_deltaFlux[6];
    flt_q4 = pnt_deltaFlux[3] - 4.* pnt_deltaFlux[4] + 6.* pnt_deltaFlux[5] - 4.* pnt_deltaFlux[6] + pnt_deltaFlux[7];



    flt_result =
              1./5.*flt_omega0*flt_q1
            +1./20.*flt_omega2*flt_q2
            +1./60.*flt_omega3*(3.*flt_q2 - 2.*flt_q3)
            +1./60.*flt_omega4*(3.*flt_q2 - 2.*flt_q3 + 3.*flt_q4)
            -1./840.*(39.*flt_q2 -14.*flt_q3 + 3.*flt_q4);

    return (flt_result);
}

float Phi_Function_W5(
        struct strct_configuration * pnt_config,
        float * pnt_flux,
        float * pnt_deltaFlux)
{
    float flt_is0, flt_is1, flt_is2;
    float flt_omegaTmp0, flt_omegaTmp1, flt_omegaTmp2;
    float flt_omega0, flt_omega2;
    float flt_q1, flt_q2;
    float flt_result;


//    Alle Glaettungsindikatoren werden eigentlich noch durch 12 dividiert. Da aber nur die Gewichtungskoeffizienten
//    von Interesse sind, wo sich die 12 wieder herauskuerzt, kann man diese weglassen.
//    Die Verzerrung durch das epsilon wird hingenommen

//    //  Berechnung der Glättungsindikatoren basierend auf FLUESSEN
//    flt_is0 = 13./12.*
//    		(pnt_flux[0]-2.*pnt_flux[1]+   pnt_flux[2])*(pnt_flux[0]-2.*pnt_flux[1]+pnt_flux[2])
//    		+ 1./4.*
//    		(pnt_flux[0]-4.*pnt_flux[1]+3.*pnt_flux[2])*(pnt_flux[0]-4.*pnt_flux[1]+3.*pnt_flux[2]);
//
//    flt_is1 = 13./12.*
//    		(pnt_flux[1]-2.*pnt_flux[2]+   pnt_flux[3])*(pnt_flux[1]-2.*pnt_flux[2]+pnt_flux[3])
//    		+ 1./4. *
//    		(pnt_flux[1]-pnt_flux[3])*(pnt_flux[1]-pnt_flux[3]);
//
//    flt_is2 = 13./12.*
//    		(pnt_flux[2]-2.*pnt_flux[3]+   pnt_flux[4])*(pnt_flux[2]-2.*pnt_flux[3]+pnt_flux[4])
//    		+ 1./4. *
//    		(3.*pnt_flux[2]-4.*pnt_flux[3]+pnt_flux[4])*(3.*pnt_flux[2]-4.*pnt_flux[3]+pnt_flux[4]);

    //  Berechnung der Glättungsindikatoren basierend auf DELTA-FLUESSEN
    flt_is0 = 13.*(pnt_deltaFlux[0] - pnt_deltaFlux[1])*(pnt_deltaFlux[0] - pnt_deltaFlux[1])
    		+ 3.*(pnt_deltaFlux[0] - 3.*pnt_deltaFlux[1])*(pnt_deltaFlux[0] - 3.*pnt_deltaFlux[1]);

    flt_is1 = 13.*(pnt_deltaFlux[1] - pnt_deltaFlux[2])*(pnt_deltaFlux[1] - pnt_deltaFlux[2])
    		+ 3.*(pnt_deltaFlux[1] + pnt_deltaFlux[2])*(pnt_deltaFlux[1] + pnt_deltaFlux[2]);

    flt_is2 = 13.*(pnt_deltaFlux[2] - pnt_deltaFlux[3])*(pnt_deltaFlux[2] - pnt_deltaFlux[3])
    		+ 3.*(3.*pnt_deltaFlux[2] - pnt_deltaFlux[3])*(3.*pnt_deltaFlux[2] - pnt_deltaFlux[3]);


//    if (flt_is0>pnt_config->flt_is_maximum)
//    	pnt_config->flt_is_maximum=flt_is0;
//    if (flt_is1>pnt_config->flt_is_maximum)
//    	pnt_config->flt_is_maximum=flt_is1;
//    if (flt_is2>pnt_config->flt_is_maximum)
//    	pnt_config->flt_is_maximum=flt_is2;
//
//    if (flt_is0<pnt_config->flt_is_minimum)
//    	pnt_config->flt_is_minimum=flt_is0;
//    if (flt_is1<pnt_config->flt_is_minimum)
//    	pnt_config->flt_is_minimum=flt_is1;
//    if (flt_is2<pnt_config->flt_is_minimum)
//    	pnt_config->flt_is_minimum=flt_is2;
//
//    pnt_config->flt_is_avrg=pnt_config->flt_is_avrg+flt_is0+flt_is1+flt_is2;
//    pnt_config->flt_is_avrg_counter+=3;

//    Berechnung der Gewichtungsfaktoren flt_omegaTmp und darauf aufbauend die resultierenden finalen Gewichtungsfaktoren flt_omega
//    flt_omegaTmp0 = pnt_config->flt_wenoOptimalerKoeffizient[0]/pow((flt_is0 + pnt_config->flt_wenoEpsilon),pnt_config->flt_wenoP);
//    flt_omegaTmp1 = pnt_config->flt_wenoOptimalerKoeffizient[1]/pow((flt_is1 + pnt_config->flt_wenoEpsilon),pnt_config->flt_wenoP);


//    Aus Optimierungsgruenden wurden die optimalen Koeffizienten direkt eingesetzt.
//    Der Nenner wurde aus selbigem Grund wie oben weggelassen (kuezrt sich fuer omega-berechnung eh raus).
    flt_omegaTmp0 = 1./((flt_is0 + pnt_config->flt_wenoEpsilon)*(flt_is0 + pnt_config->flt_wenoEpsilon));
    flt_omegaTmp1 = 6./((flt_is1 + pnt_config->flt_wenoEpsilon)*(flt_is1 + pnt_config->flt_wenoEpsilon));
    flt_omegaTmp2 = 3./((flt_is2 + pnt_config->flt_wenoEpsilon)*(flt_is2 + pnt_config->flt_wenoEpsilon));

    flt_omega0 = flt_omegaTmp0/(flt_omegaTmp0 + flt_omegaTmp1 + flt_omegaTmp2);
    flt_omega2 = flt_omegaTmp2/(flt_omegaTmp0 + flt_omegaTmp1 + flt_omegaTmp2);


//    Berechnung der kombinierten Flüsse (dient nur der Übersicht)
    flt_q1 = pnt_deltaFlux[0] - 2.* pnt_deltaFlux[1] + pnt_deltaFlux[2];
    flt_q2 = pnt_deltaFlux[1] - 2.* pnt_deltaFlux[2] + pnt_deltaFlux[3];



    flt_result =
              1./3.*flt_omega0*flt_q1
            +1./6.*(flt_omega2-0.5)*flt_q2;

    return (flt_result);
}

void CalcFluxesInXiDirection(
        struct strct_configuration * pnt_config,
        struct strct_mesh * pnt_mesh,
        struct strct_U * pnt_U_RK,
        struct strct_Flux * pnt_Flux,
        struct strct_Flux * pnt_Flux_PlusHalf,
        struct strct_Flux * pnt_Q)
{
    int i,j,k,ijk,ijk_StencilPlus,ijk_StencilMinus;
    int iPlus1jk,iPlus2jk;
    int iMinus1jk;
#if SPACEORDER==9
    int iPlus3jk,iPlus4jk,iMinus2jk,iMinus3jk;;
#endif

    int int_c;
    int int_iStencilStart,int_iStencilEnd,int_iStencil,int_iStencilPlus,int_iStencilMinus;

    float lambdaMax;

    float flt_jacobian_PlusHalf;


    float * flt_fluxPlus;
    float * flt_fluxMinus;
    float * flt_deltaFluxPlus;
    float * flt_deltaFluxMinus;

    float flt_Phi_FunctionPlus;
    float flt_Phi_FunctionMinus;
    float * flt_Phi_FunctionSum;
    float * flt_Phi_FunctionSumRightEigen;
    float * flt_Theta_Function;

    flt_fluxPlus = (float *)calloc(pnt_config->int_SpaceOrder, sizeof(float));
    flt_fluxMinus = (float *)calloc(pnt_config->int_SpaceOrder, sizeof(float));
    flt_deltaFluxPlus = (float *)calloc((pnt_config->int_SpaceOrder-1), sizeof(float));
    flt_deltaFluxMinus = (float *)calloc((pnt_config->int_SpaceOrder-1), sizeof(float));
    flt_Phi_FunctionSum = (float *)calloc(pnt_config->int_conservationEquations, sizeof(float));
    flt_Phi_FunctionSumRightEigen = (float *)calloc(pnt_config->int_conservationEquations, sizeof(float));
    flt_Theta_Function = (float *)calloc(pnt_config->int_conservationEquations, sizeof(float));

//    Berechnung der Flüsse in den Gitterpunkten
//    Die Berechnung der Gitterpunkte an den i+1/2 Stellen, die benötigt werden
//    geschieht später
    for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
    {
        for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
        {
            for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
            {
 ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

                pnt_Flux->Mass[ijk]=
                        pnt_mesh->jacobian[ijk]*(
 pnt_U_RK->rho[ijk]*pnt_U_RK->theta1[ijk]);

                pnt_Flux->xiMomentum[ijk]=
                        pnt_mesh->jacobian[ijk]*(
 pnt_U_RK->rho[ijk]*pnt_U_RK->theta1[ijk]*pnt_U_RK->u[ijk]*pnt_mesh->BC_Corrector[ijk]+
 pnt_mesh->xi_x[ijk]*pnt_U_RK->p[ijk]*
                                pnt_config->flt_Upsilon);

                pnt_Flux->etaMomentum[ijk]=
                        pnt_mesh->jacobian[ijk]*(
 pnt_U_RK->rho[ijk]*pnt_U_RK->theta1[ijk]*pnt_U_RK->v[ijk]*pnt_mesh->BC_Corrector[ijk]+
 pnt_mesh->xi_y[ijk]*pnt_U_RK->p[ijk]*
                                pnt_config->flt_Upsilon);

                pnt_Flux->zetaMomentum[ijk]=
                        pnt_mesh->jacobian[ijk]*(
 pnt_U_RK->rho[ijk]*pnt_U_RK->theta1[ijk]*pnt_U_RK->w[ijk]*pnt_mesh->BC_Corrector[ijk]+
 pnt_mesh->xi_z[ijk]*pnt_U_RK->p[ijk]*
                                pnt_config->flt_Upsilon);

                pnt_Flux->Energy[ijk]=
                        pnt_mesh->jacobian[ijk]*(
 pnt_U_RK->rho[ijk]*pnt_U_RK->theta1[ijk]*(
                                        pnt_U_RK->e[ijk]+
 pnt_U_RK->p[ijk]/pnt_U_RK->rho[ijk]*
 pnt_config->flt_Upsilon));

            }
        }
    }





    for (i=pnt_config->int_iStartReal-1; i <= pnt_config->int_iEndReal; i++)
    {
        for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
        {
            for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
            {
 ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
 iPlus1jk=(i+1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
 iPlus2jk=(i+2)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
 iMinus1jk=(i-1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

#if SPACEORDER==9
 iPlus3jk=(i+3)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
 iPlus4jk=(i+4)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
 iMinus2jk=(i-2)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
 iMinus3jk=(i-3)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
#endif

 flt_jacobian_PlusHalf=(pnt_mesh->jacobian[ijk]+pnt_mesh->jacobian[iPlus1jk])/2.0;

//                Berechne Eigenvektor für die aktuelle i,j,k Position
                CalcEigenVectorsInXiDirection(
                            pnt_config,
                            pnt_mesh,
                            pnt_U_RK,
                            ijk,
                            iPlus1jk,
                            iMinus1jk);


//                Berechne für aktuellen Stencilbereich LambdaMax
                lambdaMax=GetLambdaMaxInXiDirection(
                        pnt_config,
                        pnt_mesh,
                        pnt_U_RK,
                        i,j,k);


                for (int_c=0; int_c < pnt_config->int_conservationEquations; int_c++)
                {
//                    Definiere aktuellen Stencilbereich für die Berechnung von LF
 int_iStencilStart=i-(pnt_config->int_SpaceOrder-1)/2;
 int_iStencilEnd=i+(pnt_config->int_SpaceOrder-1)/2;

                    for (int_iStencil=int_iStencilStart; int_iStencil <= int_iStencilEnd; int_iStencil++)
                    {
                        int_iStencilPlus=int_iStencil;
 int_iStencilMinus=int_iStencilEnd-(int_iStencil-int_iStencilStart)+1;

 ijk_StencilPlus=int_iStencilPlus*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
 ijk_StencilMinus=int_iStencilMinus*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;


 flt_fluxPlus[int_iStencil-int_iStencilStart]=
                                0.5*
                                (
                                flt_leftEigenvector[int_c][0]*
 (pnt_Flux->Mass[ijk_StencilPlus]+lambdaMax*pnt_U_RK->rho[ijk_StencilPlus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilPlus])
                                +
                                flt_leftEigenvector[int_c][1]*
 (pnt_Flux->xiMomentum[ijk_StencilPlus]+lambdaMax*pnt_U_RK->rho[ijk_StencilPlus]*pnt_U_RK->u[ijk_StencilPlus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilPlus])
                                +
                                flt_leftEigenvector[int_c][2]*
 (pnt_Flux->etaMomentum[ijk_StencilPlus]+lambdaMax*pnt_U_RK->rho[ijk_StencilPlus]*pnt_U_RK->v[ijk_StencilPlus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilPlus])
                                +
                                flt_leftEigenvector[int_c][3]*
 (pnt_Flux->zetaMomentum[ijk_StencilPlus]+lambdaMax*pnt_U_RK->rho[ijk_StencilPlus]*pnt_U_RK->w[ijk_StencilPlus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilPlus])
                                +
                                flt_leftEigenvector[int_c][4]*
 (pnt_Flux->Energy[ijk_StencilPlus]+lambdaMax*pnt_U_RK->rho[ijk_StencilPlus]*pnt_U_RK->e[ijk_StencilPlus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilPlus])
                                );
 flt_fluxMinus[int_iStencil-int_iStencilStart]=
                                0.5*
                                (
                                flt_leftEigenvector[int_c][0]*
 (pnt_Flux->Mass[ijk_StencilMinus]-lambdaMax*pnt_U_RK->rho[ijk_StencilMinus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilMinus])
                                +
                                flt_leftEigenvector[int_c][1]*
 (pnt_Flux->xiMomentum[ijk_StencilMinus]-lambdaMax*pnt_U_RK->rho[ijk_StencilMinus]*pnt_U_RK->u[ijk_StencilMinus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilMinus])
                                +
                                flt_leftEigenvector[int_c][2]*
 (pnt_Flux->etaMomentum[ijk_StencilMinus]-lambdaMax*pnt_U_RK->rho[ijk_StencilMinus]*pnt_U_RK->v[ijk_StencilMinus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilMinus])
                                +
                                flt_leftEigenvector[int_c][3]*
 (pnt_Flux->zetaMomentum[ijk_StencilMinus]-lambdaMax*pnt_U_RK->rho[ijk_StencilMinus]*pnt_U_RK->w[ijk_StencilMinus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilMinus])
                                +
                                flt_leftEigenvector[int_c][4]*
 (pnt_Flux->Energy[ijk_StencilMinus]-lambdaMax*pnt_U_RK->rho[ijk_StencilMinus]*pnt_U_RK->e[ijk_StencilMinus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilMinus])
                                );
                    }

//                    Berechnung von DELTA_LF
                    for (int_iStencil=0; int_iStencil <= pnt_config->int_SpaceOrder-2; int_iStencil++)
                    {
 flt_deltaFluxPlus[int_iStencil]=flt_fluxPlus[int_iStencil+1]-flt_fluxPlus[int_iStencil];
 flt_deltaFluxMinus[int_iStencil]=flt_fluxMinus[int_iStencil]-flt_fluxMinus[int_iStencil+1];
                    }


#if SPACEORDER==9
                    //    Aufruf der Phi-Funktion fuer positive Eigenwerte
                    flt_Phi_FunctionPlus=Phi_Function_W9(
                            pnt_config,
                            flt_fluxPlus,
                            flt_deltaFluxPlus);
                    //    Aufruf der Phi-Funktion fuer negative Eigenwerte
                    flt_Phi_FunctionMinus=Phi_Function_W9(
                            pnt_config,
                            flt_fluxMinus,
                            flt_deltaFluxMinus);
#endif
#if SPACEORDER==5
					//    Aufruf der Phi-Funktion fuer positive Eigenwerte
					flt_Phi_FunctionPlus=Phi_Function_W5(
							pnt_config,
							flt_fluxPlus,
							flt_deltaFluxPlus);
					//    Aufruf der Phi-Funktion fuer negative Eigenwerte
					flt_Phi_FunctionMinus=Phi_Function_W5(
							pnt_config,
							flt_fluxMinus,
							flt_deltaFluxMinus);
#endif
 flt_Phi_FunctionSum[int_c]=(-flt_Phi_FunctionPlus+flt_Phi_FunctionMinus);

                }
                for (int_c=0; int_c < pnt_config->int_conservationEquations; int_c++)
                {
                    flt_Phi_FunctionSumRightEigen[int_c]=
 flt_Phi_FunctionSum[0]*flt_rightEigenvector[int_c][0]+
 flt_Phi_FunctionSum[1]*flt_rightEigenvector[int_c][1]+
 flt_Phi_FunctionSum[2]*flt_rightEigenvector[int_c][2]+
 flt_Phi_FunctionSum[3]*flt_rightEigenvector[int_c][3]+
 flt_Phi_FunctionSum[4]*flt_rightEigenvector[int_c][4];
                }

#if SPACEORDER==9
                //    Berechnung des Zentralterms über die Theta-Funktion
                Theta_Function_W9(
                    pnt_Flux,
                    flt_Theta_Function,
                    ijk,
                    iPlus1jk,
                    iPlus2jk,
                    iPlus3jk,
                    iPlus4jk,
                    iMinus1jk,
                    iMinus2jk,
                    iMinus3jk);
#endif
#if SPACEORDER==5
				//    Berechnung des Zentralterms über die Theta-Funktion
				Theta_Function_W5(
					pnt_Flux,
					flt_Theta_Function,
					ijk,
					iPlus1jk,
					iPlus2jk,
					iMinus1jk);
#endif

                //    Berechnung des resultierenden Fluss an der Stelle i+1/2 (entspricht i)
                //    Theta=Zentralterm
                //    Phi=WENO-Terme
 pnt_Flux_PlusHalf->Mass[ijk]=flt_Phi_FunctionSumRightEigen[0]+flt_Theta_Function[0];
 pnt_Flux_PlusHalf->xiMomentum[ijk]=flt_Phi_FunctionSumRightEigen[1]+flt_Theta_Function[1];
 pnt_Flux_PlusHalf->etaMomentum[ijk]=flt_Phi_FunctionSumRightEigen[2]+flt_Theta_Function[2];
 pnt_Flux_PlusHalf->zetaMomentum[ijk]=flt_Phi_FunctionSumRightEigen[3]+flt_Theta_Function[3];
 pnt_Flux_PlusHalf->Energy[ijk]=flt_Phi_FunctionSumRightEigen[4]+flt_Theta_Function[4];
            }
        }
    }

    for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
    {
        for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
        {
            for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
            {
 ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
 iMinus1jk=(i-1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

 pnt_Q->Mass[ijk]=pnt_Q->Mass[ijk]-(pnt_Flux_PlusHalf->Mass[ijk]-pnt_Flux_PlusHalf->Mass[iMinus1jk]);
 pnt_Q->xiMomentum[ijk]=pnt_Q->xiMomentum[ijk]-(pnt_Flux_PlusHalf->xiMomentum[ijk]-pnt_Flux_PlusHalf->xiMomentum[iMinus1jk]);
 pnt_Q->etaMomentum[ijk]=pnt_Q->etaMomentum[ijk]-(pnt_Flux_PlusHalf->etaMomentum[ijk]-pnt_Flux_PlusHalf->etaMomentum[iMinus1jk]);
 pnt_Q->zetaMomentum[ijk]=pnt_Q->zetaMomentum[ijk]-(pnt_Flux_PlusHalf->zetaMomentum[ijk]-pnt_Flux_PlusHalf->zetaMomentum[iMinus1jk]);
 pnt_Q->Energy[ijk]=pnt_Q->Energy[ijk]-(pnt_Flux_PlusHalf->Energy[ijk]-pnt_Flux_PlusHalf->Energy[iMinus1jk]);
            }
        }
    }

    free(flt_fluxPlus);
    free(flt_fluxMinus);
    free(flt_deltaFluxPlus);
    free(flt_deltaFluxMinus);
    free(flt_Phi_FunctionSum);
    free(flt_Phi_FunctionSumRightEigen);
    free(flt_Theta_Function);
}

void CalcFluxesInEtaDirection(
        struct strct_configuration * pnt_config,
        struct strct_mesh * pnt_mesh,
        struct strct_U * pnt_U_RK,
        struct strct_Flux * pnt_Flux,
        struct strct_Flux * pnt_Flux_PlusHalf,
        struct strct_Flux * pnt_Q)
{
    int i,j,k,ijk,ijk_StencilPlus,ijk_StencilMinus;
    int ijPlus1k,ijPlus2k;
    int ijMinus1k;
#if SPACEORDER==9
    int ijPlus3k,ijPlus4k,ijMinus2k,ijMinus3k;;
#endif
    int int_c;
    int int_jStencilStart,int_jStencilEnd,int_jStencil,int_jStencilPlus,int_jStencilMinus;

    float lambdaMax;

    float flt_jacobian_PlusHalf;

    float * flt_fluxPlus;
    float * flt_fluxMinus;
    float * flt_deltaFluxPlus;
    float * flt_deltaFluxMinus;

    float flt_Phi_FunctionPlus;
    float flt_Phi_FunctionMinus;
    float * flt_Phi_FunctionSum;
    float * flt_Phi_FunctionSumRightEigen;
    float * flt_Theta_Function;

    flt_fluxPlus = (float *)calloc(pnt_config->int_SpaceOrder, sizeof(float));
    flt_fluxMinus = (float *)calloc(pnt_config->int_SpaceOrder, sizeof(float));
    flt_deltaFluxPlus = (float *)calloc((pnt_config->int_SpaceOrder-1), sizeof(float));
    flt_deltaFluxMinus = (float *)calloc((pnt_config->int_SpaceOrder-1), sizeof(float));
    flt_Phi_FunctionSum = (float *)calloc(pnt_config->int_conservationEquations, sizeof(float));
    flt_Phi_FunctionSumRightEigen = (float *)calloc(pnt_config->int_conservationEquations, sizeof(float));
    flt_Theta_Function = (float *)calloc(pnt_config->int_conservationEquations, sizeof(float));

//    Berechnung der Flüsse in den Gitterpunkten
//    Die Berechnung der Gitterpunkte an den j+1/2 Stellen, die benötigt werden
//    geschieht später
    for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
    {
        for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
        {
            for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
            {
 ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

                pnt_Flux->Mass[ijk]=
                        pnt_mesh->jacobian[ijk]*(
 pnt_U_RK->rho[ijk]*pnt_U_RK->theta2[ijk]);

                pnt_Flux->xiMomentum[ijk]=
                        pnt_mesh->jacobian[ijk]*(
 pnt_U_RK->rho[ijk]*pnt_U_RK->theta2[ijk]*pnt_U_RK->u[ijk]*pnt_mesh->BC_Corrector[ijk]+
 pnt_mesh->eta_x[ijk]*pnt_U_RK->p[ijk]*
                                pnt_config->flt_Upsilon);

                pnt_Flux->etaMomentum[ijk]=
                        pnt_mesh->jacobian[ijk]*(
 pnt_U_RK->rho[ijk]*pnt_U_RK->theta2[ijk]*pnt_U_RK->v[ijk]*pnt_mesh->BC_Corrector[ijk]+
 pnt_mesh->eta_y[ijk]*pnt_U_RK->p[ijk]*
                                pnt_config->flt_Upsilon);

                pnt_Flux->zetaMomentum[ijk]=
                        pnt_mesh->jacobian[ijk]*(
 pnt_U_RK->rho[ijk]*pnt_U_RK->theta2[ijk]*pnt_U_RK->w[ijk]*pnt_mesh->BC_Corrector[ijk]+
 pnt_mesh->eta_z[ijk]*pnt_U_RK->p[ijk]*
                                pnt_config->flt_Upsilon);

                pnt_Flux->Energy[ijk]=
                        pnt_mesh->jacobian[ijk]*(
 pnt_U_RK->rho[ijk]*pnt_U_RK->theta2[ijk]*(
                                        pnt_U_RK->e[ijk]+
 pnt_U_RK->p[ijk]/pnt_U_RK->rho[ijk]*
 pnt_config->flt_Upsilon));

            }
        }
    }




    for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
    {
        for (j=pnt_config->int_jStartReal-1; j <= pnt_config->int_jEndReal; j++)
        {
            for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
            {
 ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
 ijPlus1k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+1)*pnt_config->int_kMeshPointsGhostCells+k;
 ijPlus2k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+2)*pnt_config->int_kMeshPointsGhostCells+k;
 ijMinus1k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-1)*pnt_config->int_kMeshPointsGhostCells+k;

#if SPACEORDER==9
 ijPlus3k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+3)*pnt_config->int_kMeshPointsGhostCells+k;
 ijPlus4k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+4)*pnt_config->int_kMeshPointsGhostCells+k;
 ijMinus2k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-2)*pnt_config->int_kMeshPointsGhostCells+k;
 ijMinus3k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-3)*pnt_config->int_kMeshPointsGhostCells+k;
#endif
 flt_jacobian_PlusHalf=(pnt_mesh->jacobian[ijk]+pnt_mesh->jacobian[ijPlus1k])/2.0;

//                Berechne Eigenvektor für die aktuelle i,j,k Position
                CalcEigenVectorsInEtaDirection(
                            pnt_config,
                            pnt_mesh,
                            pnt_U_RK,
                            ijk,
                            ijPlus1k,
                            ijMinus1k);


//                Berechne für aktuellen Stencilbereich LambdaMax
                lambdaMax=GetLambdaMaxInEtaDirection(
                        pnt_config,
                        pnt_mesh,
                        pnt_U_RK,
                        i,j,k);


                for (int_c=0; int_c < pnt_config->int_conservationEquations; int_c++)
                {
//                    Definiere aktuellen Stencilbereich für die Berechnung von LF
 int_jStencilStart=j-(pnt_config->int_SpaceOrder-1)/2;
 int_jStencilEnd=j+(pnt_config->int_SpaceOrder-1)/2;
                    for (int_jStencil=int_jStencilStart; int_jStencil <= int_jStencilEnd; int_jStencil++)
                    {
                        int_jStencilPlus=int_jStencil;
 int_jStencilMinus=int_jStencilEnd-(int_jStencil-int_jStencilStart)+1;

 ijk_StencilPlus=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+int_jStencilPlus*pnt_config->int_kMeshPointsGhostCells+k;
 ijk_StencilMinus=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+int_jStencilMinus*pnt_config->int_kMeshPointsGhostCells+k;


 flt_fluxPlus[int_jStencil-int_jStencilStart]=
                                0.5*
                                (
                                flt_leftEigenvector[int_c][0]*
 (pnt_Flux->Mass[ijk_StencilPlus]+lambdaMax*pnt_U_RK->rho[ijk_StencilPlus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilPlus])
                                +
                                flt_leftEigenvector[int_c][1]*
 (pnt_Flux->xiMomentum[ijk_StencilPlus]+lambdaMax*pnt_U_RK->rho[ijk_StencilPlus]*pnt_U_RK->u[ijk_StencilPlus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilPlus])
                                +
                                flt_leftEigenvector[int_c][2]*
 (pnt_Flux->etaMomentum[ijk_StencilPlus]+lambdaMax*pnt_U_RK->rho[ijk_StencilPlus]*pnt_U_RK->v[ijk_StencilPlus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilPlus])
                                +
                                flt_leftEigenvector[int_c][3]*
 (pnt_Flux->zetaMomentum[ijk_StencilPlus]+lambdaMax*pnt_U_RK->rho[ijk_StencilPlus]*pnt_U_RK->w[ijk_StencilPlus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilPlus])
                                +
                                flt_leftEigenvector[int_c][4]*
 (pnt_Flux->Energy[ijk_StencilPlus]+lambdaMax*pnt_U_RK->rho[ijk_StencilPlus]*pnt_U_RK->e[ijk_StencilPlus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilPlus])
                                );
 flt_fluxMinus[int_jStencil-int_jStencilStart]=
                                0.5*
                                (
                                flt_leftEigenvector[int_c][0]*
 (pnt_Flux->Mass[ijk_StencilMinus]-lambdaMax*pnt_U_RK->rho[ijk_StencilMinus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilMinus])
                                +
                                flt_leftEigenvector[int_c][1]*
 (pnt_Flux->xiMomentum[ijk_StencilMinus]-lambdaMax*pnt_U_RK->rho[ijk_StencilMinus]*pnt_U_RK->u[ijk_StencilMinus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilMinus])
                                +
                                flt_leftEigenvector[int_c][2]*
 (pnt_Flux->etaMomentum[ijk_StencilMinus]-lambdaMax*pnt_U_RK->rho[ijk_StencilMinus]*pnt_U_RK->v[ijk_StencilMinus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilMinus])
                                +
                                flt_leftEigenvector[int_c][3]*
 (pnt_Flux->zetaMomentum[ijk_StencilMinus]-lambdaMax*pnt_U_RK->rho[ijk_StencilMinus]*pnt_U_RK->w[ijk_StencilMinus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilMinus])
                                +
                                flt_leftEigenvector[int_c][4]*
 (pnt_Flux->Energy[ijk_StencilMinus]-lambdaMax*pnt_U_RK->rho[ijk_StencilMinus]*pnt_U_RK->e[ijk_StencilMinus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilMinus])
                                );
                    }

//                    Berechnung von DELTA_LF
                    for (int_jStencil=0; int_jStencil <= pnt_config->int_SpaceOrder-2; int_jStencil++)
                    {
 flt_deltaFluxPlus[int_jStencil]=flt_fluxPlus[int_jStencil+1]-flt_fluxPlus[int_jStencil];
 flt_deltaFluxMinus[int_jStencil]=flt_fluxMinus[int_jStencil]-flt_fluxMinus[int_jStencil+1];
                    }


#if SPACEORDER==9
                    //    Aufruf der Phi-Funktion fuer positive Eigenwerte
                    flt_Phi_FunctionPlus=Phi_Function_W9(
                            pnt_config,
                            flt_fluxPlus,
                            flt_deltaFluxPlus);
                    //    Aufruf der Phi-Funktion fuer negative Eigenwerte
                    flt_Phi_FunctionMinus=Phi_Function_W9(
                            pnt_config,
                            flt_fluxMinus,
                            flt_deltaFluxMinus);
#endif
#if SPACEORDER==5
					//    Aufruf der Phi-Funktion fuer positive Eigenwerte
					flt_Phi_FunctionPlus=Phi_Function_W5(
							pnt_config,
							flt_fluxPlus,
							flt_deltaFluxPlus);
					//    Aufruf der Phi-Funktion fuer negative Eigenwerte
					flt_Phi_FunctionMinus=Phi_Function_W5(
							pnt_config,
							flt_fluxMinus,
							flt_deltaFluxMinus);
#endif
 flt_Phi_FunctionSum[int_c]=(-flt_Phi_FunctionPlus+flt_Phi_FunctionMinus);

                }
                for (int_c=0; int_c < pnt_config->int_conservationEquations; int_c++)
                {
                    flt_Phi_FunctionSumRightEigen[int_c]=
 flt_Phi_FunctionSum[0]*flt_rightEigenvector[int_c][0]+
 flt_Phi_FunctionSum[1]*flt_rightEigenvector[int_c][1]+
 flt_Phi_FunctionSum[2]*flt_rightEigenvector[int_c][2]+
 flt_Phi_FunctionSum[3]*flt_rightEigenvector[int_c][3]+
 flt_Phi_FunctionSum[4]*flt_rightEigenvector[int_c][4];
                }

#if SPACEORDER==9
                //    Berechnung des Zentralterms über die Theta-Funktion
                Theta_Function_W9(
                    pnt_Flux,
                    flt_Theta_Function,
                    ijk,
                    ijPlus1k,
                    ijPlus2k,
                    ijPlus3k,
                    ijPlus4k,
                    ijMinus1k,
                    ijMinus2k,
                    ijMinus3k);
#endif
#if SPACEORDER==5
				//    Berechnung des Zentralterms über die Theta-Funktion
				Theta_Function_W5(
					pnt_Flux,
					flt_Theta_Function,
					ijk,
					ijPlus1k,
					ijPlus2k,
					ijMinus1k);
#endif

                //    Berechnung des resultierenden Fluss an der Stelle j+1/2 (entspricht j)
                //    Theta=Zentralterm
                //    Phi=WENO-Terme
 pnt_Flux_PlusHalf->Mass[ijk]=flt_Phi_FunctionSumRightEigen[0]+flt_Theta_Function[0];
 pnt_Flux_PlusHalf->xiMomentum[ijk]=flt_Phi_FunctionSumRightEigen[1]+flt_Theta_Function[1];
 pnt_Flux_PlusHalf->etaMomentum[ijk]=flt_Phi_FunctionSumRightEigen[2]+flt_Theta_Function[2];
 pnt_Flux_PlusHalf->zetaMomentum[ijk]=flt_Phi_FunctionSumRightEigen[3]+flt_Theta_Function[3];
 pnt_Flux_PlusHalf->Energy[ijk]=flt_Phi_FunctionSumRightEigen[4]+flt_Theta_Function[4];
            }
        }
    }

    for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
    {
        for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
        {
            for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
            {
 ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
 ijMinus1k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-1)*pnt_config->int_kMeshPointsGhostCells+k;

 pnt_Q->Mass[ijk]=pnt_Q->Mass[ijk]-(pnt_Flux_PlusHalf->Mass[ijk]-pnt_Flux_PlusHalf->Mass[ijMinus1k]);
 pnt_Q->xiMomentum[ijk]=pnt_Q->xiMomentum[ijk]-(pnt_Flux_PlusHalf->xiMomentum[ijk]-pnt_Flux_PlusHalf->xiMomentum[ijMinus1k]);
 pnt_Q->etaMomentum[ijk]=pnt_Q->etaMomentum[ijk]-(pnt_Flux_PlusHalf->etaMomentum[ijk]-pnt_Flux_PlusHalf->etaMomentum[ijMinus1k]);
 pnt_Q->zetaMomentum[ijk]=pnt_Q->zetaMomentum[ijk]-(pnt_Flux_PlusHalf->zetaMomentum[ijk]-pnt_Flux_PlusHalf->zetaMomentum[ijMinus1k]);
 pnt_Q->Energy[ijk]=pnt_Q->Energy[ijk]-(pnt_Flux_PlusHalf->Energy[ijk]-pnt_Flux_PlusHalf->Energy[ijMinus1k]);
            }
        }
    }

    free(flt_fluxPlus);
    free(flt_fluxMinus);
    free(flt_deltaFluxPlus);
    free(flt_deltaFluxMinus);
    free(flt_Phi_FunctionSum);
    free(flt_Phi_FunctionSumRightEigen);
    free(flt_Theta_Function);


//####################OUTPUT FUER FLUESSE AN WAND #########################
//	if(strcmp(pnt_config->BC_Bottom,pnt_config->BCWallViscous)==0)
//		{
//		i=16;
//	k=pnt_config->int_kStartReal;
//	printf("masse:n");
//	for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jStartReal+2; j++){
//	ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//	printf("%f ",pnt_Flux->Mass[ijk]);
//	}printf("\n");
//	printf("eta-impuls:n");
//	for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jStartReal+2; j++){
//	ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//	printf("%f ",pnt_Flux->etaMomentum[ijk]);
//	}printf("\n");
//		}
}

void CalcFluxesInZetaDirection(
        struct strct_configuration * pnt_config,
        struct strct_mesh * pnt_mesh,
        struct strct_U * pnt_U_RK,
        struct strct_Flux * pnt_Flux,
        struct strct_Flux * pnt_Flux_PlusHalf,
        struct strct_Flux * pnt_Q)
{
    int i,j,k,ijk,ijk_StencilPlus,ijk_StencilMinus;
    int ijkPlus1,ijkPlus2;
    int ijkMinus1;
#if SPACEORDER==9
    int ijkPlus3,ijkPlus4,ijkMinus2,ijkMinus3;;
#endif
    int int_c;
    int int_kStencilStart,int_kStencilEnd,int_kStencil,int_kStencilPlus,int_kStencilMinus;

    float lambdaMax;

    float flt_jacobian_PlusHalf;

    float * flt_fluxPlus;
    float * flt_fluxMinus;
    float * flt_deltaFluxPlus;
    float * flt_deltaFluxMinus;

    float flt_Phi_FunctionPlus;
    float flt_Phi_FunctionMinus;
    float * flt_Phi_FunctionSum;
    float * flt_Phi_FunctionSumRightEigen;
    float * flt_Theta_Function;

    flt_fluxPlus = (float *)calloc(pnt_config->int_SpaceOrder, sizeof(float));
    flt_fluxMinus = (float *)calloc(pnt_config->int_SpaceOrder, sizeof(float));
    flt_deltaFluxPlus = (float *)calloc((pnt_config->int_SpaceOrder-1), sizeof(float));
    flt_deltaFluxMinus = (float *)calloc((pnt_config->int_SpaceOrder-1), sizeof(float));
    flt_Phi_FunctionSum = (float *)calloc(pnt_config->int_conservationEquations, sizeof(float));
    flt_Phi_FunctionSumRightEigen = (float *)calloc(pnt_config->int_conservationEquations, sizeof(float));
    flt_Theta_Function = (float *)calloc(pnt_config->int_conservationEquations, sizeof(float));

//    Berechnung der Flüsse in den Gitterpunkten
//    Die Berechnung der Gitterpunkte an den j+1/2 Stellen, die benötigt werden
//    geschieht später
    for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
    {
        for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
        {
            for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
            {
 ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

                pnt_Flux->Mass[ijk]=
                        pnt_mesh->jacobian[ijk]*(
 pnt_U_RK->rho[ijk]*pnt_U_RK->theta3[ijk]);

                pnt_Flux->xiMomentum[ijk]=
                        pnt_mesh->jacobian[ijk]*(
 pnt_U_RK->rho[ijk]*pnt_U_RK->theta3[ijk]*pnt_U_RK->u[ijk]*pnt_mesh->BC_Corrector[ijk]+
 pnt_mesh->zeta_x[ijk]*pnt_U_RK->p[ijk]*
                                pnt_config->flt_Upsilon);

                pnt_Flux->etaMomentum[ijk]=
                        pnt_mesh->jacobian[ijk]*(
 pnt_U_RK->rho[ijk]*pnt_U_RK->theta3[ijk]*pnt_U_RK->v[ijk]*pnt_mesh->BC_Corrector[ijk]+
 pnt_mesh->zeta_y[ijk]*pnt_U_RK->p[ijk]*
                                pnt_config->flt_Upsilon);

                pnt_Flux->zetaMomentum[ijk]=
                        pnt_mesh->jacobian[ijk]*(
 pnt_U_RK->rho[ijk]*pnt_U_RK->theta3[ijk]*pnt_U_RK->w[ijk]*pnt_mesh->BC_Corrector[ijk]+
 pnt_mesh->zeta_z[ijk]*pnt_U_RK->p[ijk]*
                                pnt_config->flt_Upsilon);

                pnt_Flux->Energy[ijk]=
                        pnt_mesh->jacobian[ijk]*(
 pnt_U_RK->rho[ijk]*pnt_U_RK->theta3[ijk]*(
                                        pnt_U_RK->e[ijk]+
 pnt_U_RK->p[ijk]/pnt_U_RK->rho[ijk]*
 pnt_config->flt_Upsilon));

            }
        }
    }




    for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
    {
        for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
        {
            for (k=pnt_config->int_kStartReal-1; k <= pnt_config->int_kEndReal; k++)
            {
 ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
 ijkPlus1=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+1;
 ijkPlus2=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+2;
 ijkMinus1=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-1;

#if SPACEORDER==9
 ijkPlus3=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+3;
  ijkPlus4=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+4;
 ijkMinus2=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-2;
 ijkMinus3=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-3;
#endif
 flt_jacobian_PlusHalf=(pnt_mesh->jacobian[ijk]+pnt_mesh->jacobian[ijkPlus1])/2.0;

//                Berechne Eigenvektor für die aktuelle i,j,k Position
                CalcEigenVectorsInZetaDirection(
                            pnt_config,
                            pnt_mesh,
                            pnt_U_RK,
                            ijk,
                            ijkPlus1,
                            ijkMinus1);


//                Berechne für aktuellen Stencilbereich LambdaMax
                lambdaMax=GetLambdaMaxInZetaDirection(
                        pnt_config,
                        pnt_mesh,
                        pnt_U_RK,
                        i,j,k);


                for (int_c=0; int_c < pnt_config->int_conservationEquations; int_c++)
                {
//                    Definiere aktuellen Stencilbereich für die Berechnung von LF
 int_kStencilStart=k-(pnt_config->int_SpaceOrder-1)/2;
 int_kStencilEnd=k+(pnt_config->int_SpaceOrder-1)/2;
                    for (int_kStencil=int_kStencilStart; int_kStencil <= int_kStencilEnd; int_kStencil++)
                    {
                        int_kStencilPlus=int_kStencil;
 int_kStencilMinus=int_kStencilEnd-(int_kStencil-int_kStencilStart)+1;

 ijk_StencilPlus=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_kStencilPlus;
 ijk_StencilMinus=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_kStencilMinus;


 flt_fluxPlus[int_kStencil-int_kStencilStart]=
                                0.5*
                                (
                                flt_leftEigenvector[int_c][0]*
 (pnt_Flux->Mass[ijk_StencilPlus]+lambdaMax*pnt_U_RK->rho[ijk_StencilPlus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilPlus])
                                +
                                flt_leftEigenvector[int_c][1]*
 (pnt_Flux->xiMomentum[ijk_StencilPlus]+lambdaMax*pnt_U_RK->rho[ijk_StencilPlus]*pnt_U_RK->u[ijk_StencilPlus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilPlus])
                                +
                                flt_leftEigenvector[int_c][2]*
 (pnt_Flux->etaMomentum[ijk_StencilPlus]+lambdaMax*pnt_U_RK->rho[ijk_StencilPlus]*pnt_U_RK->v[ijk_StencilPlus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilPlus])
                                +
                                flt_leftEigenvector[int_c][3]*
 (pnt_Flux->zetaMomentum[ijk_StencilPlus]+lambdaMax*pnt_U_RK->rho[ijk_StencilPlus]*pnt_U_RK->w[ijk_StencilPlus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilPlus])
                                +
                                flt_leftEigenvector[int_c][4]*
 (pnt_Flux->Energy[ijk_StencilPlus]+lambdaMax*pnt_U_RK->rho[ijk_StencilPlus]*pnt_U_RK->e[ijk_StencilPlus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilPlus])
                                );
 flt_fluxMinus[int_kStencil-int_kStencilStart]=
                                0.5*
                                (
                                flt_leftEigenvector[int_c][0]*
 (pnt_Flux->Mass[ijk_StencilMinus]-lambdaMax*pnt_U_RK->rho[ijk_StencilMinus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilMinus])
                                +
                                flt_leftEigenvector[int_c][1]*
 (pnt_Flux->xiMomentum[ijk_StencilMinus]-lambdaMax*pnt_U_RK->rho[ijk_StencilMinus]*pnt_U_RK->u[ijk_StencilMinus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilMinus])
                                +
                                flt_leftEigenvector[int_c][2]*
 (pnt_Flux->etaMomentum[ijk_StencilMinus]-lambdaMax*pnt_U_RK->rho[ijk_StencilMinus]*pnt_U_RK->v[ijk_StencilMinus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilMinus])
                                +
                                flt_leftEigenvector[int_c][3]*
 (pnt_Flux->zetaMomentum[ijk_StencilMinus]-lambdaMax*pnt_U_RK->rho[ijk_StencilMinus]*pnt_U_RK->w[ijk_StencilMinus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilMinus])
                                +
                                flt_leftEigenvector[int_c][4]*
 (pnt_Flux->Energy[ijk_StencilMinus]-lambdaMax*pnt_U_RK->rho[ijk_StencilMinus]*pnt_U_RK->e[ijk_StencilMinus]*flt_jacobian_PlusHalf)//*pnt_mesh->jacobian[ijk_StencilMinus])
                                );
                    }

//                    Berechnung von DELTA_LF
                    for (int_kStencil=0; int_kStencil <= pnt_config->int_SpaceOrder-2; int_kStencil++)
                    {
 flt_deltaFluxPlus[int_kStencil]=flt_fluxPlus[int_kStencil+1]-flt_fluxPlus[int_kStencil];
 flt_deltaFluxMinus[int_kStencil]=flt_fluxMinus[int_kStencil]-flt_fluxMinus[int_kStencil+1];
                    }


#if SPACEORDER==9
                    //    Aufruf der Phi-Funktion fuer positive Eigenwerte
                    flt_Phi_FunctionPlus=Phi_Function_W9(
                            pnt_config,
                            flt_fluxPlus,
                            flt_deltaFluxPlus);
                    //    Aufruf der Phi-Funktion fuer negative Eigenwerte
                    flt_Phi_FunctionMinus=Phi_Function_W9(
                            pnt_config,
                            flt_fluxMinus,
                            flt_deltaFluxMinus);
#endif
#if SPACEORDER==5
					//    Aufruf der Phi-Funktion fuer positive Eigenwerte
					flt_Phi_FunctionPlus=Phi_Function_W5(
							pnt_config,
							flt_fluxPlus,
							flt_deltaFluxPlus);
					//    Aufruf der Phi-Funktion fuer negative Eigenwerte
					flt_Phi_FunctionMinus=Phi_Function_W5(
							pnt_config,
							flt_fluxMinus,
							flt_deltaFluxMinus);
#endif
 flt_Phi_FunctionSum[int_c]=(-flt_Phi_FunctionPlus+flt_Phi_FunctionMinus);

                }
                for (int_c=0; int_c < pnt_config->int_conservationEquations; int_c++)
                {
                    flt_Phi_FunctionSumRightEigen[int_c]=
 flt_Phi_FunctionSum[0]*flt_rightEigenvector[int_c][0]+
 flt_Phi_FunctionSum[1]*flt_rightEigenvector[int_c][1]+
 flt_Phi_FunctionSum[2]*flt_rightEigenvector[int_c][2]+
 flt_Phi_FunctionSum[3]*flt_rightEigenvector[int_c][3]+
 flt_Phi_FunctionSum[4]*flt_rightEigenvector[int_c][4];
                }


#if SPACEORDER==9
                //    Berechnung des Zentralterms über die Theta-Funktion
                Theta_Function_W9(
                    pnt_Flux,
                    flt_Theta_Function,
                    ijk,
                    ijkPlus1,
                    ijkPlus2,
                    ijkPlus3,
                    ijkPlus4,
                    ijkMinus1,
                    ijkMinus2,
                    ijkMinus3);
#endif
#if SPACEORDER==5
				//    Berechnung des Zentralterms über die Theta-Funktion
				Theta_Function_W5(
					pnt_Flux,
					flt_Theta_Function,
					ijk,
					ijkPlus1,
					ijkPlus2,
					ijkMinus1);
#endif

                //    Berechnung des resultierenden Fluss an der Stelle j+1/2 (entspricht j)
                //    Theta=Zentralterm
                //    Phi=WENO-Terme
 pnt_Flux_PlusHalf->Mass[ijk]=flt_Phi_FunctionSumRightEigen[0]+flt_Theta_Function[0];
 pnt_Flux_PlusHalf->xiMomentum[ijk]=flt_Phi_FunctionSumRightEigen[1]+flt_Theta_Function[1];
 pnt_Flux_PlusHalf->etaMomentum[ijk]=flt_Phi_FunctionSumRightEigen[2]+flt_Theta_Function[2];
 pnt_Flux_PlusHalf->zetaMomentum[ijk]=flt_Phi_FunctionSumRightEigen[3]+flt_Theta_Function[3];
 pnt_Flux_PlusHalf->Energy[ijk]=flt_Phi_FunctionSumRightEigen[4]+flt_Theta_Function[4];
            }
        }
    }

    for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
    {
        for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
        {
            for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
            {
 ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
 ijkMinus1=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-1;

 pnt_Q->Mass[ijk]=pnt_Q->Mass[ijk]-(pnt_Flux_PlusHalf->Mass[ijk]-pnt_Flux_PlusHalf->Mass[ijkMinus1]);
 pnt_Q->xiMomentum[ijk]=pnt_Q->xiMomentum[ijk]-(pnt_Flux_PlusHalf->xiMomentum[ijk]-pnt_Flux_PlusHalf->xiMomentum[ijkMinus1]);
 pnt_Q->etaMomentum[ijk]=pnt_Q->etaMomentum[ijk]-(pnt_Flux_PlusHalf->etaMomentum[ijk]-pnt_Flux_PlusHalf->etaMomentum[ijkMinus1]);
 pnt_Q->zetaMomentum[ijk]=pnt_Q->zetaMomentum[ijk]-(pnt_Flux_PlusHalf->zetaMomentum[ijk]-pnt_Flux_PlusHalf->zetaMomentum[ijkMinus1]);
 pnt_Q->Energy[ijk]=pnt_Q->Energy[ijk]-(pnt_Flux_PlusHalf->Energy[ijk]-pnt_Flux_PlusHalf->Energy[ijkMinus1]);
            }
        }
    }

    free(flt_fluxPlus);
    free(flt_fluxMinus);
    free(flt_deltaFluxPlus);
    free(flt_deltaFluxMinus);
    free(flt_Phi_FunctionSum);
    free(flt_Phi_FunctionSumRightEigen);
    free(flt_Theta_Function);
}

void CalcEigenVectorsInXiDirection(
        struct strct_configuration * pnt_config,
        struct strct_mesh * pnt_mesh,
        struct strct_U * pnt_U_RK,
        int ijk,
        int iPlus1jk,
        int iMinus1jk)
{
    float u,v,w,rho;
    float xi_x_hat,xi_y_hat,xi_z_hat;
//    float eta_x_hat,eta_y_hat,eta_z_hat;
//    float zeta_x_hat,zeta_y_hat,zeta_z_hat;
    float gamma,c,M;
    float sqrt_2;

    sqrt_2=sqrt(2.);

    gamma=pnt_config->flt_gammaNumber;

    u=0.5*(pnt_U_RK->u[ijk]+pnt_U_RK->u[iPlus1jk]);
    v=0.5*(pnt_U_RK->v[ijk]+pnt_U_RK->v[iPlus1jk]);
    w=0.5*(pnt_U_RK->w[ijk]+pnt_U_RK->w[iPlus1jk]);
    rho=0.5*(pnt_U_RK->rho[ijk]+pnt_U_RK->rho[iPlus1jk]);
    c=0.5*(pnt_U_RK->c[ijk]+pnt_U_RK->c[iPlus1jk]);
    M=sqrt(u*u+v*v+w*w)/c;

//    float sqrt_rho=sqrt(pnt_U_RK->rho[ijk]);
//    float sqrt_rhoP1=sqrt(pnt_U_RK->rho[iPlus1jk]);
//
//    u=(pnt_U_RK->u[ijk]*sqrt_rhoP1+pnt_U_RK->u[iPlus1jk]*sqrt_rho)
//    		/(sqrt_rhoP1+sqrt_rho);
//    v=(pnt_U_RK->v[ijk]*sqrt_rhoP1+pnt_U_RK->v[iPlus1jk]*sqrt_rho)
//    		/(sqrt_rhoP1+sqrt_rho);
//    w=(pnt_U_RK->w[ijk]*sqrt_rhoP1+pnt_U_RK->w[iPlus1jk]*sqrt_rho)
//    		/(sqrt_rhoP1+sqrt_rho);
//    c=(pnt_U_RK->c[ijk]*sqrt_rhoP1+pnt_U_RK->c[iPlus1jk]*sqrt_rho)
//        		/(sqrt_rhoP1+sqrt_rho);
//    rho=sqrt_rho*sqrt_rhoP1;
//    M=sqrt(u*u+v*v+w*w)/c;

    xi_x_hat=0.5*
            (
 pnt_mesh->xi_x[ijk]/sqrt(pow(pnt_mesh->xi_x[ijk],2)+pow(pnt_mesh->xi_y[ijk],2)+pow(pnt_mesh->xi_z[ijk],2.))+
 pnt_mesh->xi_x[iPlus1jk]/sqrt(pow(pnt_mesh->xi_x[iPlus1jk],2)+pow(pnt_mesh->xi_y[iPlus1jk],2)+pow(pnt_mesh->xi_z[iPlus1jk],2.))
            );
    xi_y_hat=0.5*
            (
 pnt_mesh->xi_y[ijk]/sqrt(pow(pnt_mesh->xi_x[ijk],2)+pow(pnt_mesh->xi_y[ijk],2)+pow(pnt_mesh->xi_z[ijk],2.))+
 pnt_mesh->xi_y[iPlus1jk]/sqrt(pow(pnt_mesh->xi_x[iPlus1jk],2)+pow(pnt_mesh->xi_y[iPlus1jk],2)+pow(pnt_mesh->xi_z[iPlus1jk],2.))
            );
    xi_z_hat=0.5*
            (
 pnt_mesh->xi_z[ijk]/sqrt(pow(pnt_mesh->xi_x[ijk],2)+pow(pnt_mesh->xi_y[ijk],2)+pow(pnt_mesh->xi_z[ijk],2.))+
 pnt_mesh->xi_z[iPlus1jk]/sqrt(pow(pnt_mesh->xi_x[iPlus1jk],2)+pow(pnt_mesh->xi_y[iPlus1jk],2)+pow(pnt_mesh->xi_z[iPlus1jk],2.))
            );

//    eta_x_hat=0.5*
//            (
//  pnt_mesh->eta_x[ijk]/sqrt(pow(pnt_mesh->eta_x[ijk],2)+pow(pnt_mesh->eta_y[ijk],2)+pow(pnt_mesh->eta_z[ijk],2.))+
//  pnt_mesh->eta_x[iPlus1jk]/sqrt(pow(pnt_mesh->eta_x[iPlus1jk],2)+pow(pnt_mesh->eta_y[iPlus1jk],2)+pow(pnt_mesh->eta_z[iPlus1jk],2.))
//            );
//    eta_y_hat=0.5*
//            (
//  pnt_mesh->eta_y[ijk]/sqrt(pow(pnt_mesh->eta_x[ijk],2)+pow(pnt_mesh->eta_y[ijk],2)+pow(pnt_mesh->eta_z[ijk],2.))+
//  pnt_mesh->eta_y[iPlus1jk]/sqrt(pow(pnt_mesh->eta_x[iPlus1jk],2)+pow(pnt_mesh->eta_y[iPlus1jk],2)+pow(pnt_mesh->eta_z[iPlus1jk],2.))
//            );
//    eta_z_hat=0.5*
//            (
//  pnt_mesh->eta_z[ijk]/sqrt(pow(pnt_mesh->eta_x[ijk],2)+pow(pnt_mesh->eta_y[ijk],2)+pow(pnt_mesh->eta_z[ijk],2.))+
//  pnt_mesh->eta_z[iPlus1jk]/sqrt(pow(pnt_mesh->eta_x[iPlus1jk],2)+pow(pnt_mesh->eta_y[iPlus1jk],2)+pow(pnt_mesh->eta_z[iPlus1jk],2.))
//            );
//
//    zeta_x_hat=0.5*
//            (
//  pnt_mesh->zeta_x[ijk]/sqrt(pow(pnt_mesh->zeta_x[ijk],2)+pow(pnt_mesh->zeta_y[ijk],2)+pow(pnt_mesh->zeta_z[ijk],2.))+
//  pnt_mesh->zeta_x[iPlus1jk]/sqrt(pow(pnt_mesh->zeta_x[iPlus1jk],2)+pow(pnt_mesh->zeta_y[iPlus1jk],2)+pow(pnt_mesh->zeta_z[iPlus1jk],2.))
//            );
//    zeta_y_hat=0.5*
//            (
//  pnt_mesh->zeta_y[ijk]/sqrt(pow(pnt_mesh->zeta_x[ijk],2)+pow(pnt_mesh->zeta_y[ijk],2)+pow(pnt_mesh->zeta_z[ijk],2.))+
//  pnt_mesh->zeta_y[iPlus1jk]/sqrt(pow(pnt_mesh->zeta_x[iPlus1jk],2)+pow(pnt_mesh->zeta_y[iPlus1jk],2)+pow(pnt_mesh->zeta_z[iPlus1jk],2.))
//            );
//    zeta_z_hat=0.5*
//            (
//  pnt_mesh->zeta_z[ijk]/sqrt(pow(pnt_mesh->zeta_x[ijk],2)+pow(pnt_mesh->zeta_y[ijk],2)+pow(pnt_mesh->zeta_z[ijk],2.))+
//  pnt_mesh->zeta_z[iPlus1jk]/sqrt(pow(pnt_mesh->zeta_x[iPlus1jk],2)+pow(pnt_mesh->zeta_y[iPlus1jk],2)+pow(pnt_mesh->zeta_z[iPlus1jk],2.))
//            );



    //Rechter Eigenvektor
    flt_rightEigenvector[0][0]=xi_x_hat;
    flt_rightEigenvector[0][1]=xi_y_hat;
    flt_rightEigenvector[0][2]=xi_z_hat;
    flt_rightEigenvector[0][3]=(1./(c*sqrt_2));
    flt_rightEigenvector[0][4]=(1./(c*sqrt_2));

    flt_rightEigenvector[1][0]=u*xi_x_hat;
    flt_rightEigenvector[1][1]=u*xi_y_hat-rho*xi_z_hat;
    flt_rightEigenvector[1][2]=u*xi_z_hat+rho*xi_y_hat;
    flt_rightEigenvector[1][3]=(1./(c*sqrt_2))*(u+xi_x_hat*c);
    flt_rightEigenvector[1][4]=(1./(c*sqrt_2))*(u-xi_x_hat*c);

    flt_rightEigenvector[2][0]=v*xi_x_hat+rho*xi_z_hat;
    flt_rightEigenvector[2][1]=v*xi_y_hat;
    flt_rightEigenvector[2][2]=v*xi_z_hat-rho*xi_x_hat;
    flt_rightEigenvector[2][3]=(1./(c*sqrt_2))*(v+xi_y_hat*c);
    flt_rightEigenvector[2][4]=(1./(c*sqrt_2))*(v-xi_y_hat*c);

    flt_rightEigenvector[3][0]=w*xi_x_hat-rho*xi_y_hat;
    flt_rightEigenvector[3][1]=w*xi_y_hat+rho*xi_x_hat;
    flt_rightEigenvector[3][2]=w*xi_z_hat;
    flt_rightEigenvector[3][3]=(1./(c*sqrt_2))*(w+xi_z_hat*c);
    flt_rightEigenvector[3][4]=(1./(c*sqrt_2))*(w-xi_z_hat*c);

 flt_rightEigenvector[4][0]=(u*u+v*v+w*w)/2.*xi_x_hat+rho*(v*xi_z_hat-w*xi_y_hat);
 flt_rightEigenvector[4][1]=(u*u+v*v+w*w)/2.*xi_y_hat+rho*(w*xi_x_hat-u*xi_z_hat);
 flt_rightEigenvector[4][2]=(u*u+v*v+w*w)/2.*xi_z_hat+rho*(u*xi_y_hat-v*xi_x_hat);
 flt_rightEigenvector[4][3]=(1./(c*sqrt_2))*((u*u+v*v+w*w)/2.+c*c/(gamma-1.)+c*(u*xi_x_hat+v*xi_y_hat+w*xi_z_hat));
 flt_rightEigenvector[4][4]=(1./(c*sqrt_2))*((u*u+v*v+w*w)/2.+c*c/(gamma-1.)-c*(u*xi_x_hat+v*xi_y_hat+w*xi_z_hat));


    //Linker Eigenvektor
 flt_leftEigenvector[0][0]=xi_x_hat*(1.0-0.5*(gamma-1.0)*(M*M))-1./rho*(v*xi_z_hat-w*xi_y_hat);
    flt_leftEigenvector[0][1]=((gamma-1.0)/(c*c))*u*xi_x_hat;
 flt_leftEigenvector[0][2]=((gamma-1.0)/(c*c))*v*xi_x_hat+(xi_z_hat)/rho;
 flt_leftEigenvector[0][3]=((gamma-1.0)/(c*c))*w*xi_x_hat-(xi_y_hat)/rho;
    flt_leftEigenvector[0][4]=-((gamma-1.0)/(c*c))*xi_x_hat;

 flt_leftEigenvector[1][0]=xi_y_hat*(1.0-0.5*(gamma-1.0)*(M*M))-1./rho*(w*xi_x_hat-u*xi_z_hat);
 flt_leftEigenvector[1][1]=((gamma-1.0)/(c*c))*u*xi_y_hat-(xi_z_hat)/rho;
    flt_leftEigenvector[1][2]=((gamma-1.0)/(c*c))*v*xi_y_hat;
 flt_leftEigenvector[1][3]=((gamma-1.0)/(c*c))*w*xi_y_hat+(xi_x_hat)/rho;
    flt_leftEigenvector[1][4]=-((gamma-1.0)/(c*c))*xi_y_hat;

 flt_leftEigenvector[2][0]=xi_z_hat*(1.0-0.5*(gamma-1.0)*(M*M))-1./rho*(u*xi_y_hat-v*xi_x_hat);
 flt_leftEigenvector[2][1]=((gamma-1.0)/(c*c))*u*xi_z_hat+(xi_y_hat)/rho;
 flt_leftEigenvector[2][2]=((gamma-1.0)/(c*c))*v*xi_z_hat-(xi_x_hat)/rho;
    flt_leftEigenvector[2][3]=((gamma-1.0)/(c*c))*w*xi_z_hat;
    flt_leftEigenvector[2][4]=-((gamma-1.0)/(c*c))*xi_z_hat;

 flt_leftEigenvector[3][0]=1./(sqrt_2*c)*((gamma-1.0)*(u*u+v*v+w*w)/2.-c*(u*xi_x_hat+v*xi_y_hat+w*xi_z_hat));
    flt_leftEigenvector[3][1]=(-(gamma-1.0)*u/c+xi_x_hat)/sqrt_2;
    flt_leftEigenvector[3][2]=(-(gamma-1.0)*v/c+xi_y_hat)/sqrt_2;
    flt_leftEigenvector[3][3]=(-(gamma-1.0)*w/c+xi_z_hat)/sqrt_2;
    flt_leftEigenvector[3][4]=(gamma-1.)/(sqrt_2*c);

 flt_leftEigenvector[4][0]=1./(sqrt_2*c)*((gamma-1.0)*(u*u+v*v+w*w)/2.+c*(u*xi_x_hat+v*xi_y_hat+w*xi_z_hat));
    flt_leftEigenvector[4][1]=(-(gamma-1.)*u/c-xi_x_hat)/sqrt_2;
    flt_leftEigenvector[4][2]=(-(gamma-1.)*v/c-xi_y_hat)/sqrt_2;
    flt_leftEigenvector[4][3]=(-(gamma-1.)*w/c-xi_z_hat)/sqrt_2;
    flt_leftEigenvector[4][4]=(gamma-1.)/(sqrt_2*c);

}

void CalcEigenVectorsInEtaDirection(
        struct strct_configuration * pnt_config,
        struct strct_mesh * pnt_mesh,
        struct strct_U * pnt_U_RK,
        int ijk,
        int ijPlus1k,
        int ijMinus1k)
{
    float u,v,w,rho;
//    float xi_x_hat,xi_y_hat,xi_z_hat;
    float eta_x_hat,eta_y_hat,eta_z_hat;
//    float zeta_x_hat,zeta_y_hat,zeta_z_hat;
    float gamma,c,M;
    float sqrt_2;

    sqrt_2=sqrt(2.);

    gamma=pnt_config->flt_gammaNumber;

    u=0.5*(pnt_U_RK->u[ijk]+pnt_U_RK->u[ijPlus1k]);
    v=0.5*(pnt_U_RK->v[ijk]+pnt_U_RK->v[ijPlus1k]);
    w=0.5*(pnt_U_RK->w[ijk]+pnt_U_RK->w[ijPlus1k]);
    rho=0.5*(pnt_U_RK->rho[ijk]+pnt_U_RK->rho[ijPlus1k]);
    c=0.5*(pnt_U_RK->c[ijk]+pnt_U_RK->c[ijPlus1k]);
    M=sqrt(u*u+v*v+w*w)/c;

//    float sqrt_rho=sqrt(pnt_U_RK->rho[ijk]);
//    float sqrt_rhoP1=sqrt(pnt_U_RK->rho[ijPlus1k]);
//
//    u=(pnt_U_RK->u[ijk]*sqrt_rhoP1+pnt_U_RK->u[ijPlus1k]*sqrt_rho)
//    		/(sqrt_rhoP1+sqrt_rho);
//    v=(pnt_U_RK->v[ijk]*sqrt_rhoP1+pnt_U_RK->v[ijPlus1k]*sqrt_rho)
//    		/(sqrt_rhoP1+sqrt_rho);
//    w=(pnt_U_RK->w[ijk]*sqrt_rhoP1+pnt_U_RK->w[ijPlus1k]*sqrt_rho)
//    		/(sqrt_rhoP1+sqrt_rho);
//    c=(pnt_U_RK->c[ijk]*sqrt_rhoP1+pnt_U_RK->c[ijPlus1k]*sqrt_rho)
//        		/(sqrt_rhoP1+sqrt_rho);
//    rho=sqrt_rho*sqrt_rhoP1;
//    M=sqrt(u*u+v*v+w*w)/c;

//    xi_x_hat=0.5*
//            (
//  pnt_mesh->xi_x[ijk]/sqrt(pow(pnt_mesh->xi_x[ijk],2)+pow(pnt_mesh->xi_y[ijk],2)+pow(pnt_mesh->xi_z[ijk],2.))+
//  pnt_mesh->xi_x[ijPlus1k]/sqrt(pow(pnt_mesh->xi_x[ijPlus1k],2)+pow(pnt_mesh->xi_y[ijPlus1k],2)+pow(pnt_mesh->xi_z[ijPlus1k],2.))
//            );
//    xi_y_hat=0.5*
//            (
//  pnt_mesh->xi_y[ijk]/sqrt(pow(pnt_mesh->xi_x[ijk],2)+pow(pnt_mesh->xi_y[ijk],2)+pow(pnt_mesh->xi_z[ijk],2.))+
//  pnt_mesh->xi_y[ijPlus1k]/sqrt(pow(pnt_mesh->xi_x[ijPlus1k],2)+pow(pnt_mesh->xi_y[ijPlus1k],2)+pow(pnt_mesh->xi_z[ijPlus1k],2.))
//            );
//    xi_z_hat=0.5*
//            (
//  pnt_mesh->xi_z[ijk]/sqrt(pow(pnt_mesh->xi_x[ijk],2)+pow(pnt_mesh->xi_y[ijk],2)+pow(pnt_mesh->xi_z[ijk],2.))+
//  pnt_mesh->xi_z[ijPlus1k]/sqrt(pow(pnt_mesh->xi_x[ijPlus1k],2)+pow(pnt_mesh->xi_y[ijPlus1k],2)+pow(pnt_mesh->xi_z[ijPlus1k],2.))
//            );

    eta_x_hat=0.5*
            (
 pnt_mesh->eta_x[ijk]/sqrt(pow(pnt_mesh->eta_x[ijk],2)+pow(pnt_mesh->eta_y[ijk],2)+pow(pnt_mesh->eta_z[ijk],2.))+
 pnt_mesh->eta_x[ijPlus1k]/sqrt(pow(pnt_mesh->eta_x[ijPlus1k],2)+pow(pnt_mesh->eta_y[ijPlus1k],2)+pow(pnt_mesh->eta_z[ijPlus1k],2.))
            );
    eta_y_hat=0.5*
            (
 pnt_mesh->eta_y[ijk]/sqrt(pow(pnt_mesh->eta_x[ijk],2)+pow(pnt_mesh->eta_y[ijk],2)+pow(pnt_mesh->eta_z[ijk],2.))+
 pnt_mesh->eta_y[ijPlus1k]/sqrt(pow(pnt_mesh->eta_x[ijPlus1k],2)+pow(pnt_mesh->eta_y[ijPlus1k],2)+pow(pnt_mesh->eta_z[ijPlus1k],2.))
            );
    eta_z_hat=0.5*
            (
 pnt_mesh->eta_z[ijk]/sqrt(pow(pnt_mesh->eta_x[ijk],2)+pow(pnt_mesh->eta_y[ijk],2)+pow(pnt_mesh->eta_z[ijk],2.))+
 pnt_mesh->eta_z[ijPlus1k]/sqrt(pow(pnt_mesh->eta_x[ijPlus1k],2)+pow(pnt_mesh->eta_y[ijPlus1k],2)+pow(pnt_mesh->eta_z[ijPlus1k],2.))
            );

//    zeta_x_hat=0.5*
//            (
//  pnt_mesh->zeta_x[ijk]/sqrt(pow(pnt_mesh->zeta_x[ijk],2)+pow(pnt_mesh->zeta_y[ijk],2)+pow(pnt_mesh->zeta_z[ijk],2.))+
//  pnt_mesh->zeta_x[ijPlus1k]/sqrt(pow(pnt_mesh->zeta_x[ijPlus1k],2)+pow(pnt_mesh->zeta_y[ijPlus1k],2)+pow(pnt_mesh->zeta_z[ijPlus1k],2.))
//            );
//    zeta_y_hat=0.5*
//            (
//  pnt_mesh->zeta_y[ijk]/sqrt(pow(pnt_mesh->zeta_x[ijk],2)+pow(pnt_mesh->zeta_y[ijk],2)+pow(pnt_mesh->zeta_z[ijk],2.))+
//  pnt_mesh->zeta_y[ijPlus1k]/sqrt(pow(pnt_mesh->zeta_x[ijPlus1k],2)+pow(pnt_mesh->zeta_y[ijPlus1k],2)+pow(pnt_mesh->zeta_z[ijPlus1k],2.))
//            );
//    zeta_z_hat=0.5*
//            (
//  pnt_mesh->zeta_z[ijk]/sqrt(pow(pnt_mesh->zeta_x[ijk],2)+pow(pnt_mesh->zeta_y[ijk],2)+pow(pnt_mesh->zeta_z[ijk],2.))+
//  pnt_mesh->zeta_z[ijPlus1k]/sqrt(pow(pnt_mesh->zeta_x[ijPlus1k],2)+pow(pnt_mesh->zeta_y[ijPlus1k],2)+pow(pnt_mesh->zeta_z[ijPlus1k],2.))
//            );



    //Rechter Eigenvektor
    flt_rightEigenvector[0][0]=eta_x_hat;
    flt_rightEigenvector[0][1]=eta_y_hat;
    flt_rightEigenvector[0][2]=eta_z_hat;
    flt_rightEigenvector[0][3]=(1./(c*sqrt_2));
    flt_rightEigenvector[0][4]=(1./(c*sqrt_2));

    flt_rightEigenvector[1][0]=u*eta_x_hat;
    flt_rightEigenvector[1][1]=u*eta_y_hat-rho*eta_z_hat;
    flt_rightEigenvector[1][2]=u*eta_z_hat+rho*eta_y_hat;
    flt_rightEigenvector[1][3]=(1./(c*sqrt_2))*(u+eta_x_hat*c);
    flt_rightEigenvector[1][4]=(1./(c*sqrt_2))*(u-eta_x_hat*c);

    flt_rightEigenvector[2][0]=v*eta_x_hat+rho*eta_z_hat;
    flt_rightEigenvector[2][1]=v*eta_y_hat;
    flt_rightEigenvector[2][2]=v*eta_z_hat-rho*eta_x_hat;
    flt_rightEigenvector[2][3]=(1./(c*sqrt_2))*(v+eta_y_hat*c);
    flt_rightEigenvector[2][4]=(1./(c*sqrt_2))*(v-eta_y_hat*c);

    flt_rightEigenvector[3][0]=w*eta_x_hat-rho*eta_y_hat;
    flt_rightEigenvector[3][1]=w*eta_y_hat+rho*eta_x_hat;
    flt_rightEigenvector[3][2]=w*eta_z_hat;
    flt_rightEigenvector[3][3]=(1./(c*sqrt_2))*(w+eta_z_hat*c);
    flt_rightEigenvector[3][4]=(1./(c*sqrt_2))*(w-eta_z_hat*c);

 flt_rightEigenvector[4][0]=(u*u+v*v+w*w)/2.*eta_x_hat+rho*(v*eta_z_hat-w*eta_y_hat);
 flt_rightEigenvector[4][1]=(u*u+v*v+w*w)/2.*eta_y_hat+rho*(w*eta_x_hat-u*eta_z_hat);
 flt_rightEigenvector[4][2]=(u*u+v*v+w*w)/2.*eta_z_hat+rho*(u*eta_y_hat-v*eta_x_hat);
 flt_rightEigenvector[4][3]=(1./(c*sqrt_2))*((u*u+v*v+w*w)/2.+c*c/(gamma-1.)+c*(u*eta_x_hat+v*eta_y_hat+w*eta_z_hat));
 flt_rightEigenvector[4][4]=(1./(c*sqrt_2))*((u*u+v*v+w*w)/2.+c*c/(gamma-1.)-c*(u*eta_x_hat+v*eta_y_hat+w*eta_z_hat));


    //Linker Eigenvektor
 flt_leftEigenvector[0][0]=eta_x_hat*(1.0-0.5*(gamma-1.0)*(M*M))-1./rho*(v*eta_z_hat-w*eta_y_hat);
    flt_leftEigenvector[0][1]=((gamma-1.0)/(c*c))*u*eta_x_hat;
 flt_leftEigenvector[0][2]=((gamma-1.0)/(c*c))*v*eta_x_hat+(eta_z_hat)/rho;
 flt_leftEigenvector[0][3]=((gamma-1.0)/(c*c))*w*eta_x_hat-(eta_y_hat)/rho;
    flt_leftEigenvector[0][4]=-((gamma-1.0)/(c*c))*eta_x_hat;

 flt_leftEigenvector[1][0]=eta_y_hat*(1.0-0.5*(gamma-1.0)*(M*M))-1./rho*(w*eta_x_hat-u*eta_z_hat);
 flt_leftEigenvector[1][1]=((gamma-1.0)/(c*c))*u*eta_y_hat-(eta_z_hat)/rho;
    flt_leftEigenvector[1][2]=((gamma-1.0)/(c*c))*v*eta_y_hat;
 flt_leftEigenvector[1][3]=((gamma-1.0)/(c*c))*w*eta_y_hat+(eta_x_hat)/rho;
    flt_leftEigenvector[1][4]=-((gamma-1.0)/(c*c))*eta_y_hat;

 flt_leftEigenvector[2][0]=eta_z_hat*(1.0-0.5*(gamma-1.0)*(M*M))-1./rho*(u*eta_y_hat-v*eta_x_hat);
 flt_leftEigenvector[2][1]=((gamma-1.0)/(c*c))*u*eta_z_hat+(eta_y_hat)/rho;
 flt_leftEigenvector[2][2]=((gamma-1.0)/(c*c))*v*eta_z_hat-(eta_x_hat)/rho;
    flt_leftEigenvector[2][3]=((gamma-1.0)/(c*c))*w*eta_z_hat;
    flt_leftEigenvector[2][4]=-((gamma-1.0)/(c*c))*eta_z_hat;

 flt_leftEigenvector[3][0]=1./(sqrt_2*c)*((gamma-1.0)*(u*u+v*v+w*w)/2.-c*(u*eta_x_hat+v*eta_y_hat+w*eta_z_hat));
    flt_leftEigenvector[3][1]=(-(gamma-1.0)*u/c+eta_x_hat)/sqrt_2;
    flt_leftEigenvector[3][2]=(-(gamma-1.0)*v/c+eta_y_hat)/sqrt_2;
    flt_leftEigenvector[3][3]=(-(gamma-1.0)*w/c+eta_z_hat)/sqrt_2;
    flt_leftEigenvector[3][4]=(gamma-1.)/(sqrt_2*c);

 flt_leftEigenvector[4][0]=1./(sqrt_2*c)*((gamma-1.0)*(u*u+v*v+w*w)/2.+c*(u*eta_x_hat+v*eta_y_hat+w*eta_z_hat));
    flt_leftEigenvector[4][1]=(-(gamma-1.)*u/c-eta_x_hat)/sqrt_2;
    flt_leftEigenvector[4][2]=(-(gamma-1.)*v/c-eta_y_hat)/sqrt_2;
    flt_leftEigenvector[4][3]=(-(gamma-1.)*w/c-eta_z_hat)/sqrt_2;
    flt_leftEigenvector[4][4]=(gamma-1.)/(sqrt_2*c);

}

void CalcEigenVectorsInZetaDirection(
        struct strct_configuration * pnt_config,
        struct strct_mesh * pnt_mesh,
        struct strct_U * pnt_U_RK,
        int ijk,
        int ijkPlus1,
        int ijkMinus1)
{
    float u,v,w,rho;
//    float xi_x_hat,xi_y_hat,xi_z_hat;
//    float eta_x_hat,eta_y_hat,eta_z_hat;
    float zeta_x_hat,zeta_y_hat,zeta_z_hat;
    float gamma,c,M;
    float sqrt_2;

    sqrt_2=sqrt(2.);

    gamma=pnt_config->flt_gammaNumber;

    u=0.5*(pnt_U_RK->u[ijk]+pnt_U_RK->u[ijkPlus1]);
    v=0.5*(pnt_U_RK->v[ijk]+pnt_U_RK->v[ijkPlus1]);
    w=0.5*(pnt_U_RK->w[ijk]+pnt_U_RK->w[ijkPlus1]);
    rho=0.5*(pnt_U_RK->rho[ijk]+pnt_U_RK->rho[ijkPlus1]);
    c=0.5*(pnt_U_RK->c[ijk]+pnt_U_RK->c[ijkPlus1]);
    M=sqrt(u*u+v*v+w*w)/c;

//    float sqrt_rho=sqrt(pnt_U_RK->rho[ijk]);
//    float sqrt_rhoP1=sqrt(pnt_U_RK->rho[ijkPlus1]);
//
//    u=(pnt_U_RK->u[ijk]*sqrt_rhoP1+pnt_U_RK->u[ijkPlus1]*sqrt_rho)
//    		/(sqrt_rhoP1+sqrt_rho);
//    v=(pnt_U_RK->v[ijk]*sqrt_rhoP1+pnt_U_RK->v[ijkPlus1]*sqrt_rho)
//    		/(sqrt_rhoP1+sqrt_rho);
//    w=(pnt_U_RK->w[ijk]*sqrt_rhoP1+pnt_U_RK->w[ijkPlus1]*sqrt_rho)
//    		/(sqrt_rhoP1+sqrt_rho);
//    c=(pnt_U_RK->c[ijk]*sqrt_rhoP1+pnt_U_RK->c[ijkPlus1]*sqrt_rho)
//        		/(sqrt_rhoP1+sqrt_rho);
//    rho=sqrt_rho*sqrt_rhoP1;
//    M=sqrt(u*u+v*v+w*w)/c;

//    xi_x_hat=0.5*
//            (
//  pnt_mesh->xi_x[ijk]/sqrt(pow(pnt_mesh->xi_x[ijk],2)+pow(pnt_mesh->xi_y[ijk],2)+pow(pnt_mesh->xi_z[ijk],2.))+
//  pnt_mesh->xi_x[ijkPlus1]/sqrt(pow(pnt_mesh->xi_x[ijkPlus1],2)+pow(pnt_mesh->xi_y[ijkPlus1],2)+pow(pnt_mesh->xi_z[ijkPlus1],2.))
//            );
//    xi_y_hat=0.5*
//            (
//  pnt_mesh->xi_y[ijk]/sqrt(pow(pnt_mesh->xi_x[ijk],2)+pow(pnt_mesh->xi_y[ijk],2)+pow(pnt_mesh->xi_z[ijk],2.))+
//  pnt_mesh->xi_y[ijkPlus1]/sqrt(pow(pnt_mesh->xi_x[ijkPlus1],2)+pow(pnt_mesh->xi_y[ijkPlus1],2)+pow(pnt_mesh->xi_z[ijkPlus1],2.))
//            );
//    xi_z_hat=0.5*
//            (
//  pnt_mesh->xi_z[ijk]/sqrt(pow(pnt_mesh->xi_x[ijk],2)+pow(pnt_mesh->xi_y[ijk],2)+pow(pnt_mesh->xi_z[ijk],2.))+
//  pnt_mesh->xi_z[ijkPlus1]/sqrt(pow(pnt_mesh->xi_x[ijkPlus1],2)+pow(pnt_mesh->xi_y[ijkPlus1],2)+pow(pnt_mesh->xi_z[ijkPlus1],2.))
//            );
//
//    eta_x_hat=0.5*
//            (
//  pnt_mesh->eta_x[ijk]/sqrt(pow(pnt_mesh->eta_x[ijk],2)+pow(pnt_mesh->eta_y[ijk],2)+pow(pnt_mesh->eta_z[ijk],2.))+
//  pnt_mesh->eta_x[ijkPlus1]/sqrt(pow(pnt_mesh->eta_x[ijkPlus1],2)+pow(pnt_mesh->eta_y[ijkPlus1],2)+pow(pnt_mesh->eta_z[ijkPlus1],2.))
//            );
//    eta_y_hat=0.5*
//            (
//  pnt_mesh->eta_y[ijk]/sqrt(pow(pnt_mesh->eta_x[ijk],2)+pow(pnt_mesh->eta_y[ijk],2)+pow(pnt_mesh->eta_z[ijk],2.))+
//  pnt_mesh->eta_y[ijkPlus1]/sqrt(pow(pnt_mesh->eta_x[ijkPlus1],2)+pow(pnt_mesh->eta_y[ijkPlus1],2)+pow(pnt_mesh->eta_z[ijkPlus1],2.))
//            );
//    eta_z_hat=0.5*
//            (
//  pnt_mesh->eta_z[ijk]/sqrt(pow(pnt_mesh->eta_x[ijk],2)+pow(pnt_mesh->eta_y[ijk],2)+pow(pnt_mesh->eta_z[ijk],2.))+
//  pnt_mesh->eta_z[ijkPlus1]/sqrt(pow(pnt_mesh->eta_x[ijkPlus1],2)+pow(pnt_mesh->eta_y[ijkPlus1],2)+pow(pnt_mesh->eta_z[ijkPlus1],2.))
//            );

    zeta_x_hat=0.5*
            (
 pnt_mesh->zeta_x[ijk]/sqrt(pow(pnt_mesh->zeta_x[ijk],2)+pow(pnt_mesh->zeta_y[ijk],2)+pow(pnt_mesh->zeta_z[ijk],2.))+
 pnt_mesh->zeta_x[ijkPlus1]/sqrt(pow(pnt_mesh->zeta_x[ijkPlus1],2)+pow(pnt_mesh->zeta_y[ijkPlus1],2)+pow(pnt_mesh->zeta_z[ijkPlus1],2.))
            );
    zeta_y_hat=0.5*
            (
 pnt_mesh->zeta_y[ijk]/sqrt(pow(pnt_mesh->zeta_x[ijk],2)+pow(pnt_mesh->zeta_y[ijk],2)+pow(pnt_mesh->zeta_z[ijk],2.))+
 pnt_mesh->zeta_y[ijkPlus1]/sqrt(pow(pnt_mesh->zeta_x[ijkPlus1],2)+pow(pnt_mesh->zeta_y[ijkPlus1],2)+pow(pnt_mesh->zeta_z[ijkPlus1],2.))
            );
    zeta_z_hat=0.5*
            (
 pnt_mesh->zeta_z[ijk]/sqrt(pow(pnt_mesh->zeta_x[ijk],2)+pow(pnt_mesh->zeta_y[ijk],2)+pow(pnt_mesh->zeta_z[ijk],2.))+
 pnt_mesh->zeta_z[ijkPlus1]/sqrt(pow(pnt_mesh->zeta_x[ijkPlus1],2)+pow(pnt_mesh->zeta_y[ijkPlus1],2)+pow(pnt_mesh->zeta_z[ijkPlus1],2.))
            );



    //Rechter Eigenvektor
    flt_rightEigenvector[0][0]=zeta_x_hat;
    flt_rightEigenvector[0][1]=zeta_y_hat;
    flt_rightEigenvector[0][2]=zeta_z_hat;
    flt_rightEigenvector[0][3]=(1./(c*sqrt_2));
    flt_rightEigenvector[0][4]=(1./(c*sqrt_2));

    flt_rightEigenvector[1][0]=u*zeta_x_hat;
    flt_rightEigenvector[1][1]=u*zeta_y_hat-rho*zeta_z_hat;
    flt_rightEigenvector[1][2]=u*zeta_z_hat+rho*zeta_y_hat;
    flt_rightEigenvector[1][3]=(1./(c*sqrt_2))*(u+zeta_x_hat*c);
    flt_rightEigenvector[1][4]=(1./(c*sqrt_2))*(u-zeta_x_hat*c);

    flt_rightEigenvector[2][0]=v*zeta_x_hat+rho*zeta_z_hat;
    flt_rightEigenvector[2][1]=v*zeta_y_hat;
    flt_rightEigenvector[2][2]=v*zeta_z_hat-rho*zeta_x_hat;
    flt_rightEigenvector[2][3]=(1./(c*sqrt_2))*(v+zeta_y_hat*c);
    flt_rightEigenvector[2][4]=(1./(c*sqrt_2))*(v-zeta_y_hat*c);

    flt_rightEigenvector[3][0]=w*zeta_x_hat-rho*zeta_y_hat;
    flt_rightEigenvector[3][1]=w*zeta_y_hat+rho*zeta_x_hat;
    flt_rightEigenvector[3][2]=w*zeta_z_hat;
    flt_rightEigenvector[3][3]=(1./(c*sqrt_2))*(w+zeta_z_hat*c);
    flt_rightEigenvector[3][4]=(1./(c*sqrt_2))*(w-zeta_z_hat*c);

 flt_rightEigenvector[4][0]=(u*u+v*v+w*w)/2.*zeta_x_hat+rho*(v*zeta_z_hat-w*zeta_y_hat);
 flt_rightEigenvector[4][1]=(u*u+v*v+w*w)/2.*zeta_y_hat+rho*(w*zeta_x_hat-u*zeta_z_hat);
 flt_rightEigenvector[4][2]=(u*u+v*v+w*w)/2.*zeta_z_hat+rho*(u*zeta_y_hat-v*zeta_x_hat);
 flt_rightEigenvector[4][3]=(1./(c*sqrt_2))*((u*u+v*v+w*w)/2.+c*c/(gamma-1.)+c*(u*zeta_x_hat+v*zeta_y_hat+w*zeta_z_hat));
 flt_rightEigenvector[4][4]=(1./(c*sqrt_2))*((u*u+v*v+w*w)/2.+c*c/(gamma-1.)-c*(u*zeta_x_hat+v*zeta_y_hat+w*zeta_z_hat));


    //Linker Eigenvektor
 flt_leftEigenvector[0][0]=zeta_x_hat*(1.0-0.5*(gamma-1.0)*(M*M))-1./rho*(v*zeta_z_hat-w*zeta_y_hat);
    flt_leftEigenvector[0][1]=((gamma-1.0)/(c*c))*u*zeta_x_hat;
 flt_leftEigenvector[0][2]=((gamma-1.0)/(c*c))*v*zeta_x_hat+(zeta_z_hat)/rho;
 flt_leftEigenvector[0][3]=((gamma-1.0)/(c*c))*w*zeta_x_hat-(zeta_y_hat)/rho;
    flt_leftEigenvector[0][4]=-((gamma-1.0)/(c*c))*zeta_x_hat;

 flt_leftEigenvector[1][0]=zeta_y_hat*(1.0-0.5*(gamma-1.0)*(M*M))-1./rho*(w*zeta_x_hat-u*zeta_z_hat);
 flt_leftEigenvector[1][1]=((gamma-1.0)/(c*c))*u*zeta_y_hat-(zeta_z_hat)/rho;
    flt_leftEigenvector[1][2]=((gamma-1.0)/(c*c))*v*zeta_y_hat;
 flt_leftEigenvector[1][3]=((gamma-1.0)/(c*c))*w*zeta_y_hat+(zeta_x_hat)/rho;
    flt_leftEigenvector[1][4]=-((gamma-1.0)/(c*c))*zeta_y_hat;

 flt_leftEigenvector[2][0]=zeta_z_hat*(1.0-0.5*(gamma-1.0)*(M*M))-1./rho*(u*zeta_y_hat-v*zeta_x_hat);
 flt_leftEigenvector[2][1]=((gamma-1.0)/(c*c))*u*zeta_z_hat+(zeta_y_hat)/rho;
 flt_leftEigenvector[2][2]=((gamma-1.0)/(c*c))*v*zeta_z_hat-(zeta_x_hat)/rho;
    flt_leftEigenvector[2][3]=((gamma-1.0)/(c*c))*w*zeta_z_hat;
    flt_leftEigenvector[2][4]=-((gamma-1.0)/(c*c))*zeta_z_hat;

 flt_leftEigenvector[3][0]=1./(sqrt_2*c)*((gamma-1.0)*(u*u+v*v+w*w)/2.-c*(u*zeta_x_hat+v*zeta_y_hat+w*zeta_z_hat));
    flt_leftEigenvector[3][1]=(-(gamma-1.0)*u/c+zeta_x_hat)/sqrt_2;
    flt_leftEigenvector[3][2]=(-(gamma-1.0)*v/c+zeta_y_hat)/sqrt_2;
    flt_leftEigenvector[3][3]=(-(gamma-1.0)*w/c+zeta_z_hat)/sqrt_2;
    flt_leftEigenvector[3][4]=(gamma-1.)/(sqrt_2*c);

 flt_leftEigenvector[4][0]=1./(sqrt_2*c)*((gamma-1.0)*(u*u+v*v+w*w)/2.+c*(u*zeta_x_hat+v*zeta_y_hat+w*zeta_z_hat));
    flt_leftEigenvector[4][1]=(-(gamma-1.)*u/c-zeta_x_hat)/sqrt_2;
    flt_leftEigenvector[4][2]=(-(gamma-1.)*v/c-zeta_y_hat)/sqrt_2;
    flt_leftEigenvector[4][3]=(-(gamma-1.)*w/c-zeta_z_hat)/sqrt_2;
    flt_leftEigenvector[4][4]=(gamma-1.)/(sqrt_2*c);

}

float GetLambdaMaxInXiDirection(
        struct strct_configuration * pnt_config,
        struct strct_mesh * pnt_mesh,
        struct strct_U * pnt_U_RK,
        int int_acutalI,
        int int_acutalJ,
        int int_acutalK)
{
    int i,j,k,ijk;
    int i_start, i_end;
    i_start=int_acutalI-(pnt_config->int_SpaceOrder-1)/2;
    i_end=int_acutalI+(pnt_config->int_SpaceOrder-1)/2+1;
    j=int_acutalJ;
    k=int_acutalK;
//    float flt_lambdaMaxPlus;
//    float flt_lambdaMaxMinus;
    float flt_lambdaMax;
    float flt_lambdaMaxActual;
    flt_lambdaMaxActual=0.0;
    for (i=i_start; i <= i_end; i++)
    {
 ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

 	 	 flt_lambdaMax=
				 fabs(pnt_U_RK->theta1[ijk])+
				 pnt_U_RK->c[ijk]*
				 sqrt(
				 pnt_mesh->xi_x[ijk]*pnt_mesh->xi_x[ijk]+
				 pnt_mesh->xi_y[ijk]*pnt_mesh->xi_y[ijk]+
				 pnt_mesh->xi_z[ijk]*pnt_mesh->xi_z[ijk]);
//        flt_lambdaMaxPlus=
//                fabs
//                (
//                pnt_U_RK->theta1[ijk]+
//                pnt_U_RK->c[ijk]*
//                sqrt(
//                pnt_mesh->xi_x[ijk]*pnt_mesh->xi_x[ijk]+
//                pnt_mesh->xi_y[ijk]*pnt_mesh->xi_y[ijk]+
//                pnt_mesh->xi_z[ijk]*pnt_mesh->xi_z[ijk])
//                );
//        flt_lambdaMaxMinus=
//                fabs
//                (
//                pnt_U_RK->theta1[ijk]-
//                pnt_U_RK->c[ijk]*
//                sqrt(
//                pnt_mesh->xi_x[ijk]*pnt_mesh->xi_x[ijk]+
//                pnt_mesh->xi_y[ijk]*pnt_mesh->xi_y[ijk]+
//                pnt_mesh->xi_z[ijk]*pnt_mesh->xi_z[ijk])
//                );
//
//        if(flt_lambdaMaxActual<flt_lambdaMaxPlus)
//        {
//            flt_lambdaMaxActual=flt_lambdaMaxPlus;
//        }
//
//        if(flt_lambdaMaxActual<flt_lambdaMaxMinus)
//        {
//            flt_lambdaMaxActual=flt_lambdaMaxMinus;
//        }
        if(flt_lambdaMaxActual<flt_lambdaMax)
        {
            flt_lambdaMaxActual=flt_lambdaMax;
        }
    }
    return(flt_lambdaMaxActual);
}

float GetLambdaMaxInEtaDirection(
        struct strct_configuration * pnt_config,
        struct strct_mesh * pnt_mesh,
        struct strct_U * pnt_U_RK,
        int int_acutalI,
        int int_acutalJ,
        int int_acutalK)
{
    int i,j,k,ijk;
    int j_start, j_end;
    j_start=int_acutalJ-(pnt_config->int_SpaceOrder-1)/2;
    j_end=int_acutalJ+(pnt_config->int_SpaceOrder-1)/2+1;
    i=int_acutalI;
    k=int_acutalK;
//    float flt_lambdaMaxPlus;
//    float flt_lambdaMaxMinus;
    float flt_lambdaMax;
    float flt_lambdaMaxActual;
    flt_lambdaMaxActual=0.0;
    for (j=j_start; j <= j_end; j++)
    {
 ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

         flt_lambdaMax=
                 fabs(pnt_U_RK->theta2[ijk])+
                 pnt_U_RK->c[ijk]*
                 sqrt(
                 pnt_mesh->eta_x[ijk]*pnt_mesh->eta_x[ijk]+
                 pnt_mesh->eta_y[ijk]*pnt_mesh->eta_y[ijk]+
                 pnt_mesh->eta_z[ijk]*pnt_mesh->eta_z[ijk]);
//        flt_lambdaMaxPlus=
//                fabs
//                (
//                pnt_U_RK->theta2[ijk]+
//                pnt_U_RK->c[ijk]*
//                sqrt(
//                pnt_mesh->eta_x[ijk]*pnt_mesh->eta_x[ijk]+
//                pnt_mesh->eta_y[ijk]*pnt_mesh->eta_y[ijk]+
//                pnt_mesh->eta_z[ijk]*pnt_mesh->eta_z[ijk])
//                );
//        flt_lambdaMaxMinus=
//                fabs
//                (
//                pnt_U_RK->theta2[ijk]-
//                pnt_U_RK->c[ijk]*
//                sqrt(
//                pnt_mesh->eta_x[ijk]*pnt_mesh->eta_x[ijk]+
//                pnt_mesh->eta_y[ijk]*pnt_mesh->eta_y[ijk]+
//                pnt_mesh->eta_z[ijk]*pnt_mesh->eta_z[ijk])
//                );
//
//        if(flt_lambdaMaxActual<flt_lambdaMaxPlus)
//        {
//            flt_lambdaMaxActual=flt_lambdaMaxPlus;
//        }
//
//        if(flt_lambdaMaxActual<flt_lambdaMaxMinus)
//        {
//            flt_lambdaMaxActual=flt_lambdaMaxMinus;
//        }
		 if(flt_lambdaMaxActual<flt_lambdaMax)
		 {
			 flt_lambdaMaxActual=flt_lambdaMax;
		 }
    }
    return(flt_lambdaMaxActual);
}

float GetLambdaMaxInZetaDirection(
        struct strct_configuration * pnt_config,
        struct strct_mesh * pnt_mesh,
        struct strct_U * pnt_U_RK,
        int int_acutalI,
        int int_acutalJ,
        int int_acutalK)
{
    int i,j,k,ijk;
    int k_start, k_end;
    k_start=int_acutalK-(pnt_config->int_SpaceOrder-1)/2;
    k_end=int_acutalK+(pnt_config->int_SpaceOrder-1)/2+1;
    i=int_acutalI;
    j=int_acutalJ;
//    float flt_lambdaMaxPlus;
//    float flt_lambdaMaxMinus;
    float flt_lambdaMax;
    float flt_lambdaMaxActual;
    flt_lambdaMaxActual=0.0;
    for (k=k_start; k <= k_end; k++)
    {
 ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

		 flt_lambdaMax=
				 fabs(pnt_U_RK->theta3[ijk])+
				 pnt_U_RK->c[ijk]*
				 sqrt(
				 pnt_mesh->zeta_x[ijk]*pnt_mesh->zeta_x[ijk]+
				 pnt_mesh->zeta_y[ijk]*pnt_mesh->zeta_y[ijk]+
				 pnt_mesh->zeta_z[ijk]*pnt_mesh->zeta_z[ijk]);
//        flt_lambdaMaxPlus=
//                fabs
//                (
//                pnt_U_RK->theta3[ijk]+
//                pnt_U_RK->c[ijk]*
//                sqrt(
//                pnt_mesh->zeta_x[ijk]*pnt_mesh->zeta_x[ijk]+
//                pnt_mesh->zeta_y[ijk]*pnt_mesh->zeta_y[ijk]+
//                pnt_mesh->zeta_z[ijk]*pnt_mesh->zeta_z[ijk])
//                );
//        flt_lambdaMaxMinus=
//                fabs
//                (
//                pnt_U_RK->theta3[ijk]-
//                pnt_U_RK->c[ijk]*
//                sqrt(
//                pnt_mesh->zeta_x[ijk]*pnt_mesh->zeta_x[ijk]+
//                pnt_mesh->zeta_y[ijk]*pnt_mesh->zeta_y[ijk]+
//                pnt_mesh->zeta_z[ijk]*pnt_mesh->zeta_z[ijk])
//                );
//
//        if(flt_lambdaMaxActual<flt_lambdaMaxPlus)
//        {
//            flt_lambdaMaxActual=flt_lambdaMaxPlus;
//        }
//
//        if(flt_lambdaMaxActual<flt_lambdaMaxMinus)
//        {
//            flt_lambdaMaxActual=flt_lambdaMaxMinus;
//        }
        if(flt_lambdaMaxActual<flt_lambdaMax)
        {
            flt_lambdaMaxActual=flt_lambdaMax;
        }
    }
    return(flt_lambdaMaxActual);
}

void Theta_Function_W9(
    struct strct_Flux * pnt_Flux,
    float * pnt_theta,
    int ijk,
    int Plus1ijk,
    int Plus2ijk,
    int Plus3ijk,
    int Plus4ijk,
    int Minus1ijk,
    int Minus2ijk,
    int Minus3ijk)
{
pnt_theta[0]=1./840.*(
            533.*(pnt_Flux->Mass[Plus1ijk]+pnt_Flux->Mass[ijk])
 -139.*(pnt_Flux->Mass[Plus2ijk]+pnt_Flux->Mass[Minus1ijk])
 +29.*(pnt_Flux->Mass[Plus3ijk]+pnt_Flux->Mass[Minus2ijk])
 -3.*(pnt_Flux->Mass[Plus4ijk]+pnt_Flux->Mass[Minus3ijk]));
pnt_theta[1]=1./840.*(
 533.*(pnt_Flux->xiMomentum[Plus1ijk]+pnt_Flux->xiMomentum[ijk])
 -139.*(pnt_Flux->xiMomentum[Plus2ijk]+pnt_Flux->xiMomentum[Minus1ijk])
 +29.*(pnt_Flux->xiMomentum[Plus3ijk]+pnt_Flux->xiMomentum[Minus2ijk])
 -3.*(pnt_Flux->xiMomentum[Plus4ijk]+pnt_Flux->xiMomentum[Minus3ijk]));
pnt_theta[2]=1./840.*(
 533.*(pnt_Flux->etaMomentum[Plus1ijk]+pnt_Flux->etaMomentum[ijk])
 -139.*(pnt_Flux->etaMomentum[Plus2ijk]+pnt_Flux->etaMomentum[Minus1ijk])
 +29.*(pnt_Flux->etaMomentum[Plus3ijk]+pnt_Flux->etaMomentum[Minus2ijk])
 -3.*(pnt_Flux->etaMomentum[Plus4ijk]+pnt_Flux->etaMomentum[Minus3ijk]));
pnt_theta[3]=1./840.*(
 533.*(pnt_Flux->zetaMomentum[Plus1ijk]+pnt_Flux->zetaMomentum[ijk])
 -139.*(pnt_Flux->zetaMomentum[Plus2ijk]+pnt_Flux->zetaMomentum[Minus1ijk])
 +29.*(pnt_Flux->zetaMomentum[Plus3ijk]+pnt_Flux->zetaMomentum[Minus2ijk])
 -3.*(pnt_Flux->zetaMomentum[Plus4ijk]+pnt_Flux->zetaMomentum[Minus3ijk]));
pnt_theta[4]=1./840.*(
 533.*(pnt_Flux->Energy[Plus1ijk]+pnt_Flux->Energy[ijk])
 -139.*(pnt_Flux->Energy[Plus2ijk]+pnt_Flux->Energy[Minus1ijk])
 +29.*(pnt_Flux->Energy[Plus3ijk]+pnt_Flux->Energy[Minus2ijk])
 -3.*(pnt_Flux->Energy[Plus4ijk]+pnt_Flux->Energy[Minus3ijk]));
}

void Theta_Function_W5(
    struct strct_Flux * pnt_Flux,
    float * pnt_theta,
    int ijk,
    int Plus1ijk,
    int Plus2ijk,
    int Minus1ijk)
{
pnt_theta[0]=1./12.*(
            -1.*(pnt_Flux->Mass[Minus1ijk]+pnt_Flux->Mass[Plus2ijk])
            +7.*(pnt_Flux->Mass[ijk]+pnt_Flux->Mass[Plus1ijk]));
pnt_theta[1]=1./12.*(
            -1.*(pnt_Flux->xiMomentum[Minus1ijk]+pnt_Flux->xiMomentum[Plus2ijk])
            +7.*(pnt_Flux->xiMomentum[ijk]+pnt_Flux->xiMomentum[Plus1ijk]));
pnt_theta[2]=1./12.*(
            -1.*(pnt_Flux->etaMomentum[Minus1ijk]+pnt_Flux->etaMomentum[Plus2ijk])
            +7.*(pnt_Flux->etaMomentum[ijk]+pnt_Flux->etaMomentum[Plus1ijk]));
pnt_theta[3]=1./12.*(
            -1.*(pnt_Flux->zetaMomentum[Minus1ijk]+pnt_Flux->zetaMomentum[Plus2ijk])
            +7.*(pnt_Flux->zetaMomentum[ijk]+pnt_Flux->zetaMomentum[Plus1ijk]));
pnt_theta[4]=1./12.*(
            -1.*(pnt_Flux->Energy[Minus1ijk]+pnt_Flux->Energy[Plus2ijk])
            +7.*(pnt_Flux->Energy[ijk]+pnt_Flux->Energy[Plus1ijk]));
}
