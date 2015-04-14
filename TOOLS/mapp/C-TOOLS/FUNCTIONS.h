#ifndef FUNCTIONS_H

#define FUNCTIONS_H

#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */

float calcDUDY(int iMeshPoints, int jMeshPoints, int kMeshPoints, int i, int j, int k, float* u,float* v,float* w,float* p,float* rho, float* x, float* y, float ReynoldsNumber)
{

	int ijk=k*jMeshPoints*iMeshPoints+0*iMeshPoints+i;
	int ijPlus1k=k*jMeshPoints*iMeshPoints+1*iMeshPoints+i;
	int ijPlus2k=k*jMeshPoints*iMeshPoints+2*iMeshPoints+i;
	
	float u_j=sqrt(u[ijk]*u[ijk]+v[ijk]*v[ijk]+w[ijk]*w[ijk]);
	float u_jPlus1=sqrt(u[ijPlus1k]*u[ijPlus1k]+v[ijPlus1k]*v[ijPlus1k]+w[ijPlus1k]*w[ijPlus1k]);
	float u_jPlus2=sqrt(u[ijPlus2k]*u[ijPlus2k]+v[ijPlus2k]*v[ijPlus2k]+w[ijPlus2k]*w[ijPlus2k]);

	//float du=-23./24.*u_j+7./8.*u_jPlus1+1./8.*u_jPlus2-1./24.*u_jPlus3; //Berechnung der Ableitung (3. Ordnung) im Zwischenpunkt der Wand

	float du=2.*u_j;
	float dy=sqrt((x[ijPlus1k]-x[ijk])*(x[ijPlus1k]-x[ijk])+(y[ijPlus1k]-y[ijk])*(y[ijPlus1k]-y[ijk]));
	
	float dy_j=0.5*sqrt((x[ijPlus1k]-x[ijk])*(x[ijPlus1k]-x[ijk])+(y[ijPlus1k]-y[ijk])*(y[ijPlus1k]-y[ijk]));
	float dy_jPlus1=dy_j+sqrt((x[ijPlus1k]-x[ijk])*(x[ijPlus1k]-x[ijk])+(y[ijPlus1k]-y[ijk])*(y[ijPlus1k]-y[ijk]));
	float dy_jPlus2=dy_j+sqrt((x[ijPlus2k]-x[ijk])*(x[ijPlus2k]-x[ijk])+(y[ijPlus2k]-y[ijk])*(y[ijPlus2k]-y[ijk]));
	
	float dudy=du/dy;
	dudy=
	-u_j*(dy_jPlus1+dy_jPlus2)/((dy_j-dy_jPlus1)*(dy_j-dy_jPlus2))
	-u_jPlus1*(dy_j+dy_jPlus2)/((dy_jPlus1-dy_j)*(dy_jPlus1-dy_jPlus2))
	-u_jPlus2*(dy_j+dy_jPlus1)/((dy_jPlus2-dy_j)*(dy_jPlus2-dy_jPlus1));
	return dudy;
}

float calcViscLength(int iMeshPoints, int jMeshPoints, int kMeshPoints, int i, int j, int k, float* u,float* v,float* w,float* p,float* rho, float* x, float* y, float ReynoldsNumber)
{
	int ij0k=k*jMeshPoints*iMeshPoints+0*iMeshPoints+i;

	float mue_w=((1.0+110.4/300.0)*pow((p[ij0k]/rho[ij0k]),1.5))/(p[ij0k]/rho[ij0k]+110.4/300.0);

	float dudy=calcDUDY(iMeshPoints, jMeshPoints, kMeshPoints, i, j, k, u, v, w, p, rho, x, y, ReynoldsNumber);

	float viscLength=mue_w/rho[ij0k]/sqrt(mue_w/rho[ij0k]*dudy)/sqrt(ReynoldsNumber);
	return viscLength;
}

float calcUTau(int iMeshPoints, int jMeshPoints, int kMeshPoints, int i, int j, int k, float* u,float* v,float* w,float* p,float* rho, float* x, float* y, float ReynoldsNumber)
{
	int ij0k=k*jMeshPoints*iMeshPoints+0*iMeshPoints+i;
	int ijk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i;

	float mue_w=((1.0+110.4/300.0)*pow((p[ij0k]/rho[ij0k]),1.5))/(p[ij0k]/rho[ij0k]+110.4/300.0);

	float viscLength=calcViscLength(iMeshPoints, jMeshPoints, kMeshPoints, i, j, k, u, v, w, p, rho, x, y, ReynoldsNumber);
	
	float uTau=mue_w/rho[ij0k]/viscLength/(ReynoldsNumber);
	return uTau;
}

float calcCF(int iMeshPoints, int jMeshPoints, int kMeshPoints, int i, int j, int k, float* u,float* v,float* w,float* p,float* rho, float* x, float* y, float ReynoldsNumber)
{
	int ij0k=k*jMeshPoints*iMeshPoints+0*iMeshPoints+i;;
	int ijk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i;

	float dudy=calcDUDY(iMeshPoints, jMeshPoints, kMeshPoints, i, j, k, u, v, w, p, rho, x, y, ReynoldsNumber);

	float mue_w=((1.0+110.4/300.0)*pow((p[ij0k]/rho[ij0k]),1.5))/(p[ij0k]/rho[ij0k]+110.4/300.0);

	float tau_w=mue_w*dudy;
	
	float vorzeichen=u[ij0k]/fabs(u[ij0k]);

	float cf=vorzeichen*2.*tau_w/ReynoldsNumber;
	printf("vorzeichen: %g cf:%g\t",vorzeichen,cf);
	
	return cf;
}





float * readDelta(int iMeshPoints)
{
	float * Delta;
	char dummy[500];
	int i,ijk;
	int status;
	Delta=calloc(6*iMeshPoints,sizeof(float));

	printf("Beginne Grenzschichtdickenverlaufs-Import\n");
	FILE * file;
	char input_gs[300];
	sprintf(input_gs,"Grenzschichtdickenverlaeufe.dat");
	file=fopen(input_gs,"r");
	fgets(dummy,499,file);
	fgets(dummy,499,file);	
	for(i=0;i<iMeshPoints;i++)
	{
		//                                            x            d99          dStern        theta      Re_theta    position_eta_d99
		fscanf(file," %f %f %f %f %f %f\n",&Delta[6*i+0],&Delta[6*i+1],&Delta[6*i+2],&Delta[6*i+3],&Delta[6*i+4],&Delta[6*i+5]);

	}
	fclose(file);	
	printf("Import fertig\n");
	return Delta;
}

float * calcDelta(int i, int iMeshPoints, int jMeshPoints, int kMeshPoints, 
	float* CoordinateX, float* CoordinateY, float* CoordinateZ, 
	float* VelocityX, float* VelocityY, float* VelocityZ)
{
	float* DeltaCalc;
	DeltaCalc=calloc(4,sizeof(float));
	
	int j,k;
	int iPlus1jk,iMinus1jk;
	int ijPlus1k,ijMinus1k;
	int ijkPlus1,ijkMinus1;

	float jacobian;

	float xi_x,eta_x,zeta_x;
	float xi_y,eta_y,zeta_y;
	float xi_z,eta_z,zeta_z;
	float x_xi,x_eta,x_zeta;
	float y_xi,y_eta,y_zeta;
	float z_xi,z_eta,z_zeta;	

	float udx,udy,udz;
	float vdx,vdy,vdz;
	float wdx,wdy,wdz;

	float u_xi,u_eta,u_zeta;
	float v_xi,v_eta,v_zeta;
	float w_xi,w_eta,w_zeta;		
	float DeltaStern=0.;
	float dn,dx,dy;

	float uPseudo,uPseudo_infty;
	float omega_x,omega_y,omega_z,omega;
	int jMax;
	float delta;
	jMax=0.6*jMeshPoints;
	k=1;
	dn=0.;
	
	for(j=1;j<=jMeshPoints;j++)
	{
		ijPlus1k=k*jMeshPoints*iMeshPoints+(j+1)*iMeshPoints+i;
		ijMinus1k=k*jMeshPoints*iMeshPoints+(j-1)*iMeshPoints+i;
		dn+=0.5*sqrt(
		(CoordinateX[ijPlus1k]-CoordinateX[ijMinus1k])*(CoordinateX[ijPlus1k]-CoordinateX[ijMinus1k])+
		(CoordinateY[ijPlus1k]-CoordinateY[ijMinus1k])*(CoordinateY[ijPlus1k]-CoordinateY[ijMinus1k]));
		if(dn>0.5)
		{
			jMax=j;
			break;
		}
	}

		jMax=0.6*jMeshPoints;
	//printf("jMax:%d\n",jMax);
	for(k=1;k<kMeshPoints-1;k++)
	{
		uPseudo_infty=0.0;
		for(j=1;j<=jMax;j++)
		{
			/****************************************************************************/
			iPlus1jk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i+1;
			iMinus1jk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i-1;
			ijPlus1k=k*jMeshPoints*iMeshPoints+(j+1)*iMeshPoints+i;
			ijMinus1k=k*jMeshPoints*iMeshPoints+(j-1)*iMeshPoints+i;
			ijkPlus1=(k+1)*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
			ijkMinus1=(k-1)*jMeshPoints*iMeshPoints+j*iMeshPoints+i;


			x_xi=(-0.5*CoordinateX[iMinus1jk]+0.5*CoordinateX[iPlus1jk]);
			y_xi=(-0.5*CoordinateY[iMinus1jk]+0.5*CoordinateY[iPlus1jk]);
			z_xi=(-0.5*CoordinateZ[iMinus1jk]+0.5*CoordinateZ[iPlus1jk]);
			x_eta=(-0.5*CoordinateX[ijMinus1k]+0.5*CoordinateX[ijPlus1k]);
			y_eta=(-0.5*CoordinateY[ijMinus1k]+0.5*CoordinateY[ijPlus1k]);
			z_eta=(-0.5*CoordinateZ[ijMinus1k]+0.5*CoordinateZ[ijPlus1k]);
			x_zeta=(-0.5*CoordinateX[ijkMinus1]+0.5*CoordinateX[ijkPlus1]);
			y_zeta=(-0.5*CoordinateY[ijkMinus1]+0.5*CoordinateY[ijkPlus1]);
			z_zeta=(-0.5*CoordinateZ[ijkMinus1]+0.5*CoordinateZ[ijkPlus1]);


			jacobian=x_xi*y_eta*z_zeta+x_eta*y_zeta*z_xi+x_zeta*y_xi*z_eta
				-x_zeta*y_eta*z_xi-x_eta*y_xi*z_zeta-x_xi*y_zeta*z_eta;


			xi_x=(y_eta*z_zeta-y_zeta*z_eta)/jacobian;
			xi_y=(x_zeta*z_eta-x_eta*z_zeta)/jacobian;
			xi_z=(x_eta*y_zeta-x_zeta*y_eta)/jacobian;
			eta_x=(y_zeta*z_xi-y_xi*z_zeta)/jacobian;
			eta_y=(x_xi*z_zeta-x_zeta*z_xi)/jacobian;
			eta_z=(x_zeta*y_xi-x_xi*y_zeta)/jacobian;
			zeta_x=(y_xi*z_eta-y_eta*z_xi)/jacobian;
			zeta_y=(x_eta*z_xi-x_xi*z_eta)/jacobian;
			zeta_z=(x_xi*y_eta-x_eta*y_xi)/jacobian;

			u_xi=(-0.5*VelocityX[iMinus1jk]+0.5*VelocityX[iPlus1jk]);
			v_xi=(-0.5*VelocityY[iMinus1jk]+0.5*VelocityY[iPlus1jk]);
			w_xi=(-0.5*VelocityZ[iMinus1jk]+0.5*VelocityZ[iPlus1jk]);
			u_eta=(-0.5*VelocityX[ijMinus1k]+0.5*VelocityX[ijPlus1k]);
			v_eta=(-0.5*VelocityY[ijMinus1k]+0.5*VelocityY[ijPlus1k]);
			w_eta=(-0.5*VelocityZ[ijMinus1k]+0.5*VelocityZ[ijPlus1k]);
			u_zeta=(-0.5*VelocityX[ijkMinus1]+0.5*VelocityX[ijkPlus1]);
			v_zeta=(-0.5*VelocityY[ijkMinus1]+0.5*VelocityY[ijkPlus1]);
			w_zeta=(-0.5*VelocityZ[ijkMinus1]+0.5*VelocityZ[ijkPlus1]);

			udx=u_xi*xi_x+u_eta*eta_x+u_zeta*zeta_x;
			udy=u_xi*xi_y+u_eta*eta_y+u_zeta*zeta_y;
			udz=u_xi*xi_z+u_eta*eta_z+u_zeta*zeta_z;
			vdx=v_xi*xi_x+v_eta*eta_x+v_zeta*zeta_x;
			vdy=v_xi*xi_y+v_eta*eta_y+v_zeta*zeta_y;
			vdz=v_xi*xi_z+v_eta*eta_z+v_zeta*zeta_z;
			wdx=w_xi*xi_x+w_eta*eta_x+w_zeta*zeta_x;
			wdy=w_xi*xi_y+w_eta*eta_y+w_zeta*zeta_y;
			wdz=w_xi*xi_z+w_eta*eta_z+w_zeta*zeta_z;

			omega_z = 0.5 * (udy-vdx);
			omega_y = 0.5 * (udz-wdx);
			omega_x = 0.5 * (vdz-wdy);
		
			omega=sqrt(omega_x*omega_x+omega_y*omega_y+omega_z*omega_z);
		
			dn=0.5*sqrt(
			(CoordinateX[ijPlus1k]-CoordinateX[ijMinus1k])*(CoordinateX[ijPlus1k]-CoordinateX[ijMinus1k])+
			(CoordinateY[ijPlus1k]-CoordinateY[ijMinus1k])*(CoordinateY[ijPlus1k]-CoordinateY[ijMinus1k]));
		
			uPseudo_infty+=omega*dn;
		}	
	
		uPseudo=0.;
		for(j=1;j<=jMax;j++)
		{
			/****************************************************************************/
			iPlus1jk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i+1;
			iMinus1jk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i-1;
			ijPlus1k=k*jMeshPoints*iMeshPoints+(j+1)*iMeshPoints+i;
			ijMinus1k=k*jMeshPoints*iMeshPoints+(j-1)*iMeshPoints+i;
			ijkPlus1=(k+1)*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
			ijkMinus1=(k-1)*jMeshPoints*iMeshPoints+j*iMeshPoints+i;


			x_xi=(-0.5*CoordinateX[iMinus1jk]+0.5*CoordinateX[iPlus1jk]);
			y_xi=(-0.5*CoordinateY[iMinus1jk]+0.5*CoordinateY[iPlus1jk]);
			z_xi=(-0.5*CoordinateZ[iMinus1jk]+0.5*CoordinateZ[iPlus1jk]);
			x_eta=(-0.5*CoordinateX[ijMinus1k]+0.5*CoordinateX[ijPlus1k]);
			y_eta=(-0.5*CoordinateY[ijMinus1k]+0.5*CoordinateY[ijPlus1k]);
			z_eta=(-0.5*CoordinateZ[ijMinus1k]+0.5*CoordinateZ[ijPlus1k]);
			x_zeta=(-0.5*CoordinateX[ijkMinus1]+0.5*CoordinateX[ijkPlus1]);
			y_zeta=(-0.5*CoordinateY[ijkMinus1]+0.5*CoordinateY[ijkPlus1]);
			z_zeta=(-0.5*CoordinateZ[ijkMinus1]+0.5*CoordinateZ[ijkPlus1]);


			jacobian=x_xi*y_eta*z_zeta+x_eta*y_zeta*z_xi+x_zeta*y_xi*z_eta
				-x_zeta*y_eta*z_xi-x_eta*y_xi*z_zeta-x_xi*y_zeta*z_eta;


			xi_x=(y_eta*z_zeta-y_zeta*z_eta)/jacobian;
			xi_y=(x_zeta*z_eta-x_eta*z_zeta)/jacobian;
			xi_z=(x_eta*y_zeta-x_zeta*y_eta)/jacobian;
			eta_x=(y_zeta*z_xi-y_xi*z_zeta)/jacobian;
			eta_y=(x_xi*z_zeta-x_zeta*z_xi)/jacobian;
			eta_z=(x_zeta*y_xi-x_xi*y_zeta)/jacobian;
			zeta_x=(y_xi*z_eta-y_eta*z_xi)/jacobian;
			zeta_y=(x_eta*z_xi-x_xi*z_eta)/jacobian;
			zeta_z=(x_xi*y_eta-x_eta*y_xi)/jacobian;

			u_xi=(-0.5*VelocityX[iMinus1jk]+0.5*VelocityX[iPlus1jk]);
			v_xi=(-0.5*VelocityY[iMinus1jk]+0.5*VelocityY[iPlus1jk]);
			w_xi=(-0.5*VelocityZ[iMinus1jk]+0.5*VelocityZ[iPlus1jk]);
			u_eta=(-0.5*VelocityX[ijMinus1k]+0.5*VelocityX[ijPlus1k]);
			v_eta=(-0.5*VelocityY[ijMinus1k]+0.5*VelocityY[ijPlus1k]);
			w_eta=(-0.5*VelocityZ[ijMinus1k]+0.5*VelocityZ[ijPlus1k]);
			u_zeta=(-0.5*VelocityX[ijkMinus1]+0.5*VelocityX[ijkPlus1]);
			v_zeta=(-0.5*VelocityY[ijkMinus1]+0.5*VelocityY[ijkPlus1]);
			w_zeta=(-0.5*VelocityZ[ijkMinus1]+0.5*VelocityZ[ijkPlus1]);

			udx=u_xi*xi_x+u_eta*eta_x+u_zeta*zeta_x;
			udy=u_xi*xi_y+u_eta*eta_y+u_zeta*zeta_y;
			udz=u_xi*xi_z+u_eta*eta_z+u_zeta*zeta_z;
			vdx=v_xi*xi_x+v_eta*eta_x+v_zeta*zeta_x;
			vdy=v_xi*xi_y+v_eta*eta_y+v_zeta*zeta_y;
			vdz=v_xi*xi_z+v_eta*eta_z+v_zeta*zeta_z;
			wdx=w_xi*xi_x+w_eta*eta_x+w_zeta*zeta_x;
			wdy=w_xi*xi_y+w_eta*eta_y+w_zeta*zeta_y;
			wdz=w_xi*xi_z+w_eta*eta_z+w_zeta*zeta_z;

			omega_z = 0.5 * (udy-vdx);
			omega_y = 0.5 * (udz-wdx);
			omega_x = 0.5 * (vdz-wdy);
		
			omega=sqrt(omega_x*omega_x+omega_y*omega_y+omega_z*omega_z);
		
			dn=0.5*sqrt(
			(CoordinateX[ijPlus1k]-CoordinateX[ijMinus1k])*(CoordinateX[ijPlus1k]-CoordinateX[ijMinus1k])+
			(CoordinateY[ijPlus1k]-CoordinateY[ijMinus1k])*(CoordinateY[ijPlus1k]-CoordinateY[ijMinus1k]));
		
			uPseudo+=omega*dn;
		
			if(uPseudo<0.98*uPseudo_infty){DeltaCalc[0]+=dn/(kMeshPoints-2);DeltaCalc[3]=j;}
			DeltaCalc[1]+=(1.-uPseudo/uPseudo_infty)*dn/(kMeshPoints-2);
			DeltaCalc[2]+=uPseudo/uPseudo_infty*(1.-uPseudo/uPseudo_infty)*dn/(kMeshPoints-2);
		}
	}

	return DeltaCalc;
}

char *replace_str(char *str, char *orig, char *rep)
{
  static char buffer[4096];
  char *p;

  if(!(p = strstr(str, orig)))  // Is 'orig' even in 'str'?
    return str;

  strncpy(buffer, str, p-str); // Copy characters from 'str' start to 'orig' st$
  buffer[p-str] = '\0';

  sprintf(buffer+(p-str), "%s%s", rep, p+strlen(orig));

  return buffer;
}

#endif /* FUNCTIONS_H */
