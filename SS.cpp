//PROJET SYSTEME SOLAIRE
//Par David Bourgeois
//libre de droit avec attribution et partage dans les mêmes conditions (CC BY-SA)

#include "header.h"

// TOUT EST EN SI!!!!!!
//variables globales
long double p[dim*5+2][np];
//p = 0 x, 1 y, 2 z, 3 vx, 4 vy, 5 vz, 6 m, 7 r, 8 xi, 9 yi, 10 zi, 11 x-1, 12 y-1, 13 z-1, 14 vx-1, 15 vy-1, 16 vz-1

int main(){
	char nom[12][np];
	unsigned short i, j, k, demitour[np];
	unsigned long long int bip[2];
	double elli[4][np][tourmax], ddt[etape][etape], m[etape][2], t, dt[3], time0, tmax, e[3];
	//elli = 0 période, 1 ex, 2 périastre, 3 apoastre

    ofstream coordonee("coo.res");

	//initialisation
	bip[0]=0;
	bip[1]=0;
	for (j= 0; j< np; j++){
		demitour[j]=0;
		for (k= 0; k< tourmax; k++){
			elli[2][j][k]=pow(10,40);
			for (i= 0; i< 4; i++){
				if(i==2){continue;}
				elli[i][j][k]=0;
			}
		}
	}

	EXTRAIRE(nom,ddt,m);
	tmax=PARAMETRAGE(dt,e);
	time0=clock();

	for( t = 0; t < tmax;){
		RKN(dt,ddt,m,e);
		t+=dt[0];
		ELLIPSE(elli, demitour, t);
		PROGRESSION(t,tmax,bip,time0,nom,elli,demitour);
	}

	coordonee.close();
	system("gnuplot -p gp.conf");
	RAPPORT(nom,elli,t,demitour,bip);

	return 42;
}

void EXTRAIRE(char nom[12][np], double ddt[etape][etape], double m[etape][2]){
	char mot[12];
	unsigned short i, j;
	double q;

	ifstream coo("cooi.res");
	for( j = 0; j < np; j++){
		for( i = 0; i < 2*dim; i++){
			//fscanf(coordonee, "%lf", &p[i][j]);
			coo>>p[i][j];
			p[i][j]*=1000;
			p[i+11][j]=p[i][j];
			if(i>2){continue;}
			p[i+8][j]=p[i][j];
		}
	}
	coo.close();

	ifstream planete("planete.res");
	for( j = 0; j < np; j++){
		planete>>mot>>p[6][j]>>p[7][j];
		for( i = 0; i < 2*dim; i++){nom[i][j]=mot[i];}
	}
	planete.close();

	ifstream butcher("butcher.res");
	for (i= 0; i< etape; i++){
		for (j= 0; j< etape; j++){
			ddt[i][j]=0;
		}
	}
	ddt[0][0] = 1.;
	for( i = 1; i < etape; i++){
		for( j = 1; j < i+1; j++){
			butcher>>ddt[i][j]>>mot[i]>>q;
			ddt[i][j]/=q;
		}
	}
	for( j = 0; j < 2; j++){
		for( i = 0; i < etape; i++){
			butcher>>m[i][j]>>mot[i]>>q;
			m[i][j]/=q;
		}
	}
	butcher.close();
}

double PARAMETRAGE(double dt[3], double e[3]){
	double tmax;

	cout<<"Donnez le temps de simulation désiré en années terreste:"<<endl;
	cin>>tmax;
	tmax*=(3600*24*365.25636042);
	e[0]=0;
	e[2]=7626.2/(tmax*tmax); //exp(-pow(10,-9)*tmax-28.); marche bien pour pour Dormand–Prince
	e[1]=0.1*e[2];
	dt[0]=1000.;
	dt[1]=e[2]*pow(10,14);
	dt[2]=tmax*pow(10,-5);

	return tmax;
}

//Calcul différentiel
void RKN(double dt[3], const double ddt[etape][etape], const double m[etape][2], double e[3]){
    const unsigned short n=dim*2;
	unsigned short i, j, k, l;
	long double p2[n+2][np], dp[n][np][etape];
	//dp =0 vx, 1 vy, 2 vz, 3 ax, 4 ay, 5 az

	//initialisation p2
	for (j= 0; j< np; j++){
		for( i = 0; i< n+2; i++){
			p2[i][j]=p[i][j];
		}
	}

	//Pas du RK
	for (k= 0; k< etape; k++){
		PFD(p2, dp, k);
		if(k==etape){break;}
		for (j= 0; j< np; j++){
			for( i = 0; i< n; i++){
				p2[i][j]=p[i][j];
				for (l= 0; l< k; l++){
					p2[i][j]+=dp[i][j][l]*ddt[k][l]*dt[0];
				}
			}
		}
	}

	//reset p2
	for (j= 0; j< np; j++){
		for( i = 0; i< n+2; i++){
			p2[i][j]=p[i][j];
		}
	}

    e[0]=0;
	//calcul de p,p2 final et erreur max
	for (j= 0; j< np; j++){
		for( i = 0; i< n; i++){
			for (k= 0; k< etape; k++){
				p[i][j]+=dt[0]*dp[i][j][k]*m[k][0];
				p2[i][j]+=dt[0]*dp[i][j][k]*m[k][1];
			}
            if(abs(p[i][j]-p2[i][j]/p[i][j])>e[0]){e[0]=abs((p[i][j]-p2[i][j])/p[i][j]);}
		}
	}
	ERREUR(dt,ddt,m,e);
}

void ERREUR(double dt[3], const double ddt[etape][etape], const double m[etape][2], double e[3]){
	unsigned short i, j;

	if(dt[0]<dt[1]){
		dt[0]=dt[1];
		return;
	}
	if(dt[0]>dt[2]){
		dt[0]=dt[2];
		return;
	}

	if(e[0]>e[2]){
		//rembobinage
		for (j= 0; j< np; j++){
			for( i = 0; i< 2*dim; i++){
				p[i][j]=p[i+3*dim+2][j];
			}
		}
		dt[0]/=2;
		RKN(dt,ddt,m,e);
	}
	if(e[0]<e[1]){dt[0]*=2;}
}

void PFD(const long double p2[2*dim+2][np], long double dp[2*dim][np][etape], const unsigned short k){
	unsigned short i, j;
	long double x, y, z, fg;

	//vitesse + initialisation acceleration 
	for( j = 0; j< np; j++){
		for( i = 0; i< dim; i++){
			dp[0+i][j][k]=p2[dim+i][j];
			dp[dim+i][j][k]=0;
		}

		//accelerations Fi/j
		for(i = 0; i< np; i++){
			if(i==j){continue;}
			x=p2[0][j]-p2[0][i];
			y=p2[1][j]-p2[1][i];
			z=p2[2][j]-p2[2][i];
			fg= -G * (p2[6][i]) / pow(NORME(x, y, z), 3);
			dp[3][j][k]+=fg*x;
			dp[4][j][k]+=fg*y;
			dp[5][j][k]+=fg*z;
		}
	}
}

void ELLIPSE(double elli[4][np][tourmax], unsigned short demitour[], const double t){
	unsigned short j, k, l;
	double x, y, z, r;

	for( j = 1; j< np; j++){
		k=short(demitour[j]/2);
		x=p[0][j]-p[0][0];
		y=p[1][j]-p[1][0];
		z=p[2][j]-p[2][0];
		r=NORME(x,y,z);

		if(elli[2][j][k]>r){elli[2][j][k]=r;}
		if(elli[3][j][k]<r){elli[3][j][k]=r;}
		elli[1][j][k]=(elli[3][j][k]-elli[2][j][k])/(elli[3][j][k]+elli[2][j][k]);
		
		if(TOUR(j,demitour)){
			elli[0][j][k]=t;
			for( l = 0; l<k-1; l++){
				elli[0][j][k]-=elli[0][j][l];
			}
		}
	}
}

bool TOUR(const unsigned short j, unsigned short demitour[]){
	if((p[0][j]-p[8][j])*(p[11][j]-p[8][j])<0){
		demitour[j]++;
		if(demitour[j]%2==0){
			return(true);
		}
	}
	return(false);
}

void PROGRESSION(const double t,const double tmax, unsigned long long int bip[2], const double time0, const char nom[12][np], const double elli[4][np][tourmax], const unsigned short demitour[]){
	short i, j;
	double time1;

    ofstream coordonee("coo.res",ios_base::app);
    bip[0]++;
    if(bip[0]==10000000000){
        bip[0]=0;
        bip[1]++;
    }

	//petit rapport pour les impatients
	if(bip[0]%100000==0){
		system("gnuplot -p gp.conf");
		time1 =((tmax/t)-1)*((clock()-time0)/1000000.);
		cout<<"##############  "<<(100*t/tmax)<<" %  ##############"<<endl;
		cout<<"Temps restant estimé= "<<time1<<" s"<<endl;
	}

    if(bip[0]%10000000==0){RAPPORT(nom,elli,t,demitour,bip);}

	//enregistrement mouvement précédent + écriture fichier
	for (j= 0; j< np; j++){
		for( i = 0; i< 2*dim; i++){
			p[i+3*dim+2][j]=p[i][j];
            if(bip[0] % 2500 == 0){coordonee << p[i][j] - p[i][0] << " ";}
		}
        if(bip[0] % 2500 == 0){coordonee << endl;}
	}

}

void RAPPORT(const char nom[12][np], const double elli[4][np][tourmax], const double t, const unsigned short demitour[], const unsigned long long int bip[2]){
	unsigned short i, j, k, l;
	double s[4][np], s2[4][np], bipbip;

	bipbip=bip[0]+pow(10,10)*bip[1];
	//initialisation
	for( j = 0; j< np; j++){
		for( i = 0; i< 4; i++){
			s[i][j]=0;
			s2[i][j]=0;
		}
	}

    ofstream rapport("rapport.txt");
	cout<<endl<<"####################################################"<<endl<<endl;
	cout<<"##############  RAPPORT D'EXECUTION  ##############"<<endl<<endl;
	cout<<"####################################################"<<endl<<endl;
	cout<<"Le programme s'est executé sur une période de "<<t<<" "<<"s"<<" et a demandé "<<bipbip<<" itérations"<<endl;
    rapport<<"####################################################"<<endl<<endl;
    rapport<<"##############  RAPPORT D'EXECUTION  ##############"<<endl<<endl;
    rapport<<"####################################################"<<endl<<endl;
    rapport<<"Le programme s'est executé sur une période de "<<t<<" "<<"s"<<" et a demandé "<<bipbip<<" itérations"<<endl;

	cout<<"NOM\tpériode\t\t\texcentricité\t\tpériastre\t\t\tapoastre"<<endl;
    rapport<<"NOM\tpériode\t\t\texcentricité\t\tpériastre\t\t\tapoastre"<<endl;
	for( j = 0; j< np; j++){
		for( i = 0; i< 4; i++){
			for( k = 0; k< demitour[j]/2; k++){
				s[i][j]+=2*elli[i][j][k]/demitour[j];
				s2[i][j]+=2*elli[i][j][k]*elli[i][j][k]/demitour[j];
			}
			s2[i][j]=pow(s2[i][j]-s[i][j]*s[i][j],0.5);
			if(i==0){
			    for( l = 0; l< 12; l++){
			        cout<<nom[l][j];
                    rapport<<nom[l][j];
			    }
			}
			cout<<"\t"<<s[i][j]<<" ± "<<s2[i][j];
            rapport<<"\t"<<s[i][j]<<" ± "<<s2[i][j];
		}
		cout<<endl;
        rapport<<endl;
	}
    cout<<endl<<"####################################################"<<endl<<endl;
    rapport.close();
}

long double NORME(double x, double y, double z){

	return pow(x*x+y*y+z*z,0.5);
}