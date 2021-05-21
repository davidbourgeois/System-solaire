// Header programme système solaire

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include "definition.h"

using namespace std;

//Constantes globales
const auto G= 6.67408 * pow(10, -11);

//prototypes fonctions
long double NORME(double x, double y, double z);
void RAPPORT(const char nom[12][np], const double elli[4][np][tourmax], double t, const unsigned short demitour[], const unsigned long long int bip[2]);
void PROGRESSION(double t, double tmax, unsigned long long int bip[2], double time0, const char nom[12][np], const double elli[4][np][tourmax], const unsigned short demitour[]);
bool TOUR(unsigned short j, unsigned short demitour[]);
void ELLIPSE(double elli[4][np][tourmax], unsigned short demitour[], double t);
void PFD(const long double p2[2*dim+2][np], long double dp[2*dim][np][etape], unsigned short k);
void ERREUR(double dt[3], const double ddt[etape][etape], const double m[etape][2], double e[3]);
void RKN(double dt[3], const double ddt[etape][etape], const double m[etape][2], double e[3]);
double PARAMETRAGE(double dt[3], double e[3]);
void EXTRAIRE(char nom[12][np], double ddt[etape][etape], double m[etape][2]);