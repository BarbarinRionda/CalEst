/*
+------------------------------------------------------------------------------+
¦       SOLVERS DE ALTO RENDIMIENTO PARA SISTEMAS DE ECUACIONES LINEALES.      ¦
+------------------------------------------------------------------------------+
*/

#include <math.h>

#include "Solvers.h"

/*
+------------------------------------------------------------------------------+
¦      SOLUCION DE SISTEMAS DE ECUACIONES POR GAUSS - JORDAN CON PIVOTEO.      ¦
+------------------------------------------------------------------------------+
*/
//******************** FUNCIONES DE OBJETO: GaussJordan ************************
//------------------------------------------------------------------------------
// Constructor de objeto.
GaussJordan::GaussJordan(unsigned num_ecua, unsigned num_term, double **mat_coef,
	double **term_ind, double tolerancia)
{
	register unsigned i, j;

	// Asigna valores:
	numecua = num_ecua;
	numterm = num_term;
	toler = tolerancia;
	coeftes = 0;
	termind = 0;
	
	// Asigna memoria para los coeficientes:
	coeftes = new double* [numecua];
	for(i=0; i<numecua; i++) coeftes[i] = new double [numecua];

	// Asigna memoria para los términos independientes:
	termind = new double* [numterm];
	for(i=0; i<numterm; i++) termind[i] = new double [numecua];

	// Asigna valores a la matriz de coeficientes:
	for(i=0; i<numecua; i++)
		for(j=0; j<numecua; j++)
			coeftes[i][j] = mat_coef[i][j];

	// Asigna valores a los vectores de términos independientes:
	for(i=0; i<numterm; i++)
		for(j=0; j<numecua; j++)
			termind[i][j] = term_ind[i][j];
}

//------------------------------------------------------------------------------
// Destructor de objeto.
GaussJordan::~GaussJordan(void)
{
	register unsigned i;
	
	if(coeftes) {
		for(i=0; i<numecua; i++) delete [] coeftes[i];
		delete [] coeftes;
	}
}

//------------------------------------------------------------------------------
// Resuelve el sistema de ecuaciones.
double **GaussJordan::Solucion(void)
{
	// Regresa el puntero al vector de soluciones si no ocurrio ningún problema.
	// Regresa 0 si no se pudo resolver el sistema.

	register unsigned i, j, k;
	double pivote, c;

	for(i=0; i<numecua; i++) {
		// Verifica que el pivote sea diferente de cero:
		if(fabs(coeftes[i][i]) <= toler)
			if(!Pivoteo(i)) return 0;

		pivote = coeftes[i][i];

		// Divide el renglón i-ésimo entre el pivote:
		for(j=0; j<numecua; j++) coeftes[i][j] /= pivote;

		for(j=0; j<numterm; j++) termind[j][i] /= pivote;

		// Hace ceros en la columna del pivote:
		for(j=0; j<numecua; j++) {
			if(i != j) {
				c = coeftes[j][i];

				for(k=0; k<numecua; k++) 
					coeftes[j][k] -= c*coeftes[i][k];

				for(k=0; k<numterm; k++)
					termind[k][j] -= c*termind[k][i];
			}
		}
	}

	return termind;
}

//------------------------------------------------------------------------------
// Intercambia renglones hasta encontrar un pivote != 0.
bool GaussJordan::Pivoteo(unsigned linea)
{
	// Regresa TRUE si se pivoteo con éxito.
	// Regresa FALSE si ya no se pudo pivotear.

	register unsigned i, j;
	double t, *temp;

	for(i=linea+1; i<numecua; i++) {
		if(fabs(coeftes[i][linea]) > toler) {
			temp = coeftes[linea];
			coeftes[linea] = coeftes[i];
			coeftes[i] = temp;

			for(j=0; j<numterm; j++) {
				t = termind[j][linea];
				termind[j][linea] = termind[j][i];
				termind[j][i] = t;
			}

			return true;
		}
	}

	return false;
}

/*
+------------------------------------------------------------------------------+
¦   SOLUCION DE SISTEMAS DE ECUACIONES POR FACTORIZACION DE CROUT CON PIVOTEO. ¦
+------------------------------------------------------------------------------+
*/
//******************** FUNCIONES DE OBJETO: FactCrout **************************
//------------------------------------------------------------------------------
// Constructor del objeto.
FactCrout::FactCrout(unsigned num_ecua, unsigned num_term, double **mat_coef, 
	double **term_ind, double tolerancia)
{
	register unsigned i, j;

	// Asigna valores:
	numecua = num_ecua;
	numterm = num_term;
	toler = tolerancia;
	coeftes = 0;
	termind = 0;

	// Asigna memoria para los coeficientes:
	coeftes = new double* [numecua];
	for(i=0; i<numecua; i++) coeftes[i] = new double [numecua];

	// Asigna memoria para los términos independientes:
	termind = new double* [numterm];
	for(i=0; i<numterm; i++) termind[i] = new double [numecua];

	// Asigna valores a la matriz de coeficientes:
	for(i=0; i<numecua; i++)
		for(j=0; j<numecua; j++)
			coeftes[i][j] = mat_coef[i][j];

	// Asigna valores a los vectores de términos independientes:
	for(i=0; i<numterm; i++)
		for(j=0; j<numecua; j++)
			termind[i][j] = term_ind[i][j];
}

//------------------------------------------------------------------------------
// Destructor del objeto.
FactCrout::~FactCrout(void)
{
	register unsigned i;
	
	if(coeftes) {
		for(i=0; i<numecua; i++) delete [] coeftes[i];
		delete [] coeftes;
	}
}

//------------------------------------------------------------------------------
// Resuelve el sistema de ecuaciones.
double **FactCrout::Solucion(void)
{
	// Regresa el puntero al vector de soluciones si no ocurrio ningún problema.
	// Regresa 0 si no se pudo resolver el sistema.

	register unsigned i, j, k;
	double sum, *suma;

	if(!(suma = new double [numterm])) return 0;

	for(i=0; i<numecua-1; i++) {
		if(fabs(coeftes[0][0]) <= toler)
			if(!Pivoteo(0)) return 0;

		coeftes[0][i+1] /= coeftes[0][0];

		for(j=1; j<=i; j++) {
			sum = 0.0;

			for(k=0; k<=j-1; k++)
				sum += coeftes[j][k]*coeftes[k][i+1];
	
			if(!coeftes[j][j]) return 0;

			coeftes[j][i+1] = (coeftes[j][i+1]-sum)/coeftes[j][j];
		}
		
		for(j=1; j<=i; j++) {
			sum = 0.0;

			for(k=0; k<=j-1; k++)
				sum += coeftes[k][j]*coeftes[i+1][k];
			
			coeftes[i+1][j] -= sum;
		}
				
		sum = 0.0;

		for(j=0; j<=i; j++)
			sum += coeftes[i+1][j]*coeftes[j][i+1];

		coeftes[i+1][i+1] -= sum;
	}
	
	if(fabs(coeftes[0][0]) <= toler)
		if(!Pivoteo(0)) return 0;

	for(i=0; i<numterm; i++)
		termind[i][0] /= coeftes[0][0];    ///////////////////
	
	for(i=1; i<numecua; i++) {
		for(j=0; j<numterm; j++) suma[j] = 0.0;
		
		for(j=0; j<i; j++)
			for(k=0; k<numterm; k++)
				suma[k] += coeftes[i][j]*termind[k][j]; //////////////
		
		if(coeftes[i][i] <= toler)
			if(!Pivoteo(i)) return 0;

		for(j=0; j<numterm; j++)
			termind[j][i] = (termind[j][i]-suma[j])/coeftes[i][i]; ////////////
	}

	for(i=(int)(numecua-2); (int)i>=0; (int)i--) {
		for(j=0; j<numterm; j++) suma[j] = 0.0;

		for(j=(int)i+1; j<=numecua-1; j++)
			for(k=0; k<numterm; k++)
				suma[k] += coeftes[(int)i][j]*termind[k][j];     ////////////////

		for(j=0; j<numterm; j++)
			termind[j][i] -= suma[j];               ///////////////////
	}

	return termind;
}


//------------------------------------------------------------------------------
// Intercambia renglones hasta encontrar un pivote != 0.
bool FactCrout::Pivoteo(unsigned linea)
{
	// Regresa TRUE si se pivoteo con éxito.
	// Regresa FALSE si ya no se pudo pivotear.

	register unsigned i, j;
	double t, *temp;

	for(i=linea+1; i<numecua; i++) {
		if(fabs(coeftes[i][linea]) > toler) {
			temp = coeftes[linea];
			coeftes[linea] = coeftes[i];
			coeftes[i] = temp;

			for(j=0; j<numterm; j++) {
				t = termind[j][linea];
				termind[j][linea] = termind[j][i];
				termind[j][i] = t;
			}

			return true;
		}
	}

	return false;
}