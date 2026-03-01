/*
+------------------------------------------------------------------------------+
¦                         S I C B A S A   d e   C. V.                          ¦
¦                                                                              ¦
¦                              PROYECTO "CALEST"                               ¦
¦                                                                              ¦
¦                 C A L C U L O   D E   E S T R U C T U R A S                  ¦
¦                                                                              ¦
¦    Programa de análisis estructural por el método matricial de rigideces,    ¦
¦   obtiene la deformación que sufre una estructura como consecuencia de la    ¦
¦                              aplicación de cargas.                           ¦
¦          .                                                                   ¦
¦                                                                              ¦
¦                  Realizó: Armando Barbarín Rionda Ramírez.                   ¦
¦                                                                              ¦
¦ (2000)                                                           Versión 1.0 ¦
+------------------------------------------------------------------------------+
*/

#include <fstream.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "AlgLineal.h"
#include "CalEst.h"
#include "Solvers.h"

// Definición de funciones varias:
void Remplaza(char *cad, char quita, char pon);

/*
+------------------------------------------------------------------------------+
¦                     F U N C I Ó N   P R I N C I P A L                        ¦
+------------------------------------------------------------------------------+
*/
//------------------------------------------------------------------------------
int main(int argv, char *argc[])
{
	// Regresa 0 si no hubo ningún error.
	// Regresa -1 si hubo algún error al cargar datos.
	// Regresa -2 si hubo algún error al realizar el cálculo.
	// Regresa -3 si hubo algún error al escribir archivo de resultados.
	// Regresa -4 si hubo algún error al escribir archivo de reporte HTML.

	if(argv < 2) {
		cout << "Use: CALEST nombre.ext\n";
		return -1;
	}
	
	// Construye objeto:
	CalEst est;
	
	if(est.Set(argc[1])) {
		cout << "Error al cargar datos del problema.\n";

		// Verifica si se ocacionó un error:
		if(est.error.estado) {
			cout << "\nError: " << est.error.codigo << "\n";
			if(est.error.mensaje[0]) cout << "Mensaje: " << est.error.mensaje << "\n";
			cout << "\n\n";
		}

		return -1;
	}

/*	if(est.Set(3, 3, 2, 1, 2, 0, 0, 2, true, "Ejemplo de armadura", "Cm", "Kg"))
		return -1;

	// Características de las barras:
	est.barra[0].Set(0, 2, 1, 1);
	est.barra[1].Set(0, 1, 1, 0);
	est.barra[2].Set(1, 2, 1, 0);

	// Coordenados de los nodos:
	est.nodo[0].Set(0, 0);
	est.nodo[1].Set(0, 600);
	est.nodo[2].Set(800, 0);

	// Restricciones:
	est.apoyo[0].Set(0, true, 0, false, 0);
	est.apoyo[1].Set(1, true, 0, true, 0);

	// Materiales:
	est.mate[0].Set(2.1e6, 4.516e1, 0, 1);
	est.mate[1].Set(2.1e6, 9.032e1, 0, 1);

	// Casos de Carga:
	if(est.CC[0].Set(1, 0, 0, false, est.nbarras, est.ngdln_max, est.ngdle_max,
		"Fuerzas en nodos")) return -1;
	
	// Cargas en nodos:
	est.CC[0].fnodo[0].Set(2, 0, -10000);
*/
	if(est.Calcula(CAL_TODO, 1e-8, true, 1)) {
		cout << "Proceso de cálculo abortado.\n";

		// Verifica si se ocaciono un error:
		if(est.error.estado) {
			cout << "\nError: " << est.error.codigo << "\n";
			if(est.error.mensaje[0]) cout << "Mensaje: " << est.error.mensaje << "\n";
			cout << "\n\n";
		}
		
		return -2;
	}

	if(est.EscribeRES("results.txt")) {
		cout << "Escritura de archivo de resultados fallida.\n";

		// Verifica si se ocaciono un error:
		if(est.error.estado) {
			cout << "\nError: " << est.error.codigo << "\n";
			if(est.error.mensaje[0]) cout << "Mensaje: " << est.error.mensaje << "\n";
			cout << "\n\n";
		}

		return -3;
	}

	return 0;
}

/*
+------------------------------------------------------------------------------+
¦                  F U N C I O N E S   D E   O B J E T O S                     ¦
+------------------------------------------------------------------------------+
*/
//********************** FUNCIONES DE OBJETO: CalEst ***************************
//------------------------------------------------------------------------------
// Constructor del objeto.
CalEst::CalEst(void)
{
	// Inicializa valores:
	mate  = 0;
	nodo  = 0;
	barra = 0;
	apoyo = 0;
	elas  = 0;
	incl  = 0;
	CC    = 0;
}

//------------------------------------------------------------------------------
// Destructor del objeto.
CalEst::~CalEst(void)
{
	Destruye();
}

//------------------------------------------------------------------------------
// Afecta el vector de la matriz ensamblada con los apoyos elásticos.
void CalEst::ApoyosElasticos(void)
{
	if(!nelas) return;
	
	register unsigned i, j, cont, renglon;

	for(i=0; i<nelas; i++) {
		cont = 0;
		renglon = (elas[i].nodo+1)*ngdln_max-(ngdln_max-0);

		for(j=0; j<renglon; j++) if(!vapoyo[j]) cont++;

		for(j=0; j<ngdln_max; j++)
			if(!vapoyo[renglon+j])
				matens[LocElemVec(cont+j, cont+j)] += elas[i].k[j];
	}
}

//------------------------------------------------------------------------------
// Afecta las matrices elementales que están sujetas a un apoyo inclinado.
int CalEst::ApoyosInclinados(unsigned bar)
{
	// Regresa 0 si no hubo ningún error.
	// Regresa -1 si la barra es de tipo incorrecto (tipo != 1).
	// Regresa -2 si la barra está sostenida por apoyos inclinado - inclidado.

	if(!nincl) return 0;

	register unsigned i;
	double *vec=0, L, ang1, ang2;

	if(!(vec = new double [10])) {
		error.MensError("MEM000");
		sprintf(error.mensaje, "No se pudo crear matriz de barra Art-Art con apoyo \
inclinado: %u", (bar+1));
		return -1;
	}

	for(i=0; i<nincl; i++) {
		if(barra[bar].nodoini == incl[i].nodo) {
			if(barra[bar].tipo != 1) {
				error.MensError("AIN001");
				sprintf(error.mensaje, "Barra: %u", (bar+1));
				return -1;
			}
			else if(barra[bar].nodofin == incl[i].nodo) {
				error.MensError("AIN002");
				sprintf(error.mensaje, "Barra: %u", (bar+1));
				return -2;
			}

			LonAng(bar, L, ang1);
			ang2 = GAR(incl[i].xincl);

			vec[0] = pow(sin(ang2)*sin(ang1)+cos(ang2)*cos(ang1), 2);
			vec[1] = sin(ang2)*cos(ang2)*(pow(sin(ang1),2)-pow(cos(ang1), 2))+
				sin(ang1)*cos(ang1)*(pow(cos(ang2), 2)-pow(sin(ang2), 2));
			vec[2] = pow(sin(ang2)*cos(ang1)-cos(ang2)*sin(ang1), 2);
			vec[3] =-cos(ang2)*pow(cos(ang1), 2)-sin(ang2)*sin(ang1)*cos(ang1);
			vec[4] = sin(ang2)*pow(cos(ang1), 2)-cos(ang2)*sin(ang1)*cos(ang1);
			vec[5] = pow(cos(ang1), 2);
			vec[6] =-sin(ang2)*pow(sin(ang1), 2)-cos(ang2)*sin(ang1)*cos(ang1);
			vec[7] =-cos(ang2)*pow(sin(ang1), 2)+sin(ang2)*sin(ang1)*cos(ang1);
			vec[8] = sin(ang1)*cos(ang1);
			vec[9] = pow(sin(ang1), 2);

			// Borra matriz anterior:
			if(barra[bar].matriz) 
				delete [] barra[bar].matriz;

			// Asigna el puntero de la nueva matriz elemental:
			if(!(barra[bar].matriz = Multiplica(mate[barra[bar].mat].prop[0]*mate[barra[bar].
				mat].prop[1]/L, 10, vec))) {
				error.MensError("MEM000");
				sprintf(error.mensaje, "No se pudo intercambiar matrices de apoyo inclinado.");
				return -1;
			}

			break;
		}
		else if(barra[bar].nodofin == incl[i].nodo) {
			if(barra[bar].tipo != 1) {
				error.MensError("AIN001");
				sprintf(error.mensaje, "Barra: %u", (bar+1));
				return -1;
			}
			else if(barra[bar].nodoini == incl[i].nodo) {
				error.MensError("AIN002");
				sprintf(error.mensaje, "Barra: %u", (bar+1));
				return -2;
			}

			LonAng(bar, L, ang1);
			ang2 = GAR(incl[i].xincl);

			vec[0] = pow(cos(ang1), 2);
			vec[1] = sin(ang1)*cos(ang1);
			vec[2] = pow(sin(ang1), 2);
			vec[3] =-cos(ang2)*pow(cos(ang1), 2)-sin(ang2)*sin(ang1)*cos(ang1);
			vec[4] =-sin(ang2)*pow(sin(ang1), 2)-cos(ang2)*sin(ang1)*cos(ang1);
			vec[5] = pow(sin(ang2)*sin(ang1)+cos(ang2)*cos(ang1), 2);
			vec[6] = sin(ang2)*pow(cos(ang1),2)-cos(ang2)*sin(ang1)*cos(ang1);
			vec[7] =-cos(ang2)*pow(sin(ang1), 2)+sin(ang2)*sin(ang1)*cos(ang1);
			vec[8] = sin(ang2)*cos(ang2)*(pow(sin(ang1), 2)-pow(cos(ang1), 2))-
				sin(ang1)*cos(ang1)*(pow(sin(ang2), 2)-pow(cos(ang2), 2));
			vec[9] = pow(sin(ang2)*cos(ang1)-cos(ang2)*sin(ang1), 2);
						
			// Borra matriz anterior:
			if(barra[bar].matriz) 
				delete [] barra[bar].matriz;

			// Asigna el puntero de la nueva matriz elemental:
			if(!(barra[bar].matriz = Multiplica(mate[barra[bar].mat].prop[0]*mate[barra[bar].
				mat].prop[1]/L, 10, vec))) {
				error.MensError("MEM000");
				sprintf(error.mensaje, "No se pudo intercambiar matrices de apoyo inclinado.");
				return -1;
			}

			break;
		}
	}

	if(vec) delete [] vec;

	return 0;
}

//------------------------------------------------------------------------------
// Asigna memoria para la resolución de la estructura.
int CalEst::AsignaMemoria1(void)
{
	// Regresa 0 si no se encontró ningún error.
	// Regresa -1 si se no se pudo asignar memoria suficiente.
	
	register unsigned i;

	// Asigna memoria para el vector de NGDLN máximos:
	if(!(vngdln_max = new unsigned [nnodos])) return -1;

	// Inicializa el vector de NGDLN máximos:
	for(i=0; i<nnodos; i++) vngdln_max[i] = ndims;

	// Asigna memoria para el vector de restricciones nodales:
	if(!(vapoyo = new int [ngdlt])) return -1;

	// Inicializa el vector de restricciones en nodos:
	for(i=0; i<ngdlt; i++) vapoyo[i] = -1;

	return 0;
}

//------------------------------------------------------------------------------
// Asigna memoria para cada solicitación de carga.
int CalEst::AsignaMemoria2(unsigned caso)
{
	// Regresa 0 si no se encontró ningún error.
	// Regresa -1 si se no se pudo asignar memoria suficiente.

	register unsigned i;

	// Asigna memoria para el vector de cargas:
	if(!(CC[caso].vcarga = new double [ngdlt])) return -1;
		
	// Inicializa el vector de cargas a cero.
	for(i=0; i<ngdlt; i++) CC[caso].vcarga[i] = 0.0;

	// Asigna memoria para el vector de desplazamientos:
	if(!(CC[caso].des = new double [ngdlt])) return -1;
	
	// Inicializa el vector de desplazamientos:
	for(i=0; i<ngdlt; i++) CC[caso].des[i] = 0.0;

	// Asigna memoria para los vectores de fuerzas resultantes en barras:
	if(!(CC[caso].vfbarra = new double* [nbarras])) return -1;
	for(i=0; i<nbarras; i++) CC[caso].vfbarra[i] = 0;

	// Asigna memoria para los vectores de fuerzas resultantes en nodos:
	if(!(CC[caso].vfnodo = new double* [nnodos])) return -1;
	for(i=0; i<nnodos; i++) CC[caso].vfnodo[i] = 0;

	// Asigna memoria para el vector de esfuerzos:
	if(!(CC[caso].esf = new double [nbarras])) return -1;

	return 0;
}

//------------------------------------------------------------------------------
// Calcula el vector de la matriz elemental de una barra Articulada - Articulada.
double *CalEst::BarraArtArt(unsigned bar)
{
	// Regresa 0 si no hubo memoria paea asignar, sino regresa el puntero.

	double *vec=0, L, ang, *vecres=0;

	switch(ndims) {
		case 2:
			// Matriz para la barra Art-Art en 2D:
			if(!(vec = new double [10])) {
				error.MensError("MEM000");
				sprintf(error.mensaje, "No se pudo crear matriz de barra: %u", (bar+1));
				return 0;
			}

			LonAng(bar, L, ang);

			vec[0] = pow(cos(ang), 2);
			vec[1] = sin(ang)*cos(ang);
			vec[2] = pow(sin(ang), 2);
			vec[3] =-pow(cos(ang), 2);
			vec[4] =-sin(ang)*cos(ang);
			vec[5] = pow(cos(ang), 2);
			vec[6] =-sin(ang)*cos(ang);
			vec[7] =-pow(sin(ang), 2);
			vec[8] = sin(ang)*cos(ang);
			vec[9] = pow(sin(ang), 2);
		break;
		case 3:
			// Matriz para la barra Art-Art en 3D:
			double ang1;
			
			if(!(vec = new double [21])) {
				error.MensError("MEM000");
				sprintf(error.mensaje, "No se pudo crear matriz de barra: %u", (bar+1));
				return 0;
			}

			LonAng(bar, L, ang, ang1);

			vec[0]  = vec[9]  = pow(cos(ang), 2)*pow(sin(ang1), 2);
			vec[1]  = vec[13] = cos(ang)*sin(ang)*pow(sin(ang1), 2);
			vec[2]  = vec[14] = pow(sin(ang), 2)*pow(sin(ang1), 2);
			vec[3]  = vec[18] = cos(ang)*sin(ang1)*cos(ang1);
			vec[4]  = vec[19] = sin(ang)*cos(ang1)*sin(ang1);
			vec[5]  = vec[20] = pow(cos(ang1), 2);
			vec[6]  =-vec[0];
			vec[7]  = vec[10] =-vec[1];
			vec[8]  = vec[15] =-vec[3];
			vec[11] =-vec[2];
			vec[12] = vec[16] =-vec[4];
			vec[17] =-vec[5];
		break;
	}

	if(!(vecres = Multiplica(mate[barra[bar].mat].prop[0]*mate[barra[bar].mat].prop[1]/L,
		(ndims==2) ? 10:21, vec))) {
		error.MensError("MEM000");
		sprintf(error.mensaje, "No se pudo asignar matriz de barra articulada.");
		return 0;
	}

	if(vec) delete [] vec;

	return vecres;
}

//------------------------------------------------------------------------------
// Calcula el vector de la matriz elemental de una barra Articulada - Empotrada.
double *CalEst::BarraArtEmp(unsigned bar)
{
	// Regresa 0 si no hubo memoria paea asignar, sino regresa el puntero.

	double *vec=0, L, ang, E, A, I;

	if(!(vec = new double [15])) {
		error.MensError("MEM000");
		sprintf(error.mensaje, "No se pudo crear matriz de barra Art-Emp: %u", (bar+1));
		return 0;
	}

	LonAng(bar, L, ang);

	E = mate[barra[bar].mat].prop[0];
	A = mate[barra[bar].mat].prop[1];
	I = mate[barra[bar].mat].prop[2];

	vec[0]  = vec[5] = E*A/L*pow(cos(ang), 2)+3*E*I/pow(L, 3)*pow(sin(ang), 2);
	vec[1]  = vec[8] = (E*(A*pow(L, 2)-3*I)*sin(2*ang))/2*pow(L, 3);
	vec[2]  = vec[9] = 3*E*I/pow(L, 3)*pow(cos(ang), 2)+E*A/L*pow(sin(ang), 2);
	vec[3]  =-vec[0];
	vec[4]  = vec[6] =-vec[1];
	vec[7]  =-vec[2];
	vec[10] =-3*E*I/pow(L, 2)*sin(ang);
	vec[11] = 3*E*I/pow(L, 2)*cos(ang);
	vec[12] =-vec[10];
	vec[13] =-vec[11];
	vec[14] = 3*E*I/L;
	
	return vec;
}

//------------------------------------------------------------------------------
// Calcula el vector de la matriz elemental de una barra Empotrada - Articulada.
double *CalEst::BarraEmpArt(unsigned bar)
{
	// Regresa 0 si no hubo memoria paea asignar, sino regresa el puntero.

	double *vec=0, L, ang, E, A, I;

	if(!(vec = new double [15])) {
		error.MensError("MEM000");
		sprintf(error.mensaje, "No se pudo crear matriz de barra Emp-Art: %u", (bar+1));
		return 0;
	}

	LonAng(bar, L, ang);

	E = mate[barra[bar].mat].prop[0];
	A = mate[barra[bar].mat].prop[1];
	I = mate[barra[bar].mat].prop[2];

	vec[0]  = vec[9]  = E*A/L*pow(cos(ang), 2)+3*E*I/pow(L, 3)*pow(sin(ang), 2);
	vec[1]  = vec[13] = (E*A/L-3*E*I/pow(L, 3))*sin(ang)*cos(ang);
	vec[2]  = vec[14] = E*A/L*pow(sin(ang), 2)+3*E*I/pow(L, 3)*pow(cos(ang), 2);
	vec[3]  =-3*E*I/pow(L, 2)*sin(ang);
	vec[4]  = 3*E*I/pow(L, 2)*cos(ang);
	vec[5]  = 3*E*I/L;
	vec[6]  =-vec[0];
	vec[7]  = vec[10] =-vec[1];
	vec[8]  =-vec[3];
	vec[11] =-vec[2];
	vec[12] =-vec[4];
	
	return vec;
}

//------------------------------------------------------------------------------
// Calcula el vector de la matriz elemental de una barra Empotrada - Empotrada.
double *CalEst::BarraEmpEmp(unsigned bar)
{
	// Regresa 0 si no hubo memoria paea asignar, sino regresa el puntero.

	double *vec=0, L, ang, E, A, I;

	if(!(vec = new double [21])) {
		error.MensError("MEM001");
		sprintf(error.mensaje, "No se pudo crear matriz de barra Emp-Emp: %u", (bar+1));
		return 0;
	}

	LonAng(bar, L, ang);

	E = mate[barra[bar].mat].prop[0];
	A = mate[barra[bar].mat].prop[1];
	I = mate[barra[bar].mat].prop[2];

	vec[0]  = vec[9]  = E*A/L*pow(cos(ang), 2)+12*E*I/pow(L, 3)*pow(sin(ang), 2);
	vec[1]  = vec[13] = (E*A/L-12*E*I/pow(L, 3))*sin(ang)*cos(ang);
	vec[2]  = vec[14] = E*A/L*pow(sin(ang), 2)+12*E*I/pow(L, 3)*pow(cos(ang), 2);
	vec[3]  = vec[15] =-6*E*I/pow(L, 2)*sin(ang);
	vec[4]  = vec[16] = 6*E*I/pow(L, 2)*cos(ang);
	vec[5]  = vec[20] = 4*E*I/L;
	vec[6]  =-vec[0];
	vec[7]  = vec[10] =-vec[1];
	vec[8]  = vec[18] =-vec[3];
	vec[11] =-vec[2];
	vec[12] = vec[19] =-vec[4];
	vec[17] = vec[5]/2;

	return vec;
}

//------------------------------------------------------------------------------
// Resuelve la estructura.
int CalEst::Calcula(unsigned cal, double toler, bool solsim, unsigned nsolver)
{
	// Regresa 0 ni no hubo ningún error.
	// Regresa -1 si no se pudo asignar memoria suficiente..
	// Regresa -2 si no se pudieron generar las matrices elementales.
	// Regresa -3 si no se pudo ensamblar el vector de cargas.
	// Regresa -4 si no se pudo generar la matriz ensamblada.
	// Regresa -5 si no se pudo resolver el sistema de ecuaciones.
	// Regresa -6 si no se pudiéron obtener las reacciones en barras.
	// Regresa -7 si no se pudiéron obtener las reacciones en nodos.

	register unsigned i;

	//cout.precision(3);

	// Asigna memoria para guardar las matrices de cálculo:
	if(AsignaMemoria1()) {
		error.MensError("MEM003");
		return -1;
	}

	// Cálculo de matrices elementales:
	if((cal & CAL_MATELEM) == CAL_MATELEM)
		for(i=0; i<nbarras; i++)
			if(MatrizElemental(i)) return -2;

	// Cálculo de la matriz ensamblada:
	if((cal & CAL_MATENSA) == CAL_MATENSA) {
		if(VectorRestricciones()) return -6;

		for(i=0; i<ncargas; i++) {
			if(AsignaMemoria2(i)) {
				error.MensError("MEM004");
				sprintf(error.mensaje, "Caso de carga: %u", i+1);
				return -1;
			}

			if(Cargas(i)) return -3;
		}
		
		if(MatrizEnsamblada()) return -4;
	
		// Afecta matriz ensamblada con rigideces de resorte:
		ApoyosElasticos();
	}

	// Resuelve el sistema de ec. y calcula los desplazamientos:
	if((cal & CAL_DESPLAZ) == CAL_DESPLAZ) {
		if(dimdis) {
			// Resuelve el sistema de ecuaciones de forma simultánea:
			if(solsim) {
				if(Resuelve(nsolver, toler)) return -5;
			}
			// Resuelve el sistema de ecuaciones de forma independiente:
			else {
				for(i=0; i<ncargas; i++)
					if(Resuelve(i, nsolver, toler)) return -5;
			}
		}

		// Calcula los desplazamientos:
		for(i=0; i<ncargas; i++)
			Desplazamientos(i);
	}

	// Calcula las reacciones en las barras:
	if((cal & CAL_REACBAR) == CAL_REACBAR)
		for(i=0; i<ncargas; i++)
			if(ReacBarras(i)) return -6;

	// Calcula las reacciones en los nodos:
	if((cal & CAL_REACNOD) == CAL_REACNOD)
		for(i=0; i<ncargas; i++)
			if(ReacNodos(i)) return -7;
	
	// Calcula los esfuerzos en las barras:
	if((cal & CAL_ESFUERZ) == CAL_ESFUERZ)
		for(i=0; i<ncargas; i++)
			Esfuerzos(i);

	return 0;
}

//------------------------------------------------------------------------------
// Construlle el vector de cargas resultante de la acción de fuerzas externas.
int CalEst::Cargas(unsigned caso)
{
	// Regresa 0 ni no hubo ningún error.
	// Regresa -1 si no hubo memoria suficiente.

	register unsigned i, j, renglon;
	double L, ang, *fep=0;

	// Reserva memoria para las FEP:
	if(!(fep = new double [ngdle_max])) {
		error.MensError("MEM000");
		sprintf(error.mensaje, "Al crear vector de FEP. Caso de carga: %u", caso+1);
		return -1;
	}

	for(i=0; i<ngdle_max; i++)
		fep[i] = 0.0;

	// Ensambla en el vector de cargas, las cargas puntuales en nodos:
	for(i=0; i<CC[caso].nfnodos; i++) {
		for(j=0; j<ngdln_max; j++) {
			renglon = (CC[caso].fnodo[i].nodo+1)*ngdln_max-(ngdln_max-j);
			CC[caso].vcarga[renglon] += CC[caso].fnodo[i].f[j];
		}
	}

	// Cargas puntuales dentro de las barras::
	if(ndims == 2) {
		// 2D:
		for(i=0; i<CC[caso].nfbarras; i++) {
			LonAng(CC[caso].fbarra[i].barra, L, ang);

			switch(barra[CC[caso].fbarra[i].barra].tipo) {
				case 1:
					CC[caso].vfep[CC[caso].fbarra[i].barra][0] = CC[caso].
						fbarra[i].f*(L-CC[caso].fbarra[i].d)/L*sin(ang);
					CC[caso].vfep[CC[caso].fbarra[i].barra][1]-= CC[caso].
						fbarra[i].f*(L-CC[caso].fbarra[i].d)/L*cos(ang);
					CC[caso].vfep[CC[caso].fbarra[i].barra][3] = CC[caso].
						fbarra[i].f*(CC[caso].fbarra[i].d)/L*sin(ang);
					CC[caso].vfep[CC[caso].fbarra[i].barra][4]-= CC[caso].
						fbarra[i].f*(CC[caso].fbarra[i].d)/L*cos(ang);
				break;
				case 2:
					// Componentes de FEP en sistema local:	
					fep[1] =-CC[caso].fbarra[i].f*pow((L-CC[caso].fbarra[i].d)
						/L, 2)*(1+2*CC[caso].fbarra[i].d/L);
					fep[2] =-CC[caso].fbarra[i].f*CC[caso].fbarra[i].d*pow((L-
						CC[caso].fbarra[i].d), 2)/pow(L, 2);
					fep[4] =-CC[caso].fbarra[i].f*pow(CC[caso].fbarra[i].d/L, 2)
						*(1+2*(L-CC[caso].fbarra[i].d)/L);
					fep[5] = CC[caso].fbarra[i].f*pow(CC[caso].fbarra[i].d, 2)
						*(L-CC[caso].fbarra[i].d)/pow(L, 2);

					// Trasladando al sistema global:
					CC[caso].vfep[CC[caso].fbarra[i].barra][0] +=-fep[1]*sin(ang);
					CC[caso].vfep[CC[caso].fbarra[i].barra][1] += fep[1]*cos(ang);
					CC[caso].vfep[CC[caso].fbarra[i].barra][2] += fep[2];
					CC[caso].vfep[CC[caso].fbarra[i].barra][3] +=-fep[4]*sin(ang);
					CC[caso].vfep[CC[caso].fbarra[i].barra][4] += fep[4]*cos(ang);
					CC[caso].vfep[CC[caso].fbarra[i].barra][5] += fep[5];				
				break;
				case 3:
					// Componentes de FEP en sistema local:	
					fep[4] =-(CC[caso].fbarra[i].f*pow(CC[caso].fbarra[i].d, 2)/
						(2*pow(L, 3)))*(L-CC[caso].fbarra[i].d+2*L);
					fep[1] =-(fep[4]+CC[caso].fbarra[i].f);
					fep[2] =-CC[caso].fbarra[i].f*(L-CC[caso].fbarra[i].d)*(pow(L,
						2)-pow(L-CC[caso].fbarra[i].d, 2))/(2*pow(L, 2));
				
					// Trasladando al sistema global:
					CC[caso].vfep[CC[caso].fbarra[i].barra][0] +=-fep[1]*sin(ang);
					CC[caso].vfep[CC[caso].fbarra[i].barra][1] += fep[1]*cos(ang);
					CC[caso].vfep[CC[caso].fbarra[i].barra][2] += fep[2];
					CC[caso].vfep[CC[caso].fbarra[i].barra][3] +=-fep[4]*sin(ang);
					CC[caso].vfep[CC[caso].fbarra[i].barra][4] += fep[4]*cos(ang);
				break;
				case 4:
					// Componentes de FEP en sistema local:	
					fep[1] =-(CC[caso].fbarra[i].f*pow(CC[caso].fbarra[i].d, 2)/
						(2*pow(L, 3)))*(L-CC[caso].fbarra[i].d+2*L);
					fep[4] =-(fep[4]+CC[caso].fbarra[i].f);
					fep[5] =-CC[caso].fbarra[i].f*(L-CC[caso].fbarra[i].d)*(pow(L,
						2)-pow(L-CC[caso].fbarra[i].d, 2))/(2*pow(L, 2));

					// Trasladando al sistema global:
					CC[caso].vfep[CC[caso].fbarra[i].barra][0] +=-fep[1]*sin(ang);
					CC[caso].vfep[CC[caso].fbarra[i].barra][1] += fep[1]*cos(ang);
					CC[caso].vfep[CC[caso].fbarra[i].barra][3] +=-fep[4]*sin(ang);
					CC[caso].vfep[CC[caso].fbarra[i].barra][4] += fep[4]*cos(ang);
					CC[caso].vfep[CC[caso].fbarra[i].barra][5] += fep[5];	
				break;
			}
		}
	
		// Cargas distribuidas rectangulares:
		for(i=0; i<CC[caso].nfrect; i++) {
			LonAng(CC[caso].frect[i].barra, L, ang);

			switch(barra[CC[caso].frect[i].barra].tipo) {
				case 1:
					double reac;
				
					reac = CC[caso].frect[i].f*L/2;
	
					CC[caso].vfep[CC[caso].frect[i].barra][0] = reac*sin(ang);
					CC[caso].vfep[CC[caso].frect[i].barra][1]-= reac*cos(ang);
					CC[caso].vfep[CC[caso].frect[i].barra][3] = reac*sin(ang);
					CC[caso].vfep[CC[caso].frect[i].barra][4]-= reac*cos(ang);
				break;
				case 2:
					// Componentes de FEP en sistema local:
					fep[1] =-CC[caso].frect[i].f*L/2;
			
					// Trasladando al sistema global:
					CC[caso].vfep[CC[caso].frect[i].barra][0] +=-fep[1]*sin(ang);
					CC[caso].vfep[CC[caso].frect[i].barra][1] += fep[1]*cos(ang);
					CC[caso].vfep[CC[caso].frect[i].barra][2] +=-CC[caso].frect[i].
						f*pow(L, 2)/12;
					CC[caso].vfep[CC[caso].frect[i].barra][3] +=-fep[1]*sin(ang);
					CC[caso].vfep[CC[caso].frect[i].barra][4] += fep[1]*cos(ang);
					CC[caso].vfep[CC[caso].frect[i].barra][5] += CC[caso].frect[i].
						f*pow(L, 2)/12;
				break;
				case 3:
					// Componentes de FEP en sistema local:
					fep[1] =-(5*CC[caso].frect[i].f*L)/8;
					fep[4] =-(3*CC[caso].frect[i].f*L)/8;
		
					// Trasladando al sistema global:
					CC[caso].vfep[CC[caso].frect[i].barra][0] +=-fep[1]*sin(ang);
					CC[caso].vfep[CC[caso].frect[i].barra][1] += fep[1]*cos(ang);
					CC[caso].vfep[CC[caso].frect[i].barra][2] +=-CC[caso].frect[i].
						f*pow(L, 2)/8;
					CC[caso].vfep[CC[caso].frect[i].barra][3] +=-fep[4]*sin(ang);
					CC[caso].vfep[CC[caso].frect[i].barra][4] += fep[4]*cos(ang);
				break;
				case 4:
					// Componentes de FEP en sistema local:
					fep[1] =-(3*CC[caso].frect[i].f*L)/8;
					fep[4] =-(5*CC[caso].frect[i].f*L)/8;
		
					// Trasladando al sistema global:
					CC[caso].vfep[CC[caso].frect[i].barra][0] +=-fep[1]*sin(ang);
					CC[caso].vfep[CC[caso].frect[i].barra][1] += fep[1]*cos(ang);
					CC[caso].vfep[CC[caso].frect[i].barra][3] +=-fep[4]*sin(ang);
					CC[caso].vfep[CC[caso].frect[i].barra][4] += fep[4]*cos(ang);
					CC[caso].vfep[CC[caso].frect[i].barra][5] +=-CC[caso].frect[i].
						f*pow(L, 2)/8;
				break;
			}
		}

		// Cargas por peso propio:
		if(CC[caso].PesoPropio) {
			double w;

			for(i=0; i<nbarras; i++) {
				LonAng(i, L, ang);
			
				// Obtiene el peso distribuido de la barra:
				w = mate[barra[i].mat].prop[1]*mate[barra[i].mat].prop[3];

				switch(barra[i].tipo) {
					case 1:
						// Cargas puntuales debidas a PP:
						CC[caso].vfep[i][1] += w*L/2;
						CC[caso].vfep[i][4] += w*L/2;
					break;
					case 2:
						// Componentes globales del peso:
						CC[caso].vfep[i][1] += w*L/2;
						CC[caso].vfep[i][2] += w*pow(L, 2)*cos(ang)/12;
						CC[caso].vfep[i][4] += w*L/2;
						CC[caso].vfep[i][5] +=-w*pow(L, 2)*cos(ang)/12;
					break;
					case 3:
						// Componentes globales del peso:
						CC[caso].vfep[i][1] += 5*w*L/8;
						CC[caso].vfep[i][2] += w*pow(L, 2)*cos(ang)/8;
						CC[caso].vfep[i][4] += 3*w*L/8;
					break;
					case 4:
						// Componentes globales del peso:
						CC[caso].vfep[i][1] += 3*w*L/8;
						CC[caso].vfep[i][4] += 5*w*L/8;
						CC[caso].vfep[i][5] += w*pow(L, 2)*cos(ang)/8;
					break;
				}
			}
		}
	}
	else {
		// 3D:
		// Cargas por peso propio:
		if(CC[caso].PesoPropio) {
			double w, ang1;

			for(i=0; i<nbarras; i++) {
				LonAng(i, L, ang, ang1);
			
				// Obtiene el peso distribuido de la barra:
				w = mate[barra[i].mat].prop[1]*mate[barra[i].mat].prop[3];
	
				// Cargas puntuales debidas a PP:
				CC[caso].vfep[i][2] += w*L/2;
				CC[caso].vfep[i][5] += w*L/2;
			}
		}
	}

	// Ensamblando FEP en el vector de cargas:
	for(i=0; i<nbarras; i++) {
		for(j=0; j<ngdln_max; j++) {
			renglon = (barra[i].nodoini+1)*ngdln_max-(ngdln_max-j);
			CC[caso].vcarga[renglon] -= CC[caso].vfep[i][j];
			renglon = (barra[i].nodofin+1)*ngdln_max-(ngdln_max-j);
			CC[caso].vcarga[renglon] -= CC[caso].vfep[i][ngdln_max+j];
		}
	}

	// Desensambla apoyos y grados de libertad inexistentes:
	unsigned nelems, nelem, cont=0;
	double *vecres=0;

	nelems = ngdlt;
		
	for(i=0; i<ngdlt; i++) {
		if(vapoyo[i]) {
			nelem = i-cont;

			if(!(vecres = EliminaElem(CC[caso].vcarga, nelems, nelem))) {
				error.MensError("MEM000");
				sprintf(error.mensaje, "Al desensamblar vector de cargas. Caso de \
carga: %u", caso+1);
				return -1;
			}
			
			if(CC[caso].vcarga) delete [] CC[caso].vcarga;
			CC[caso].vcarga = vecres;
			cont++;
		}
	}
	
	if(fep) delete [] fep;

	return 0;
}

//------------------------------------------------------------------------------
// Completa el vector de desplazamientos con los conocidos y los resultantes.
void CalEst::Desplazamientos(unsigned caso)
{
	register unsigned i, j, nodo, drest, cont=0;

	for(i=0; i<ngdlt; i++) {
		switch(vapoyo[i]) {
			case 0:
				CC[caso].des[i] = CC[caso].desres[cont];
				cont++;
			break;
			case 1:
				nodo = (unsigned)(i/ngdln_max);
				
				for(j=0; j<napoyos; j++) {
					if(apoyo[j].nodo == nodo) {
						drest = i-nodo*ngdln_max;
						CC[caso].des[i] = apoyo[j].d[drest];
						j = napoyos;
					}
				}
			break;
		}
	}
}

//------------------------------------------------------------------------------
// Destruye el objeto.
void CalEst::Destruye(void)
{
	Res_Calculo :: Destruye();

	// Libera memoria:
	if(mate) {
   		delete [] mate;
		mate = 0;
		nmats = 0;
	}

	if(nodo) {
   		delete [] nodo;
		nodo = 0;
		nnodos = 0;
	}

	if(barra) {
   		delete [] barra;
		barra = 0;
		nbarras = 0;
	}
	
	if(apoyo) {
   		delete [] apoyo;
		apoyo = 0;
		napoyos = 0;
	}

	if(elas) {
   		delete [] elas;
		elas = 0;
		nelas = 0;
	}

	if(incl) {
   		delete [] incl;
		incl = 0;
		nincl = 0;
	}
	
	if(CC) {
   		delete [] CC;
		CC = 0;
		ncargas = 0;
	}

	nprop = 0;
	ndims = 0;
}

//------------------------------------------------------------------------------
// Escribe los resultados en un archivo con formato ASCII.
int CalEst::EscribeRES(char *nomarch)
{
	// Regresa 0 si no hubo ningún error.
	// Regresa -1 si no se puede crear el archivo.

	register unsigned i, j, k;

	// Inicia la escritura de resultados:
	ofstream arch(nomarch);

	if(!arch) return -1;

    arch << "\t\t\t\t\t\tSICBASA Motor de Cálculo v. 1.0\n\n";
	arch << "Este reporte del CALEST 1.0 es para uso interno de SICBASA cualquier interpretación\n\
incorrecta del mismo es responsabilidad del usuario.\n\n";

	// Escribe el título:
	arch << "Título del problema; " << titulo << "\n\n";
	
	// Escribe las unidades utilizadas de longitud y fuerza:
	arch << "UNIDADES:\n";
	arch << "Longitud: " << uLongitud << "\n";
	arch << "Fuerza:   " << uFuerza << "\n\n";
	
	// Escribe los parámetros generales del problema:
	arch << "CARACTERÍSTICA GENERALES DEL PROBLEMA:\n";
	arch << "No. de nodos:      " << nnodos << "\n";
	arch << "No. de barras:     " << nbarras << "\n";
	arch << "No. de apoyos:     " << napoyos << "\n";
	arch << "No. de C de Carga: " << ncargas << "\n";
	arch << "No. de GDL max:    " << ngdln_max << "\n";
	arch << "No. de mats:       " << nmats << "\n";
	arch << "No. de dims.:      " << ndims << "\n";
	arch << "No. de ap. elás.:  " << nelas << "\n";
	arch << "No. de ap. incl.:  " << nincl << "\n";
	arch << "Ind. esc. de mat.: " << ((matrices) ? "Encendido" : "Apagado") << "\n\n";
	
	// Escribe las coordenadas de los nodos:
	arch << "\nCARACTERÍSTICAS DE LOS NODOS:\n";
	arch << "No.\tX\tY" << ((ndims == 3) ? "\tZ\n" : "\n");
	
	for(i=0; i<nnodos; i++) {
		arch << (i+1) << "\t";
		for(j=0; j<ndims; j++) arch << nodo[i].coord[j] << "\t";
		arch << "\n";
	}

	// Escribe las características de las barras:
	arch << "CARACTERÍSTICAS DE LAS BARRAS:\n";
	arch << "No.\tNo. mat.\tNodo ini.\tNodo fin.\tTipo\n";

	for(i=0; i<nbarras; i++) {
		arch << (i+1) << "\t";
		arch << (barra[i].mat+1) << "\t\t\t";
		arch << (barra[i].nodoini+1) << "\t\t\t";
		arch << (barra[i].nodofin+1) << "\t\t\t";
		arch << barra[i].tipo << "\n";
	}

	// Escribe las restricciones:
	arch << "\nCARACTERÍSTICAS DE LOS APOYOS RÍGIDOS:\n";
	arch << "No.\tNodo\tDr. rest.\tDesp.\n";

	for(i=0; i<napoyos; i++) {
		arch << (i+1) << "\t";
		arch << (apoyo[i].nodo+1) << "\t";
		for(j=0; j<ngdln_max; j++) arch << apoyo[i].drest[j] << "\t";
		for(j=0; j<ngdln_max; j++) arch << apoyo[i].d[j] << "\t";
		arch << "\n";
	}

	// Escribe las propiedades de los materiales:
	arch << "\nCARACTERÍSTICAS DE LOS MATERIALES:\n";
	arch << "No.\tE\tA\tI\tPeso esp\n";

	for(i=0; i<nmats; i++) {
		arch << (i+1) << "\t";
		for(j=0; j<nprop; j++) arch << mate[i].prop[j] << "\t";
		arch  << "\n";
	}

	// Escribe los apoyos elásticos:
	arch << "\nCARACTERÍSTICAS DE LOS APOYOS ELÁSTICOS:\n";
	arch << "No.\tNodo\tX\tY\tGiro\n";

	for(i=0; i<nelas; i++) {
		arch << (i+1) << "\t";
		arch << (elas[i].nodo+1) << "\t";
		for(j=0; j<ngdln_max; j++) arch << elas[i].k[j] << "\t";
		arch << "\n";
	}

	// Escribe los apoyos inclinados:
	arch << "\nCARACTERÍSTICAS DE LOS APOYOS INCLINADOS:\n";
	arch << "No.\tNodo\tÁngulo\n";

	for(i=0; i<nincl; i++) {
		arch << (i+1) << "\t";
		arch << (incl[i].nodo+1) << "\t";
		arch << incl[i].xincl << "\n";
	}

	// Escribe los casos de carga:
	for(i=0; i<ncargas; i++) {
		// Escribe la linea de cargas:
		arch << "\nCASO DE CARGA: " << (i+1) << "\n";
		arch << "Título del caso: " << CC[i].titulo << "\n\n";

		arch << "CARACTERÍSTICAS GENERALES DEL CASO DE CARGA:\n";
		arch << "No. de cargas punt. en nodos:   " << CC[i].nfnodos << "\n";
		arch << "No. de cargas punt. en barras:  " << CC[i].nfbarras << "\n";
		arch << "No. de cargas unif. repartidas: " << CC[i].nfrect << "\n";
		arch << "Peso propio:                    " << CC[i].PesoPropio << "\n\n";

		// Escribe las fuerzas en los nodos:
		arch << "FUERZAS EN PUNTUALES EN NODOS:\n";
		arch << "No.\tNodo\tX\tY\tM\n";

		for(j=0; j<CC[i].nfnodos; j++) {
			arch << (j+1) << "\t";
			arch << (CC[i].fnodo[j].nodo+1) << "\t\t";
			for(k=0; k<ngdln_max; k++) arch << CC[i].fnodo[j].f[k] << "\t";
			arch << "\n";
		}

		// Escribe las fuerzas en las barras:
		arch << "\nFUERZAS EN PUNTUALES EN BARRAS:\n";
		arch << "No.\tBarra\tCarga\tDistancia\n";

		for(j=0; j<CC[i].nfbarras; j++) {
			arch << (j+1) << "\t";
			arch << (CC[i].fbarra[j].barra+1) << "\t\t";
			arch << CC[i].fbarra[j].f << "\t\t";
			arch << CC[i].fbarra[j].d << "\n";
		}

		// Escribe las fuerzas rectangulares:
		arch << "\nFUERZAS EN RECTANGULARES EN BARRAS:\n";
		arch << "No.\tBarra\tCarga\n";

		for(j=0; j<CC[i].nfrect; j++) {
			arch << (j+1) << "\t";
			arch << (CC[i].frect[j].barra+1) << "\t";
			arch << CC[i].frect[j].f << "\n";
		}
	}

	// RESULTADOS:

	arch << "\n\n******************** RESULTADOS DEL CÁLCULO ********************\n\n";

	// Escribe las matrices (si se desea):
	if(matrices) {
		// Escribe matrices elementales:
		arch << "MATRICES ELEMENTALES DE RIGIDEZ:\n";

		for(i=0; i<nbarras; i++) {
			arch << "\nMatriz de barra: " << (i+1) << "\n";

			for(j=0; j<barra[i].ngdle; j++) {
				for(k=0; k<barra[i].ngdle; k++)
					arch << barra[i].matriz[LocElemVec(j, k)] << "\t";
				arch << "\n";
			}

			// Escribe vector de FEP:
			arch << "\nVector de FEP:\n";

			for(j=0; j<ncargas; j++) {
				for(k=0; k<barra[i].ngdle; k++) arch << CC[j].vfep[i][k] 
					<< "\t";
				arch << "\n";
			}
		}

		// Escribe matriz ensamblada disminuida:
		arch << "\n\nMATRIZ DE RIGIDEZ GLOBAL DISMINUIDA:\n";

		for(i=0; i<dimdis; i++) {
			for(j=0; j<dimdis; j++) arch << matens[LocElemVec(i, j)] << "\t";
		
			arch << "\n";
		}
	}

	// Resultados por cada caso de carga:
	for(i=0; i<ncargas; i++) {
		arch << "\n\nRESULTADOS PARA EL CASO DE CARGA: " << (i+1) << "\n";

		// Escribe el vector de cargas ensamblado por cada caso:
		arch << "\nVector de cargas ensamblado por cada caso:\n";

		for(j=0; j<dimdis; j++)
			arch << CC[i].vcarga[j] << "\t";

		arch << "\n";
		
		// Escribe el vector de desplazamientos resultantes para cada caso:
		arch << "\nVector de desplazamientos resultantes:\n";

		for(j=0; j<dimdis; j++) 
			arch << CC[i].desres[j] << "\t";

		arch << "\n";

		// Escribe el vector de desplazamientos de nodos (para cada caso):
		arch << "\nDesplazamientos nodales por caso:\n";

		for(j=0; j<nnodos; j++) {
			arch << "Nodo " << (j+1) << ":\t";

			for(k=0; k<ngdln_max; k++) 
				arch << CC[i].des[(j+1)*ngdln_max-(ngdln_max-k)] << "\t";
			
			arch << "\n";
		}

		arch << "\n";

		// Escribe las fuerzas en nodos de las barras (para cada caso):
		arch << "\nFuerzas en nodos de barras, por caso:\n";

		for(j=0; j<nbarras; j++) {
			arch << "Barra " << (j+1) << ":\t";

			for(k=0; k<barra[j].ngdle; k++)
				arch << CC[i].vfbarra[j][k] << "\t";

			arch << "\n";
		}

		arch << "\n";

		// Escribe las reacciones (para cada caso):
		arch << "\nReacciones, por caso:\n";

		for(j=0; j<napoyos; j++) {
			arch << "Apoyo " << (j+1) << ":\t";

			for(k=0; k<vngdln_max[apoyo[j].nodo]; k++) {
				arch << -CC[i].vfnodo[apoyo[j].nodo][k] << "\t";
			}
		
			arch << "\n";
		}

		arch << "\n";

		// Escribe los esfuerzos (para cada caso):
		arch << "\nEsfuerzos, por caso:\n";
		arch << "Mínimo: " << CC[i].esf_min << "\n";
		arch << "Máximo: " << CC[i].esf_max << "\n\n";

		for(j=0; j<nbarras; j++) {
			arch << "Barra " << (j+1) << ":\t";
			arch << CC[i].esf[j] << "\n";
		}

		arch << "\n";
	}

	// Cierra el archivo:
	arch.close();

	return 0;
}

//------------------------------------------------------------------------------
// Asigna memoria para cargar la geometría y propiedades de la estructura.
int CalEst::Inicializa(void)
{
	// Regresa 0 si no se encontró ningún error.
	// Regresa -1 si se no se pudo asignar memoria suficiente.

	register unsigned i;

	// Asignación de memoria para los arreglos de barras:
	if(!(nodo  = new Dat_Nodos [nnodos])) return -1;
	if(!(barra = new Dat_Barras [nbarras])) return -1;
	if(!(apoyo = new Dat_Apoyos [napoyos])) return -1;
	if(!(mate  = new Dat_Materiales [nmats])) return -1;
	if(!(elas  = new Dat_Apoyos_Elasticos [nelas])) return -1;
	if(!(incl  = new Dat_Apoyos_Inclinados [nincl])) return -1;
	if(!(CC    = new Dat_CasosCarga [ncargas])) return -1;

	// Asigna memoria para las coordenadas de los nodos:
	for(i=0; i<nnodos; i++)
		if(nodo[i].Inicializa(ndims)) return -1;

	// Asigna memoria para las propiedades de los materiales:
	for(i=0; i<nmats; i++)
		if(mate[i].Inicializa(nprop)) return -1;

	// Asigna memoria para los apoyos:
	for(i=0; i<napoyos; i++)
		if(apoyo[i].Inicializa(ngdln_max)) return -1;

	// Asigna memoria para los apoyos elásticos:
	for(i=0; i<nelas; i++)
		if(elas[i].Inicializa(ngdln_max)) return -1;

	return 0;
}

//------------------------------------------------------------------------------
// Lee los datos de un archivo ASCII en formato MECA.
int CalEst::LeeDAT(char *nomarch)
{
	// Regresa 0 si no hubo ningún error.
	// Regresa -1 si no existe el archivo.
	// Regresa -2 si no es archivo de datos.
	// Regresa -3 si no hubo memoria suficiente para cargar datos.
	// Regresa -4 si la Versión no coincide.
	// Regresa -5 si en un tipo de problema que no puede resolver CalEst.

	register unsigned i, j, k;
	unsigned temp, basura, ngdln, tipo;
	char mensaje[256];

	// Inicia la lectura de datos:
	ifstream arch(nomarch);

	if(!arch) {
		error.MensError("LEC001");
		return -1;
	}

     // Verifica si es un archivo de datos del MECA:
	arch >> mensaje;

	if(strcmp(mensaje, "PreMECA")) {
		error.MensError("LEC002");
		arch.close();
		return -2;
	}

	// Verifica si es una versión válida:
	arch >> mensaje;

	if(strcmp(mensaje, "1.0")) {
		error.MensError("LEC003");
		sprintf(error.mensaje, "Archivo de versión %s", mensaje);
		arch.close();
		return -4;
	}

	// Lee el título:
	arch >> titulo;
	Remplaza(titulo, '_', ' ');

	// Lee las unidades utilizadas de longitud y fuerza:
	arch >> uLongitud;
	arch >> uFuerza;
	Remplaza(uLongitud, '_', ' ');
	Remplaza(uFuerza, '_', ' ');
	
	// Lee la línea de parámetros:
	arch >> nnodos;
	arch >> nbarras;
	arch >> napoyos;
	arch >> ncargas;
	arch >> tipo;
	arch >> basura;
	arch >> nmats;
	arch >> basura;
	arch >> temp;
	matrices = (temp) ? true : false;
	arch >> ndims;
	arch >> nelas;
	arch >> nincl;
	arch >> basura;

	nprop = (ndims == 2 || tipo == 1) ? 4 : 7;

	// Verifica si se puede resolver el problema.
	if(ndims == 3 && tipo != 1) {
		error.MensError("LEC004");
		arch.close();
		return -5;
	}

	// Grados de libertad máximos por nodo (según dimensión):
	if(ndims == 2 && tipo == 1) ngdln = 2;
	else if((ndims == 2 && tipo != 1) || (ndims == 3 && tipo ==1)) ngdln = 3;
	else ngdln = 6;
	ngdln_max = (ndims == 2 || tipo == 1) ? 3 : 6;
	ngdle_max = ngdln_max*2;
	ngdlt = ngdln_max*nnodos;

   	// Asigna memoria para los datos del archivo:
	if(Inicializa()) {
		arch.close();
		error.MensError("MEM001");
		return -3;
	}

	// Lee las características de las barras:
	for(i=0; i<nbarras; i++) {
		arch >> barra[i].mat;
		barra[i].mat--;
		arch >> barra[i].nodoini;
		barra[i].nodoini--;
		arch >> barra[i].nodofin;
		barra[i].nodofin--;

		if(tipo == 3) arch >> barra[i].tipo;
		else barra[i].tipo = tipo;
		
		switch(barra[i].tipo) {
			case 1:	// Articulado - Articulado
				barra[i].ngdlni = barra[i].ngdlnf = (ndims == 2) ? 2 : 3;
			break;
			case 2:	// Empotrado - Empotrado
				barra[i].ngdlni = barra[i].ngdlnf = (ndims == 2) ? 3 : 6;
			break;
			case 3:	// Empotrado - Articulado
				barra[i].ngdlni = (ndims == 2) ? 3 : 6;
				barra[i].ngdlnf = (ndims == 2) ? 2 : 3;
			break;
			case 4:	// Articulado - Empotrado
				barra[i].ngdlni = (ndims == 2) ? 2 : 3;
				barra[i].ngdlnf = (ndims == 2) ? 3 : 6;
			break;
			case 5: // Empotrado - Libre
				barra[i].ngdlni = (ndims == 2) ? 3 : 6;
				barra[i].ngdlnf = 0;
			break;
			case 6: // Libre - Empotrado
				barra[i].ngdlni = 0;
				barra[i].ngdlnf = (ndims == 2) ? 3 : 6;
			break;
		}				

		barra[i].ngdle = barra[i].ngdlni+barra[i].ngdlnf;
	}

	// Lee las coordenadas de los nodos:
	for(i=0; i<nnodos; i++)
		for(j=0; j<ndims; j++) arch >> nodo[i].coord[j];

	// Lee las restricciones (apoyos):
	for(i=0; i<napoyos; i++) {
		arch >> apoyo[i].nodo;
		apoyo[i].nodo--;

		for(j=0; j<ngdln; j++) {
			arch >> temp;
			apoyo[i].drest[j] = (temp) ? true : false;
		}

		for(j=0; j<ngdln; j++) arch >> apoyo[i].d[j];

		if(ndims == 2 && tipo == 1) {
			apoyo[i].drest[2] = false;
			apoyo[i].d[2] = 0.0;
		}
	}

	// Lee las propiedades de los materiales:
	for(i=0; i<nmats; i++)
		for(j=0; j<nprop; j++) arch >> mate[i].prop[j];

	// Lee los apoyos elásticos:
	for(i=0; i<nelas; i++) {
		arch >> elas[i].nodo;
		elas[i].nodo--;

		for(j=0; j<ngdln; j++) arch >> elas[i].k[j];
		if(ndims == 2 && tipo == 1) elas[i].k[2] = 0.0;
	}

	// Lee los apoyos inclinados:
	for(i=0; i<nincl; i++) {
		arch >> incl[i].nodo;
		incl[i].nodo--;
		arch >> incl[i].xincl;
	}

	// Lee los casos de carga:
	for(i=0; i<ncargas; i++) {
		// Lee la linea de cargas:
		arch >> CC[i].titulo;
		Remplaza(CC[i].titulo, '_', ' ');
		arch >> CC[i].nfnodos;
		arch >> CC[i].nfbarras;
		arch >> CC[i].nfrect;
		arch >> basura;
		arch >> basura;
		arch >> temp;
		CC[i].PesoPropio = (temp) ? true : false;
		
		CC[i].nnodos = nnodos;
		CC[i].nbarras = nbarras;

		// Asigna memoria para los Casos de Carga:
		if(CC[i].Inicializa(ngdln_max, ngdle_max)) {
			error.MensError("MEM002");
			return -3;
		}

		// Lee las fuerzas en los nodos:
		for(j=0; j<CC[i].nfnodos; j++) {
			arch >> CC[i].fnodo[j].nodo;
			CC[i].fnodo[j].nodo--;

			for(k=0; k<ngdln; k++) arch >> CC[i].fnodo[j].f[k];
			if(ndims == 2 && tipo == 1) CC[i].fnodo[j].f[2] = 0.0;
		}

		// Lee las fuerzas en las barras:
		for(j=0; j<CC[i].nfbarras; j++) {
			arch >> CC[i].fbarra[j].barra;
			CC[i].fbarra[j].barra--;
			arch >> CC[i].fbarra[j].f;
			arch >> CC[i].fbarra[j].d;
		}

		// Lee las fuerzas rectangulares:
		for(j=0; j<CC[i].nfrect; j++) {
			arch >> CC[i].frect[j].barra;
			CC[i].frect[j].barra--;
			arch >> CC[i].frect[j].f;
		}
	}

	// Cierra el archivo:
	arch.close();

	return 0;
}

//------------------------------------------------------------------------------
// Obtiene la longitud y el ángulo de una barra en 2D.
void CalEst::LonAng(unsigned bar, double &longitud, double &angulo)
{
	double x1, y1, x2, y2;

	x1 = nodo[barra[bar].nodoini].coord[0];
	y1 = nodo[barra[bar].nodoini].coord[1];
	x2 = nodo[barra[bar].nodofin].coord[0];
	y2 = nodo[barra[bar].nodofin].coord[1];

	// Obtiene la longitud de la barra:
	longitud = sqrt(pow((x2-x1), 2)+pow((y2-y1), 2));

	// Obtiene el ángulo de la barra:
	angulo = atan2(y2-y1, x2-x1);
}

//------------------------------------------------------------------------------
// Obtiene la longitud y los ángulos de una barra en 3D.
void CalEst::LonAng(unsigned bar, double &longitud, double &angulo, double &angulo1)
{
	double x1, y1, z1, x2, y2, z2, Lxy;

	x1 = nodo[barra[bar].nodoini].coord[0];
	y1 = nodo[barra[bar].nodoini].coord[1];
	z1 = nodo[barra[bar].nodoini].coord[2];
	x2 = nodo[barra[bar].nodofin].coord[0];
	y2 = nodo[barra[bar].nodofin].coord[1];
	z2 = nodo[barra[bar].nodofin].coord[2];

	// Obtiene la longitud de la barra:
	longitud = sqrt(pow(x2-x1, 2)+pow(y2-y1, 2)+pow(z2-z1, 2));
	Lxy = sqrt(pow(x2-x1, 2)+pow(y2-y1, 2));

	// Obtiene los ángulos de la barra:
	angulo = atan2(y2-y1, x2-x1);
	angulo1 = atan2(Lxy, z2-z1);
}

//------------------------------------------------------------------------------
// Calcula los esfuerzos en las barras.
void CalEst::Esfuerzos(unsigned caso)
{
	register unsigned i;
	double L, ang1;
	
	// Cálculo de esfuerzos:
	if(ndims == 2) {
		for(i=0; i<nbarras; i++) {
			LonAng(i, L, ang1);

			CC[caso].esf[i] = -(
				CC[caso].vfbarra[i][0]*cos(ang1)+
				CC[caso].vfbarra[i][1]*sin(ang1)
				)/mate[barra[i].mat].prop[1];
		}
	}
	else {
		double ang2;

		for(i=0; i<nbarras; i++) {
			LonAng(i, L, ang1, ang2);

			CC[caso].esf[i] = -(
				CC[caso].vfbarra[i][0]*cos(ang1)*cos(ang2)+
				CC[caso].vfbarra[i][1]*sin(ang1)*cos(ang2)+
				CC[caso].vfbarra[i][2]*sin(ang2)
				)/mate[barra[i].mat].prop[1];
		}
	}

	// Encuentra el esfuerzo máximo y mínimo:
	CC[caso].esf_min =  1e10;
	CC[caso].esf_max = -1e10;

	for(i=0; i<nbarras; i++) {
		CC[caso].esf_min = MIN(CC[caso].esf_min, CC[caso].esf[i]);
		CC[caso].esf_max = MAX(CC[caso].esf_max, CC[caso].esf[i]);
	}
}

//------------------------------------------------------------------------------
// Selecciona y calcula el tipo de matriz elemental dependiendo del tipo de barra.
int CalEst::MatrizElemental(unsigned bar)
{
	// Regresa 0 si no hubo ningún error.
	// Regresa -1 si hubo algún error.
	
	// Cálculo de las matrices elementales según el tipo de sujeción de la barra:
	switch(barra[bar].tipo) {
		case 1:
			if(!(barra[bar].matriz = BarraArtArt(bar))) return -1;
		break;
		case 2:
			if(!(barra[bar].matriz = BarraEmpEmp(bar))) return -1;
		break;
		case 3:
			if(!(barra[bar].matriz = BarraEmpArt(bar))) return -1;
		break;
		case 4:
			if(!(barra[bar].matriz = BarraArtEmp(bar))) return -1;
		break;
	}

	if(ApoyosInclinados(bar)) return -1;
	
	return 0;
}

//------------------------------------------------------------------------------
// Ensambla la matriz global en un vector.
int CalEst::MatrizEnsamblada(void)
{
	// Regresa 0 si no hubo ningún error.
	// Regresa -1 si no hubo memoria para asignar.

	register unsigned i, j, k, l, r, c, ren, col;

	// Reserva memoria para ensamblaje:
	c = DimVec(dimdis);

	if(!(matens = new double [c])) {
		error.MensError("MEM000", "Al generar matriz ensamblada.");
		return -1;
	}

	// Inicializa la matriz ensamblada a cero:
	for(i=0; i<c; i++) matens[i] = 0.0;

	// Inserta elementos de las matrices elementales a la matriz global:
	for(i=0; i<nbarras; i++) {
		// Cuadrante I de la matriz elemental:
		for(j=0; j<barra[i].ngdlni; j++) {
			for(k=0; k<barra[i].ngdlni; k++) {
				r = (barra[i].nodoini+1)*ngdln_max-(ngdln_max-j);
				c = (barra[i].nodoini+1)*ngdln_max-(ngdln_max-k);
							
				if(!vapoyo[r] && !vapoyo[c] && c<=r) {
					ren = col = 0;
					for(l=0; l<=r; l++) if(!vapoyo[l]) ren++;
					for(l=0; l<=c; l++) if(!vapoyo[l]) col++;

					matens[LocElemVec(ren-1, col-1)] += barra[i].matriz
						[LocElemVec(j, k)];
				}
			}
		}

		// Cuadrante II de la matriz elemental:
		for(j=0; j<barra[i].ngdlni; j++) {
			for(k=0; k<barra[i].ngdlnf; k++) {
				r = (barra[i].nodoini+1)*ngdln_max-(ngdln_max-j);
				c = (barra[i].nodofin+1)*ngdln_max-(ngdln_max-k);

				if(!vapoyo[r] && !vapoyo[c] && c<=r) {
					ren = col = 0;
					for(l=0; l<=r; l++) if(!vapoyo[l]) ren++;
					for(l=0; l<=c; l++) if(!vapoyo[l]) col++;

					matens[LocElemVec(ren-1, col-1)] += barra[i].matriz
						[LocElemVec(j, barra[i].ngdlni+k)];
				}
			}
		}

		// Cuadrante III de la matriz elemental:
		for(j=0; j<barra[i].ngdlnf; j++) {
			for(k=0; k<barra[i].ngdlni; k++) {
				r = (barra[i].nodofin+1)*ngdln_max-(ngdln_max-j);
				c = (barra[i].nodoini+1)*ngdln_max-(ngdln_max-k);

				if(!vapoyo[r] && !vapoyo[c] && c<=r) {
					ren = col = 0;
					for(l=0; l<=r; l++) if(!vapoyo[l]) ren++;
					for(l=0; l<=c; l++) if(!vapoyo[l]) col++;

					matens[LocElemVec(ren-1, col-1)] += barra[i].matriz
						[LocElemVec(barra[i].ngdlni+j, k)];
				}
			}
		}

		// Cuadrante IV de la matriz elemental:
		for(j=0; j<barra[i].ngdlnf; j++) {
			for(k=0; k<barra[i].ngdlnf; k++) {
				r = (barra[i].nodofin+1)*ngdln_max-(ngdln_max-j);
				c = (barra[i].nodofin+1)*ngdln_max-(ngdln_max-k);

				if(!vapoyo[r] && !vapoyo[c] && c<=r) {
					ren = col = 0;
					for(l=0; l<=r; l++) if(!vapoyo[l]) ren++;
					for(l=0; l<=c; l++) if(!vapoyo[l]) col++;

					matens[LocElemVec(ren-1, col-1)] += barra[i].matriz
						[LocElemVec(barra[i].ngdlni+j, barra[i].ngdlni+k)];
				}
			}
		}
	}

/*	ofstream arch("matriz.txt");

	for(i=0; i<DimVec(dimdis); i++)
		arch << matens[i] << "\t";

	arch.close();*/

	return 0;
}

//------------------------------------------------------------------------------
// Reacciones en barras.
int CalEst::ReacBarras(unsigned caso)
{
	// Regresa 0 ni no hubo ningún error.
	// Regresa -1 si no hubo memoria suficiente.
	// Regresa -2 si no se pudieron obtener las fuerzas en las barras.

	register unsigned i, j, renglon;
	double *vdb=0, *fep=0;

	// Asigna memoria para el vector se desplazamientos por barra:
	if(!(vdb = new double [ngdle_max])) {
		error.MensError("MEM000");
		sprintf(error.mensaje, "Obtención de vectores de reacciones en barras.\n\
Caso de carga: %u",	caso+1);
		return -1;
	}

	// Asigna memoria para el vector de FEP reducido:
	if(!(fep = new double [ngdle_max])) {
		error.MensError("MEM000");
		sprintf(error.mensaje, "Obtención de vectores de reaccioness en barras.\n\
Caso de carga: %u", caso+1);
		return -1;
	}

	for(i=0; i<nbarras; i++) {
		for(j=0; j<ngdle_max; j++) vdb[j] = 0.0;

		/* Forma el vector de desplazamientos por barra (VDB) y el vector de FEP
			reducido: */
		for(j=0; j<barra[i].ngdlni; j++) {
			renglon = (barra[i].nodoini+1)*ngdln_max-(ngdln_max-j);

			vdb[j] = CC[caso].des[renglon];
			fep[j] = CC[caso].vfep[i][j];
		}
		
		for(j=0; j<barra[i].ngdlnf; j++) {
			renglon = (barra[i].nodofin+1)*ngdln_max-(ngdln_max-j);

			vdb[barra[i].ngdlni+j] = CC[caso].des[renglon];
			fep[barra[i].ngdlni+j] = CC[caso].vfep[i][ngdln_max+j];
		}

 		// Multiplica la matriz elemental por el vector de desplazamientos:
		double **matriz=0, *vecres;
		
		matriz = VecMatSim(barra[i].matriz, barra[i].ngdle);

		if(!(CC[caso].vfbarra[i] = Multiplica(matriz, barra[i].ngdle, vdb, barra[i].ngdle))) {
			error.MensError("FUE000");	
			sprintf(error.mensaje, "Caso de carga: %u", caso+1);
			return -2;
		}

		// Libera memoria:
		if(matriz) {
			for(j=0; j<barra[i].ngdle; j++)
				delete [] matriz[j];

			delete [] matriz;
		}

		// Suma las FEP:
		if(!(vecres = Suma(CC[caso].vfbarra[i], barra[i].ngdle,
			fep, barra[i].ngdle))) {
			error.MensError("FUE000");
			sprintf(error.mensaje, "Caso de carga: %u", caso+1);
			return -2;
		}
		else {
			if(CC[caso].vfbarra[i]) delete [] CC[caso].vfbarra[i];
			CC[caso].vfbarra[i] = vecres;
		}
	}

	if(vdb) delete [] vdb;
	if(fep) delete [] fep;

	return 0;
}

//------------------------------------------------------------------------------
// Reacciones en nodos.
int CalEst::ReacNodos(unsigned caso)
{
	// Regresa 0 ni no hubo ningún error.
	// Regresa -1 si no hubo memoria suficiente.

	register unsigned i, j;

	// Asigna memoria para
	for(i=0; i<nnodos; i++)
		if(!(CC[caso].vfnodo[i] = new double [vngdln_max[i]])) {
			error.MensError("MEM000");
			sprintf(error.mensaje, "Obtención de vectores de reacciones en nodos.\n\
Caso de carga: %u", caso+1);
			return -1;
		}

	// Inicializa el vector de fuerzas resultantes en nodos:
	for(i=0; i<nnodos; i++)
		for(j=0; j<vngdln_max[i]; j++)
			CC[caso].vfnodo[i][j] = 0.0;

	// Suma las reacciones de las barras a cada nodo conectado:
	for(i=0; i<nbarras; i++) {
		for(j=0; j<barra[i].ngdlni; j++)
			CC[caso].vfnodo[barra[i].nodoini][j] -= CC[caso].vfbarra[i][j];
		for(j=0; j<barra[i].ngdlnf; j++)
			CC[caso].vfnodo[barra[i].nodofin][j] -= CC[caso].vfbarra[i][barra[i].
				ngdlni+j];
	}

	return 0;
}

//------------------------------------------------------------------------------
// Resolución del sistema de ecuaciones por diversos metodos de manera simultánea.
int CalEst::Resuelve(unsigned nsolver, double toler)
{
	// Regresa 0 si no hubo ningún error.
	// Regresa -1 si no se pudo asignar memoria para resolver el sistema.
	// Regresa -2 si no se pudo resolver el sistema.
	// Regresa -3 si no se pudo resolver el sistema de forma simultánea.

	register unsigned i;
	double **matriz=0, **matres=0, **termind=0;
	
	// Convierte el vector de la matriz ensamblada a matriz simétrica:
	if(!(matriz = VecMatSim(matens, dimdis))) {
		error.MensError("MEM000");
		sprintf(error.mensaje, "No se pudo crear la matriz a resolver.");
		return -1;
	}

	if(!(termind = new double* [ncargas])) {
		error.MensError("MEM000");
		sprintf(error.mensaje, "No se pudo resolver el sistema de forma simultánea.");
		return -3;
	}

	// Copia los punteros del vector de cargas a una matriz de vectores de carga:
	for(i=0; i<ncargas; i++)
		termind[i] = CC[i].vcarga;

	switch(nsolver) {
		case 1:
			{
				FactCrout FC(dimdis, ncargas, matriz, termind, toler);

				if(!(matres = FC.Solucion())) {
					error.MensError("SOL000");
					sprintf(error.mensaje, "Solver: Factorización de Crout.\n\
Solución simultánea");
					return -2;
				}
			}
		break;
		case 2:
			{
				GaussJordan GJ(dimdis, ncargas, matriz, termind, toler);

				if(!(matres = GJ.Solucion())) {
					error.MensError("SOL000");
					sprintf(error.mensaje, "Solver: Eliminación Gaussiana.\n\
Solución simultánea");
					return -2;
				}
			}
		break;
	}

	// Copia la matriz resultante al los vectores resultantes:
	for(i=0; i<ncargas; i++)
		CC[i].desres = matres[i];

	// Libera memoria:
	if(matriz) {
		for(i=0; i<dimdis; i++)
			delete [] matriz[i];

		delete [] matriz;
	}

	if(matres) delete [] matres;

	if(termind) delete [] termind;
	
	return 0;
}

//------------------------------------------------------------------------------
// Resolución del sistema de ecuaciones por diversos metodos de manera independiente.
int CalEst::Resuelve(unsigned caso, unsigned nsolver, double toler)
{
	// Regresa 0 si no hubo ningún error.
	// Regresa -1 si no se pudo asignar memoria para resolver el sistema.
	// Regresa -2 si no se pudo resolver el sistema.
	// Regresa -3 si no se pudo resolver el sistema de forma simultánea.

	register unsigned i;
	double **matriz=0, **matres=0, **termind=0;
	
	// Convierte el vector de la matriz ensamblada a matriz simétrica:
	if(!(matriz = VecMatSim(matens, dimdis))) {
		error.MensError("MEM000");
		sprintf(error.mensaje, "No se pudo crear la matriz a resolver.");
		return -1;
	}

	if(!(termind = new double* [1])) {
		error.MensError("MEM000");
		sprintf(error.mensaje, "No se pudo resolver el sistema de forma independiente.");
		return -3;
	}

	// Copia los punteros del vector de cargas a una matriz de vectores de carga:
	termind[0] = CC[caso].vcarga;

	switch(nsolver) {
		case 1:
			{
				FactCrout FC(dimdis, 1, matriz, termind, toler);

				if(!(matres = FC.Solucion())) {
					error.MensError("SOL000");
					sprintf(error.mensaje, "Solver: Factorización de Crout.\n\
Caso de carga: %u", caso+1);
					return -2;
				}
			}
		break;
		case 2:
			{
				GaussJordan GJ(dimdis, 1, matriz, termind, toler);

				if(!(matres = GJ.Solucion())) {
					error.MensError("SOL000");
					sprintf(error.mensaje, "Solver: Eliminación Gaussiana.\n\
Caso de carga: %u", caso+1);
					return -2;
				}
			}
		break;
	}

	// Copia la matriz resultante al los vectores resultantes:
	CC[caso].desres = matres[0];

	// Libera memoria:
	if(matriz) {
		for(i=0; i<dimdis; i++)
			delete [] matriz[i];

		delete [] matriz;
	}

	if(matres) delete [] matres;

	if(termind) delete [] termind;
	
	return 0;
}

//------------------------------------------------------------------------------
// Establece el nombre del archivo de donde se obtendrán los datos.
int CalEst::Set(char *nomarch)
{
	// Inicializa valores:
	mate  = 0;
	nodo  = 0;
	barra = 0;
	apoyo = 0;
	elas  = 0;
	incl  = 0;
	CC    = 0;

	if(LeeDAT(nomarch)) return -1;
	
	return 0;
}

//------------------------------------------------------------------------------
// Establece los valores de las variables para inicializar el objeto.
int CalEst::Set(
	unsigned nNodos,
	unsigned nBarras, 
	unsigned nApoyos,
	unsigned nCargas,
	unsigned nMats,
	unsigned nElas,
	unsigned nIncl,
	unsigned nDim,
	bool Matrices,
	char *Titulo,
	char *ULongitud,
	char *UFuerza)
{
	// Inicializa valores:
	nnodos  = nNodos;
	nbarras = nBarras;
	napoyos = nApoyos;
	ncargas = nCargas;
	nmats   = nMats;
	nelas   = nElas;
	nincl   = nIncl;
	matrices = Matrices;
	ndims    = nDim;
	if(Titulo) strcpy(titulo, Titulo);
	if(ULongitud) strcpy(uLongitud, ULongitud);
	if(UFuerza) strcpy(uFuerza, UFuerza);

	//nprop = (ndims == 2 || tipo == 1) ? 4 : 7;
	nprop = 4;

	mate  = 0;
	nodo  = 0;
	barra = 0;
	apoyo = 0;
	elas  = 0;
	incl  = 0;
	CC    = 0;

	// Grados de libertad máximos por nodo (según dimensión):
	//ngdln_max = (ndims == 2 || tipo == 1) ? 3 : 6;
	ngdln_max = 3;
	ngdle_max = ngdln_max*2;
	ngdlt = ngdln_max*nnodos;

	if(Inicializa()) return -1;

	return 0;
}

//------------------------------------------------------------------------------
// Encuentra las restricciones a ensamblar.
int CalEst::VectorRestricciones(void)
{
	// Regresa 0 si no hubo ningún error.
	// Regresa -1 si no hay restricciones libres a resolver.

	register unsigned i, j, k, r;
	bool bandera = false;

	// Encuentra el NGDLN máximo para cada nodo:
	for(i=0; i<nbarras; i++) {
		switch(barra[i].tipo) {
			case 2:
				vngdln_max[barra[i].nodoini] = MAX(vngdln_max[barra[i].nodoini],
					ngdln_max);
				vngdln_max[barra[i].nodofin] = MAX(vngdln_max[barra[i].nodofin],
					ngdln_max);
			break;
			case 3:
				vngdln_max[barra[i].nodoini] = MAX(vngdln_max[barra[i].nodoini],
					ngdln_max);
			break;
			case 4:
				vngdln_max[barra[i].nodofin] = MAX(vngdln_max[barra[i].nodofin],
					ngdln_max);
			break;
		}
	}

	for(i=0; i<nnodos; i++) {
		for(j=0; j<napoyos; j++) {
			if(apoyo[j].nodo == i) {
				for(k=0; k<vngdln_max[apoyo[j].nodo]; k++) {
					r = (apoyo[j].nodo+1)*ngdln_max-(ngdln_max-k);
					vapoyo[r] = apoyo[j].drest[k];
				}

				bandera = true;
			}
		}

		if(!bandera) {
			for(j=0; (int)j<vngdln_max[i]; j++) {
				r = (i+1)*ngdln_max-(ngdln_max-j);
				vapoyo[r] = 0;
			}
		}
		else bandera = false;
	}

	// Busca la dimensión que tendrá la matriz ensamblada disminuida:
	dimdis = 0;

	for(i=0; i<ngdlt; i++)
		if(!vapoyo[i]) dimdis++;

	return 0;
}

//******************** FUNCIONES DE OBJETO: Dat_Apoyos *************************
//------------------------------------------------------------------------------
// Destructor de objeto.
Dat_Apoyos::~Dat_Apoyos(void)
{
	Destruye();
}

//------------------------------------------------------------------------------
// Destruye objeto.
void Dat_Apoyos::Destruye(void)
{
	// Libera memoria:
	if(drest) {
		delete [] drest;
		drest = 0;
	}

	if(d) {
		delete [] d;
		d = 0;
	}
}

//------------------------------------------------------------------------------
// Inicializa el objeto.
int Dat_Apoyos::Inicializa(unsigned ngdln_max)
{
	// Regresa 0 si no hubo ningún error.
	// Regresa -1 si no se pudo asignar memoria.
	
	if(!(drest = new bool [ngdln_max])) return -1;
	if(!(d     = new double [ngdln_max])) return -1;

	return 0;
}

//------------------------------------------------------------------------------
// Establece los valores de las variables para inicializar el objeto.
void Dat_Apoyos::Set(unsigned nNodo, bool drX, double dX, bool drY, double dY,
	bool drM, double dM)
{
	nodo = nNodo;
	drest[0] = drX;
	d[0] = dX;
	drest[1] = drY;
	d[1] = dY;
	drest[2] = drM;
	d[2] = dM;
}

//**************** FUNCIONES DE OBJETO: Dat_Apoyos_Elasticos *******************
//------------------------------------------------------------------------------
// Destructor de objeto.
Dat_Apoyos_Elasticos::~Dat_Apoyos_Elasticos(void)
{
	Destruye();
}

//------------------------------------------------------------------------------
// Destruye objeto.
void Dat_Apoyos_Elasticos::Destruye(void)
{
	// Limpia memoria:
	if(k) {
		delete [] k;
		k = 0;
	}
}

//------------------------------------------------------------------------------
// Inicializa el objeto.
int Dat_Apoyos_Elasticos::Inicializa(unsigned ngdln_max)
{
	// Regresa 0 si no hubo ningún error.
	// Regresa -1 si no se pudo asignar memoria.

	if(!(k = new double [ngdln_max])) return -1;

	return 0;
}

//------------------------------------------------------------------------------
// Establece los valores de las variables para inicializar el objeto.
void Dat_Apoyos_Elasticos::Set(unsigned nNodo, double kX, double kY, double kM)
{
	nodo = nNodo;
	k[0] = kX;
	k[1] = kY;
	k[2] = kM;
}

//**************** FUNCIONES DE OBJETO: Dat_Apoyos_Inclinados ******************
//------------------------------------------------------------------------------
// Establece los valores de las variables para inicializar el objeto.
void Dat_Apoyos_Inclinados::Set(unsigned nNodo, double xIncl)
{
	nodo = nNodo;
	xincl = xIncl;
}

//******************* FUNCIONES DE OBJETO: Dat_Materiales **********************
//------------------------------------------------------------------------------
// Destructor de objeto.
Dat_Materiales::~Dat_Materiales(void)
{
	Destruye();
}

//------------------------------------------------------------------------------
// Destruye objeto.
void Dat_Materiales::Destruye(void)
{
	// Libera memoria:
	if(prop) {
		delete [] prop;
		prop = 0;
	}
}

//------------------------------------------------------------------------------
// Inicializa el objeto.
int Dat_Materiales::Inicializa(unsigned nprop)
{
	// Regresa 0 si no hubo ningún error.
	// Regresa -1 si no se pudo asignar memoria.

	if(!(prop = new double [nprop])) return -1;

	return 0;
}

//------------------------------------------------------------------------------
// Establece los valores de las variables para inicializar el objeto.
void Dat_Materiales::Set(double E, double A, double I, double PE)
{
	prop[0] = E;
	prop[1] = A;
	prop[2] = I;
	prop[3] = PE;
}

//********************* FUNCIONES DE OBJETO: Dat_Nodos *************************
//------------------------------------------------------------------------------
// Destructor de objeto.
Dat_Nodos::~Dat_Nodos(void)
{
	Destruye();
}

//------------------------------------------------------------------------------
// Destruye objeto.
void Dat_Nodos::Destruye(void)
{
	// Libera memoria:
	if(coord) {
		delete [] coord;
		coord = 0;
	}
}

//------------------------------------------------------------------------------
// Inicializa el objeto.
int Dat_Nodos::Inicializa(unsigned ndim)
{
	// Regresa 0 si no hubo ningún error.
	// Regresa -1 si no se pudo asignar memoria.

	if(!(coord = new double [ndim])) return -1;

	return 0;
}

//------------------------------------------------------------------------------
// Establece los valores de las variables para inicializar el objeto en 2D.
void Dat_Nodos::Set(double cX, double cY)
{
	coord[0] = cX;
	coord[1] = cY;
}

//------------------------------------------------------------------------------
// Establece los valores de las variables para inicializar el objeto en 3D.
void Dat_Nodos::Set(double cX, double cY, double cZ)
{
	coord[0] = cX;
	coord[1] = cY;
	coord[2] = cZ;
}

//******************** FUNCIONES DE OBJETO: Dat_Barras *************************
//------------------------------------------------------------------------------
// Destructor de objeto.
Dat_Barras::~Dat_Barras(void)
{
	Destruye();
}

//------------------------------------------------------------------------------
// Destruye objeto
void Dat_Barras::Destruye(void)
{
	// Libera memoria:
	if(matriz) {
		delete [] matriz;
		matriz = 0;
	}
}

//------------------------------------------------------------------------------
// Establece los valores de las variables para inicializar el objeto.
void Dat_Barras::Set(unsigned NodoIni, unsigned NodoFin, unsigned Tipo, unsigned 
	nMat, unsigned nDim)
{
	nodoini = NodoIni;
	nodofin = NodoFin;
	tipo = Tipo;
	mat = nMat;

	switch(tipo) {
		case 1:
			ngdlni = ngdlnf = (nDim == 2) ? 2 : 3;
		break;
		case 2:
			ngdlni = ngdlnf = (nDim == 2) ? 3 : 6;
		break;
		case 3:
			ngdlni = (nDim == 2) ? 3 : 6;
			ngdlnf = (nDim == 2) ? 2 : 3;
		break;
		case 4:
			ngdlni = (nDim == 2) ? 2 : 3;
			ngdlnf = (nDim == 2) ? 3 : 6;
		break;
	}				

	ngdle = ngdlni+ngdlnf;
}

//****************** FUNCIONES DE OBJETO: Dat_CasosCarga ***********************
//------------------------------------------------------------------------------
// Constructor del objeto.
Dat_CasosCarga::Dat_CasosCarga(void)
{
	// Inicializa valores:
	fnodo  = 0;
	fbarra = 0;
	frect  = 0;
}

//------------------------------------------------------------------------------
// Destructor del objeto.
Dat_CasosCarga::~Dat_CasosCarga(void)
{
	Destruye();
}

//------------------------------------------------------------------------------
// Destruye el objeto.
void Dat_CasosCarga::Destruye(void)
{
	// Libera memoria:
	if(fnodo) {
		delete [] fnodo;
		fnodo = 0;
	}

	if(fbarra) {
		delete [] fbarra;
		fbarra = 0;
	}

	if(frect) {
		delete [] frect;
		frect = 0;
	}
}

//------------------------------------------------------------------------------
// Establece los valores de las variables para inicializar el objeto.
int Dat_CasosCarga::Set(
		unsigned nFNodos,
		unsigned nFBarras,
		unsigned nFRect,
		bool PPropio,
		unsigned nNodos,
		unsigned nBarras,
		unsigned ngdln_max,
		unsigned ngdle_max,
		char *Titulo
	)
{
	// Regresa 0 si no hubo ningún error.
	// Regresa -1 si no se pudo asignar memoria.
	
	// Inicializa valores:
	nfnodos  = nFNodos;
	nfbarras = nFBarras;
	nfrect   = nFRect;
	PesoPropio = PPropio;
	nnodos   = nNodos;
	nbarras  = nBarras;
	if(Titulo) strcpy(titulo, Titulo);

	fnodo  = 0;
	fbarra = 0;
	frect  = 0;

	if(Inicializa(ngdln_max, ngdle_max)) return -1;

	return 0;
}

//------------------------------------------------------------------------------
// Inicializa el objeto.
int Dat_CasosCarga::Inicializa(unsigned ngdln_max, unsigned ngdle_max)
{
	// Regresa 0 si no se encontró ningún error.
	// Regresa -1 si se no se pudo asignar memoria suficiente.

	register unsigned i, j;

	// Asigna memoria para las fuerzas en los casos de carga:
	if(!(fnodo  = new Dat_CNodos [nfnodos])) return -1;
	if(!(fbarra = new Dat_CBarras [nfbarras])) return -1;
	if(!(frect  = new Dat_CRectangulares [nfrect])) return -1;
	
	// Asigna memoria a las fuerzas de los nodos:
	for(i=0; i<nfnodos; i++)
		if(fnodo[i].Inicializa(ngdln_max)) return -1;

	// Asigna memoria para las fuerzas de empotramiento perfecto (FEP):
	if(!(vfep = new double* [nbarras])) return -1;

	for(i=0; i<nbarras; i++)
		if(!(vfep[i] = new double [ngdle_max])) return -1;

	// Inicializa las FEP a cero:
	for(i=0; i<nbarras; i++)
		for(j=0; j<ngdle_max; j++)
			vfep[i][j] = 0.0;

	return 0;
}

//********************* FUNCIONES DE OBJETO: Dat_CBarras ***********************
//------------------------------------------------------------------------------
// Establece los valores de las variables para inicializar el objeto.
void Dat_CBarras::Set(unsigned nBarra, double D, double F)
{
	barra = nBarra;
	d = D;
	f = F;
}

//********************* FUNCIONES DE OBJETO: Dat_CNodos ************************
//------------------------------------------------------------------------------
// Destructor de objeto.
Dat_CNodos::~Dat_CNodos(void)
{
	Destruye();
}

//------------------------------------------------------------------------------
// Destruye objeto.
void Dat_CNodos::Destruye(void)
{
	// Limpia memoria:
	if(f) {
		delete [] f;
		f = 0;
	}
}

//------------------------------------------------------------------------------
// Inicializa el objeto.
int Dat_CNodos::Inicializa(unsigned ngdln_max)
{
	// Regresa 0 si no hubo ningún error.
	// Regresa -1 si no se pudo asignar memoria.

	if(!(f = new double [ngdln_max])) return -1;

	return 0;
}

//------------------------------------------------------------------------------
// Establece los valores de las variables para inicializar el objeto.
void Dat_CNodos::Set(unsigned nNodo, double fX, double fY, double fM)
{
	nodo = nNodo;
	f[0] = fX;
	f[1] = fY;
	f[2] = fM;
}

//******************* FUNCIONES DE OBJETO: Dat_CRectangulares ******************
//------------------------------------------------------------------------------
// Establece los valores de las variables para inicializar el objeto.
void Dat_CRectangulares::Set(unsigned nBarra, double F)
{
	barra = nBarra;
	f = F;
}

//******************* FUNCIONES DE OBJETO: Res_CalculoCC ***********************
//------------------------------------------------------------------------------
// Constructor del objeto.
Res_CalculoCC::Res_CalculoCC(void)
{
	// Inicializa valores:
	desres  = 0;
	des     = 0;
	esf		= 0;
	vcarga  = 0;
	vfep    = 0;
	vfbarra = 0;
	vfnodo  = 0;
}

//------------------------------------------------------------------------------
// Destructor del objeto.
Res_CalculoCC::~Res_CalculoCC(void)
{
	Destruye();
}

//------------------------------------------------------------------------------
// Destruye el objeto.
void Res_CalculoCC::Destruye(void)
{
	register unsigned i;

	// Libera memoria:
	if(desres) {
		delete [] desres;
		desres = 0;
	}

	if(des) {
		delete [] des;
		des = 0;
	}
	
	if(esf) {
		delete [] esf;
		esf = 0;
	}

	if(vcarga) {
		delete [] vcarga;
		vcarga = 0;
	}
	
	if(vfep) {
		for(i=0; i<nbarras; i++) delete [] vfep[i];
		delete [] vfep;
		vfep = 0;
	}

	if(vfbarra) {
		for(i=0; i<nbarras; i++) delete [] vfbarra[i];
		delete [] vfbarra;
		vfbarra = 0;
	}

	if(vfnodo) {
		for(i=0; i<nnodos; i++) delete [] vfnodo[i];
		delete [] vfnodo;
		vfnodo = 0;
	}
}

//******************** FUNCIONES DE OBJETO: Res_Calculo ************************
//------------------------------------------------------------------------------
// Constructor del objeto.
Res_Calculo::Res_Calculo(void)
{
	// Inicializa valores:
	vngdln_max = 0;
	vapoyo     = 0;
	matens     = 0;
}

//------------------------------------------------------------------------------
// Destructor del objeto.
Res_Calculo::~Res_Calculo(void)
{
	Destruye();
}

//------------------------------------------------------------------------------
// Destruye el objeto.
void Res_Calculo::Destruye(void)
{
	// Libera memoria:
	if(vngdln_max) {
   		delete [] vngdln_max;
		vngdln_max = 0;
	}
	
	if(vapoyo) {
   		delete [] vapoyo;
		vapoyo = 0;
	}
	
	if(matens) {
   		delete [] matens;
		matens = 0;
   }
}

/*
+------------------------------------------------------------------------------+
¦                       F U N C I O N E S   V A R I A S                        ¦
+------------------------------------------------------------------------------+
*/
//------------------------------------------------------------------------------
// Reemplaza un caracter por otro.
void Remplaza(char *cad, char quita, char pon)
{
	register unsigned i;
    unsigned longitud = strlen(cad);

	for(i=0; i<longitud; i++)
		if(cad[i] == quita) cad[i] = pon;
}
