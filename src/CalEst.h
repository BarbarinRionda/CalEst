/*
+------------------------------------------------------------------------------+
¦     DEFINICIÓN DE LOS OBJETOS NECESARIOS PARA EL CALCULO DE ESTRUCTURAS      ¦
+------------------------------------------------------------------------------+
*/

#ifndef _CALEST

#define _CALEST

#include "EdoError.h"

//°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
// Definiciones:
#define PI	3.14159265358979323846
#define GAR(g) (g*PI/180)
#define RAG(r) (r*180/PI)
#define MAX(a, b)  ((a > b) ? a : b)
#define MIN(a, b)  ((a < b) ? a : b)

//°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
// Definiciones de indicadores de cálculo:
#define CAL_MATELEM		0x00000001	// Indicador de cálculo de matrices elementales.
#define CAL_MATENSA		0x00000011	// Indicador de cálculo de matriz ensamblada.
#define CAL_DESPLAZ		0x00000111	// Indicador de cálculo de desplazamientos.
#define CAL_REACBAR		0x00001111	// Indicador de cálculo de reacciones en barras.
#define CAL_REACNOD		0x00011111	// Indicador de cálculo de reacciones en nodos.
#define CAL_ESFUERZ		0x00101111	// Indicador de cálculo de esfuerzos.
#define CAL_TODO		0x11111111	// Indicador de cálculo total.

//>>>>>>>>>>>>>>>>>>>>>> OBJETOS GEOMÉTRICOS Y DE PROPIEDADES <<<<<<<<<<<<<<<<<<
//****************************** Objeto: Dat_Nodos *****************************
//°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
// Objeto: datos de nodos:
class Dat_Nodos {
public:
	double *coord;					// Coordenadas del nodo.

	Dat_Nodos(void) { coord = 0; }
	~Dat_Nodos(void);
	void Destruye(void);
	int Inicializa(unsigned ndim);
	void Set(
		double cX,					// Coordenada X.
		double cY					// Coordenada Y.
	);
	void Set(
		double cX,					// Coordenada X.
		double cY,					// Coordenada Y.
		double cZ					// Coordenada Z.
	);
};

//****************************** Objeto: Dat_Barras ****************************
//°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
// Objeto: datos de barras:
class Dat_Barras {
public:
	unsigned mat;					// Tipo de material de la barra.
	unsigned tipo;					// Tipo de barra (1=Art-Art, 2=Emp-Emp, 3=Emp-Art, 4=Art-Emp)
	unsigned nodoini;				// Nodo inicial.
	unsigned nodofin;				// Nodo final.
	unsigned ngdlni;				// NGDL en el nodo inicial.
	unsigned ngdlnf;				// NGDL en el nodo final.
	unsigned ngdle;					// NGDL en el elemento = ngdlni+ngdlnf.
	double *matriz;					// Matriz de rigidez de la barra.
	
	Dat_Barras(void) { }
	~Dat_Barras(void);
	void Destruye(void);
	void Set(
		unsigned NodoIni,			// Nodo inicial.
		unsigned NodoFin,			// Nodo final.
		unsigned Tipo,				// Tipo de barra (1=Art-Art, 2=Emp-Emp, 3=Emp-Art, 4=Art-Emp)
		unsigned nMat=0,			// Número de material.
		unsigned nDim=2				// Número de dimensiones.
	);
};

//**************************** Objeto: Dat_Materiales **************************
//°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
// Objeto: datos de materiales:
class Dat_Materiales {
public:
	double *prop;					// Propiedades del material.

	Dat_Materiales(void) { prop = 0; }
	~Dat_Materiales(void);
	void Destruye(void);
	int Inicializa(unsigned nprop);
	void Set(
		double E,					// Módulo de elasticidad.
		double A,					// Área de la sección.
		double I,					// Momento de inercia.
		double PE					// Peso específico.
	);
};

//****************************** Objeto: Dat_Apoyos ****************************
//°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
// Objeto: datos de restricciones (apoyos):
class Dat_Apoyos {
public:
	unsigned nodo;					// Número de nodo.
	bool *drest;					// Especifica la dirección restringida.
	double *d;						// Valor del desplazamiento de la restricción.

	Dat_Apoyos(void) { drest = 0; d = 0; }
	~Dat_Apoyos(void);
	void Destruye(void);
	int Inicializa(unsigned ngdln_max);
	void Set(
		unsigned nNodo,				// Número de nodo.
		bool     drX,				// Especifica la dirección X restringida.
		double   dX,				// Valor del desplazamiento de la restricción en X.
		bool     drY,				// Especifica la dirección Y restringida.
		double   dY,				// Valor del desplazamiento de la restricción en Y.
		bool     drM=false,			// Especifica la dirección Z restringida.
		double   dM=0				// Valor del desplazamiento de la restricción en Z.
	);
};

//************************ Objeto: Dat_Apoyos_Elasticos ************************
//°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
// Objeto: apoyos elásticos:
class Dat_Apoyos_Elasticos {
public:
	unsigned nodo;					// Nodo que contiene el apoyo elástico.
	double *k;						// Rigideces del resorte.

	Dat_Apoyos_Elasticos(void) { k = 0; }
	~Dat_Apoyos_Elasticos(void);
	void Destruye(void);
	int Inicializa(unsigned ngdln_max);
	void Set(
		unsigned nNodo,				// Número de nodo.
		double   kX,				// Rigidez del resorte en la dirección X.
		double   kY,				// Rigidez del resorte en la dirección Y.
		double   kM=0				// Rigidez del resorte del momento (res. helicoidal).
	);
};

//************************ Objeto: Dat_Apoyos_Inclinados ***********************
//°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
// Objeto: apoyos inclinados:
class Dat_Apoyos_Inclinados {
public:
	unsigned nodo;					// Nodo que contiene el apoyo inclinado.
	double xincl;					// Inclinación del apoyo (en grados).

	Dat_Apoyos_Inclinados(void) { }
	void Set(
		unsigned nNodo,				// Nodo que contiene el apoyo inclinado.
		double   xIncl				// Inclinación del apoyo (en grados).
	);
};

//>>>>>>>>>>>>>>>>>>>>>>>>>>>> OBJETOS DE CASOS DE CARGA <<<<<<<<<<<<<<<<<<<<<<<
//******************************* Objeto: Dat_CNodos ***************************
//°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
// Objeto: cargas en los nodos:
class Dat_CNodos {
public:
	unsigned nodo;					// Número de nodo.
	double	 *f;					// Datos de las fuerzas.

	Dat_CNodos(void) { f = 0; }
	~Dat_CNodos(void);
	void Destruye(void);
	int Inicializa(unsigned ngdln_max);
	void Set(
		unsigned nNodo,				// Número de nodo.
		double   fX,				// Carga en la dirección X.
		double   fY,				// Carga en la dirección X.
		double   fM=0				// Momento.
	);
};

//***************************** Objeto: Dat_CBarras ****************************
//°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
// Objeto: cargas concentradas dentro de las barras:
class Dat_CBarras {
public:
	unsigned barra;					// Número de barra.
	double	 d;						// Distancia del nodo inicial.
	double	 f;						// Carga perpendicular a la barra.

	Dat_CBarras(void) { };
	void Set(
		unsigned nBarra,			// Número de barra.
		double   D,					// Distancia del nodo inicial.
		double   F					// Carga perpendicular a la barra.
	);
};

//*************************** Objeto: Dat_CRectangulares ***********************
//°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
// Objeto: cargas distribuidas rectangulares:
class Dat_CRectangulares {
public:
	unsigned barra;					// Número de barra.
	double	 f;						// Carga perpendicular a la barra.

	Dat_CRectangulares(void) { }
	void Set(
		unsigned nBarra,			// Número de barra.
		double   F					// Carga perpendicular a la barra.
	);
};

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>> OBJETOS DE RESULTADOS <<<<<<<<<<<<<<<<<<<<<<<<<<
//****************************** Objeto: Res_Calculo ***************************
//°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
// Objeto: resultados del cálculo según geometría y propiedades:
class Res_Calculo {
public:
	unsigned ngdln_max;				// No. máximo de grados de libertad por nodo.
	unsigned ngdle_max;				// No. máximo de grados de libertad por elemento.
	unsigned ngdlt;					// No. de grados de libertad totales.
	unsigned *vngdln_max;			// Vector de NGDLN máximo por nodo.
	int		 *vapoyo;				// Vector de restriciones nodales. (-1=res. no existente, 0=libre de mov., 1=sujeto).
	unsigned dimdis;				// Dimensión de la matriz ensamblada disminuida.
	double	 *matens;				// Vector de matriz ensamblada global (disminuida).
	
	Res_Calculo(void);
	~Res_Calculo(void);
	void Destruye(void);
};

//**************************** Objeto: Res_CalculoCC ***************************
//°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
// Objeto: resultados del cálculo por cada caso de carga:
class Res_CalculoCC {
public:
	unsigned nbarras;
	unsigned nnodos;

	double *desres;					// Desplazamientos resultantes del cálculo.
	double *des;					// Desplazamientos.
	double *esf;					// Vector de esfuerzos elementales.
	double *vcarga;					// Vector de cargas externas.
	double **vfnodo;				// Vector de fuerzas resultantes en nodos.
	double **vfbarra;				// Vector de fuerzas resultantes en barras.
	double **vfep;					// Matriz de vectores de FEP.
	double esf_min;					// Esfuerzo mínimo.
	double esf_max;					// Esfuerzo máximo.

	Res_CalculoCC(void);
	~Res_CalculoCC(void);
	void Destruye(void);
};

//>>>>>>>>>>>>>>>>>>>>>>>>>>> OBJETO DE CASOS DE CARGA <<<<<<<<<<<<<<<<<<<<<<<<<
//*************************** Objeto: Dat_CasosCarga ***************************
//°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
// Objeto: casos de carga:
class Dat_CasosCarga: public Res_CalculoCC {
public:
	char titulo[256];				// Título del caso de carga.
	unsigned nfnodos;				// Número de cargas nodales.
	unsigned nfbarras;				// Número de cargas en las barras.
	unsigned nfrect;				// Número de cargas uniformemente distribuidas.
	bool PesoPropio;				// Indica si se debe tomar en cuenta el peso propio.
	Dat_CNodos *fnodo;				// Cargas en los nodos.
	Dat_CBarras *fbarra;    		// Cargas en las barras.
	Dat_CRectangulares *frect;		// Cargas uniformemente repartidas.

	Dat_CasosCarga(void);
	~Dat_CasosCarga(void);
	void Destruye(void);
	int Inicializa(unsigned ngdln_max, unsigned ngdle_max);
	int Set(
		unsigned nFNodos,			// Número de cargas nodales.
		unsigned nFBarras,			// Número de cargas en las barras.
		unsigned nFRect,			// Número de cargas uniformemente distribuidas.
		bool     PPropio,			// Indica si se debe tomar en cuenta el peso propio.
		unsigned nNodos,			// Número de nodos.
		unsigned nBarras,			// Número de barras.
		unsigned ngdln_max,			// Número de grados de libertad nodales máximos.
		unsigned ngdle_max,			// Número de grados de libertad elementales máximos.
		char     *Titulo=0			// Título del caso de carga.
	);
};

//>>>>>>>>>>>>>>>>>>>>>>>>>>>> OBJETO DE ESTRUCTURA <<<<<<<<<<<<<<<<<<<<<<<<<<<<
//******************************** Objeto: CalEst ******************************
//°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
class CalEst: public Res_Calculo {
public:
	char	 titulo[256];			// Título del problema.
    char	 uLongitud[256];		// Unidades de longitud.
    char	 uFuerza[256];			// Unidades de fuerza.
	unsigned nnodos;				// No. de puntos nodales.
	unsigned nbarras;				// No. de barras.
	unsigned napoyos;				// No. de restricciones (apoyos).
	unsigned ncargas;				// No. de casos de carga a la que está sujeta la estructura.
	unsigned nmats;					// No. de materiales.
	unsigned nprop;					// No. de propiedades de los materiales.
	unsigned ndims;					// No. de dimensiones.
	unsigned nelas;					// No. de apoyos elásticos.
	unsigned nincl;					// No. de apoyos inclinados.
    bool	 matrices;				// Indicador de salida de matrices.
	Dat_Materiales *mate;			// Propiedades de los materiales.
	Dat_Nodos *nodo;				// Datos de los nodos.
	Dat_Barras *barra;				// Datos de las barras.
	Dat_Apoyos *apoyo;				// Datos de las restricciones (apoyos).
	Dat_Apoyos_Elasticos *elas;		// Datos de los apoyos elásticos.
	Dat_Apoyos_Inclinados *incl;	// Datos de los apoyos inclinados.
	Dat_CasosCarga *CC;    			// Datos de los casos de carga.
	
private:
	void   ApoyosElasticos(void);
	int    ApoyosInclinados(unsigned bar);
	int    AsignaMemoria1(void);
	int    AsignaMemoria2(unsigned caso);
	double *BarraArtArt(unsigned bar);
	double *BarraArtEmp(unsigned bar);
	double *BarraEmpArt(unsigned bar);
	double *BarraEmpEmp(unsigned bar);
	int    Cargas(unsigned caso);
	void   Desplazamientos(unsigned caso);
	void   Esfuerzos(unsigned caso);
	int    Inicializa(void);
	int    LeeDAT(char *nomarch);
	void   LonAng(unsigned bar, double &longitud, double &angulo);
	void   LonAng(unsigned bar, double &longitud, double &angulo, double &angulo1);
	int    MatrizElemental(unsigned bar);
	int    MatrizEnsamblada(void);
	int    ReacBarras(unsigned caso);
	int    ReacNodos(unsigned caso);
	int	   Resuelve(unsigned nsolver, double toler);
	int    Resuelve(unsigned caso, unsigned nsolver, double toler);
	int    VectorRestricciones(void);

public:
	EdoError error;

	CalEst(void);
	~CalEst(void);
	int Calcula(unsigned cal, double toler, bool solsim=true, unsigned nsolver=1);
	void Destruye(void);
	int EscribeRES(char *nomarch);
	int Set(char *nomarch);
	int Set(
		unsigned nNodos,			// No. de puntos nodales.
		unsigned nBarras,			// No. de barras.
		unsigned nApoyos,			// No. de restricciones (apoyos).
		unsigned nCargas,			// No. de casos de carga a la que está sujeta la estructura.
		unsigned nMats,				// No. de materiales.
		unsigned nElas,				// No. de apoyos elásticos.
		unsigned nIncl,				// No. de apoyos inclinados.
		unsigned nDim=2,			// No. de dimensiones.
		bool     Matrices=false,	// Indicador de salida de matrices.
		char     *Titulo=0,			// Título del problema.
		char     *ULongitud=0,		// Unidades de longitud.
		char     *UFuerza=0			// Unidades de fuerza.
	);
};

#endif