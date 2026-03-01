/*
+------------------------------------------------------------------------------+
¦         F U N C I O N E S   P A R A   A L G E B R A   L I N E A L            ¦
+------------------------------------------------------------------------------+
*/

//------------------------------------------------------------------------------
/* Función para encontrar la dimensión de una matriz simetrica, dada la dimensión
de el vector que guarda su triangular.*/
unsigned DimMat(unsigned dimv)
{
	register unsigned i, res;

	res = dimv;
	
	for(i=0; i<dimv; i++) {
		res -= i;
		if(!res) break;
	}

	return i;
}

//------------------------------------------------------------------------------
/* Función para encontrar la dimensión de un vector que contenga la matriz triangular
dada la dimensión de la matriz cuadrada simetrica.*/
unsigned DimVec(unsigned dimm)
{
	register unsigned i, res=0;
	
	for(i=dimm; i>0; i--) res += i;

	return res;
}

//------------------------------------------------------------------------------
// Localiza la posición de un elemento de una matriz simétrica en un vector,
unsigned LocElemVec(unsigned ren, unsigned col)
{
	// Regresa la posición en el vector.

	register unsigned i, j, cont=0;

	if(col > ren) {
		i= ren;
		ren = col;
		col = i;
	}

	for(i=0; i<ren; i++)
		for(j=0; j<=i; j++)
			cont++;

	cont += col;

	return cont;
}

//------------------------------------------------------------------------------
// Convierte una matriz simétrica a un vector.
double *MatSimVec(double **mat, unsigned dimv)
{
	// Regresa la dirección del vector.
	// Regresa 0 si no se pudo asignar memoria.

	register unsigned i, j, cont=0, dimm;
	double *vec;

	if(!(vec = new double [dimv])) return 0;

	dimm = DimMat(dimv);

	// Pasa la información debajo de la diagonal principal a un vector:
	for(i=0; i<dimm; i++) {
		for(j=0; j<=i; j++) {
			vec[cont] = mat[i][j];
			cont++;
		}
	}

	return vec;
}

//------------------------------------------------------------------------------
// Función para multiplicar un escalar por un vector.
double *Multiplica(double fac1, unsigned dimv, double *fac2)
{
	// Regresa 0 si no se pudo efectuar la multiplicación.
	// Regresa el puntero del vector.

	register unsigned i;
	double *res;

	if(!(res = new double [dimv])) return 0;

	for(i=0; i<dimv; i++)
		res[i] = fac1*fac2[i];

	return res;
}

//------------------------------------------------------------------------------
// Función para multiplicar una matriz cuadrada por un vector.
double *Multiplica(double **mat, unsigned dimm, double *vec, unsigned dimv)
{
	// Regresa 0 si no se pudo realizar la multiplicación.
	// Regresa el puntero del vector.

	if(dimm != dimv) return 0;

	register unsigned i, j;
	double *vecres;

	// Asigna memoria para el vector resultante:
	if(!(vecres = new double [dimv])) return 0;

	// Inicializa el vector resultante:
	for(i=0; i<dimm; i++) vecres[i] = 0.0;

	// Realiza la multiplicación:
	for(i=0; i<dimm; i++)
		for(j=0; j<dimm; j++)
			vecres[i] += mat[i][j]*vec[j];

	return vecres;
}

//------------------------------------------------------------------------------
// Función para sumar dos vectores:
double *Suma(double *vec1, unsigned r1, double *vec2, unsigned r2)
{
	// Regresa 0 si no se pudo realizar la suma.
	// Regresa el puntero del vector resultante de la suma.

	if(r1 != r2) return 0;

	register unsigned i;
	double *vec3;

	// Asigna memoria para el vector resultante:
	if(!(vec3 = new double [r1])) return 0;

	// Realiza la suma:
	for(i=0; i<r1; i++)
		vec3[i] = vec1[i]+vec2[i];

	return vec3;
}

//------------------------------------------------------------------------------
// Convierte un vector en una matriz simétrica.
double **VecMatSim(double *vec, unsigned dimm)
{
	// Regresa la dirección de la matriz.
	// Regresa 0 si no se pudo asignar memoria.

	register unsigned i, j, cont=0;
	double **mat;

	if(!(mat = new double* [dimm])) return 0;
	for(i=0; i<dimm; i++) if(!(mat[i] = new double [dimm])) return 0;

	// Inserta el vector debajo de la diagonal principal y completa matriz:
	for(i=0; i<dimm; i++)
		for(j=0; j<=i; j++) {
			mat[i][j] = mat[j][i] = vec[cont];
			cont++;
		}
	
	return mat;
}

//------------------------------------------------------------------------------
// Elimina un elemento de un vector.
double *EliminaElem(double *vector, unsigned &nelems, unsigned nelem)
{
	// Regresa 0 si no hubo memoria suficiente para asignar.
	// Regresa la direccion del puntero.
	
	register unsigned i, j=0;
	double *vec;

	// Asigna memoria a el vector disminuido:
	if(!(vec = new double [nelems-1])) return 0;

	// Copia los elementos no eliminados:
	if(nelem >= 0) {
		for(i=0; i<nelems; i++) {
			if(i != nelem) {
				vec[j] = vector[i];
				j++;
			}
		}
		
		nelems--;
	}

	return vec;
}