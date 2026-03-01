/*
+------------------------------------------------------------------------------+
¦                DECLARACIÓN DE FUNCIONES DE ALGEBRA LINEAL                    ¦
+------------------------------------------------------------------------------+
*/

unsigned DimMat(unsigned dimv);
unsigned DimVec(unsigned dimm);
unsigned LocElemVec(unsigned ren, unsigned col);
double *MatSimVec(double **mat, unsigned dimv);
double *Multiplica(double fac1, unsigned dimv, double *fac2);
double *Multiplica(double **mat, unsigned dimm, double *vec, unsigned dimv);
double *Suma(double *vec1, unsigned r1, double *vec2, unsigned r2);
double **VecMatSim(double *vec, unsigned dimm);
double *EliminaElem(double *vector, unsigned &nelems, unsigned nelem);
