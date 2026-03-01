/*
+------------------------------------------------------------------------------+
¦               DEFINICIÓN DEL OBJETO DE MANIPULACIÓN DE ERRORES               ¦
+------------------------------------------------------------------------------+
*/

#ifndef EDOERROR
#define EDOERROR

//******************************* Objeto: EdoError *****************************
//°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
class EdoError {
public:
	EdoError(void);

	char codigo[7];							// Código de error (3 letras, 3 números)
	bool estado;							// Indica si se ha ocacionado un error.
	char mensaje[256];						// Mensaje adicional.

	void MensError(const char cod[7]);
	void MensError(const char cod1[4], unsigned cod2);
	void MensError(const char cod[7], char *mens);
};

#endif
