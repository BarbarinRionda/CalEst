/*
+------------------------------------------------------------------------------+
¦       F U N C I O N E S   P A R A   M A N E J O   D E   E R R O R E S        ¦
+------------------------------------------------------------------------------+
*/

#include <string.h>
#include <stdio.h>

#include "EdoError.h"

//********************* FUNCIONES DE OBJETO: EdoError **************************
//------------------------------------------------------------------------------
// Constructor de objeto.
EdoError::EdoError(void)
{
	estado = false;
	mensaje[0] = '\0';
}

//------------------------------------------------------------------------------
// Establece el mensaje de error a partir de una cadena.
void EdoError::MensError(const char cod[7])
{
	estado = true;
	
	strcpy(codigo, cod);
}

//------------------------------------------------------------------------------
// Establece el mensaje de error a partir de una cadena y un entero.
void EdoError::MensError(const char cod1[4], unsigned cod2)
{
	estado = true;

	sprintf(codigo, "%s%u", cod1, cod2);
}

//------------------------------------------------------------------------------
// Establece el mensaje de error a partir de una cadena y un mensaje adicional.
void EdoError::MensError(const char cod[7], char *mens)
{
	estado = true;
	
	strcpy(codigo, cod);
	strcpy(mensaje, mens);
}
