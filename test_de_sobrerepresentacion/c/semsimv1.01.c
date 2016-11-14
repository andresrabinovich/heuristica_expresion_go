#include <R.h>
#include <Rdefines.h>

double getListElement(char **nombres, const char *str, int len_nombres){
	int i, elmt = 0;
	for (i = 0; i < len_nombres; i++){
		if(strcmp(nombres[i], str) == 0) {
			elmt = i;
			break;
		}
	}
	return elmt;
}

SEXP semsim(SEXP anotaciones, SEXP sic){
	int categoria1, categoria2, k, l, m, len_anotaciones, len_nombres, len, len2, iic;
	double mica = 0;
	double *p_semsimM, *ic;
	SEXP semsimM, dimnames, names, names_categorias;
	len_anotaciones = Rf_length(anotaciones);

	ic = REAL(sic);

	//Traemos los nombres de la matriz de ic	
	names = getAttrib(sic, R_NamesSymbol);
	len_nombres = length(names);
	char *nombres[len_nombres];		
	for(k = 0; k < len_nombres; ++k){
		nombres[k] = R_alloc(strlen(CHAR(STRING_ELT(names, k))), sizeof(char));
		strcpy(nombres[k], CHAR(STRING_ELT(names, k))); 
	}

	//Traemos los nombres de la lista de anotaciones para armar la matriz de semsim
	names_categorias = getAttrib(anotaciones, R_NamesSymbol);
	
	//Generamos la matriz de semsim y la asociamos a un puntero para modificarla
	PROTECT(semsimM = allocMatrix(REALSXP, len_anotaciones, len_anotaciones));
	p_semsimM = REAL(semsimM);
	
	/*
	//Guardamos toda la lista en memoria para que sea más rápido el acceso
	int *p_anotaciones[len_anotaciones];
	for(categoria1 = 0; categoria1 < len_anotaciones; ++categoria1){	
		//Guardo la primer categoría en un vector porque vamos a estar accediendo mucho
		len = Rf_length(VECTOR_ELT(anotaciones, categoria1));			
		int padres[len];		
		for(k = 0; k < len; ++k){
			//padres[k] = R_alloc(strlen(CHAR(STRING_ELT(VECTOR_ELT(anotaciones, categoria1), k))), sizeof(char));
			//strcpy(padres[k], CHAR(STRING_ELT(VECTOR_ELT(anotaciones, categoria1), k))); 
			//padres[k] = (int*)R_alloc(1, sizeof(int));
			padres[k] = INTEGER(VECTOR_ELT(anotaciones, categoria1))[k];
			printf("%d\n", padres[k]);
			
		}
		p_anotaciones[categoria1] = padres;
	}
*/
		
	//Categoria por categoria vamos buscando los padres en común
	for(categoria1 = 0; categoria1 < len_anotaciones; ++categoria1){
		

		//Guardo la primer categoría en un vector porque vamos a estar accediendo mucho
		len = Rf_length(VECTOR_ELT(anotaciones, categoria1));			
		//char *padres[len];		
		//for(k = 0; k < len; ++k){
		//	padres[k] = R_alloc(strlen(CHAR(STRING_ELT(VECTOR_ELT(anotaciones, categoria1), k))), sizeof(char));
		//	strcpy(padres[k], CHAR(STRING_ELT(VECTOR_ELT(anotaciones, categoria1), k))); 
		//}

		//Recorro la segunda categoria
		for(categoria2 = categoria1; categoria2 < len_anotaciones; ++categoria2){
			
			len2 = Rf_length(VECTOR_ELT(anotaciones, categoria2));	
			mica = 0;

			//Vamos buscando matches entre las dos categorías		
			for(l = 0; l < len; ++l){
				for(m = 0; m < len2; ++m){
					if(INTEGER(VECTOR_ELT(anotaciones, categoria1))[l] == INTEGER(VECTOR_ELT(anotaciones, categoria2))[m]){
						iic = INTEGER(VECTOR_ELT(anotaciones, categoria1))[l];
						if(mica < ic[iic]) mica = ic[iic];
					}
				}
			}
			p_semsimM[categoria1*len_anotaciones + categoria2] = mica;
			p_semsimM[categoria2*len_anotaciones + categoria1] = mica;
			//Rprintf("Match: %s vs. %s (%g)\n", CHAR(STRING_ELT(getAttrib(anotaciones, R_NamesSymbol), categoria1)), CHAR(STRING_ELT(getAttrib(anotaciones, R_NamesSymbol), categoria2)), mica);
		}
	}


	PROTECT(dimnames = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(dimnames, 0, names_categorias);
	SET_VECTOR_ELT(dimnames, 1, names_categorias);
	dimnamesgets(semsimM, dimnames);
	UNPROTECT(2);
	return(semsimM);
}
