#include <R.h>
#include <Rdefines.h>

double getListElement(char **nombres, const char *str, int len_nombres){
	int i, elmt = 0;
	for (i = 0; i < len_nombres; i++){
		//Rprintf("%s\n", str);
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
	
	//Guardamos toda la lista en memoria para que sea más rápido el acceso
	char **p_anotaciones[len_anotaciones];
	for(categoria1 = 0; categoria1 < len_anotaciones; ++categoria1){	
		len = Rf_length(VECTOR_ELT(anotaciones, categoria1));			
		p_anotaciones[categoria1] = malloc(len*sizeof(char*));		
		for(k = 0; k < len; ++k){
			p_anotaciones[categoria1][k] = R_alloc(strlen(CHAR(STRING_ELT(VECTOR_ELT(anotaciones, categoria1), k))), sizeof(char));
			strcpy(p_anotaciones[categoria1][k], CHAR(STRING_ELT(VECTOR_ELT(anotaciones, categoria1), k))); 
		}
	}
len_anotaciones = 2;

	//Categoria por categoria vamos buscando los padres en común
	for(categoria1 = 0; categoria1 < len_anotaciones; ++categoria1){

		len = Rf_length(VECTOR_ELT(anotaciones, categoria1));

		//Recorro la segunda categoria
		for(categoria2 = categoria1; categoria2 < len_anotaciones; ++categoria2){
			len2 = Rf_length(VECTOR_ELT(anotaciones, categoria2));	
			mica = 0;

			//Vamos buscando matches entre las dos categorías		
			for(l = 0; l < len; ++l){
				for(m = 0; m < len2; ++m){
					//Rprintf("%s vs %s\n", p_anotaciones[categoria1][l], p_anotaciones[categoria2][m]);
					if(strcmp(p_anotaciones[categoria1][l], p_anotaciones[categoria2][m]) == 0){
						iic = getListElement(nombres, p_anotaciones[categoria1][l], len_nombres);
					//	iic = getListElement(nombres, "UNO", len_nombres);
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

	//Limpiamos
	//for(categoria1 = 0; categoria1 < len_anotaciones; ++categoria1){	
	//	free(p_anotaciones[categoria1]);		
	//}

	UNPROTECT(2);
	return(semsimM);
}
