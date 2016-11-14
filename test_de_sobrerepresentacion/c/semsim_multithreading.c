#include <R.h>
#include <Rdefines.h>
#include <pthread.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

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
int min(int a, int b){
	if(a < b) return a;
	return b;
}
struct s_argumentos {
	SEXP anotaciones; 
	int len_nombres;	
	int len_anotaciones;
	double *p_semsimM;
	double *ic;
	char **nombres;
};

struct s{
	struct s_argumentos *argumentos;
	int desde;
	int hasta;
};

/* this function is run by the second thread */
void *mica_multithread(void *v_argumentos)
{
	struct s *argumentos = v_argumentos;
	int categoria1, categoria2, k, l, m, len, len2, iic;
	double mica;
	//Categoria por categoria vamos buscando los padres en común
	for(categoria1 = argumentos->desde; categoria1 < argumentos->hasta; ++categoria1){
		
		//Guardo la primer categoría en un vector porque vamos a estar accediendo mucho
		len = Rf_length(VECTOR_ELT(argumentos->argumentos->anotaciones, categoria1));			
		char *padres[len];		
		for(k = 0; k < len; ++k){
			padres[k] = R_alloc(strlen(CHAR(STRING_ELT(VECTOR_ELT(argumentos->argumentos->anotaciones, categoria1), k))), sizeof(char));
			strcpy(padres[k], CHAR(STRING_ELT(VECTOR_ELT(argumentos->argumentos->anotaciones, categoria1), k))); 
		}

		//Recorro la segunda categoria
		for(categoria2 = categoria1; categoria2 < argumentos->argumentos->len_anotaciones; ++categoria2){
			len2 = Rf_length(VECTOR_ELT(argumentos->argumentos->anotaciones, categoria2));	
			mica = 0;

			//Vamos buscando matches entre las dos categorías		
			for(l = 0; l < len; ++l){
				for(m = 0; m < len2; ++m){
					if(strcmp(padres[l], CHAR(STRING_ELT(VECTOR_ELT(argumentos->argumentos->anotaciones, categoria2), m))) == 0){
						iic = getListElement(argumentos->argumentos->nombres, padres[l], argumentos->argumentos->len_nombres);
						if(mica < argumentos->argumentos->ic[iic]) mica = argumentos->argumentos->ic[iic];
					}
				}
			}
			printf("%d, %d, %g\n", categoria1, categoria2, mica);
			argumentos->argumentos->p_semsimM[categoria1*argumentos->argumentos->len_anotaciones + categoria2] = mica;
			argumentos->argumentos->p_semsimM[categoria2*argumentos->argumentos->len_anotaciones + categoria1] = mica;
		}
	}
	return NULL;

}

SEXP semsim(SEXP anotaciones, SEXP sic){
	int categoria1, k, len_anotaciones, len_nombres, len;
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
		//Guardo la primer categoría en un vector porque vamos a estar accediendo mucho
		len = Rf_length(VECTOR_ELT(anotaciones, categoria1));			
		char *padres[len];		
		for(k = 0; k < len; ++k){
			padres[k] = R_alloc(strlen(CHAR(STRING_ELT(VECTOR_ELT(anotaciones, categoria1), k))), sizeof(char));
			strcpy(padres[k], CHAR(STRING_ELT(VECTOR_ELT(anotaciones, categoria1), k))); 
		}
		p_anotaciones[categoria1] = padres;
	}


	struct s_argumentos argumentos;
	argumentos.anotaciones = anotaciones; 
	argumentos.len_nombres = len_nombres;	
	argumentos.len_anotaciones = len_anotaciones;
	argumentos.nombres = nombres;
	argumentos.p_semsimM = p_semsimM;
	argumentos.ic = ic;
	
	int i_threads = 10;
	struct s args[i_threads];
	int carga = ceil((float)len_anotaciones/(float)i_threads);

	args[0].argumentos = &argumentos;
	args[0].desde = 0;
	args[0].hasta = carga;
	printf("%d, %d\n", args[0].desde, args[0].hasta);
	for(k = 1; k < i_threads; k++){
		args[k].argumentos = &argumentos;
		args[k].desde = args[k-1].hasta;
		args[k].hasta = min(args[k-1].hasta + carga, len_anotaciones);
		printf("%d, %d\n", args[k].desde, args[k].hasta);
	}

	pthread_t thread[i_threads];
	
	for(k = 0; k < i_threads; k++){
		pthread_create(&thread[k], NULL, mica_multithread, &args[k]);
	}
	for(k = 0; k < i_threads; k++){
		pthread_join(thread[k], NULL);
	}

	PROTECT(dimnames = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(dimnames, 0, names_categorias);
	SET_VECTOR_ELT(dimnames, 1, names_categorias);
	dimnamesgets(semsimM, dimnames);
	UNPROTECT(2);
	return(semsimM);
}
