#include <stdio.h>
# include<R.h>
# include<Rmath.h>
# include<stdio.h>
# include<stdlib.h>



double  *read_variables(FILE *fi, int NB_VARIABLES,int NB_INDIVIDUALS)
{
	double *data=calloc(NB_VARIABLES*NB_INDIVIDUALS, sizeof(double));
	int x=0;
	while ((!feof(fi)) && (!ferror(fi)) && x <NB_VARIABLES){
		int i;
		if (x<NB_VARIABLES){
			if( x>NB_VARIABLES )
				printf("error : too many genes : %i - %i\n", x, NB_VARIABLES );
			for (i=0;i<NB_INDIVIDUALS;i++)
				fscanf (fi, "%le\t", &data[x*NB_VARIABLES+ i]);

		}
		x++;

	}
	printf("x = %d\t %d\t %d\n",x,NB_VARIABLES,NB_INDIVIDUALS);
	return(data);

};
// fonction updb.lasso du script de David
double updb(int *valr, int *vals, double *matYZ, double *matZZ, int *ncolZZ, int *nrowYZ, double *mathb, double *vectPsi, double *lambda){
    int i, j, k, r, s, nbf, m;
    r = *valr;
    s = *vals;
    nbf = *ncolZZ;
    m = *nrowYZ;

    double solm, solp, sol;
  //  solm = 0.0;//l'initialisation de cette variable est inutile
    //solp = 0.0;// ici aussi
    sol = 0.0;
    double reel = *lambda;
    // calcul intermediaire de sum(hb[r,-s]*Czz[-s,s])
    double somA = 0.0;
    for (i=0;i<nbf;i++) {
        if (i!=(s)) {
			somA += mathb[(i*m)+(r)] * matZZ[((s)*nbf)+i];
        }
    }

    solm = (matYZ[(s)*m+(r)] - (reel * vectPsi[r]/2) - somA) / matZZ[(s)*nbf+(s)];

    solp = (matYZ[(s)*m+(r)] + (reel * vectPsi[r]/2) - somA) / matZZ[(s)*nbf+(s)];

    double *hb0= (double *) calloc((nbf*m),sizeof(double));
    if(hb0==NULL){
		printf("Erreur allocation memoire\n");
		exit(1);
    }
    // creation de la matrice hb0= a hb
    for (j=0;j<m;j++){
        for (k=0;k<nbf;k++){
			hb0[k*m+j] = mathb[k*m+j];
        }
    }
    // initialisation de hb0[r,s]
    hb0[(s)*m+(r)] = solp;
    // calcul intermediaire de sum(hb0[r,]*Czz[,s])
    double somB = 0.0;
    for (i=0;i<nbf;i++) {
        somB += hb0[(i*m)+(r)] * matZZ[((s)*nbf)+i];
    }

    // creation d'une variation signsolp correspondant a la fonction sign(solp) de R
    double dfp = 0.0;
    int signsolp = 0;
    if (solp>0){
        signsolp = 1;
	}
	else if (solp==0){
		signsolp = 0;
	}
	else {signsolp = -1;
	}


    dfp = -2*(matYZ[(s)*m+(r)]/vectPsi[r]) + 2*(somB/vectPsi[r]) + reel*signsolp;

    // nouvelle initialisation de hb0[r,s]
    hb0[(s)*m+(r)] = solm;
    // calcul intermediaire de sum(hb0[r,]*Czz[,s]) avec une nouvelle initialisation de hb0
    double somC = 0.0;
    for (i=0;i<nbf;i++) {
        somC += hb0[(i*m)+(r)] * matZZ[((s)*nbf)+i];
    }

    // creation d'une variation signsolm correspondant a la fonction sign(solm) de R
    double dfm = 0.0;
    int signsolm = 0;
    if (solm>0){
        signsolm = 1;
	}
	else if (solm==0){
		signsolm = 0;
	}
	else {signsolm = -1;}

    dfm = -2*(matYZ[(s)*m+(r)]/vectPsi[r]) + 2*(somC/vectPsi[r]) + reel * signsolm;


    //correspond a la fonction sol0 du script de David
    hb0[(s)*m+(r)]=0;
    // calcul intermediaire de sum(hb0[r,]*Czz[,s]) avec une nouvelle initialisation de hb0
    double somD = 0.0;
    for (i=0;i<nbf;i++) {
        somD += hb0[(i*m)+(r)] * matZZ[((s)*nbf)+i];
    }
    double resp = 0.0;
    double resm = 0.0;
    int sign =0;
    resp = -2*(matYZ[(s)*m+(r)]/vectPsi[r]) + 2*(somD/vectPsi[r]) + reel;
    resm = -2*(matYZ[(s)*m+(r)]/vectPsi[r]) + 2*(somD/vectPsi[r]) - reel;

    if ((resp>=0) & (resm>=0)) {sign = 0;}
    else if ((resp>=0) & (resm<0)) {sign = 1;}
	else if ((resp<0) & (resm<0)) {sign = 0;}
	else sign = 1;


    if (sign==1) {sol = 0;}
	else if (fabs(dfm)<(10e-9))
	{sol = solm;}
	else if (fabs(dfp)<(10e-9)) {sol = solp;}

    free(hb0);
    return(sol);

}

// fonction sql du script de David
// cette fonction retourne un vecteur. Le vecteur est recupere du produit d'une matrice et de sa transposee,
// puis de la selection du triangle inferieur.
double *fctsql(double *matM, int *nrowM, int *ncolM){
    int i, j, nrow, ncol, k;
    nrow = *nrowM;
    ncol = *ncolM;

	double *M = calloc(nrow*nrow,sizeof(double));
	for (i=0;i<nrow;i++){
		for (j=0;j<nrow;j++){
			M[j*nrow+i]=0;
		}
	}
    double *tmatM = calloc(ncol*nrow,sizeof(double));
	for (i=0;i<nrow;i++){
		for (j=0;j<ncol;j++){
			tmatM[i*ncol+j]= matM[j*nrow +i];

		}
	}

    double s;
    for(i=0;i<nrow;i++){
        for(k=0;k<nrow;k++){
			s = 0;
            for(j=0;j<ncol;j++){
				s = s + matM[i+(nrow*j)] * tmatM[j+(ncol*k)];
				M[i+(nrow*k)] = s;
            }
        }
    }
    free(tmatM);

    double *Vres = calloc((nrow*(nrow-1)/2),sizeof(double));
	for(k=0;k<(nrow*(nrow-1)/2);k++){
		Vres[k]=0;
	}

    k = 0;
    for (j=0;j<nrow;j++){
        for (i=j+1;i<nrow;i++){
            Vres[k]=M[j*nrow+i];
            k += 1;
        }

    }
	free(M);//je remet les free au bon endroit
   // free(Vres);
    return(Vres);

	//ici tu sortais de la fonction avec le return avant de liberer la memoire
    //free(M);
    //free(Vres);
}

// fonction appelee diag qui calcule le vecteur hPsi de la fonction cycliccd du script de David dans R
double *diag(double *matS, double *mathb, double *matCyz, int *nrowS, int *ncolCyz){
    int m, nbf, i, j,k;
    m = *nrowS;
    nbf = *ncolCyz;

	// matrice provisoire pour stocker le produit matriciel
	double *M = calloc(m*m,sizeof(double));
	for (i=0;i<m;i++){
		for (j=0;j<m;j++){
			M[j*m+i]=0;
		}
	}
    // matrice provisoire pour stocker la transpose
    double *tmatCyz = calloc(m*nbf,sizeof(double));

	// ok si nbf<m
	for (i=0;i<m;i++){
		for (j=0;j<nbf;j++){
			tmatCyz[i*nbf+j]= matCyz[j*m +i];

		}
	}
    // matrice provisoire pour stocker la soustraction des matrices
    double *diff = calloc(m*m,sizeof(double));
	for (i=0;i<m;i++){
		for (j=0;j<m;j++){
			diff[j*m+i]=0;
		}
	}

    // produit de la matrice hb et du transpose de Cyz
    double s;
    for(i=0;i<m;i++){
		for(k=0;k<m;k++){
			s = 0;
            for(j=0;j<nbf;j++){
				s = s + mathb[i+(m*j)] * tmatCyz[j+(nbf*k)];
				M[i+(m*k)] = s;
            }
        }
	}


    // soustraction S-hb*t(Cyz)
    for (i=0;i<m;i++){
        for (j=0;j<m;j++){
            diff[j*m+i] = matS[j*m+i] - M[j*m+i];
        }
    }
    double *vectres=calloc(m,sizeof(double));

    for (k=0;k<m;k++){
        vectres[k] = diff[k*m+k];
    }

    // libere la memoire
    free(tmatCyz);
    free(M);
    free(diff);
	//free(vectres);//Je remet ici aussi
    return(vectres);

}


void cycliccd(double *matS, double *matYZ, double *matZZ, int *nrowYZ, int *ncolZZ, double *vectPsi, int *matzeros, double *lambda, double *minerr, int *maxiter, double *resPsi, double *resB){

    int r, s, i,j,k,m,nbf, maxit, stop, iter;
    m = *nrowYZ;
    nbf = *ncolZZ;
    maxit = *maxiter;
    stop = 0;
    iter = 0;

    double resupdb = 0.0;
    double miner = *minerr;
    double mean = 0.0;


	double *hb = calloc(m*nbf,sizeof(double));
	double *hB = calloc(m*nbf,sizeof(double));
	for (i=0;i<m;i++){
		for (j=0;j<nbf;j++){
			hb[j*m+i]=0;
			hB[j*m+i]=0;
		}
	}
    double *hPsi = calloc(m,sizeof(double));
	for (k=0;k<m;k++){
		hPsi[k] = vectPsi[k];
	}
   // double *temp ;//= calloc(m,sizeof(double));
    double *sqlhb = calloc(m*(m-1)/2,sizeof(double));
  //  double *sqlhbtemp;// =calloc(m*(m-1)/2,sizeof(double));//Attention tu declare 2 fois la meme variable
    double *sqlhB = calloc(m*(m-1)/2,sizeof(double));
   // double *sqlhBtemp = calloc(m*(m-1)/2,sizeof(double));//Attention tu declare 2 fois la meme variable
    double *meantemp = calloc(m*(m-1)/2,sizeof(double));
    double meansom = 0.0;

    while(stop!=1){

        iter += 1;
		for (i=0;i<m;i++){
			for (j=0;j<nbf;j++){
				hB[j*m+i] = hb[j*m+i];
			}
		}
		for (r=0;r<m;r++){
			for (s=0;s<nbf;s++){
				if (matzeros[s*m+r]==1) {resupdb = 0.0;}
				else {resupdb = updb(&r,&s,matYZ,matZZ,&nbf,&m,hb,hPsi,lambda);}

				hb[s*m+r] = resupdb;
			}
		}

        // calcul de hPsi (avec la fonction diag)
        double *temp = diag(matS, hb, matYZ, &m, &nbf);

		for (k=0;k<m;k++){
			hPsi[k] = temp[k];
		}

        // calcul de mean(abs(sql(hb)-sql(hB)))
        double *sqlhbtemp = fctsql(hb, &m, &nbf);

		for (k=0;k<(m*(m-1)/2);k++){
			sqlhb[k] = sqlhbtemp[k];
		}

         double *sqlhBtemp = fctsql(hB, &m, &nbf);

		for (k=0;k<(m*(m-1)/2);k++){
			sqlhB[k] = sqlhBtemp[k];
		}

        meansom = 0.0;
        for (k=0;k<m*(m-1)/2;k++) {
            meantemp[k] = fabs(sqlhb[k] - sqlhB[k]);
            meansom = meansom + meantemp[k];
        }

        mean = meansom / (m*(m-1)/2);
		//Rprintf("meansom %f\n",meansom);
		//Rprintf("mean %f\n",mean);


        if ((iter>maxit)|(mean<miner)) {
            stop = 1;}
        else stop = 0;


		free(temp);
		free(sqlhbtemp);
		free(sqlhBtemp);
	}

	for (k=0;k<m;k++){
		resPsi[k] = hPsi[k];
	}
	for (i=0;i<m;i++){
		for (j=0;j<nbf;j++){
            resB[j*m+i] = hb[j*m+i];
		}
	}

    // libere la memoire
    free(hb);
    free(hB);
    free(hPsi);
   // free(temp);
    free(sqlhb);
    //free(sqlhbtemp);
    free(sqlhB);
   // free(sqlhBtemp);
    free(meantemp);
}
