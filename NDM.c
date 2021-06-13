//PROGRAMA MEDIDA POLARIZACION FOTON CON MATRICES DE PAULI

#include <stdio.h>
#include "complex.h"
#include "gsl_rng.h"
 
#define tam 2
#define TAM 4

#define error 0.0001

#define omega 1
#define T 3.142	//Depende de omega

#define fotones 10000

#define deltat 0.0001
#define pasos 10

#define N 50000 



////////////////////////- FUNCIONES VARIABLES COMPLEJAS -//////////////////////////
///////////////////////////////////////////////////////////////////////////////////
fcomplex PdtoEscalar(fcomplex v1[tam], fcomplex v2[tam]);
int Bra (fcomplex ket[tam], fcomplex bra[tam]);
int CProyector (fcomplex ket[tam], fcomplex bra[tam], fcomplex proyector[tam][tam]);
int ActOperador (fcomplex Operador[tam][tam], fcomplex ket[tam], fcomplex resultado[tam]);
fcomplex ValorEsp (fcomplex Operador[tam][tam], fcomplex bra[tam], fcomplex ket[tam]);
int Vadd(fcomplex V1[tam], fcomplex V2[tam], fcomplex SUMA[tam]);
int Vsub(fcomplex V1[tam], fcomplex V2[tam], fcomplex RESTA[tam]);
int Madd(fcomplex M1[tam][tam], fcomplex M2[tam][tam], fcomplex SUMA[tam][tam]);
int Msub(fcomplex M1[tam][tam], fcomplex M2[tam][tam], fcomplex RESTA[tam][tam]);
int Mmul (fcomplex M1[tam][tam], fcomplex M2[tam][tam], fcomplex PDTO[tam][tam]);

fcomplex PdtoEscalarSist(fcomplex v1[TAM], fcomplex v2[TAM]);
int BraSist (fcomplex ket[TAM], fcomplex bra[TAM]);
int CProyectorSist (fcomplex ket[TAM], fcomplex bra[TAM], fcomplex proyector[TAM][TAM]);
int ActOperadorSist (fcomplex Operador[TAM][TAM], fcomplex ket[TAM], fcomplex resultado[TAM]);
fcomplex ValorEspSist (fcomplex Operador[TAM][TAM], fcomplex bra[TAM], fcomplex ket[TAM]);
int VaddSist(fcomplex V1[TAM], fcomplex V2[TAM], fcomplex SUMA[TAM]);
int VsubSist(fcomplex V1[TAM], fcomplex V2[TAM], fcomplex RESTA[TAM]);
int VPdtoExt(fcomplex V1[tam], fcomplex V2[tam], fcomplex PdtoExt[TAM]);
int MaddSist(fcomplex M1[TAM][TAM], fcomplex M2[TAM][TAM], fcomplex SUMA[TAM][TAM]);
int MsubSist(fcomplex M1[TAM][TAM], fcomplex M2[TAM][TAM], fcomplex RESTA[TAM][TAM]);
int MmulSist (fcomplex M1[TAM][TAM], fcomplex M2[TAM][TAM], fcomplex PDTO[TAM][TAM]);
int MPdtoExt(fcomplex M1[tam][tam], fcomplex M2[tam][tam], fcomplex PDTOEXT[TAM][TAM]);
//////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////- FUNCIONES PROGRAMA -/////////////////////////////////
//Usar solo si no es la identidad
int OpSist(int eje1, int eje2, fcomplex Px[tam][tam], fcomplex Py[tam][tam], fcomplex Pz[tam][tam], fcomplex I[tam][tam], fcomplex OpSist[TAM][TAM]);

gsl_rng *t;

int main ()
{
	//VARIABLES DE FOTONES, INDEPENDIENTES
	fcomplex foton1[tam], foton2[tam]; //Est inicial, |0> |1>
	fcomplex Px[tam][tam], Py[tam][tam], Pz[tam][tam], I[tam][tam]; //Polariz respecto eje, matriz Pauli; identidad
	
	
	fcomplex x_up[tam],x_down[tam]; //Autovector (+1) y (-1) de Px
	fcomplex y_up[tam],y_down[tam];
	fcomplex z_up[tam],z_down[tam];	

	
	//VARIABLES DE SISTEMA, CONJUNTO
	fcomplex sistema[TAM];
	fcomplex ObservableX1[TAM][TAM], ObservableZ1[TAM][TAM], ObservableX2[TAM][TAM], ObservableZ2[TAM][TAM];
	
	//Fotones entrelazados
	fcomplex psimas[TAM], psimenos[TAM], phimas[TAM], phimenos[TAM];
	
	fcomplex entrelazado[TAM];
	fcomplex separable[TAM];
	
	fcomplex aux1[TAM], aux2[TAM];
	fcomplex Maux1[TAM][TAM],Maux2[TAM][TAM];
	
	//OTRAS
	
	
	int i,j,k,l, control;	//Bucles, llamar a funciones que devuelven entero que no interesa
	
	int eje1,eje2;	//Medida sobre 1er foton, 2o foton
	int oscilacion;
	int num_pasos;
	
	double norma;	//Renormalizar estados
	
	fcomplex pdto, valorX1, valorZ1, valorX2, valorZ2; //PdtoEscalar para norma, ValorEsp en medida	
	
	double promedioX1[N], promedioZ1[N], promedioX2[N], promedioZ2[N];
	
	double prob1, prob2;
	double aleatorio; //Numero sorteo
	
	int up, down; //Autovalores de operadores fotones
	fcomplex X_up[tam][tam], X_down[tam][tam]; 	//Proyectores Px


	fcomplex X1P1[TAM][TAM], X1P2[TAM][TAM];	//Proyectores asociados a autoestados X1
		
	
	double time;
	fcomplex H[TAM][TAM];	//Hamiltoniano del sistema
	fcomplex Hsist[TAM];	//Hsist= H|sist>
	fcomplex final[TAM];
	
	//VARIABLES NUM ALEATORIOS
	int seed=5164982;
	extern gsl_rng *t;
	
	
	//FICHEROS
	FILE *f1;
	
	f1=fopen("NDM.txt","w"); //Nombre del fichero donde guarda los datos
	
	
	//INICIALIZAR GENERADOR ALEATORIO
            t=gsl_rng_alloc(gsl_rng_taus); //Inicializar num aleat y tipo de generador
            gsl_rng_set(t,seed); //Pasar semilla para iniciar generador
	
//-------------------------------------------------------------------------------------------
//EST INICIAL FOTON 1, |0>
//EST INICIAL FOTON 2, |1>
	
	//OPCION PROGRAMA
	foton1[0]=Complex(1.,0.);
	foton1[1]=Complex(0.,0.);

	foton2[0]=Complex(0.,0.);
	foton2[1]=Complex(1.,0.);

	
//FOTONES ENTRELAZADOS. BASE DE BELL
	//Phi+=1/sqrt(2) * (00+11)
	control=VPdtoExt(foton1, foton1, aux1);
	control=VPdtoExt(foton2, foton2, aux2);
	control=VaddSist(aux1, aux2, phimas);
	
			
	pdto=PdtoEscalarSist(phimas, phimas);
	norma=Cabs(pdto);
	if(norma!=1)
		for(i=0;i<TAM;i++)
			phimas[i]=RCmul((1/sqrt(norma)), phimas[i]);
			
	//Phi-=1/sqrt(2) * (00-11)
	control=VPdtoExt(foton1, foton1, aux1);
	control=VPdtoExt(foton2, foton2, aux2);
	control=VsubSist(aux1, aux2, phimenos);
			
	pdto=PdtoEscalarSist(phimenos, phimenos);
	norma=Cabs(pdto);
	if(norma!=1)
		for(i=0;i<TAM;i++)
			phimenos[i]=RCmul((1/sqrt(norma)), phimenos[i]);
			
	//Psi+=1/sqrt(2) * (01+10)
	control=VPdtoExt(foton1, foton2, aux1);
	control=VPdtoExt(foton2, foton1, aux2);
	control=VaddSist(aux1, aux2, psimas);
	
			
	pdto=PdtoEscalarSist(psimas, psimas);
	norma=Cabs(pdto);
	if(norma!=1)
		for(i=0;i<TAM;i++)
			psimas[i]=RCmul((1/sqrt(norma)), psimas[i]);
			
	//Psi-=1/sqrt(2) * (01-10)
	control=VPdtoExt(foton1, foton2, aux1);
	control=VPdtoExt(foton2, foton1, aux2);
	control=VsubSist(aux1, aux2, psimenos);
			
	pdto=PdtoEscalarSist(psimenos, psimenos);
	norma=Cabs(pdto);
	if(norma!=1)
		for(i=0;i<TAM;i++)
			psimenos[i]=RCmul((1/sqrt(norma)), psimenos[i]);						
		
		
		
			
//FOTONES entrelazado=1/sqrt(3) * (00-10-11)
	control=VPdtoExt(foton1, foton1, aux1);
	control=VPdtoExt(foton2, foton1, aux2);
	control=VsubSist(aux1, aux2, entrelazado);
	control=VPdtoExt(foton2, foton2, aux1);
	control=VsubSist(entrelazado, aux1, entrelazado);
			
	pdto=PdtoEscalarSist(entrelazado, entrelazado);
	norma=Cabs(pdto);
	if(norma!=1)
		for(i=0;i<TAM;i++)
			entrelazado[i]=RCmul((1/sqrt(norma)), entrelazado[i]);
			
// FOTONES separable=1/2 * (00+01+10+11)= (0>+1>)1 x (0>+1>)2
	control=VaddSist(psimas,phimas,separable);
	
	pdto=PdtoEscalarSist(separable, separable);
	norma=Cabs(pdto);
	if(norma!=1)
		for(i=0;i<TAM;i++)
			separable[i]=RCmul((1/sqrt(norma)), separable[i]);								
					
//--------------------------------------------------------------------------------------------
//OBSERVABLES, AUTOVECTORES Y PROYECTORES
			
	//PAULI X, 1
	Px[0][0]=Complex(0.,0.);
	Px[0][1]=Complex(1.,0.);
	Px[1][0]=Complex(1.,0.);
	Px[1][1]=Complex(0.,0.);
	
		//AUTOVECTORES
		x_up[0]=Complex(1./sqrt(2.),0.);
		x_up[1]=Complex(1./sqrt(2.),0.);
		x_down[0]=Complex(-1./sqrt(2.),0.);
		x_down[1]=Complex(1./sqrt(2.),0.);
		
	
	//PAULI Y, 2
	Py[0][0]=Complex(0.,0.);
	Py[0][1]=Complex(0.,-1.);
	Py[1][0]=Complex(0.,1.);
	Py[1][1]=Complex(0.,0.);
	
		//AUTOVECTORES
		y_up[0]=Complex(0.,-1./sqrt(2.));
		y_up[1]=Complex(1./sqrt(2.),0.);
		y_down[0]=Complex(0.,1./sqrt(2.));
		y_down[1]=Complex(1./sqrt(2.),0.);
		
	
	//PAULI Z, 3
	Pz[0][0]=Complex(1.,0.);
	Pz[0][1]=Complex(0.,0.);
	Pz[1][0]=Complex(0.,0.);
	Pz[1][1]=Complex(-1.,0.);
	
	
		//AUTOVECTORES
		z_up[0]=Complex(1.,0.);
		z_up[1]=Complex(0.,0.);
		z_down[0]=Complex(0.,0.);
		z_down[1]=Complex(1.,0.);
		
	up=1;
	down=-1;	
	
	//IDENTIDAD I, 0
	//AUTOVALOR 1, TODOS VECTORES SON AUTOESTADOS PORQUE VAN A ELLOS MISMOS
	I[0][0]=Complex(1.,0.);
	I[0][1]=Complex(0.,0.);
	I[1][0]=Complex(0.,0.);
	I[1][1]=Complex(1.,0.);
	
	
//-------------------------------------------------------------------------------------------
//Inicializacion de promedios

	
	for(i=0;i<N;i++)
	{
		promedioX1[i]=0.;
		promedioZ1[i]=0.;
		promedioX2[i]=0.;
		promedioZ2[i]=0.;
	}

			
//-----------------------------------------------------------------------------------------------	
//VARIABLES CONSTANTES DE LOS BUCLES
	
//OBSERVABLES------------------------------------------------------------------------------------	
	//eje1, eje2
	control=OpSist(1, 0, Px, Py, Pz, I, ObservableX1);	//Definir Operador del sistema, OpSist

	control=CProyector(x_up,x_up,X_up);
	control=CProyector(x_down,x_down,X_down);
	control=MPdtoExt(X_up,I,X1P1);
	control=MPdtoExt(X_down,I,X1P2);	

	
	control=OpSist(3, 0, Px, Py, Pz, I, ObservableZ1);	//Definir Operador del sistema, OpSist
	control=OpSist(0, 1, Px, Py, Pz, I, ObservableX2);	//Definir Operador del sistema, OpSist
	control=OpSist(0, 3, Px, Py, Pz, I, ObservableZ2);	//Definir Operador del sistema, OpSist
	
	
//DEFINE EL ESTADO INICIAL DEL HAMILTONIANO------------------------------------------------------	
	//HAMILTONIANO ENTRELAZADO Z1 X Z2 	
//	control=MPdtoExt(Pz,Pz,H); 
//		for(i=0;i<TAM;i++)
//			for(j=0;j<TAM;j++)
//				H[i][j]=RCmul(omega,H[i][j]);
				
	//HAMILTONIANO SEPARABLE Z1 + Z2 	
	control=MPdtoExt(Pz,I,Maux1);
	control=MPdtoExt(I,Pz,Maux2);
	control=MaddSist(Maux1,Maux2,H); 
		for(i=0;i<TAM;i++)
			for(j=0;j<TAM;j++)
				H[i][j]=RCmul(omega,H[i][j]);
							
//---------------------------------------------------------------------------------------------
//EV. TEMPORAL + MEDIDA
	num_pasos=T/(pasos*deltat);

for(l=0;l<fotones;l++)
{
	printf("%i\n",l);
	time=0.;
	
	for(i=0;i<TAM;i++)
		sistema[i]=separable[i];		//DEFINE EL ESTADO INICIAL DEL SISTEMA
	
	
	oscilacion=1;
	
	
	for(k=0;k<N;k++) //En cada ciclo se avanza 0.001s
	{		
	   	valorX1=ValorEspSist (ObservableX1, sistema, sistema);	//Calculo valor esperado del observable, valor
		valorZ1=ValorEspSist (ObservableZ1, sistema, sistema);
		valorX2=ValorEspSist (ObservableX2, sistema, sistema);
		valorZ2=ValorEspSist (ObservableZ2, sistema, sistema);
		
		promedioX1[k]+=valorX1.r;
		promedioZ1[k]+=valorZ1.r;
		promedioX2[k]+=valorX2.r;
		promedioZ2[k]+=valorZ2.r;
					
	   	if(k==(oscilacion*num_pasos/4)+1000) //ELIGE TIEMPO MEDIDA DEL SISTEMA. oscilacion*num_pasos=T
	   	{
			
			prob1=Cabs(ValorEspSist (X1P1, sistema, sistema));
			prob2=Cabs(ValorEspSist (X1P2, sistema, sistema));
		
		
		//COLAPSO
			aleatorio=gsl_rng_uniform(t);
		
			if(aleatorio<=prob1) //Colapso a autoestado1 dado por el proyector X1P1
			{
				control=ActOperadorSist (X1P1, sistema, sistema);
				for(i=0;i<TAM;i++)
					sistema[i]=RCmul((1/sqrt(prob1)), sistema[i]);		
			}
		
			else  //Colapso a autoestado2
			{
				control=ActOperadorSist (X1P2, sistema, sistema);
				for(i=0;i<TAM;i++)
					sistema[i]=RCmul((1/sqrt(prob2)), sistema[i]);		
			}
					
			oscilacion++;

		}
		
		pdto=PdtoEscalarSist(sistema, sistema);
		norma=Cabs(pdto);
				
		//EV. TEMPORAL
		for(j=0;j<pasos;j++)		
		{
			control=ActOperadorSist(H,sistema,Hsist);
			for(i=0;i<TAM;i++)
				Hsist[i]=Cmul(Complex(0.,deltat),Hsist[i]);
			
			control=VsubSist(sistema, Hsist, final);
			
			pdto=PdtoEscalarSist(final, final);
			norma=Cabs(pdto);
			for(i=0;i<TAM;i++)
				sistema[i]=RCmul((1/sqrt(norma)), final[i]);	
		}
		
		time=time+pasos*deltat;
		
	}			
			
}

	time=0.;
	for(k=0;k<N;k++)
	{
		time=time+pasos*deltat;
		promedioX1[k]=promedioX1[k]/fotones;
		promedioZ1[k]=promedioZ1[k]/fotones;
		promedioX2[k]=promedioX2[k]/fotones;
		promedioZ2[k]=promedioZ2[k]/fotones;
		fprintf(f1,"%lf\t%lf\t%lf\t%lf\t%lf\n", time, promedioX1[k], promedioZ1[k], promedioX2[k],promedioZ2[k]);	//Tiempo, medida, valorX1, valorZ1, valorX2, valorZ2
	}
	
	fclose(f1);
		
	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////- FUNCIONES PROGRAMA -////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

int OpSist(int eje1, int eje2, fcomplex Px[tam][tam], fcomplex Py[tam][tam], fcomplex Pz[tam][tam], fcomplex I[tam][tam], fcomplex OpSist[TAM][TAM])
//Op sobre sistema de dos Op sobre fotones individuales
{
	int control;
	
	if(eje1==0)
	{
		if(eje2==0) control=MPdtoExt(I, I, OpSist);
		else if(eje2==1) control=MPdtoExt(I, Px, OpSist);
		else if(eje2==2) control=MPdtoExt(I, Py, OpSist);
		else control=MPdtoExt(I, Pz, OpSist);
	}
	else if(eje1==1)
	{
		if(eje2==0) control=MPdtoExt(Px, I, OpSist);
		else if(eje2==1) control=MPdtoExt(Px, Px, OpSist);
		else if(eje2==2) control=MPdtoExt(Px, Py, OpSist);
		else control=MPdtoExt(Px, Pz, OpSist);
	}
	else if(eje1==2)
	{
		if(eje2==0) control=MPdtoExt(Py, I, OpSist);
		else if(eje2==1) control=MPdtoExt(Py, Px, OpSist);
		else if(eje2==2) control=MPdtoExt(Py, Py, OpSist);
		else control=MPdtoExt(Py, Pz, OpSist);
	}
	else
	{
		if(eje2==0) control=MPdtoExt(Pz, I, OpSist);
		else if(eje2==1) control=MPdtoExt(Pz, Px, OpSist);
		else if(eje2==2) control=MPdtoExt(Pz, Py, OpSist);
		else control=MPdtoExt(Pz, Pz, OpSist);
	}	
	
	return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////- FUNCIONES VARIABLES COMPLEJAS -//////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
fcomplex PdtoEscalar(fcomplex v1[tam], fcomplex v2[tam])
//Devuelve el producto escalar de dos vectores
{
	int i;
	fcomplex aux;
	
	aux=Complex(0.,0.);
	for(i=0;i<tam;i++)
		aux=Cadd(aux,Cmul(Conjg(v1[i]),v2[i]));
	
	return aux;	
}


int Bra (fcomplex ket[tam], fcomplex bra[tam])
//Devuelve el vector conjugado del ket
{
	int i;
	
	for(i=0;i<tam;i++)
		bra[i]=Conjg(ket[i]);
		
	return 0;		
}


int CProyector (fcomplex ket[tam], fcomplex bra[tam], fcomplex proyector[tam][tam])	
//Devuelve matriz proyector de dos vectores, bra y ket
{
	
	int i,j;
	
	int braaux;
	fcomplex aux[tam];
	
	braaux=Bra(bra,aux);
	
	for(i=0;i<tam;i++)
		for(j=0;j<tam;j++)
		{	
			
			proyector[i][j]=Cmul(ket[i],aux[j]);
		}
		
	return	0;
}

int ActOperador (fcomplex Operador[tam][tam], fcomplex ket[tam], fcomplex resultado[tam])
//Devuelve el vector del producto de una matriz y otro vector
{
	int i,j;
	fcomplex aux;
	
	for(i=0;i<tam;i++)
	{
		aux=Complex(0.,0.);
		for(j=0;j<tam;j++)
		{
			aux=Cadd(aux,Cmul(Operador[i][j],ket[j]));
		}
		resultado[i]=aux;
	}
	
	return 0;
}

fcomplex ValorEsp (fcomplex Operador[tam][tam], fcomplex ket[tam], fcomplex bra[tam])
//Devuelve el escalar de vector-matriz-vector
{
	fcomplex valor, pto;

	int control;	
	double norma;
	
	fcomplex dcha[tam];
	
	control=ActOperador (Operador, ket, dcha);
	
	valor=PdtoEscalar(bra, dcha); //El conjugado del bra se calcula dentro de PdtoEscalar
	
	pto=PdtoEscalar(bra,ket);
	norma=Cabs(pto);
	
	valor=RCmul(1./norma,valor);
	
	return valor;
}

////////////////////////////////////////////////////////////////////////////////////////////
//OPERACIONES CON VECTORES

int Vadd(fcomplex V1[tam], fcomplex V2[tam], fcomplex SUMA[tam])
//Devuelve V1+V2
{
	int j;
	
	for(j=0;j<tam;j++)
		SUMA[j]=Cadd(V1[j],V2[j]);
			
	return 0;
}

int Vsub(fcomplex V1[tam], fcomplex V2[tam], fcomplex RESTA[tam])
//Devuelve V1-V2
{
	int j;
	
	for(j=0;j<tam;j++)
		RESTA[j]=Csub(V1[j],V2[j]);
			
	return 0;
}


//OPERACIONES CON MATRICES

int Madd(fcomplex M1[tam][tam], fcomplex M2[tam][tam], fcomplex SUMA[tam][tam])
//Devuelve M1+M2
{
	int i,j;
	
	for(i=0;i<tam;i++)
		for(j=0;j<tam;j++)
			SUMA[i][j]=Cadd(M1[i][j],M2[i][j]);
			
	return 0;
}

int Msub(fcomplex M1[tam][tam], fcomplex M2[tam][tam], fcomplex RESTA[tam][tam])
//Devuelve M1-M2
{
	int i,j;
	
	for(i=0;i<tam;i++)
		for(j=0;j<tam;j++)
			RESTA[i][j]=Csub(M1[i][j],M2[i][j]);
			
	return 0;
}

int Mmul (fcomplex M1[tam][tam], fcomplex M2[tam][tam], fcomplex PDTO[tam][tam])
//Devuelve M1*M2
{
	int i,j,k;
	
	
	
	for(i=0;i<tam;i++)
		for(j=0;j<tam;j++)
		{
			PDTO[i][j]=Complex(0.,0.);
			for(k=0;k<tam;k++)
				PDTO[i][j]=Cadd(PDTO[i][j], Cmul(M1[i][k],M2[k][j]));
		}	
	return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////
fcomplex PdtoEscalarSist(fcomplex v1[TAM], fcomplex v2[TAM])
//Devuelve el producto escalar de dos vectores
{
	int i;
	fcomplex aux;
	
	aux=Complex(0.,0.);
	for(i=0;i<TAM;i++)
		aux=Cadd(aux,Cmul(Conjg(v1[i]),v2[i]));
	
	return aux;	
}


int BraSist (fcomplex ket[TAM], fcomplex bra[TAM])
//Devuelve el vector conjugado del ket
{
	int i;
	
	for(i=0;i<TAM;i++)
		bra[i]=Conjg(ket[i]);
		
	return 0;		
}


int CProyectorSist (fcomplex ket[TAM], fcomplex bra[TAM], fcomplex proyector[TAM][TAM])	
//Devuelve matriz proyector de dos vectores, bra y ket
{
	
	int i,j;
	
	int braaux;
	fcomplex aux[TAM];
	
	braaux=BraSist(bra,aux);
	
	for(i=0;i<TAM;i++)
		for(j=0;j<TAM;j++)
		{	
			
			proyector[i][j]=Cmul(ket[i],aux[j]);
		}
		
	return	0;
}

int ActOperadorSist (fcomplex Operador[TAM][TAM], fcomplex ket[TAM], fcomplex resultado[TAM])
//Devuelve el vector del producto de una matriz y otro vector
{
	int i,j;
	fcomplex aux;
	
	for(i=0;i<TAM;i++)
	{
		aux=Complex(0.,0.);
		for(j=0;j<TAM;j++)
		{
			aux=Cadd(aux,Cmul(Operador[i][j],ket[j]));
		}
		resultado[i]=aux;
	}
	
	return 0;
}

fcomplex ValorEspSist (fcomplex Operador[TAM][TAM], fcomplex ket[TAM], fcomplex bra[TAM])
//Devuelve el escalar de vector-matriz-vector
{
	fcomplex valor,pto;

	int control;
	double norma;
	fcomplex dcha[TAM];
			
	control=ActOperadorSist (Operador, ket, dcha);
	
	valor=PdtoEscalarSist(bra, dcha); //El conjugado del bra se calcula dentro de PdtoEscalar7
		
	return valor;
}

////////////////////////////////////////////////////////////////////////////////////////////
//OPERACIONES CON VECTORES

int VaddSist(fcomplex V1[TAM], fcomplex V2[TAM], fcomplex SUMA[TAM])
//Devuelve V1+V2
{
	int j;
	
	for(j=0;j<TAM;j++)
		SUMA[j]=Cadd(V1[j],V2[j]);
			
	return 0;
}

int VsubSist(fcomplex V1[TAM], fcomplex V2[TAM], fcomplex RESTA[TAM])
//Devuelve V1-V2
{
	int j;
	
	for(j=0;j<TAM;j++)
		RESTA[j]=Csub(V1[j],V2[j]);
			
	return 0;
}

int VPdtoExt(fcomplex V1[TAM], fcomplex V2[TAM], fcomplex PdtoExt[TAM])
{
	int i,j;
	
	for(i=0;i<tam;i++)
		for(j=0;j<tam;j++)
			PdtoExt[tam*i+j]=Cmul(V1[i],V2[j]); 
	
	return 0;			
}

//OPERACIONES CON MATRICES

int MaddSist(fcomplex M1[TAM][TAM], fcomplex M2[TAM][TAM], fcomplex SUMA[TAM][TAM])
//Devuelve M1+M2
{
	int i,j;
	
	for(i=0;i<TAM;i++)
		for(j=0;j<TAM;j++)
			SUMA[i][j]=Cadd(M1[i][j],M2[i][j]);
			
	return 0;
}

int MsubSist(fcomplex M1[TAM][TAM], fcomplex M2[TAM][TAM], fcomplex RESTA[TAM][TAM])
//Devuelve M1-M2
{
	int i,j;
	
	for(i=0;i<TAM;i++)
		for(j=0;j<TAM;j++)
			RESTA[i][j]=Csub(M1[i][j],M2[i][j]);
			
	return 0;
}

int MmulSist (fcomplex M1[TAM][TAM], fcomplex M2[TAM][TAM], fcomplex PDTO[TAM][TAM])
//Devuelve M1*M2
{
	int i,j,k;
	
	
	
	for(i=0;i<TAM;i++)
		for(j=0;j<TAM;j++)
		{
			PDTO[i][j]=Complex(0.,0.);
			for(k=0;k<TAM;k++)
				PDTO[i][j]=Cadd(PDTO[i][j], Cmul(M1[i][k],M2[k][j]));
		}	
	return 0;
}


int MPdtoExt(fcomplex M1[tam][tam], fcomplex M2[tam][tam], fcomplex PdtoExt[TAM][TAM])
{
	int i,j,k,l;

	for(i=0;i<tam;i++)
		for(j=0;j<tam;j++)
			for(k=0;k<tam;k++)
				for(l=0;l<tam;l++)
					PdtoExt[tam*i+k][tam*j+l]=Cmul(M1[i][j],M2[k][l]);			
	
	return 0;	
}
