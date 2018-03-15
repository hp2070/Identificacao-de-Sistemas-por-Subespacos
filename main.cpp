#include <iostream>
#include <armadillo>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <time.h>



using namespace std;
using namespace arma;

typedef long clock_t;
double clkini=0,clkfim;
#define CLOCKS_PER_SEC 1000 // NÃO TENHO CERTEZA DO CLOCK INTERNO

void plot(const char* diretorio,mat x){
    FILE *arq; // ABRINDO ARQUIVO
    arq = fopen(diretorio, "w");

    int nrowX=0, ncolX=0;

    nrowX = x.n_rows;
    ncolX = x.n_cols;

    printf("%s : {%d,%d} \n",diretorio,nrowX,ncolX);
    for (int i=0; i<nrowX; i++){
        for(int j=0; j<ncolX; j++){
            fprintf(arq,"%.16f ",x(i,j));
        }
        fprintf(arq,"\n\n\n\n");
    }
    fclose(arq);

}
mat load(const char* diretorio,int nrow, int ncol){

    FILE *arq; // ABRINDO ARQUIVO
    arq = fopen(diretorio, "rt");

    mat prov(nrow,ncol); // MATRIZ PROVISÓRIA

    int range = 16; // TAMANHO DE CADA NUMERO
    char linha[ncol*(range+1)];

    int i = 0; // VARIÁVEIS AUXILIARES
    double num = 0.0000000000000000;
    double expo = 0.0000000000000000;



    while (i<nrow)
  {
      fgets(linha, ncol*(range+2), arq);
      num = 0;
	  for(int k=0; k<ncol; k++){

            num = linha[(range*k)+3] -   48.0000000000000000;
            num += (linha[(range*k)+5] - 48.0000000000000000)/10.0000000000000000;
            num += (linha[(range*k)+6] - 48.00000000000000)/100.0000000000000000;
            num += (linha[(range*k)+7] - 48.0000000000000000)/1000.0000000000000000;
            num += (linha[(range*k)+8] - 48.0000000000000000)/10000.0000000000000000;
            num += (linha[(range*k)+9] - 48.0000000000000000)/100000.0000000000000000;
            num += (linha[(range*k)+10] -48.0000000000000000)/1000000.0000000000000000;
            num += (linha[(range*k)+11] -48.0000000000000000)/10000000.0000000000000000;

            expo = (linha[(range*k)+14] - 48.0000000000000000) * 10.0000000000000000;
            expo += (linha[(range*k)+15 ] - 48.0000000000000000);

            if(linha[(range*k)+13] == '-') expo = -1.0000000000000000*expo;
            if(linha[(range*k)+2] == '-') num = -1.0000000000000000*num;

            num = num*pow(10.0000000000000000,expo);
            prov(i,k) = num;
           // printf("%.7f || %.1f\n",num,expo);
	  }
      i++;
  }
  fclose(arq);
  return prov;
}
mat blkhank(mat y,unsigned i,unsigned j){
    int l = 0,nd = 0;
    mat H;

    l = y.n_rows;
    nd = y.n_cols;

    if(l > nd){ y = y.t(); l = y.n_rows; nd = y.n_cols;}
    if(i<0){printf("ERRO - I nao pode ser negativo");return 0;}
    if(j<0){printf("ERRO - J nao pode ser negativo");return 0;}
    if(j> nd-i+1){printf("ERRO - Matriz de Hankel muito Grande"); return 0;}



    H.zeros(2*l*i,j);

    unsigned cont = 0;
    unsigned cont2 = 0;
    for(unsigned k=1; k<(i*2)+1; k++){
        cont =0;
         for(unsigned nent=(k-1)*l; nent<k*l; nent++){
                cont2 = 0;
                for(unsigned elem=k-1; elem <j+k-2; elem++){
                    H(nent,cont2) =  y(cont,elem);
                    cont2 = cont2 + 1;
                }
               cont = cont + 1;
         }
    }
    return H;
}
mat copyrange(mat prim,unsigned fromi,unsigned toi,unsigned fromj, unsigned toj,bool cut){
    mat prov;
    prov.set_size(toi,toj);
    //printf("{%d:%d} {%d:%d} \n",fromi,toi,fromj,toj);
    for(unsigned i=fromi; i<toi; i++){
        for(unsigned j=fromj; j<toj; j++){
            prov(i,j) = prim(i,j);
        }
    }
    if(cut == true){
        if(fromi > 0){
        prov.shed_rows(0,fromi-1);
        }
        if(fromj > 0){
        prov.shed_cols(0,fromj-1);
        }
    }

    return prov;
}
double tempoExecucao(){
    clkfim = (double) clock()/CLOCKS_PER_SEC;
    double result = (clkfim-clkini)*1000;
    printf(" (%.3f ms)\n",result);
    clkini = (double)clock()/CLOCKS_PER_SEC;
    return result;
}
void det_stat(mat y,mat u,unsigned i,unsigned n,mat &A,mat &B,mat &C, mat &D, unsigned debug){

    mat u_idnormal,Y,R,Q,prov;
    mat y_idnormal,U;
    mat Rf,Rp,Rp2,Ru,Rfp,Rfp2,Rpp,Rpp2,Ob,Obm,erroTamanho;
    mat Udec,Vdec,U1dec;
    mat gam, gamm,gam_inv,gamm_inv;
    mat Xi, Xip;
    mat Lhs,Rhs,sol;
    vec Sdec;
    double somatime;

    double j  = 0.0000000000000000;
    int mi2 = 0;
    int m = 0, nu=0, l=0, ny=0;

    m = u.n_rows; nu = u.n_cols;
    l = y.n_rows; ny = y.n_cols;

    clkini = (double)clock()/CLOCKS_PER_SEC;

    if(debug == 1){
        plot("u_id.txt",u);
        plot("y_id.txt",y);
    }

    j = nu-2*i+1;
    u_idnormal = u/sqrt(j);      ///NORMALIZAÇÃO DAS MATRIZES DE IDENTIFICAÇÃO
    y_idnormal = y/sqrt(j);

    if(debug == 1){
        printf("\nJ = %f\n",j);
        printf("\nI = %d\n",i);
        plot("u_idnormal.txt",u_idnormal);
        plot("y_idnormal.txt",y_idnormal);
    }

    printf("Computing Hankel Matrices Y and U...");
    Y = blkhank(y_idnormal,i,j);                   /// GERANDO MATRIZES DE HANKEL
    U = blkhank(u_idnormal,i,j);
    if(debug == 1){
        plot("Y.txt",Y);
        plot("U.txt",U);
    }
    somatime += tempoExecucao();

    printf("QR Decomposition...");
    qr_econ(Q,R,join_cols(U,Y).t());               /// DECOMPOSIÇÃO QR
    R = R.t();

    if(debug == 1) plot("RBig.txt",R);

    R = copyrange(R,0,2*i*(m+l),0,2*i*(m+l),true); ///SEGMENTO IMPORTANTE DE R
    mi2  = 2*m*i;
    if(debug == 1){
         plot("RSecc.txt",R);
         printf("\nMI2 = %d \n",mi2);
    }
    somatime += tempoExecucao();

    printf("Computing Rf,Rp,Ru,Rfp,Rpp...");
    Rf = copyrange(R,(2*m+l)*i,2*(m+l)*i,0,R.n_cols,true);  /// Rf

    Rp = copyrange(R,0,m*i,0,R.n_cols,true);
    Rp2 = copyrange(R,(2*m*i),(2*m+l)*i,0,R.n_cols,true);   /// Rp
    Rp = join_cols(Rp,Rp2);

    Ru = copyrange(R,m*i,2*m*i,0,mi2,true);                 /// Ru

    Rfp = copyrange(Rf,0,Rf.n_rows,0,mi2,true);             /// Rfp
    Rfp2 = Rfp*pinv(Ru,0.00000000000000001,"std");
    Rfp2 = Rfp2*Ru;
    Rfp = Rfp - Rfp2;
    Rfp2.clear();
    Rfp2 = copyrange(Rf,0,Rf.n_rows,mi2,2*(m+l)*i,true);
    Rfp = join_rows(Rfp,Rfp2);

    Rpp = copyrange(Rp,0,Rp.n_rows,0,mi2,true);            /// Rpp
    Rpp2 = Rpp*pinv(Ru,0.00000000000000001,"std");
    Rpp2 = Rpp2*Ru;
    Rpp = Rpp - Rpp2;
    Rpp2.clear();
    Rpp2 = copyrange(Rp,0,Rp.n_rows,mi2,2*(m+l)*i,true);
    Rpp = join_rows(Rpp,Rpp2);

    if(debug == 1){
    plot("Ru.txt",Ru);
    plot("Rp.txt",Rp);
    plot("Rf.txt",Rf);
    plot("Rfp.txt",Rfp);
    plot("Rpp.txt",Rpp);
    }
    somatime += tempoExecucao();

    printf("Computing oblique projection Ob...");
    erroTamanho = copyrange(Rpp,0,Rpp.n_rows,(2*m+l)*i-2*l,((2*m+l)*i)+1,true);
    erroTamanho = norm(erroTamanho,"fro");

    if(debug == 1 )plot("erroTamanho.txt",erroTamanho);

    if(erroTamanho(0,0) < 1e-10){ // EVITAR ERRO DE TAMANHO
        Ob = Rpp.t();
        Ob = pinv(Ob,0.00000000000000001,"std");
        Ob = Ob.t();
        Ob = Rfp*Ob;
        Ob = Ob*Rp;             /// CALCULO DA PRIMEIRA PROJEÇÃO OBLIQUA

    }
    else{

        Ob = Rpp;                                   // Ob = Rpp
        Ob = pinv(Ob,0.00000000000000001,"std");    // Ob = 1/Rpp;
        Ob = Rfp*Ob;                                // Ob = Rfp/Rpp;
        Ob = Ob*Rp;                                    // Ob = (Rfp/Rpp)*Rp;
    }
    if(debug == 1) plot("Ob.txt",Ob);
    somatime += tempoExecucao();

    printf("Computing SVD...");
    svd(Udec,Sdec,Vdec,Ob);

    if(debug == 1){
    plot("Udec.txt",Udec);
    plot("Sdec.txt",Sdec);
    plot("Vdec.txt",Vdec);
    }


    U1dec = copyrange(Udec,0,Udec.n_rows,0,n,true);
    if(debug == 1) plot("U1dec.txt",U1dec);

    somatime += tempoExecucao();

    printf("Computing Gam and Gamm...");
    gam = U1dec*diagmat(sqrt(copyrange(Sdec,0,n,0,1,true)));     /// DETERMINANDO GAM E GAMM
    gamm = copyrange(U1dec,0,l*(i-1),0,U1dec.n_cols,true)*diagmat(sqrt(copyrange(Sdec,0,n,0,1,true)));
    gam_inv = pinv(gam,0.00000000000000001,"std");
    gamm_inv = pinv(gamm,0.00000000000000001,"std");

    if(debug == 1){
    plot("gam.txt",gam);
    plot("gamm.txt",gamm);
    plot("gam_inv.txt",gam_inv);
    plot("gamm_inv.txt",gamm_inv);
    }
    somatime += tempoExecucao();

    Rf.clear();Rp.clear();Rp2.clear();Ru.clear();Rfp.clear();Rfp2.clear();Rpp.clear();Rpp2.clear();

    printf("Computing Rf2,Rp2,Ru2,Rfp2,Rpp2...");
    Rf = copyrange(R,(2*m+l)*i+l,2*(m+l)*i,0,R.n_cols,true);

    Rp = copyrange(R,0,m*(i+1),0,R.n_cols,true);                    /// Rp DA SEGUNDA PROJEÇÃO OBLIQUA
    Rp2 = copyrange(R,(2*m*i),(2*m+l)*i+l,0,R.n_cols,true);
    Rp = join_cols(Rp,Rp2);

    Ru = copyrange(R,m*i+m,2*m*i,0,mi2,true);                      /// Ru DA SEGUNDA PROJEÇÃO OBLIQUA



    Rfp = copyrange(Rf,0,Rf.n_rows,0,mi2,true);                    /// Rfp DA SEGUNDA PROJEÇÃO OBLIQUA
    Rfp2 = Rfp*pinv(Ru,0.00000000000000001,"std");
    Rfp2 = Rfp2*Ru;
    Rfp = Rfp - Rfp2;
    Rfp2.clear();
    Rfp2 = copyrange(Rf,0,Rf.n_rows,mi2,2*(m+l)*i,true);
    Rfp = join_rows(Rfp,Rfp2);

    Rpp = copyrange(Rp,0,Rp.n_rows,0,mi2,true);                    /// Rpp DA SEGUNDA PROJEÇÃO OBLIQUA
    Rpp2 = Rpp*pinv(Ru,0.00000000000000001,"std");
    Rpp2 = Rpp2*Ru;
    Rpp = Rpp - Rpp2;
    Rpp2.clear();
    Rpp2 = copyrange(Rp,0,Rp.n_rows,mi2,2*(m+l)*i,true);
    Rpp = join_rows(Rpp,Rpp2);

    if(debug == 1){
    plot("Ru2.txt",Ru);
    plot("Rp2.txt",Rp);
    plot("Rf2.txt",Rf);
    plot("Rfp2.txt",Rfp);
    plot("Rpp2.txt",Rpp);
     }

    somatime += tempoExecucao();

    printf("Computing Oblique projection Obm...");
    erroTamanho = copyrange(Rpp,0,Rpp.n_rows,(2*m+l)*i-2*l,((2*m+l)*i)+1,true);
    erroTamanho = norm(erroTamanho,"fro");

    if(debug == 1 )plot("erroTamanho.txt",erroTamanho);

    if(erroTamanho(0,0) < 1e-10){ // EVITAR ERRO DE TAMANHO
        Obm = Rpp.t();
        Obm = pinv(Obm,0.00000000000000001,"std");
        Obm = Obm.t();
        Obm = Rfp*Obm;
        Obm = Obm*Rp;             /// CALCULO DA PRIMEIRA PROJEÇÃO OBLIQUA

    }
    else{

        Obm = Rpp;                                   // Ob = Rpp
        Obm = pinv(Obm,0.00000000000000001,"std");    // Ob = 1/Rpp;
        Obm = Rfp*Obm;                                // Ob = Rfp/Rpp;
        Obm = Obm*Rp;                                    // Ob = (Rfp/Rpp)*Rp;
    }

    if(debug == 1) plot("Obm.txt",Obm);
    somatime += tempoExecucao();

    printf("Computing Xi and Xip...");
    Xi = gam_inv * Ob;
    Xip = gamm_inv * Obm;

    if(debug == 1){
    plot("Xi.txt",Xi);
    plot("Xip.txt",Xip);
    }
    somatime += tempoExecucao();

    Rhs = join_cols(Xi,copyrange(R,m*i,m*(i+1),0,R.n_cols,true)); ///   SISTEMA DE EQUAÇÕES X*A = B
    Lhs = join_cols(Xip,copyrange(R,(2*m+l)*i,(2*m+l)*i+l,0,R.n_cols,true));

    printf("Solving Rhs/Lhs ...");
    sol = solve(Rhs.t(),Lhs.t());
    sol = sol.t();

    if(debug == 1){
    plot("Rhs.txt",Rhs);
    plot("Lhs.txt",Lhs);
    plot("sol.txt",sol);
    }
    somatime += tempoExecucao();

    printf("Extracting Matrices A,B,C,D...");
    A = copyrange(sol,0,n,0,n,true);    ///EXTRAINDO MATRIZ A DA SOLUÇÃO
    B = copyrange(sol,0,n,n,n+m,true);  ///EXTRAINDO MATRIZ B DA SOLUÇÃO
    C = copyrange(sol,n,n+l,0,n,true);  ///EXTRAINDO MATRIZ C DA SOLUÇÃO
    D = copyrange(sol,n,n+l,n,n+m,true);///EXTRAINDO MATRIZ D DA SOLUÇÃO

    printf("\nElapsed Time : (%.3f ms)\n",somatime);


}
void subid(mat y, mat u, unsigned i, unsigned n, mat &A, mat &B, mat &C, mat &D, mat &K, unsigned debug){

    mat u_idnormal,Y,R,prov,P,Q,L1,L2,M,X,totm,N;
    mat y_idnormal,U;
    mat Rf,Rp,Rp2,Ru,Rfp,Rfp2,Rpp,Rpp2,Ob,erroTamanho,WOW,WOW2;
    mat Udec,Vdec,U1dec;
    mat gam,gamm,gam_inv,gamm_inv;
    mat Xi, Xip;
    mat Lhs,Rhs,sol,res,sol_bd;
    vec Sdec;
    double somatime;

    double j  = 0.0000000000000000;
    int mi2 = 0;
    int m = 0, nu=0, l=0, ny=0;

    m = u.n_rows; nu = u.n_cols;
    l = y.n_rows; ny = y.n_cols;

    if(debug == 1){
        plot("u_id2.txt",u);
        plot("y_id2.txt",y);
    }

    j = nu-2*i+1;
    u_idnormal = u/sqrt(j);      ///NORMALIZAÇÃO DAS MATRIZES DE IDENTIFICAÇÃO
    y_idnormal = y/sqrt(j);

    if(debug == 1){
        printf("\nJ = %f\n",j);
        printf("\nI = %d\n",i);
        plot("u_idnormal.txt",u_idnormal);
        plot("y_idnormal.txt",y_idnormal);
    }

    clkini = (double)clock()/CLOCKS_PER_SEC;

    printf("Computando Matriz de Henkel Y e U...");
    Y = blkhank(y_idnormal,i,j);                   /// GERANDO MATRIZES DE HANKEL
    U = blkhank(u_idnormal,i,j);
    if(debug == 1){
        plot("Y.txt",Y);
        plot("U.txt",U);
    }
    somatime += tempoExecucao();

    printf("Decomposicao QR...");
    qr_econ(Q,R,join_cols(U,Y).t());               /// DECOMPOSIÇÃO QR
    R = R.t();

    if(debug == 1) plot("RBig.txt",R);

    R = copyrange(R,0,2*i*(m+l),0,2*i*(m+l),true); ///SEGMENTO IMPORTANTE DE R
    mi2  = 2*m*i;
    if(debug == 1){
         plot("R.txt",R);
         printf("\nMI2 = %d \n",mi2);
    }
    somatime += tempoExecucao();

    printf("Computando Rf,Rp,Ru,Rfp,Rpp...");
    Rf = copyrange(R,(2*m+l)*i,2*(m+l)*i,0,R.n_cols,true);  /// Rf

    Rp = copyrange(R,0,m*i,0,R.n_cols,true);
    Rp2 = copyrange(R,(2*m*i),(2*m+l)*i,0,R.n_cols,true);   /// Rp
    Rp = join_cols(Rp,Rp2);

    Ru = copyrange(R,m*i,2*m*i,0,mi2,true);                 /// Ru

    Rfp = copyrange(Rf,0,Rf.n_rows,0,mi2,true);             /// Rfp
    Rfp2 = Rfp*pinv(Ru,0.00000000000000001,"std");
    Rfp2 = Rfp2*Ru;
    Rfp = Rfp - Rfp2;
    Rfp2.clear();
    Rfp2 = copyrange(Rf,0,Rf.n_rows,mi2,2*(m+l)*i,true);
    Rfp = join_rows(Rfp,Rfp2);

    Rpp = copyrange(Rp,0,Rp.n_rows,0,mi2,true);            /// Rpp
    Rpp2 = Rpp*pinv(Ru,0.00000000000000001,"std");
    Rpp2 = Rpp2*Ru;
    Rpp = Rpp - Rpp2;
    Rpp2.clear();
    Rpp2 = copyrange(Rp,0,Rp.n_rows,mi2,2*(m+l)*i,true);
    Rpp = join_rows(Rpp,Rpp2);

    if(debug == 1){
    plot("Ru.txt",Ru);
    plot("Rp.txt",Rp);
    plot("Rf.txt",Rf);
    plot("Rfp.txt",Rfp);
    plot("Rpp.txt",Rpp);
    }
    somatime += tempoExecucao();

    printf("Computando Projecao Obliqua Ob...");
    erroTamanho = copyrange(Rpp,0,Rpp.n_rows,(2*m+l)*i-2*l,((2*m+l)*i)+1,true);
    erroTamanho = norm(erroTamanho,"fro");
    if(debug == 1 )plot("erroTamanho.txt",erroTamanho);
    if(erroTamanho(0,0) < 1e-10){ // EVITAR ERRO DE TAMANHO
        Ob = Rpp.t();
        Ob = pinv(Ob,0.00000000000000001,"std");
        Ob = Ob.t();
        Ob = Rfp*Ob;
        Ob = Ob*Rp;             /// CALCULO DA PRIMEIRA PROJEÇÃO OBLIQUA

    }
    else{
        Ob = Rpp;                                   // Ob = Rpp
        Ob = pinv(Ob,0.00000000000000001,"std");    // Ob = 1/Rpp;
        Ob = Rfp*Ob;                                // Ob = Rfp/Rpp;
        Ob = Ob*Rp;                                    // Ob = (Rfp/Rpp)*Rp;
    }
    if(debug == 1) plot("Ob.txt",Ob);
    somatime += tempoExecucao();

    printf("Computando WOW...");
    WOW = copyrange(Ob,0,Ob.n_rows,0,mi2,true)*pinv(Ru,0.00000000000000001,"std"); //(Ob(:,1:mi2)/Ru)
    WOW = WOW*Ru; //(Ob(:,1:mi2)/Ru)*Ru
    WOW = copyrange(Ob,0,Ob.n_rows,0,mi2,true) - WOW; //Ob(:,1:mi2) - (Ob(:,1:mi2)/Ru)*Ru
    WOW2 = copyrange(Ob,0,Ob.n_rows,mi2,2*(m+l)*i,true);
    WOW = join_rows(WOW,WOW2);
    if(debug == 1)  plot("WOW.txt",WOW);
    somatime += tempoExecucao();

    printf("Computando SVD...");
    svd(Udec,Sdec,Vdec,WOW);
    if(debug == 1){
    plot("Udec.txt",Udec);
    plot("Sdec.txt",Sdec);
    plot("Vdec.txt",Vdec);
    }
    U1dec = copyrange(Udec,0,Udec.n_rows,0,n,true);
    if(debug == 1) plot("U1dec.txt",U1dec);
    somatime += tempoExecucao();

    printf("Computando Gam e Gamm...");
    gam = U1dec*diagmat(sqrt(copyrange(Sdec,0,n,0,1,true)));     /// DETERMINANDO GAM E GAMM
    gamm = copyrange(gam,0,l*(i-1),0,gam.n_cols,true);
    gam_inv = pinv(gam,0.00000000000000001,"std");
    gamm_inv = pinv(gamm,0.00000000000000001,"std");
    if(debug == 1){
    plot("gam.txt",gam);
    plot("gamm.txt",gamm);
    plot("gam_inv.txt",gam_inv);
    plot("gamm_inv.txt",gamm_inv);
    }
    somatime += tempoExecucao();

    printf("Solucao Rhs/Lhs ...");
    Rhs = join_rows(gam_inv*copyrange(R,(2*m+l)*i,2*(m+l)*i,0,(2*m+l)*i,true),zeros(n,l));
    Rhs = join_cols(Rhs,copyrange(R,m*i,2*m*i,0,(2*m+l)*i+l,true));
    Lhs = join_cols(gamm_inv*copyrange(R,(2*m+l)*i+l,2*(m+l)*i,0,(2*m+l)*i+l,true),copyrange(R,(2*m+l)*i,(2*m+l)*i+l,0,(2*m+l)*i+l,true));
    sol = solve(Rhs.t(),Lhs.t());
    sol = sol.t();
    if(debug == 1){
    plot("Rhs.txt",Rhs);
    plot("Lhs.txt",Lhs);
    plot("sol.txt",sol);
    }
    somatime += tempoExecucao();

    printf("Extraindo A , C e Residuais...");
    A = copyrange(sol,0,n,0,n,true);
    C = copyrange(sol,n,n+l,0,n,true);
    if(debug == 1){
    plot("A.txt",A);
    plot("C.txt",C);
    }
    res = Lhs - (sol*Rhs);

    plot("res.txt",res);

    gam.clear();
    gam = C;
    mat provis;
    unsigned aux = 0;
    gam.resize(l*i,n);
    provis.set_size(l,n);
    somatime += tempoExecucao();

    printf("Produtos de Kronecker ...");
    for(unsigned k=2; k<=i; k++){
    provis = copyrange(gam,(k-2)*l,(k-1)*l,0,gam.n_cols,true);
    provis = provis*A;
    aux = 0;
    for(unsigned linha=(k-1)*l; linha<k*l; linha++){
        for(unsigned coluna=0; coluna<gam.n_cols; coluna++){
            gam(linha,coluna) = provis(aux,coluna);
        }
        aux++;
    }
    provis.clear();
    }
    somatime += tempoExecucao();

    printf("Computando Gam e Gamm...");
    gamm.clear();
    gam_inv.clear();
    gamm_inv.clear();
    gamm = copyrange(gam,0,l*(i-1),0,gam.n_cols,true);
    gam_inv = pinv(gam,0.00000000000000001,"std");
    gamm_inv = pinv(gamm,0.00000000000000001,"std");
    if(debug == 1){
    plot("gamm2.txt",gamm);
    plot("gam_inv2.txt",gam_inv);
    plot("gamm_inv2.txt",gamm_inv);
    }
    somatime += tempoExecucao();

    printf("Solucao Rhs/Lhs ...");
    Rhs = join_rows(gam_inv*copyrange(R,(2*m+l)*i,2*(m+l)*i,0,(2*m+l)*i,true),zeros(n,l));
    Rhs = join_cols(Rhs,copyrange(R,m*i,2*m*i,0,(2*m+l)*i+l,true));
    Lhs = join_cols(gamm_inv*copyrange(R,(2*m+l)*i+l,2*(m+l)*i,0,(2*m+l)*i+l,true),copyrange(R,(2*m+l)*i,(2*m+l)*i+l,0,(2*m+l)*i+l,true));
    if(debug == 1){
    plot("Rhs2.txt",Rhs);
    plot("Lhs2.txt",Lhs);
    }
    somatime += tempoExecucao();

    printf("Computando P,Q,L1,L2,M,X ...");
    P = Lhs - join_cols(A,C)*copyrange(Rhs,0,n,0,Rhs.n_cols,true);
    P = copyrange(P,0,P.n_rows,0,(2*m*i),true);
    Q = copyrange(R,m*i,2*m*i,0,2*m*i,true);
    if(debug == 1){
    plot("P.txt",P);
    plot("Q.txt",Q);
    }
    L1 = A*gam_inv;
    L2 = C*gam_inv;
    M = join_rows(zeros(n,l),gamm_inv);
    X = join_cols(join_rows(eye(l,l),zeros(l,n)),join_rows(zeros(l*(i-1),l),gamm));
    if(debug == 1){
    plot("L1.txt",L1);
    plot("L2.txt",L2);
    plot("M.txt",M);
    plot("X.txt",X);
    }
    somatime += tempoExecucao();

    printf("Computando totm ...");
    totm.clear();
    for(unsigned k=1; k<=i; k++){
        N = join_cols(join_rows(copyrange(M,0,M.n_rows,(k-1)*l,l*i,true)-copyrange(L1,0,L1.n_rows,(k-1)*l,l*i,true),zeros(n,(k-1)*l)),join_rows(-1*copyrange(L2,0,L2.n_rows,(k-1)*l,l*i,true),zeros(l,(k-1)*l)));

        if(k==1){
            provis.clear();
            provis = eye(l,l) + copyrange(N,n,n+l,0,l,true);
            aux = 0;
            for(unsigned linha=n; linha<n+l; linha++){
                for(unsigned coluna=0; coluna<l; coluna++){
                    N(linha,coluna) = provis(aux,coluna);
                }
            aux++;
            }
        }
        N = N*X;

        if(k==1){
        totm = kron(copyrange(Q,(k-1)*m,k*m,0,Q.n_cols,true).t(),N);
        }else{
        totm = totm + kron(copyrange(Q,(k-1)*m,k*m,0,Q.n_cols,true).t(),N);
        }
    }
    if(debug == 1){
    plot("N.txt",N);
    plot("totm.txt",totm);
    }
    somatime += tempoExecucao();

    printf("Computando P ...");
    unsigned aux2 = 0;
    aux =0;
    provis.clear();
    provis.set_size(P.n_cols*P.n_rows,1);
    if(debug == 1){
    plot("P.txt",P);
    }


    for(aux = 0; aux < P.n_cols; aux++){
        for(unsigned k = 0; k < P.n_rows; k++){
            provis(aux2,0) = P(k,aux);
            aux2++;
        }
    }

    P = provis;
    if(debug == 1){
    plot("P.txt",P);
    }
    somatime += tempoExecucao();

    printf("Solucao totm/P ...");
    sol.clear();
    sol = solve(totm,P);
    sol.reshape(n+l,m);
    sol_bd = sol;
    somatime += tempoExecucao();

    printf("Extraindo B e D ...");
    D = copyrange(sol_bd,0,l,0,sol_bd.n_cols,true);
    B = copyrange(sol_bd,l,l+n,0,sol_bd.n_cols,true);
    somatime += tempoExecucao();

    printf("\nElapsed Time : (%.3f ms)\n",somatime);
}



int main()
{
    mat u,u_id,u_val;
    mat y,y_id,y_val;

    mat A,B,C,D,K;

    int n_id = 0;                    // QUANTIDADE DE ELEMENTOS PARA IDENTIFICAÇÃO
    int n_val = 0;                   // QUANTIDADE DE ELEMENTOS PARA VALIDAÇÃO
    int i = 10;
    unsigned debug;
    unsigned nordem =0;

    char entradas[50],saidas[50];
    unsigned numsaidas,numentradas;
    unsigned stop;

    printf("                       _________________________________________________________________________\n");
    printf("                      |Subspace Identification   29.10.2017                                     |\n");
    printf("                      |Prof. Dr. Santos Demetrio Borjas                - DMAT UFRN              |\n");
    printf("                      |Guilherme Afonso Pillon de C. A. Pessoa         - DEE  UFRN              |\n");
    printf("                      |_________________________________________________________________________| \n\n");


    printf("Input Matrices :");
    scanf("%s",&entradas);


    printf("Output Matrices :");
    scanf("%s",&saidas);

    printf("Number of Inputs :");
    scanf("%d",&numentradas);
    //numentradas = 7;

    printf("Number of Outputs :");
    scanf("%d",&numsaidas);
    //numsaidas = 6;

    printf("Qte. Valores para ID :");
    scanf("%d",&n_id);
    //n_id = 2000;

    printf("Qte. Valores para VAL :");
    scanf("%d",&n_val);
    //n_val = 1000;

    printf("System Order :");
    scanf("%d",&nordem);
    //nordem = 11;


    printf("Generate Debug Files [0][1] :");
    scanf("%d",&debug);


    u = load(entradas,n_id+n_val,numentradas);
    y = load(saidas,n_id+n_val,numsaidas);

    int m  = u.n_rows;
    int nu = u.n_cols;

    int l  = y.n_rows;
    int ny = y.n_cols;


    printf("\n\nNumber of Inputs: %d \n",min(m,nu));
    printf("Number of Outputs  : %d \n",min(l,ny));
    printf("Identification Samples  : %d \n",n_id);
    printf("Validation Samples  : %d \n",n_val);
    printf("Number of row blocks : %d \n\n",i);



    if(nu < m){ u = u.t(); m = u.n_rows; nu = u.n_cols;}
    if(ny < l){ y = y.t(); l = y.n_rows; ny = y.n_cols;}

    u_id = copyrange(u,0,min(m,nu),0,n_id,true);  ///IDENTIFICAÇÃO 0:n_id
    y_id = copyrange(y,0,min(l,ny),0,n_id,true);  ///IDENTIFICAÇÃO 0:n_id
    u_val = copyrange(u,0,min(m,nu),n_id,n_val+n_id,true); ///VALIDAÇÃO n_id:n_val
    y_val = copyrange(y,0,min(l,ny),n_id,n_val+n_id,true); ///VALIDAÇÃO n_id:n_val




    //det_stat(y_id,u_id,i,nordem,A,B,C,D,debug);
    subid(y_id,u_id,i,nordem,A,B,C,D,K,debug);

    plot("A.txt",A);
    plot("B.txt",B);
    plot("C.txt",C);
    plot("D.txt",D);

    printf("Identification complete successfully ! \n\n\n");
    scanf("%d",&stop);

    return 0;
}

