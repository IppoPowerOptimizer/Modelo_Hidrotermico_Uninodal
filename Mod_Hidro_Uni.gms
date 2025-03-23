
*==============================================================
****************************************************************
*********MODELO HIDROTERMICO UNINODAL***************************
****************************************************************
*==============================================================

*==============================================================
****************************************************************
*************HAYRO ANTHONY PUMALLOCLLA CASTILLA*****************
****************************************************************
*==============================================================

*Declaracion de conjuntos
Sets

gt      Generadores térmicos    /gt1/
gh      Generadores hidro       /gh1/
t       periodos                /t1*t3/
pargt   Parámetros de hidro     /Pmin, Pmax, c/
pargh   Parámetros de termicas  /Pmin, Pmax, ci/
;

table
Datahidro
        Pmin    Pmax    c       
gh1      0      450     0.8        
;

Table
Dataterm
       Pmin     Pmax    ci
gt1    0        600     5
;

Parameters

period(t) Horas por bloque
/
t1      8
t2      10
t3      6
/

dem(t) Demanda por cada bloque
/
t1      350
t2      700
t3      500
/

meta(gh) consumo de agua diario 
/
gh1     1.5
/
;

Variables

C_gt        Costo de operacion
pt(gt,t)    Potencia generada por la central termica por bloque horario
ph(gh,t)    Potencia generada por la central hidro por bloque horario
;

Scalar

c_rac       Coste de racionamiento
/8000/
;

Positive Variable
raci(t)     racionamiento
;
pt.lo(gt,t) = Dataterm(gt, 'Pmin');
pt.up(gt,t) = Dataterm(gt, 'Pmax');
ph.lo(gh,t) = Datahidro(gh,'Pmin');
ph.up(gh,t) = Datahidro(gh,'Pmax');

Equations   

FO          Funcion objetivo
BE(t)       Balance energetico
BH(gh)      Balance de volumen de agua
;


FO..    C_gt =e= sum(t, sum(gt, Dataterm(gt, 'ci') * pt(gt,t)) * period(t)) + sum(t, raci(t) * c_rac* period(t)) ;
BE(t).. (sum(gt, pt(gt, t)) + sum(gh, ph(gh, t)) + raci(t)) * period(t)  =e= dem(t)* period(t) ;
BH(gh).. Sum(t, Datahidro(gh, 'c') * ph(gh,t) * period(t)) * 3600 =l= meta(gh) * power(10, 6);

Model MHU /all/;

Solve MHU using lp minimizing C_gt ;
