clf();         // close current figure
clear          // clear all pasta variables
xdel(winsid()) // close all windows
A = [  -0.2032588          0  0  0  0  0
        0.2032588  -0.052646  0  0  0  0
                0          0  0  0  0  0
                0          0  0  0  0  0
                0          0  0  0  0  0
                0          0  0  0  0  0];
B = [0.2126574  0
             0  0
             0  0
             0  0
             0  0
             0  0];
C = [1  0  0  0  0  0
     0  1  0  0  0  0
     0  0  1  0  0  0
     0  0  0  1  0  0
     0  0  0  0  1  0
     0  0  0  0  0  1];
D =[0  0
    0  0
    0  0
    0  0
    0  0
    0  0];
    

// Reducing the model
// X=[u, v, r]', U=[tau_u, tau_r]', Y=[u r]'
//los puntos de interrupción para proporcionar si desea un tratamiento lineal por partes.
ap=A(1:3,1:3);  
bp=B(1:3,:); 
cp=[C(1,1:3); 
    C(3,1:3)];
dp=[D(1,1:2);
    D(3,1:2)];
// Controllability and Observability
// Cc=[B, AB, A^2 B,..., A^(n-1) B]
Cc = cont_mat(ap, bp)
rankCc=rank(Cc)
//
// O=[C; CA; CA^2;...; CA^(n-1) ]
O = obsv_mat(ap, cp)
rankO=rank(O)
// verify if the rank of Cc is n, dimension of a
// verify if the rank of O is n, dimension of a

/*             Plot singular values of LTI the model                      */
G = syslin('c', ap, bp, cp, dp);
poles=spec(ap)
tr  = trzeros(G)

w = logspace(-3,3);
sv = svplot(G,w);

 
scf(1);
plot2d("ln", w, 20*log(sv')/log(10))
xgrid(12)
xtitle("Singular values plot","Frequency (rad/s)", "Amplitude 
(dB)");
Z=ss2tf(G)
scf(0);
plot2d("ln", w, 20*log(sv')/log(10))
xgrid(12)
xtitle("Singular values plot","Frequency (rad/s)", "Amplitude (dB)");


/*                                 Scaling                                 */
d2r=%pi/180;
r2d=180/%pi;

//su = diag( [1/200] );   // scaling input
                                    // tau_u_max=300 [N] tau_r_max=100 [Nm]
//sx = diag([1/1]);  // scaling state 

//sy = diag([1/1]);        // scaling output
su = diag( [1/200,1/(200*0.8)] );   // scaling input
                                    // tau_u_max=300 [N] tau_r_max=100 [Nm]
sx = diag([1/1,1/0.5,1/(50*d2r)]);  // scaling state 

sy = diag([1/1,1/(50*d2r)]);        // scaling output


ap_s = sx*ap*inv(sx)
bp_s = sx*bp*inv(su)
cp_s = sy*cp*inv(sx)
dp_s = sy*dp*inv(su)

//g=minreal(g)
G = syslin('c', ap_s, bp_s, cp_s, dp_s);
w = logspace(-3,3);
sv = svplot(G,w);

scf(1);
plot2d("ln", w, 20*log(sv')/log(10))
xgrid(12)
xtitle("Singular values plot","Frequency (rad/s)", "Amplitude (dB)");
Z=ss2tf(G)

scf(0);
plot2d("ln", w, 20*log(sv')/log(10))
xgrid(12)
xtitle("Singular values plot","Frequency (rad/s)", "Amplitude (dB)");
/*                                 Scaling                                 */
d2r=%pi/180;
r2d=180/%pi;

su = diag( [1/(200*0.8)] );   // scaling input
                                    // tau_u_max=300 [N] tau_r_max=100 [Nm]
sx = diag([1/1,1/(50*d2r)]);  // scaling state 

sy = diag([1/(50*d2r)]);        // scaling output

ap_s = sx*ap*inv(sx)
bp_s = sx*bp*inv(su)
cp_s = sy*cp*inv(sx)
dp_s = sy*dp*inv(su)

//g=minreal(g)
G = syslin('c', ap_s, bp_s, cp_s, dp_s);
w = logspace(-3,3);
sv = svplot(G,w);

scf(1);
plot2d("ln", w, 20*log(sv')/log(10))
xgrid(12)
xtitle("Singular values plot","Frequency (rad/s)", "Amplitude (dB)");


ms=1.7;// 0.3;%1.5;    % guarantee overshot Mp < 6dB = 20*log10(2) 
wbs=0.23;//0.05;%0.23;
ee=1e-3;//1e-4
ki=1; // used to give more accurate adjustment to the cut-off frequency wbs
      // by default set it to 1
//           --------     WT Data    ------------
mt=1.3;//1.00;    % guarantee overshot Mp < 2dB = 20*log10(1.26)
wbt=4.1;//9.1;%4.1;
ee=1e-3;//1e-4

//           --------     WS     ------------

s=poly(0,'s');
ws1=(s/ms+wbs)/(s+wbs*ee),
ws2=ws1;
ws=[ws1,0;0,ws2]
//Ws=syslin('c',ws)
Ws=blockdiag(ws1)

//           --------     WT     ------------

s=poly(0,'s');
wt1=(s+wbt/mt)/(ee*s+wbt),
wt2=wt1;
wt=[wt1,0;0,wt2]
//Wt=syslin('c',wt)
Wt=blockdiag(wt1)


//           --------     WR     ------------
s=poly(0,'s');
wr1=s/s,
wr2=wr1;
wr=[wr1,0;0,wr2]
wr=blockdiag(wr1)


// ------------------ Plot weighting functions

svs = svplot(Ws,w);
svt = svplot(Wt,w);
scf(2);
plot2d("ln", w, [-20*log(svs')/log(10) -20*log(svt')/log(10)])
xgrid(12)
xtitle("Singular values plot inv(Ws) and inv(Wt)","Frequency (rad/s)", "Amplitude (dB)");


[P,r]=augment(G,'ST');
//[P,r]=augment(g,'SRT');
//P = blockdiag(Ws,Wt,eye(G))*P1;
//P=minreal(P);

//trick to tackle when "D12 is not full rank"
//P.C(1,4)=0.0001;
P.D(1,2)=0.00001;
P.D(2,2)=1;

//r=[2,2]
romin=0.0001
romax=2000
nmax=100;
//[K,ro]=h_inf(P,r,romin,romax,nmax)
// alternatives
//[AK,BK,CK,DK,(RCOND)] = hinf(P.A,P.B,P.C,P.D,2,2,4)
K = ccontrg(P, r, 1.3) // this is good for me and for this system

// -------------- Analysis of the Feeedback Control System

[Se,Re,Te]=sensi(G,K) // S=(I+GK)^-1, T=I-S=GK(I+GK)^-1

// ------------------ Plot weighting functions

svS = svplot(Se,w);
svT = svplot(Te,w);
scf(3);
plot2d("ln", w,[-20*log(svs')/log(10) 20*log(svS')/log(10)],[-1 -1 2 2],leg="$\overline{\sigma}(W_S^{-1})$@$\underline{\sigma}(W_S^{-1})$@$\overline{\sigma}(S)$@$\underline{\sigma}(S)$")
xtitle("","Frequency (rad/s)", "Amplitude (dB)");
xgrid(12)
//set(gca(),"auto_clear","off")
xtitle("","Frequency (rad/s)", "Amplitude (dB)");


scf(4);
plot2d("ln", w,[-20*log(svt')/log(10) 20*log(svT')/log(10)],[-1 -1 2 2],leg="$\overline{\sigma}(W_T^{-1})$@$\underline{\sigma}(W_T^{-1})$@$\overline{\sigma}(T)$@$\underline{\sigma}(T)$")
xtitle("","Frequency (rad/s)", "Amplitude (dB)");
xgrid(12)
xtitle("","Frequency (rad/s)", "Amplitude (dB)");

// --------------- Open loop and Closed loop analysis --------------

sysOL=G*K
sysCL=G*K*inv(eye(1,1)+G*K)
// eigenvalues in LHP (Left Half Plane), stable according to RH
spec(sysCL.A)

sv1= svplot(sysOL,w);
sv2= svplot(sysCL,w);
scf(5);
plot2d("ln", w,[20*log(sv1')/log(10) 20*log(sv2')/log(10)],[2 3 5 6],leg="$\overline{\sigma}(GK)$@$\underline{\sigma}(GK)$@$\overline{\sigma}(T)$@$\underline{\sigma}(T)$")
xtitle("","Frequency (rad/s)", "Amplitude (dB)");
xgrid(12)
//set(gca(),"auto_clear","off")
xtitle("","Frequency (rad/s)", "Amplitude (dB)");
