{ Fill in the following sections (removing comment marks ! if necessary),
  and delete those that are unused.}
TITLE 'Spatial_Solow_Model'     { the problem identification }
COORDINATES cartesian1  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  u (threshold=0.01)
  v (threshold=0.01)
SELECT         { method controls }

!ERRLIM =1e-5
!REGRID   = off
!Ngrid=1000

DEFINITIONS    { parameter definitions }

theta=0.7
eps=1e-3
epsg=eps^(2*(1-theta)/(2-theta))

A=2*GAMMAF(1+1/(1-theta))*GAMMAF(1/2+2/(1-theta))/GAMMAF(1/2+1/(1-theta))/GAMMAF(1+2/(1-theta)) ! Asymptotic value of spike height
l=1  ! Domain length
n=1 ! Number of spike
ll=l/n; ! distance between spikes
eta=(2/(1-theta)*A^(1-theta)*(cosh(ll)/sinh(ll))^theta)^(2/(2-theta)) ! Parameter defined in the paper
k0=eps^(2/(2-theta))*A^(1-theta)*(cosh(ll)/sinh(ll))^(2/(2-theta)) !  Asymptotic value of Capital 
y=(1-theta)*sqrt(eta)/2*x/eps^( (2-2*theta)/(2-theta) ) ! Inner space variable
u_s=A*(cosh( y ) )^(2/(theta-1))  ; !Asymptotic solution of Labor
v_s=eps^(2/(2-theta))*(sqrt(eta)/sinh(ll) )*cosh((x-ll))   +6.02*eps^2*cosh(x-l)-2*eps^2/(1-theta)*ln( 2*cosh(y)  )  +eps^(2/(2-theta))*sqrt(eta)*x; ! Asymptotic solution of Captial

INITIAL VALUES 


u=exp(-50*x^2)
v=2*1e-4


EQUATIONS        { PDE's, one for each variable }
 u:0.01*dt(u)=eps^2*div(grad(u))-dx(u*dx(v))+eps^2*u*(1-u)
 v:dt(v)=div(grad(v))-v+(abs(u)^(1-theta))*(abs(v)^theta)

! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  REGION 1       { For each material region }
 start(0) line to (l) 
     !point periodic(x-l)


 
 
 TIME 0 TO 300000  { if time dependent }
MONITORS         { show progress }
for cycle=10
  elevation (u,u_s)  from (0) to (l) 
  elevation(v,v_s)  from (0) to (l)
PLOTS            { save result displays }

for time =300000
transfer(u,v) file="sstheta0d7_0d001_bd.dat"

END
