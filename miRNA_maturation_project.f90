!!! developer: Raunak Adhikary 
!!! project name: miRNA maturtion project 
!!! authors: Raunak Adhikary, Dipjyoti Das
!!! This code generates mRNA, protein, pre-miRNA, and mature-miRNA numbers corresponding 
!!! to the bound and unbound state of the miRNA gene at any time 

program gillespie_miRNA_maturation
implicit none 

!!! define the variables -------------------------------------------

double precision, dimension(13) :: c
double precision, dimension(11) :: a
integer :: react_count,mu,nu,seed
double precision :: r,p,s,m,g_b,T,T2,TINT,TPRINT,summ,R1,R2,R2A0,A0,tao,steady_t
double precision :: r_u,r_b,p_u,p_b,s_u,s_b,m_u,m_b
double precision :: run_t1,run_t2

!Character(len=40):: path

call cpu_time(run_t1)     

!!!! define initial values and reaction constants -------------------------------------------

T=0.0d0; 	T2=550000.0d0;               	!!! T = time that will evolve, T2 = stopping time 
TINT=0.125d0; 	steady_t = 50000.0d0		!!! TINT = interval time, steady_t = starting time of steady-state window


!!!! define initial values mRNA (r), protein (p), pre-miRNA (s), mature miRNA (m) and miRNA gene state (g_b) ------------
!!! r_u = initial value of mRNA in 'unbound' state, r_b = initial value of mRNA in 'bound' state. 
!!! Similar for other variables. 

r_u=0.0d0; r_b=0.0d0; r=0.0d0; p_u=0.0d0; p_b=0.0d0; p=0.0d0		 
s_u=0.0d0; s_b=0.0d0; s=0.0d0; m_u=0.0d0; m_b=0.0d0; m=0.0d0; 

g_b=0.0d0					!!!g_b = 0 => miRNA-gene is 'unbound', g_b = 1 => miRNA-gene is 'bound'

react_count = 0					!!! counts the number of reaction in time interval [0,T2]

!!!! define reaction constants --------------------------------------------------------------

 c(1)= 0.100d0                               	!!!c(1)= (k_r) mRNA transcription rate 
 c(2)= 0.0004d0                               	!!!c(2)= (g_r) mRNA degradation rate
 c(3)= 0.017d0                                	!!!c(3)= (k_p) protein translation rate
 c(4)= 0.00002d0                              	!!!c(4)= (g_p) protein degradation rate
 c(5)= 0.001d0                              	!!!c(5)= (k_b) miRNA gene protein binding rate
 c(6)= 0.01d0                               	!!!c(6)= (k_u) miRNA gene protein unbinding rate
 c(7)= 0.5d0                                  	!!!c(7)= (k_s_0) miRNA gene transcription rate in 'unbound' state
 c(8)= 0.05d0                                 	!!!c(8)= (k_s) miRNA gene transcription rate in 'bound' state
 c(9)= 0.00017d0                             	!!!c(9)= (g_s) pre-miRNA degradation rate
 c(10)= 0.04d0                               	!!!c(10)= (gamma_s) mRNA-pre-miRNA co-degradation rate
 c(11)= 1.0d0                              	!!!c(11)= (k_mat) pre-miRNA maturation rate 
 c(12)= 0.000017d0                            	!!!c(12)= (g_m) mature miRNA degradation rate
 c(13)= 0.04d0                               	!!!c(13)= (gamma_m) mRNA-mature-miRNA co-degradation rate

!!!! assign value of T to TPRINT --------------------------------------------------------------

TPRINT=T

!!! open a data file and write the values ------------------------------------------------------

!path='k_mat=1.0/k_r=0.34500/'
 
!open(2,File=trim(adjustl(path))//'parameter_values.txt', status='unknown',position='append')
open(2,File='parameter_values.txt', status='unknown',position='append')

!open(3,File=trim(adjustl(path))//'output.dat', status='unknown',position='append')
open(3,File='output.dat', status='unknown',position='append')

write(2,5) c(1),c(2),c(3),c(4),c(5),c(6),c(7),c(8),c(9),c(10),c(11),c(12),c(13),steady_t,T2,TINT
5 format("##","k_r=",1X,F15.10,1X,"; g_r=",1X,F15.10,1X,"; k_p=",1X,F15.10,1X,"; g_p=",1X,F15.10,1X, &
"; k_b=",1X,F15.10,1X,"; k_u=",1X,F15.10,1X,"; k_s_0=",1X,F15.10,1X,"; k_s=",1X,F15.10,1X, &
"; g_s=",1X,F15.10,1X,"; gamma_s =",1X,F15.10,1X,"; k_mat =",1X,F15.10,1X,"; g_m =",1X,F15.10,1X, &
"; gamma_m =",1X,F15.10,'; steady_t = ',1X,F12.4,1X, '; T2 = ',1X,F12.4,1X, 'TINT = ',1X,F12.7,1X,'; run = 5 ')
flush(2)

write(2,*)'## (1)TPRINT   (2)r_u   (3)r_b   (4)r   (5)p_u   (6)p_b   (7)p   (8)s_u   (9)s_b   (10)s   &
		(11)m_u   (12)m_b   (13)m    (14)g_b    (15)reaction_count'
flush(2)

!!!! call random number generator -------------------------------------------------------------------

   call SYSTEM_CLOCK(COUNT=seed)          !!! SYSTEM_CLOCK is a inbuild fortran subroutine to call a random number generator 
                                          !!!    and      'count=seed'
                                          !!! gives the seed for the random number generator srand
   call srand(seed)

!!!! find values of a(i) as functions of c(i) and state of the system ----------------------------------------------------

10 a(1) = c(1)                               	!!! a(1) = k_r
   a(2) = c(2)*r                             	!!! a(2) = r*g_r
   a(3) = c(3)*r                             	!!! a(3) = r*k_p
   a(4) = c(4)*p                             	!!! a(4) = p*g_p

   if (g_b.EQ.0.0d0) then 			!!! g_b = 0 means 'unbound' and g_b = 1 means 'bound'
           a(5) = c(5)*p                    	!!! a(5) = p*k_b
   else                                     	!!! or
           a(5) = c(6)                      	!!! a(6) = k_u          
   end if 

   if (g_b.EQ.0.0d0) then			!!! when g_b = 0 i.e. 'unbound' miRNA transcription rate = k_s_0
           a(6) = c(7)                      	!!! a(6) = k_s_0
   else                                     	!!! or g_b = 1 i.e. 'bound' miRNA transcription rate = k_s
           a(6) = c(8)                      	!!! a(6) = k_s           
   end if 
   
   a(7) = c(9)*s                            	!!! a(7) = s*g_s
   a(8) = c(10)*r*s                         	!!! a(8) = r*s*gamma_s
   a(9) = c(11)*s                           	!!! a(9) = s*k_mat
   a(10) = c(12)*m                          	!!! a(10) = m*g_m
   a(11) = c(13)*r*m                        	!!! a(11) = r*m*gamma_m

   a0 = a(1)+a(2)+a(3)+a(4)+a(5)+a(6)+a(7)+a(8)+a(9)+a(10)+a(11)
   
!!!! calculate R1 & R2 by random number generator ---------------------------------------------------

20  R1 = rand(0)
    R2 = rand(0)
    
    if (R1.LT.1E-30.OR.R2.LT.1E-30) then
       go to 20
    end if

!!!! calculate tau and increment time T by tau --------------------------------------------------------------------
   
   tao = dlog(1/R1)/a0

21 T = T + dlog(1/R1)/a0
    
22 if (T.LT.TPRINT) then
   go to 25
   end if
   
23 if (TPRINT.GE.steady_t) then

	write(3,9) TPRINT,r_u,r_b,r,p_u,p_b,p,s_u,s_b,s,m_u,m_b,m,g_b
9 	format(F20.10,1X,F30.10,1X,F30.10,1X,F30.10,1X,F30.10,1X,F30.10,1X,F30.10,1X,F30.10,1X,F30.10,1X,F30.10,1X, &
	F30.10,1X,F30.10,1X,F30.10,1X,F30.10)
	flush(3)

   end if

!!! increase TPRINT --------------------------------------------------------------------------------------------------

 if (TPRINT.GE.T2) go to 46
 
  TPRINT=TPRINT + TINT
  go to 22

!!! calculate mu to find which reaction will occur ------------------------------------------ 

25 R2A0=R2*A0
   summ=0

26 do nu=1,11
      mu=nu
      summ=summ+ a(nu)
      if(summ.GE.R2A0) go to 30
   end do

!!! to calculate the state of the system --------------------------------------------

30 go to(31,32,33,34,35,36,37,38,39,40,41), mu

31 r=r+1.0d0                                  	!!!--------- Reaction: mRNA coding gene -> mRNA (k_r)  ----------
   react_count= react_count + 1
   if (g_b.EQ.0.0d0) then
        r_u = r; p_u = p; s_u = s; m_u = m
        r_b = 0.0d0; p_b = 0.0d0; s_b = 0.0d0; m_b = 0.0d0
   end if
   if (g_b.EQ.1.0d0) then
        r_b = r; p_b = p; s_b = s; m_b = m
        r_u = 0.0d0; p_u = 0.0d0; s_u = 0.0d0; m_u = 0.0d0
   end if
   go to 45
   
32 r=r-1.0d0                                  	!!!-------- Reaction: mRNA -> degraded mRNA (g_r)  ----------------
   react_count= react_count + 1
   if (g_b.EQ.0.0d0) then
        r_u = r; p_u = p; s_u = s; m_u = m
        r_b = 0.0d0; p_b = 0.0d0; s_b = 0.0d0; m_b = 0.0d0
   end if
   if (g_b.EQ.1.0d0) then
        r_b = r; p_b = p; s_b = s; m_b = m
        r_u = 0.0d0; p_u = 0.0d0; s_u = 0.0d0; m_u = 0.0d0
   end if
   go to 45

33 p=p+1.0d0                                   	!!!-------- Reaction: mRNA -> mRNA + protein (k_p)  ----------------
   react_count= react_count + 1
   if (g_b.EQ.0.0d0) then
        r_u = r; p_u = p; s_u = s; m_u = m
        r_b = 0.0d0; p_b = 0.0d0; s_b = 0.0d0; m_b = 0.0d0
   end if
   if (g_b.EQ.1.0d0) then
        r_b = r; p_b = p; s_b = s; m_b = m
        r_u = 0.0d0; p_u = 0.0d0; s_u = 0.0d0; m_u = 0.0d0
   end if
   go to 45
   
34 p=p-1.0d0                                    !!!-------- Reaction: protein -> degraded protein (g_p)  -------------
   react_count= react_count + 1
   if (g_b.EQ.0.0d0) then
        r_u = r; p_u = p; s_u = s; m_u = m
        r_b = 0.0d0; p_b = 0.0d0; s_b = 0.0d0; m_b = 0.0d0
   end if
   if (g_b.EQ.1.0d0) then
        r_b = r; p_b = p; s_b = s; m_b = m
        r_u = 0.0d0; p_u = 0.0d0; s_u = 0.0d0; m_u = 0.0d0
   end if
   go to 45
     
35 if (g_b.EQ.0.0d0) then
   
   p=p-1.0d0
   g_b=g_b+1.0d0 				!!!---- Reaction: protein + miRNA gene -> miRNA gene protein complex (k_b) ----
   react_count=react_count + 1
   
   r_b = r; p_b = p; s_b = s; m_b = m
   r_u = 0.0d0; p_u = 0.0d0; s_u = 0.0d0; m_u = 0.0d0
   go to 45

   else
 
   p=p+1.0d0                                	!!!---- Reaction: miRNA gene protein complex -> protein + miRNA gene (k_u) ----
   g_b=g_b-1.0d0
   react_count = react_count + 1
   
   r_u = r; p_u = p; s_u = s; m_u = m
   r_b = 0.0d0; p_b = 0.0d0; s_b = 0.0d0; m_b = 0.0d0
   go to 45

   end if

36 if (g_b.EQ.0.0d0) then

   s=s+1.0d0                                  	!!!------- Reaction: miRNA-gene -> s (k_s_0)  ------------------
   react_count=react_count + 1
   
   r_u = r; p_u = p; s_u = s; m_u = m
   r_b = 0.0d0; p_b = 0.0d0; s_b = 0.0d0; m_b = 0.0d0
   go to 45

   else
 
   s=s+1.0d0                               	!!!------- Reaction: miRNA-gene -> s (k_s)  --------------------
   react_count=react_count + 1
   
   r_b = r; p_b = p; s_b = s; m_b = m
   r_u = 0.0d0; p_u = 0.0d0; s_u = 0.0d0; m_u = 0.0d0
   go to 45 
 
   end if
   
37 s=s-1.0d0                                   	!!!------- Reaction: miRNA -> degraded miRNA (g_s) ---------------
   react_count=react_count +1
   if (g_b.EQ.0.0d0) then
        r_u = r; p_u = p; s_u = s; m_u = m
   	r_b = 0.0d0; p_b = 0.0d0; s_b = 0.0d0; m_b = 0.0d0
   end if
   if (g_b.EQ.1.0d0) then
        r_b = r; p_b = p; s_b = s; m_b = m
   	r_u = 0.0d0; p_u = 0.0d0; s_u = 0.0d0; m_u = 0.0d0
   end if
   go to 45
   
38 r=r-1.0d0
   s=s-1.0d0                                   	!!!------- Reaction: mRNA + pre-miRNA -> decay (gamma_s)  --------
   react_count=react_count +1
   if (g_b.EQ.0.0d0) then
        r_u = r; p_u = p; s_u = s; m_u = m
   	r_b = 0.0d0; p_b = 0.0d0; s_b = 0.0d0; m_b = 0.0d0
   end if
   if (g_b.EQ.1.0d0) then
        r_b = r; p_b = p; s_b = s; m_b = m
   	r_u = 0.0d0; p_u = 0.0d0; s_u = 0.0d0; m_u = 0.0d0
   end if
   go to 45
   
39 s=s-1.0d0
   m=m+1.0d0
   react_count=react_count +1			!!!------- Reaction: pre-miRNA -> mature miRNA (k_mat)  --------------
   if (g_b.EQ.0.0d0) then
        r_u = r; p_u = p; s_u = s; m_u = m
   	r_b = 0.0d0; p_b = 0.0d0; s_b = 0.0d0; m_b = 0.0d0
   end if
   if (g_b.EQ.1.0d0) then
        r_b = r; p_b = p; s_b = s; m_b = m
   	r_u = 0.0d0; p_u = 0.0d0; s_u = 0.0d0; m_u = 0.0d0
   end if
   go to 45
   
40 m=m-1.0d0
   react_count=react_count +1			!!!------ Reaction: mature-miRNA -> decay (g_m)  --------------
   if (g_b.EQ.0.0d0) then
        r_u = r; p_u = p; s_u = s; m_u = m
   	r_b = 0.0d0; p_b = 0.0d0; s_b = 0.0d0; m_b = 0.0d0
   end if
   if (g_b.EQ.1.0d0) then
        r_b = r; p_b = p; s_b = s; m_b = m
   	r_u = 0.0d0; p_u = 0.0d0; s_u = 0.0d0; m_u = 0.0d0
   end if
   go to 45
   
41 r=r-1.0d0
   m=m-1.0d0
   react_count=react_count +1			!!!------ Reaction: mRNA + mature-miRNA -> decay (gamma_m)  ---------
   if (g_b.EQ.0.0d0) then
        r_u = r; p_u = p; s_u = s; m_u = m
   	r_b = 0.0d0; p_b = 0.0d0; s_b = 0.0d0; m_b = 0.0d0
   end if
   if (g_b.EQ.1.0d0) then
        r_b = r; p_b = p; s_b = s; m_b = m
   	r_u = 0.0d0; p_u = 0.0d0; s_u = 0.0d0; m_u = 0.0d0
   end if
   go to 45
   
45 if (T.LT.T2) then 
   go to 10                                  
   end if

46 call cpu_time(run_t2)
write(2,*) '# ---------------------END OF SIMULATION run = 1 -------------------- time taken =',(run_t2-run_t1)
write(2,*)

end program gillespie_miRNA_maturation
