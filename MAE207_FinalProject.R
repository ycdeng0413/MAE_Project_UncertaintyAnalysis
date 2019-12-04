 # we firstly apply modified thompson tau approach to remove the outliers of the time dataset of each balloon
P <- 0.95
# for balloon1
Time_balloon1 <-c(4.92,5.14,4.80,4.97,5.09,5.15,4.78,4.89,4.98,4.86)
N_balloon1 <-length(Time_balloon1)
t_balloon1 <-qt(p=(1+P)/2,df= N_balloon1-2)
tau_balloon1 <- t_balloon1*(N_balloon1-1)/(sqrt(N_balloon1)*sqrt(N_balloon1-2+ t_balloon1^2))
deltaMax_balloon1 <- tau_balloon1*sd(Time_balloon1)
delta_balloon1 <- abs(Time_balloon1-mean(Time_balloon1))
quartz(w=5.5,h=5.5)
plot(delta_balloon1,ylim=c(0,1.1*max(c(delta_balloon1, deltaMax_balloon1))),main='Outlier of time of Balloon1')
abline(h= deltaMax_balloon1,col="red")
# then we subset the outlier for balloon 1
Time_balloon1<-subset(Time_balloon1, delta_balloon1 <deltaMax_balloon1)
# for balloon2
Time_balloon2 <-c(5.47,5.67,5.81,5.74,5.69,5.62,5.53,5.51,5.58,5.49)
N_balloon2 <-length(Time_balloon2)
t_balloon2 <-qt(p=(1+P)/2,df= N_balloon2-2)
tau_balloon2 <- t_balloon2*(N_balloon2-1)/(sqrt(N_balloon2)*sqrt(N_balloon2-2+ t_balloon2^2))
deltaMax_balloon2 <- tau_balloon2*sd(Time_balloon2)
delta_balloon2 <- abs(Time_balloon2-mean(Time_balloon2))
quartz(w=5.5,h=5.5)
plot(delta_balloon2,ylim=c(0,1.1*max(c(delta_balloon2, deltaMax_balloon2))),main='Outlier of time Balloon2')
abline(h= deltaMax_balloon2,col="red")
# then we subset the outlier for balloon 2
Time_balloon2<-subset(Time_balloon2, delta_balloon2 <deltaMax_balloon2)
# for balloon3
Time_balloon3 <-c(5.78,5.46,5.64,5.71,6.27,5.77,5.52,6.02,5.55,5.64)
N_balloon3 <-length(Time_balloon3)
t_balloon3 <-qt(p=(1+P)/2,df= N_balloon3-2)
tau_balloon3 <- t_balloon3*(N_balloon3-1)/(sqrt(N_balloon3)*sqrt(N_balloon3-2+ t_balloon3^2))
deltaMax_balloon3 <- tau_balloon3*sd(Time_balloon3)
delta_balloon3 <- abs(Time_balloon3-mean(Time_balloon3))
quartz(w=5.5,h=5.5)
plot(delta_balloon3,ylim=c(0,1.1*max(c(delta_balloon3, deltaMax_balloon3))),main='Outlier of time Balloon3')
abline(h= deltaMax_balloon3,col="red")
# then we subset the outlier for balloon 3
Time_balloon3<-subset(Time_balloon3, delta_balloon3 <deltaMax_balloon3)

# now we calculate the mean and standard deviation for the time of each balloon
t_balloon1<-mean(Time_balloon1) #unit:s
s_balloon1<-sd(Time_balloon1)
t_balloon2<-mean(Time_balloon2) #unit:s
s_balloon2<-sd(Time_balloon2)
t_balloon3<-mean(Time_balloon3) #unit:s
s_balloon3<-sd(Time_balloon3)

# standard uncertainty of time
u_t_balloon1<- s_balloon1
u_t_balloon2<- s_balloon2
u_t_balloon3<- s_balloon3

# standard uncertainty of height
# the s_height is 0.025 cm and the beta_height is 0.01cm
s_height<- 0.025/100 #unit:m
b_height<- 0.01/2/100 #unit:m
u_height<- sqrt(s_height^2+b_height^2)

# standard uncertainty of diameter of each balloon
# each s_d_balloon is 0.025 cm and each beta_d_balloon is 0.01 cm
s_d_balloon1<- 0.025/100 #unit:m
b_d_balloon1 <- 0.01/2/100 #unit:m
u_d_balloon1 <- sqrt(s_d_balloon1 ^2+ b_d_balloon1 ^2)
u_d_balloon2 <- u_d_balloon1
u_d_balloon3 <- u_d_balloon1

# standard uncertainty of mass of each balloon
# each s_m_balloon is 0 g and each beta_d_balloon is 0.02g
s_m_balloon1<- 0
s_m_balloon2<- 0
s_m_balloon3<- 0
b_m_balloon1<- 0.02/2/1000 #unit:kg
b_m_balloon2<-  0.02/2/1000 #unit:kg
b_m_balloon3<-  0.02/2/1000 #unit:kg
u_m_balloon1 <- sqrt(s_m_balloon1 ^2+ b_m_balloon1 ^2)
u_m_balloon2 <- u_m_balloon1
u_m_balloon3 <- u_m_balloon1

#####Uncertainty Interval of total force#####
### Taylor Series ###
print('Using Taylor Series Approach')
g<-9.8 #gravity acceleration
constant<-0.3039
# for balloon 1, standard uncertainty of the force is
d_exp<- (25/100) #unit:m
h_exp<- (440.5/100) #unit:m
t_exp<-t_balloon1
m_exp<-2/1000 #unit:kg
theta_d<- -constant*2*d_exp*(h_exp^2)/(t_exp^2)
theta_h<- -constant*2*(d_exp^2)*(h_exp)/(t_exp^2)
theta_t<- constant*2*(d_exp^2)*(h_exp^2)/(t_exp^3)
theta_m<-9.8
u_f_total_balloon1<-sqrt(theta_d^2*u_d_balloon1^2+theta_h^2*u_height^2+theta_t^2*u_t_balloon1^2+theta_m^2*u_m_balloon1^2) 
U_f_total_balloon1<-2* u_f_total_balloon1
f_total_balloon1<-m_exp*g-constant*(d_exp^2)*(h_exp^2)/(t_exp^2)
print(paste0('The 95% Uncertainty Interval of balloon 1 (d=25cm): ', f_total_balloon1,'±', U_f_total_balloon1,' (N)'))

# for balloon 2, standard uncertainty of the force is
d_exp<- (28/100) #unit:m
h_exp<- (440.5/100) #unit:m
t_exp<-t_balloon1
m_exp<-2.5/1000 #unit:kg
theta_d<- -constant*2*d_exp*(h_exp^2)/(t_exp^2)
theta_h<- -constant*2*(d_exp^2)*(h_exp)/(t_exp^2)
theta_t<- constant*2*(d_exp^2)*(h_exp^2)/(t_exp^3)
theta_m<-9.8
u_f_total_balloon2<-sqrt(theta_d^2*u_d_balloon2^2+theta_h^2*u_height^2+theta_t^2*u_t_balloon2^2+theta_m^2*u_m_balloon2^2) 
U_f_total_balloon2<-2* u_f_total_balloon2
f_total_balloon2<-m_exp*g-constant*(d_exp^2)*(h_exp^2)/(t_exp^2)
print(paste0('The 95% Uncertainty Interval of balloon 2 (d=28cm): ', f_total_balloon2,'±', U_f_total_balloon2,' (N)'))

# for balloon 3, standard uncertainty of the force is
d_exp<- (30/100) #unit:m
h_exp<- (440.5/100) #unit:m
t_exp<-t_balloon1
m_exp<-3/1000 #unit:kg
theta_d<- -constant*2*d_exp*(h_exp^2)/(t_exp^2)
theta_h<- -constant*2*(d_exp^2)*(h_exp)/(t_exp^2)
theta_t<- constant*2*(d_exp^2)*(h_exp^2)/(t_exp^3)
theta_m<-9.8
u_f_total_balloon3<-sqrt(theta_d^2*u_d_balloon3^2+theta_h^2*u_height^2+theta_t^2*u_t_balloon3^2+theta_m^2*u_m_balloon3^2) 
U_f_total_balloon3<-2* u_f_total_balloon3
f_total_balloon3<-m_exp*g-constant*(d_exp^2)*(h_exp^2)/(t_exp^2)
print(paste0('The 95% Uncertainty Interval of balloon 3 (d=30cm): ', f_total_balloon3,'±', U_f_total_balloon3,' (N)'))

# we now calculate the relative uncertainty
print(paste0('The 95% Relative Uncertainty of balloon 1 (d=25cm): ', U_f_total_balloon1/f_total_balloon1*100,' (%)'))
print(paste0('The 95% Relative Uncertainty of balloon 2 (d=28cm): ', U_f_total_balloon2/f_total_balloon2*100,' (%)'))
print(paste0('The 95% Relative Uncertainty of balloon 3 (d=30cm): ', U_f_total_balloon3/f_total_balloon3*100,' (%)'))

### Monte Carlo Method ###
print('---------------------------------------------------------------------------------------------')
print('Using Monte Carlo Approach')
# for balloon 1, standard uncertainty of the force is
set.seed(100)
M<-1e6
P<-0.95
h_mc<- rnorm(n=M,mean=440.5/100,sd= u_height)
m_mc<- rnorm(n=M,mean=2/1000,sd= u_m_balloon1)
d_mc<- rnorm(n=M,mean=25/100,sd= u_d_balloon1)
t_mc<- rnorm(n=M,mean=t_balloon1,sd= u_t_balloon1)
MC_f_total_balloon1<- m_mc*g-constant*(h_mc^2)*(d_mc^2)/(t_mc)^2
Pmin<- (1-P)/2
Pmax<- (1+P)/2
MC_f_total_balloon1_UI<-quantile(MC_f_total_balloon1,probs=c(Pmin,Pmax))
bar1<- (as.numeric(MC_f_total_balloon1_UI[1])+as.numeric(MC_f_total_balloon1_UI[2]))/2
u1<- bar1-as.numeric(MC_f_total_balloon1_UI[1])
print(paste0('The 95% Uncertainty Interval of balloon 1 (d=25cm): ',bar,'±', u,' (N)'))

# for balloon 2, standard uncertainty of the force is
M<-1e6
P<-0.95
h_mc<- rnorm(n=M,mean=440.5/100,sd= u_height)
m_mc<- rnorm(n=M,mean=2.5/1000,sd= u_m_balloon1)
d_mc<- rnorm(n=M,mean=28/100,sd= u_d_balloon1)
t_mc<- rnorm(n=M,mean=t_balloon1,sd= u_t_balloon1)
MC_f_total_balloon2<- m_mc*g-constant*(h_mc^2)*(d_mc^2)/(t_mc)^2
Pmin<- (1-P)/2
Pmax<- (1+P)/2
MC_f_total_balloon2_UI<-quantile(MC_f_total_balloon2,probs=c(Pmin,Pmax))
bar2<- (as.numeric(MC_f_total_balloon2_UI[1])+as.numeric(MC_f_total_balloon2_UI[2]))/2
u2<- bar2-as.numeric(MC_f_total_balloon2_UI[1])
print(paste0('The 95% Uncertainty Interval of balloon 2 (d=28cm): ',bar,'±', u,' (N)'))

# for balloon 3, standard uncertainty of the force is
M<-1e6
P<-0.95
h_mc<- rnorm(n=M,mean=440.5/100,sd= u_height)
m_mc<- rnorm(n=M,mean=3/1000,sd= u_m_balloon1)
d_mc<- rnorm(n=M,mean=30/100,sd= u_d_balloon1)
t_mc<- rnorm(n=M,mean=t_balloon1,sd= u_t_balloon1)
MC_f_total_balloon3<- m_mc*g-constant*(h_mc^2)*(d_mc^2)/(t_mc)^2
Pmin<- (1-P)/2
Pmax<- (1+P)/2
MC_f_total_balloon3_UI<-quantile(MC_f_total_balloon3,probs=c(Pmin,Pmax))
bar3<- (as.numeric(MC_f_total_balloon3_UI[1])+as.numeric(MC_f_total_balloon3_UI[2]))/2
u3<- bar3-as.numeric(MC_f_total_balloon3_UI[1])
print(paste0('The 95% Uncertainty Interval of balloon 3 (d=30cm): ',bar,'±', u,' (N)'))

# we now calculate the relative uncertainty
print(paste0('The 95% Relative Uncertainty of balloon 1 (d=25cm): ', u1/bar1*100,' (%)'))
print(paste0('The 95% Relative Uncertainty of balloon 2 (d=28cm): ', u2/bar2*100,' (%)'))
print(paste0('The 95% Relative Uncertainty of balloon 3 (d=30cm): ', u3/bar3*100,' (%)'))
#####Balance Check with the law of conservation#####
print('---------------------------------------------------------------------------------------------')
print('Balance Check')
# for balloon 1
Theoretical_balloon1<- 0.002*(2*4.405)/(t_balloon1)^2
abs_balloon1<- abs(Theoretical_balloon1-f_total_balloon1)
print(paste0('for balloon1, abs(F_theorectical-F_experiment), ',abs_balloon1))
print(paste0('for balloon1, expanded uncertainty, ', U_f_total_balloon1))
if (abs_balloon1<U_f_total_balloon1){print('The total force of the balloon 1 is conserved')}else{print('The total force of the balloon 1 is not conserved')}
print('')
# for balloon 2
Theoretical_balloon2<- 0.0025*(2*4.405)/(t_balloon2)^2
abs_balloon2<- abs(Theoretical_balloon2-f_total_balloon2)
print(paste0('for balloon2, abs(F_theorectical-F_experiment), ',abs_balloon2))
print(paste0('for balloon2, expanded uncertainty, ', U_f_total_balloon2))
if (abs_balloon2<U_f_total_balloon2){print('The total force of the balloon 2 is conserved')}else{print('The total force of the balloon 2 is not conserved')}
print('')
# for balloon 3
Theoretical_balloon3<- 0.003*(2*4.405)/(t_balloon3)^2
abs_balloon3<- abs(Theoretical_balloon3-f_total_balloon3)
print(paste0('for balloon3, abs(F_theorectical-F_experiment), ',abs_balloon3))
print(paste0('for balloon3, expanded uncertainty, ', U_f_total_balloon3))
if (abs_balloon3<U_f_total_balloon3){print('The total force of the balloon 3 is conserved')}else{print('The total force of the balloon 3 is not conserved')}
