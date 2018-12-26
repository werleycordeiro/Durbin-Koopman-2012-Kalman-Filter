#--------------------------------------------------------------------#
# Illustrations 2.2.5 & 2.4.3  - Durbin & Koopman (2012) P.16 & P.23 #
#--------------------------------------------------------------------# 

# Dados
dados <- read.csv("https://www.dropbox.com/s/ynqjo7k54kdihyx/nile.txt?dl=1")
dadosplot <- dados
dados <- ts(dados,start=c(1871))
len<-length(dados)

# Valores dados no livro
sigma.eps<-15099
sigma.eta<-1469.1

# Filtro de Kalman
kalman1 <- function(dados,sigma.eps,sigma.eta){

# Definindo os vetores
at <-matrix(0,len,1)
Pt <-matrix(0,len,1)
att<-matrix(0,len,1)
Ptt<-matrix(0,len,1)
vt <-matrix(0,len,1)
Ft <-matrix(0,len,1)
K  <-matrix(0,len,1)

# Informação a priori (informações que foram repassadas no enunciado do problema.)
at[1]<- 0
Pt[1]<-1e7
Ft[1]<- Pt[1] + sigma.eps
K[1] <- Pt[1] / Ft[1]

# Informação amostral (isto é, os dados estão entrando no sistema.)
vt[1]  <- dados[1]-at[1]
att[1] <- at[1] + K[1] * vt[1]
Ptt[1] <- Pt[1] * (1-K[1])

	for(i in 1:len)
	{
	#Previsão
	at[i+1]<-at[i] + K[i] * vt[i]
	Pt[i+1]<- Pt[i] * (1-K[i]) + sigma.eta
	
	#Atualização 
	vt[i+1]  <- dados[i+1]-at[i+1]
	Ft[i+1]  <- Pt[i+1] + sigma.eps
	K[i+1]   <- Pt[i+1] / Ft[i+1]
	att[i+1] <- at[i+1] + K[i+1] * vt[i+1]
	Ptt[i+1] <- Pt[i+1] * (1-K[i+1])
	}
return(list(Ptt=Ptt,att=att,vt=vt,Ft=Ft,K=K,at=at,Pt=Pt))
}

results1 <-kalman1(dados,sigma.eps,sigma.eta)

# Intervalo de confiança
upper <- results1$att[1:100] + (1.645*(sqrt(results1$Ptt[1:100]))) 
lower <- results1$att[1:100] - (1.645*(sqrt(results1$Ptt[1:100])))

#----------#
# Fig. 2.1 #
#----------#

par(mfrow=c(2,2))
plot(dadosplot[,1])
lines(results1$att[1:100],col="blue")
lines(upper,col="gray")
lines(lower,col="gray")
ts.plot(results1$Pt[2:101])
ts.plot(results1$vt[2:100])
abline(a = 0, b = 0)
ts.plot(results1$Ft[2:100])

#------------#
# Suavizador #
#------------#

# Realocando os resultados anteriores
P <- results1$Pt
a <- results1$at
F <- results1$Ft 
v <- results1$vt 
K <- results1$K

#definindo os vetores
L <- matrix(0,(len-1),1)
r <- matrix(0,len+1,1) 
alpha <- matrix(0,len,1)
N <- matrix(0,len+1,1)
V <- matrix(0,len,1)

# Computar L, r, alpha, V, N
for(i in 1:len){L[i]<-(1-K[i])}

r[101]<-0 # Há diferença entre a indexação do livro e a do código. Aqui, r_{t} e r_{t-1} estão no mesmo vetor, pois isso r[101]<-0 

for(i in 1:100){
r[101-i] <- (v[101-i] / F[101-i]) + L[101-i] * r[102-i]
}

for(i in 1:100){
alpha[101-i] <- a[101-i] + P[101-i] * r[101-i]
}

N[101]<-0
# N[101] <- 1 / F[100] # N[100] = 0
for(i in 1:100){
N[101-i] <- (1 / F[101-i]) + ((L[101-i])^2) * N[102-i]
}

for(i in 1:100){
V[101-i] <- P[101-i] - ((P[101-i])^2) * N[101-i]
}

# Intervalo de confiança 90%
upper <- alpha[1:100] + (1.645*(sqrt(V[1:100]))) 
lower <- alpha[1:100] - (1.645*(sqrt(V[1:100])))

#----------#
# Fig. 2.2 #
#----------#
par(mfrow=c(2,2))
plot(dadosplot[,1])
lines(alpha[2:100],col='blue')
lines(upper,col='gray')
lines(lower,col='gray')
ts.plot(V)
ts.plot(r)
abline(a = 0, b = 0, h = 9)
ts.plot(N)

