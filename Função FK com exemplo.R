#Vamos desenvolver uma função que implementa o filtro de Kalman para a série univariada y_{t} 
#quando F_{t} e H_{t} não são necessariamente identidades.

# -param pega as variâncias; isto é, param = (Q,R) no qual R é um escalar é Q é uma matrix diagonal 
# -y é os dados.
# -optim é um Booleano: se T, então um dos quatro critérios (determinado pelo maxlik e
# 	outofsample) é retornado; caso contrário, o critério bem como o vetor de estado e as variâncias são retornados. 
# -outofsample e maxlik são booleanos que determinam o critério de otimização (
# 	OLS ou ML, in- ou out-of-sample, respectivamente).
# -xi10 e P10 são os valores iniciais para iniciar a recursão do FK,
# 	as dimensões respectivas devem corresponder a Q, isto é, length(param)-1.
# -Fm e H correspondem as matrizes do sistema.
# -if time_varying<-T então a matriz do sistema H é preenchida com os dados
# 	(modelo com parâmetros variando no tempo)
#
KF_gen<-function(param,y,opti,outofsample,maxlik,xi10,P10,Fm,H,time_varying,x_data)
{
	len<-length(y)
	if (length(param)>2)
	{
			Q<-diag(param[1:(length(param)-1)]^2)
		} else
		{
			Q<-as.matrix(param[1]^2)
		}
		R<-param[length(param)]^2
	Pttm1<-array(dim=c(dim(Q),len+1))
		Pttm1[,,1:(len+1)]<-P10
	Ptt<-array(dim=c(dim(Q),len))
		xttm1<-xi10
	logl<-0.
		xittm1<-matrix(nrow=len,ncol=(length(param)-1))
		xittm1[1,]<-xi10
		xitt<-xittm1
# Se time_varying==T, então preenchemos os dados dentro de H (SSM com coeficientes variando no tempo)
	if (time_varying)
	{
		if (is.null(x_data))
		{
# Para um modelo autoregressivo, preenchemos os valores passados de y. 
					H<-c(y[1:dim(Q)[2]])
# Precisamos dos primeiros y[1:p] em H: portanto a primeira equação será para t = p+1
			anf<-dim(Q)[2]+1
		} else
		{
# Para um modelo de regressão, preencha os dados explicativos.
					H<-x_data[1,]
			anf<-1
		}
	} else
	{
		anf<-1
	}
# Recursão do filtro: inicia em i = dim(Q)[2] para um modelo de série temporal
# 	e em i=1 para um modelo de regressão.

		for (i in anf:len) #i<-1 H<-c(1,0) xitt[,2]
		{
# Kalman-Gain
			He<-(H%*%(Pttm1[,,i]%*%H))[1,1]+R
			epshatoutofsample<-y[i]-(H%*%xttm1)[1,1]
			xittm1[i,]<-xttm1
			xtt<-xttm1+Pttm1[,,i]%*%H*epshatoutofsample/He
			epshatinsample<-y[i]-(H%*%xtt)[1,1]
			xitt[i,]<-xtt
			xttm1<-Fm%*%xtt
			Ptt[,,i]<-Pttm1[,,i]-((Pttm1[,,i]%*%H)%*%(H%*%Pttm1[,,i]))/He
			Pttm1[,,i+1]<-Fm%*%Ptt[,,i]%*%t(Fm)+Q
	if (time_varying)
	{
		if (is.null(x_data))
		{
# Para um modelo autoregressivo, preenchemos os valores passados de y em H.
				H<-c(y[i-dim(Q)[2]+1:(dim(Q)[2])])
		} else
		{
# Para um modelo de regressão, preenchemos os dados explicativos em H.
				H<-x_data[min(len,i+1),]
		}
	}
# Aqui nós especificamos o critério de otimização: OLS ou
# ML, in-sample ou out-of-sample.
			if (outofsample)
			{
				if (maxlik)
				{
						logl<-logl+log(He)+epshatoutofsample^2/He
				} else
				{
						logl<-logl+epshatoutofsample^2
				}
			} else
			{
				if (maxlik)
				{
						logl<-logl+log(He)+epshatinsample^2/He
				} else
				{
						logl<-logl+epshatinsample^2
				}
			}
		}
		if (opti)
		{
			return(logl/len)
		} else
	{
			return(list(logl=logl/len,xitt=xitt,xittm1=xittm1,Ptt=Ptt,Pttm1=Pttm1))
	}
}
#
#--------------------------------------------------------------------------------------------------------------------------------
#
#Vamos desenvolver um exemplo de regressão simulada baseado em 
#	y_{t} = b_{0t} + b_{1t}x_{1t} + b_{2t}x_{2t} + b_{3t}x_{3t} + erro_{t}
#	no qual x_{it}, i = 1, 2, 3 são três processos de passeio aleatório independentes e os coeficientes de regressão 
#	b_{it}, i = 0,...,3 também são passeios aleatórias.

len<-100
set.seed(1)
ndim<-4
x_data<-matrix(ncol=ndim,nrow=len)
x_data[,1]<-rep(1,len)
for (i in 1:3)
{
x_data[,i+1]<-cumsum(rnorm(len))
}
# Os coeficientes da regressão são passeios aleatórios também.
beta<-cbind(cumsum(rnorm(len)),cumsum(rnorm(len)),cumsum(rnorm(len)),cumsum(rnorm(len)))
y<-apply(x_data*beta,1,sum)+rnorm(len)


# Estime os coeficientes desconhecidos por OLS (função 'lm') e retornar as estimativas.

lm_obj<-lm(y~x_data-1)
summary(lm_obj)$coef


# Aplique a função KF_gen acima para determinar os coeficientes variando no tempo.

#Otimize Q e R = sigma^2. Use a priori difusa P10<−diag (rep (1000000, ndim)), 
#	defina um valor de inicialização arbitrário para as variáveis desconhecidas 
#	param <−rep (var (y), ndim + 1), defina opti <-T, fora da amostra <- T, maxlik <-T 
#	(otimização da máxima verossimilhança fora da amostra); use xi_{1|0} = 0 e F = Id.

# inicialização arbitrátia
param<-rep(var(y),ndim+1)
xi10<-rep(0,ndim)
# priori difusa
P10<-diag(rep(1000000,ndim))
Fm<-diag(rep(1,ndim))
# otimização numérica
time_varying<-T
opti<-T
outofsample<-T
maxlik<-T
objopt<-nlminb(start=param,objective=KF_gen,
y=y,opti=opti,outofsample=outofsample,maxlik=maxlik,xi10=xi10,P10=P10,Fm=Fm,
H=H,time_varying=time_varying,x_data=x_data)
param<-objopt$par


# Mostrar o vetor de parâmetros otimizado parma^2 correspondente às variâncias
# 	em Q e R; estimar os coeficientes de regressão resultante e fazer um gráfico do vetor de estado ótimo.

param^2
opti<-F
# Determinar o vetor de estado para o parâmetro 'param' ótimo
obj<-KF_gen(param,y,opti,outofsample,maxlik,xi10,P10,Fm,H,time_varying,x_data)
ymin<-min(apply(obj$xitt,1,min))
ymax<-max(apply(obj$xitt,1,max))
par(mfrow=c(2,1))
ts.plot(obj$xitt[,1],ylim=c(ymin,ymax),col="blue",main=paste("Coef. estimados (linha sólida) vs. coef. observados",sep=""),xlab="",ylab="")
colo<-c("blue","red","green","brown")
for (i in 2:ndim)
lines(obj$xitt[,i],col=colo[i])
for (i in 1:ndim)
lines(beta[,i],col=colo[i],lty=2)
ymin<-min(apply(beta,1,min))
ymax<-max(apply(beta,1,max))
ts.plot(rep(summary(lm_obj)$coef [1,1],len),ylim=c(ymin,ymax),col="blue",main=paste("OLS (coeficientes fixos) vs. coef. observados",sep=""),xlab="",ylab="")
for (i in 2:ndim)
lines(rep(summary(lm_obj)$coef [i,1],len),col=colo[i])
for (i in 1:ndim)
lines(beta[,i],col=colo[i],lty=2)
dev.off()


