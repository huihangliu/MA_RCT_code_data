rm(list=ls()) # clear variables
library("snowfall") # upload the snowfall package
#############################################################
# the main function defined and used by myself with Case 1
myfun=function(count){
  ##----------package required for the main function----------##
  library("ncvreg") 
  library("quadprog")
  ##----------functions required for the main function----------##
  # nearest neighbor pairing
  nnp=function(Y,delta,X,n){
    til_Y=matrix(0,n,1)
    p=dim(X)[2]
    nt=sum(delta)
    nc=n-nt
    
    # treated group
    Ytobs=matrix(Y[delta>0],nt,1)
    Xtobs=X[delta>0,]
    # control group
    Ycobs=matrix(Y[delta==0],nc,1)
    Xcobs=X[delta==0,]
    
    i=1
    while(i<=n){
      if (delta[i]==1){
        ac=matrix(rep(X[i,],nc),nc,p,byrow=TRUE)-Xcobs
        bc=ac%*%t(ac)
        cc=diag(bc)
        indc=which(cc==min(cc),arr.ind=TRUE)[1]
        til_Y[i]=Y[i]-Ycobs[indc]
      }else{
        at=matrix(rep(X[i,],nt),nt,p,byrow=TRUE)-Xtobs
        bt=at%*%t(at)
        ct=diag(bt)
        indt=which(ct==min(ct),arr.ind=TRUE)[1]
        til_Y[i]=Ytobs[indt]-Y[i]
      }
      i=i+1
    }
    return(til_Y)
  }
  
  # logistic regression model
  KLloss=function(alpha){
    temp=Xpi%*%alpha
    W=exp(temp)/(1+exp(temp))
    KLloss=-(sum(delta*log(W)+(1-delta)*log(1-W)))
  }
  
  # TECV method and TEEM method
  TECV_TEEM=function(Y,X,delta,n,Mn){
    round=100
    val_TECV=matrix(0,round,Mn)
    val_TEEM=matrix(0,round,Mn)
    r=1
    while(r<=round){
      nn=sample(n)
      nn_train=nn[seq(1,n/2,1)]
      nn_test=nn[seq(n/2+1,n,1)]
      
      I1=diag(n/2)
      Y1=matrix(Y[nn_train],n/2,1)
      X1=X[nn_train,]
      delta1=matrix(delta[nn_train],n/2,1)
      Delta1=diag(as.list(delta1))
      
      Y2=matrix(Y[nn_test],n/2,1)
      X2=X[nn_test,]
      delta2=matrix(delta[nn_test],n/2,1)
      
      XMn2=X2[,seq(1,Mn,1)]
      til_Y2=nnp(Y2,delta2,XMn2,n/2)
      m=1
      while(m<=Mn){
        Xm1=t(t(X1[,seq(1,m,1)]))
        Xm2=t(t(X2[,seq(1,m,1)]))
        t_invXXmcc1=solve(t(Xm1)%*%Delta1%*%Xm1)
        t_hmucc1=Xm1%*%t_invXXmcc1%*%t(Xm1)%*%Delta1%*%Y1
        t_hecc1=Delta1%*%(Y1-t_hmucc1)
        t_hsigma_21=(sum(delta1))^(-1)%*%t(t_hecc1)%*%t_hecc1
        t_hmucc2=Xm2%*%t_invXXmcc1%*%t(Xm1)%*%Delta1%*%Y1
        
        c_invXXmcc1=solve(t(Xm1)%*%(I1-Delta1)%*%Xm1)
        c_hmucc1=Xm1%*%c_invXXmcc1%*%t(Xm1)%*%(I1-Delta1)%*%Y1
        c_hecc1=(I1-Delta1)%*%(Y1-c_hmucc1)
        c_hsigma_21=(n/2-sum(delta1))^(-1)%*%t(c_hecc1)%*%c_hecc1
        c_hmucc2=Xm2%*%c_invXXmcc1%*%t(Xm1)%*%(I1-Delta1)%*%Y1               
        
        hsigma_21=as.numeric(t_hsigma_21+c_hsigma_21)
        
        hmucc2=t_hmucc2-c_hmucc2
        hecc=til_Y2-hmucc2
        val_TECV[r,m]=t(hecc)%*%hecc
        val_TEEM[r,m]=prod(exp(-0.5/hsigma_21*(hecc^2)))/sqrt(hsigma_21)

        m=m+1
      }
      val_TEEM[r,]=val_TEEM[r,]/sum(val_TEEM[r,])
      r=r+1
    }
    Val_TECV=colMeans(val_TECV)
    model_TECV=which(Val_TECV==min(Val_TECV),arr.ind=TRUE)[1]
    ww_TEEM=colMeans(val_TEEM)
    return(c(model_TECV,ww_TEEM))
  }
  ##---------- settings for the parameters ----------##
  R=matrix(seq(0.1,0.9,0.1),1,9) # considered population R^2
  n=800
  G=5
  lenR=length(R)

  
  mRiskG_MMA=matrix(0,1,lenR) 
  mRiskI_MMA=matrix(0,1,lenR)
  mRiskTE_MMA=matrix(0,1,lenR)
  mRiskCATE_JMA=matrix(0,1,lenR)
  
  mRiskAIC=matrix(0,1,lenR)
  mRiskSAIC=matrix(0,1,lenR)
  mRiskBIC=matrix(0,1,lenR)
  mRiskSBIC=matrix(0,1,lenR)
  mRiskTECV=matrix(0,1,lenR)
  mRiskTEEM=matrix(0,1,lenR)
  mRiskEW=matrix(0,1,lenR)
  
  # Design-1
  J=10^3
  JJ=matrix(seq(1,J,1),J,1)
  ct=1;t_theta=ct*(1/JJ);t_beta_2=2;t_beta_3=0;
  cc=0.75;c_theta=cc*(1/JJ);c_beta_2=1;c_beta_3=0;
  alphatr=matrix(c(0,0.8,0),3,1)
  varmu=2.0412
  
  # Design-2 所用参数与 matlab-Design-1-probit的相同
  # ct=1;t_theta=ct*(1/JJ);t_beta=1
  # cc=0.75;c_theta=cc*(1/JJ);c_beta=0.5
  # alphatr=matrix(c(0,0.25),2,1)
  # varmu=0.5420
  
  # Design-3 所用参数与 matlab-Design-1-logistic的相同
  # ct=1;t_theta=ct*(1/JJ);t_beta=1
  # cc=0.75;c_theta=cc*(1/JJ);c_beta=0.5
  # alphatr=matrix(c(0,0.4),2,1)
  # varmu=0.5420
  
  # Design-4 所用参数与 matlab-Design-2-logistic的相同
  # ct=1;t_theta=ct*(1/JJ);t_beta=2
  # cc=0.75;c_theta=cc*(1/JJ);c_beta=1
  # alphatr=matrix(c(0,0.5),2,1)
  # varmu=2.0412
  
  # # Design-5
  # J=10^3
  # JJ=matrix(seq(1,J,1),J,1)
  # ct=0.9;t_theta=ct*(1/JJ);t_beta_2=0.3;t_beta_3=-0.3;
  # cc=0.9;c_theta=cc*(1/JJ);c_beta_2=0.3;c_beta_3=0;
  # alphatr=matrix(c(0,0.7,1),3,1)
  # varmu=0.18
  
  # Design-6
  # J=10^3
  # JJ=matrix(seq(1,J,1),J,1)
  # ct=0.9;t_theta=ct*(1/JJ);t_beta_2=0.6;t_beta_3=-0.3;
  # cc=0.9;c_theta=cc*(1/JJ);c_beta_2=0.3;c_beta_3=0;
  # alphatr=matrix(c(0,0.8,1),3,1)
  # varmu=0.36
  
  # Design-7
  # J=10^3
  # JJ=matrix(seq(1,J,1),J,1)
  # ct=1.2;t_theta=ct*(1/JJ);t_beta_2=0;t_beta_3=-0.5;
  # cc=0.9;c_theta=cc*(1/JJ);c_beta_2=-0.5;c_beta_3=0;
  # alphatr=matrix(c(0,0.8,1),3,1)
  # varmu=1
  
  # Design-8
  # J=10^3
  # JJ=matrix(seq(1,J,1),J,1)
  # ct=1.2;t_theta=ct*(1/JJ);t_beta_2=0;t_beta_3=-0.4;
  # cc=0.9;c_theta=cc*(1/JJ);c_beta_2=-0.4;c_beta_3=0;
  # alphatr=matrix(c(0,0.8,0.4),3,1)
  # varmu=0.7
  
  
  Sigma=sqrt((1-R)/R*0.5*varmu)

  I=diag(n)

  RiskG_MMA=matrix(0,G,lenR) #the MSE of G-MMA
  RiskI_MMA=matrix(0,G,lenR) #the MSE of I-MMA
  RiskTE_MMA=matrix(0,G,lenR) #the MSE of Zhao-TEMA
  RiskCATE_JMA=matrix(0,G,lenR) #the MSE of CATE-JMA
  
  RiskAIC=matrix(0,G,lenR) #the MSE of AIC
  RiskSAIC=matrix(0,G,lenR) #the MSE of SAIC
  RiskBIC=matrix(0,G,lenR) #the MSE of BIC
  RiskSBIC=matrix(0,G,lenR) #the MSE of SBIC
  RiskTECV=matrix(0,G,lenR) #the MSE of TECV
  RiskTEEM=matrix(0,G,lenR) #the MSE of TECV
  RiskEW=matrix(0,G,lenR) #the MSE of equal weight
    
  Mn=floor(1.5*n^(1/3)) #the number of candidate models
  ww0=matrix(rep(1,Mn)/Mn,Mn,1) #the initial value for calculating weights
  s=matrix(seq(1,Mn,1),Mn,1) #the dimension of each candidate model
  AMn=cbind(matrix(1,Mn,1),diag(rep(1,Mn)),-1*diag(rep(1,Mn)))
  bMn=rbind(1,matrix(0,Mn,1),-1*matrix(1,Mn,1))
    
  t_tr=matrix(0,1,Mn)
  c_tr=matrix(0,1,Mn)
    
  hmucd=matrix(0,n,Mn)
  hecd=matrix(0,n,Mn)
    
  t_hmucc=matrix(0,n,Mn)
  t_hecc=matrix(0,n,Mn)    
  c_hmucc=matrix(0,n,Mn) 
  c_hecc=matrix(0,n,Mn)
    
  hmucc=matrix(0,n,Mn)
  hecc=matrix(0,n,Mn)
    
  hmuhpi=matrix(0,n,Mn)
  tmuhpi=matrix(0,n,Mn)
  tehpi=matrix(0,n,Mn)
  
  AIC=matrix(0,Mn,1)
  BIC=matrix(0,Mn,1)
    
  X=matrix(0,n,J)
  X[,1]=matrix(1,n,1)
    
  g=1
  while (g<=G){
    X[,seq(2,J,1)]=matrix(rnorm(n*(J-1)),n,J-1)
    t_mu=X%*%t_theta+t_beta_2*matrix(X[,2],n,1)^2+t_beta_3*matrix(X[,3],n,1)^2
    c_mu=X%*%c_theta+c_beta_2*matrix(X[,2],n,1)^2+c_beta_3*matrix(X[,3],n,1)^2
    mu=t_mu-c_mu #the objective CATE
      
      #------generating the treatment indicator variable------#
      Xpi=X[,c(1,2,3)]
      temp=Xpi%*%alphatr
      piX=exp(temp)/(1+exp(temp))
      u=matrix(runif(n),n,1)
      delta=(u<piX)*1
      Delta=diag(as.list(delta))
      
      nt=sum(delta) #sample size of treated group
      nc=n-nt #sample size of control group
      XMn=X[,seq(1,Mn,1)] # the largest model
      
      #------estimating the selection probability function------#
      #--- Case-1 model Mn + SCAD ---#
      cvfit=cv.ncvreg(XMn[,seq(2,Mn,1)],delta,family="binomial",penalty="SCAD")
      summary(cvfit)
      halpha=matrix(coef(cvfit),Mn,1)
      htemp=XMn%*%halpha
      
      #--- Case-2 MLE ---#
      # Re_alpha=optim(par=c(0,0,0),fn=KLloss)
      # halpha=matrix(as.vector(unlist(Re_alpha[1])),3,1)
      # htemp=Xpi%*%halpha
      #------estimating the selection probability function------#
      
      hpiX=exp(htemp)/(1+exp(htemp))
      Deltahpi=diag(as.list(hpiX))
      invDeltahpi=diag(as.list(1/hpiX))
      invIDeltahpi=diag(as.list(1/(1-hpiX)))
      
      r=1
      while (r<=lenR){
        sigma=Sigma[r]
        t_e=matrix(rnorm(n,0,sigma),n,1)
        c_e=matrix(rnorm(n,0,sigma),n,1)
        t_Y=t_mu+t_e
        c_Y=c_mu+c_e
        Y=Delta%*%t_Y+(I-Delta)%*%c_Y
        
        tc_Y=t_Y-c_Y #G-MMA
        til_Y=nnp(Y,delta,XMn,n) #Zhao-TEMA            
        Zhpi=invDeltahpi%*%Delta%*%Y-invIDeltahpi%*%(I-Delta)%*%Y #CATE-JMA-1
        t_Z=Delta%*%Y;c_Z=(I-Delta)%*%Y #CATE-JMA-2
        
        m=1
        while(m<=Mn){
          #G-MMA
          Xm=t(t(X[,seq(1,m,1)]))
          invXXm=solve(t(Xm)%*%Xm)
          Pmcd=Xm%*%invXXm%*%t(Xm)
          hmucd[,m]=Pmcd%*%tc_Y;
          hecd[,m]=tc_Y-hmucd[,m]
          
          #I-MMA
          t_invXXmcc=solve(t(Xm)%*%Delta%*%Xm)
          t_Pmcc=Xm%*%t_invXXmcc%*%t(Xm)%*%Delta  
          t_hmucc[,m]=t_Pmcc%*%Y
          t_hecc[,m]=Delta%*%(Y-t_hmucc[,m])
          
          c_invXXmcc=solve(t(Xm)%*%(I-Delta)%*%Xm)
          c_Pmcc=Xm%*%c_invXXmcc%*%t(Xm)%*%(I-Delta)                
          c_hmucc[,m]=c_Pmcc%*%Y
          c_hecc[,m]=(I-Delta)%*%(Y-c_hmucc[,m])
          
          #TE-MMA
          hmucc[,m]=t_hmucc[,m]-c_hmucc[,m]
          hecc[,m]=til_Y-hmucc[,m]
          t_tr[m]=sum(diag(t_Pmcc))
          c_tr[m]=sum(diag(c_Pmcc))
          
          # CATE-JMA
          Dmcd=diag((1-diag(Pmcd))^(-1))
          tPmcd=Dmcd%*%(Pmcd-I)+I
          tmuhpi[,m]=tPmcd%*%Zhpi
          tehpi[,m]=Zhpi-tPmcd%*%Zhpi              
          hmuhpi[,m]=Pmcd%*%Zhpi
          
          #IC
          hsigma_2m=n^(-1)%*%t(hecc[,m])%*%hecc[,m]
          AIC[m]=n*log(hsigma_2m)+2*s[m]
          BIC[m]=n*log(hsigma_2m)+log(n)*s[m]
          
          m=m+1
        }
        
        #G-MMA
        hsigmacd_Mn=(n-s[Mn])^(-1)*t(hecd[,Mn])%*%hecd[,Mn]
        a1cd=t(hecd)%*%hecd
        a2cd=s%*%hsigmacd_Mn
        scd_MMA=solve.QP(Dmat = a1cd,dvec=-a2cd,Amat = AMn,bvec = bMn,meq = 1)
        wwcd_MMA=t(t(as.vector(unlist(scd_MMA[1]))))
        wwcd_MMA=wwcd_MMA*(wwcd_MMA>0)
        wwcd_MMA=wwcd_MMA/sum(wwcd_MMA)             
        hmucd_MMA=hmucd%*%wwcd_MMA
        #RiskG_MMA[g,r]=sum((mu-hmucd_MMA)^2)
        RiskG_MMA[g,r]=mean((mu-hmucd_MMA)^2)
        
        # I-MMA
        t_hsigma_Mn=(nt-s[Mn])^(-1)%*%t(t_hecc[,Mn])%*%t_hecc[,Mn]
        t_a1cc=t(t_hecc)%*%t_hecc
        t_a2cc=s%*%t_hsigma_Mn
        stwwcc=solve.QP(Dmat = t_a1cc,dvec=-t_a2cc,Amat = AMn,bvec = bMn,meq = 1)
        t_wwcc=t(t(as.vector(unlist(stwwcc[1]))))
        t_wwcc=t_wwcc*(t_wwcc>0)
        t_wwcc=t_wwcc/sum(t_wwcc)
        t_hmucc_MMA=t_hmucc%*%t_wwcc
        
        c_hsigma_Mn=(nc-s[Mn])^(-1)%*%t(c_hecc[,Mn])%*%c_hecc[,Mn]
        c_a1cc=t(c_hecc)%*%c_hecc
        c_a2cc=s%*%c_hsigma_Mn
        scwwcc=solve.QP(Dmat = c_a1cc,dvec=-c_a2cc,Amat = AMn,bvec = bMn,meq = 1)
        c_wwcc=t(t(as.vector(unlist(scwwcc[1]))))
        c_wwcc=c_wwcc*(c_wwcc>0)
        c_wwcc=c_wwcc/sum(c_wwcc)
        c_hmucc_MMA=c_hmucc%*%c_wwcc
        
        hmucc_MMA=t_hmucc_MMA-c_hmucc_MMA
        #RiskI_MMA[g,r]=sum((mu-hmucc_MMA)^2)
        RiskI_MMA[g,r]=mean((mu-hmucc_MMA)^2)
        
        # TE-MMA
        a1cc=t(hecc)%*%hecc
        a2cc=t_hsigma_Mn%*%t_tr+c_hsigma_Mn%*%c_tr
        swwcc=solve.QP(Dmat = a1cc,dvec=-a2cc,Amat = AMn,bvec = bMn,meq = 1)
        wwcc=t(t(as.vector(unlist(swwcc[1]))))
        wwcc=wwcc*(wwcc>0)
        wwcc=wwcc/sum(wwcc)
        hmuTE_MMA=hmucc%*%wwcc 
        #RiskTE_MMA[g,r]=sum((mu-hmuTE_MMA)^2)
        RiskTE_MMA[g,r]=mean((mu-hmuTE_MMA)^2)
        
        # CATE-JMA
        ta1hpi=t(tehpi)%*%tehpi
        ta2hpi=matrix(0,Mn,1)
        swwCATE_JMA=solve.QP(Dmat = ta1hpi,dvec=ta2hpi,Amat = AMn,bvec = bMn,meq = 1)
        wwCATE_JMA=t(t(as.vector(unlist(swwCATE_JMA[1]))))
        wwCATE_JMA=wwCATE_JMA*(wwCATE_JMA>0);
        wwCATE_JMA=wwCATE_JMA/sum(wwCATE_JMA)
        hmuhpiCATE_JMA=hmuhpi%*%wwCATE_JMA
        #RiskCATE_JMA[g,r]=sum((mu-hmuhpiCATE_JMA)^2)
        RiskCATE_JMA[g,r]=mean((mu-hmuhpiCATE_JMA)^2)
        
        #IC
        indAIC=which(AIC==min(AIC),arr.ind=TRUE)[1]
        hmuAIC=hmucc[,indAIC]
        #RiskAIC[g,r]=sum((mu-hmuAIC)^2)
        RiskAIC[g,r]=mean((mu-hmuAIC)^2)
        mAIC=AIC-AIC[indAIC]
        wwAIC=exp(-0.5*mAIC)/sum(exp(-0.5*mAIC))
        hmuSAIC=hmucc%*%wwAIC
        #RiskSAIC[g,r]=sum((mu-hmuSAIC)^2)
        RiskSAIC[g,r]=mean((mu-hmuSAIC)^2)
        
        indBIC=which(BIC==min(BIC),arr.ind=TRUE)[1]
        hmuBIC=hmucc[,indBIC]
        #RiskBIC[g,r]=sum((mu-hmuBIC)^2)
        RiskBIC[g,r]=mean((mu-hmuBIC)^2)
        mBIC=BIC-BIC[indBIC]
        wwBIC=exp(-0.5*mBIC)/sum(exp(-0.5*mBIC))
        hmuSBIC=hmucc%*%wwBIC
        #RiskSBIC[g,r]=sum((mu-hmuSBIC)^2)
        RiskSBIC[g,r]=mean((mu-hmuSBIC)^2)
        
        #TECV and TEEM
        Re_TE=TECV_TEEM(Y,X,delta,n,Mn)
        indTECV=Re_TE[1]
        hmuTECV=hmucc[,indTECV]
        #RiskTECV[g,r]=sum((mu-hmuTECV)^2)
        RiskTECV[g,r]=mean((mu-hmuTECV)^2)
        ww_TEEM=matrix(Re_TE[seq(2,Mn+1,1)],Mn,1)
        hmuTEEM=hmucc%*%ww_TEEM
        #RiskTEEM[g,r]=sum((mu-hmuTEEM)^2)
        RiskTEEM[g,r]=mean((mu-hmuTEEM)^2)
                       
        #EW
        wwEW=matrix(1,Mn,1)/Mn
        hmuEW=hmuhpi%*%wwEW
        #RiskEW[g,r]=sum((mu-hmuEW)^2)
        RiskEW[g,r]=mean((mu-hmuEW)^2)
        
        r=r+1
      }
      g=g+1
    }
    
    mRiskG_MMA=colMeans(RiskG_MMA)
    mRiskI_MMA=colMeans(RiskI_MMA)
    mRiskTE_MMA=colMeans(RiskTE_MMA)
    mRiskCATE_JMA=colMeans(RiskCATE_JMA)
    
    mRiskAIC=colMeans(RiskAIC)
    mRiskSAIC=colMeans(RiskSAIC)
    mRiskBIC=colMeans(RiskBIC)
    mRiskSBIC=colMeans(RiskSBIC)
    mRiskTECV=colMeans(RiskTECV)
    mRiskTEEM=colMeans(RiskTEEM)
    mRiskEW=colMeans(RiskEW)
    
    ReRisk=cbind(mRiskG_MMA,mRiskAIC,mRiskBIC,mRiskSAIC,mRiskSBIC,mRiskTECV,mRiskTEEM,mRiskI_MMA,mRiskTE_MMA,mRiskEW,mRiskCATE_JMA)
}

##-----------------test--------------------##
sfInit(parallel = TRUE, cpus = 20)
countt=20
ReRisk=sfLapply(seq(1,countt,1), myfun)

reRisk=matrix(0,1,99)
count=1
while(count<=countt){
  reRisk=reRisk+matrix(unlist(ReRisk[count]),1,99)
  count=count+1
}
cReRisk=reRisk/countt
sfStop()

rRiskAIC=cReRisk[,seq(9+1,2*9,1)]/cReRisk[,seq(1,9,1)]
rRiskBIC=cReRisk[,seq(2*9+1,3*9,1)]/cReRisk[,seq(1,9,1)]
rRiskSAIC=cReRisk[,seq(3*9+1,4*9,1)]/cReRisk[,seq(1,9,1)]
rRiskSBIC=cReRisk[,seq(4*9+1,5*9,1)]/cReRisk[,seq(1,9,1)]
rRiskTECV=cReRisk[,seq(5*9+1,6*9,1)]/cReRisk[,seq(1,9,1)]
rRiskTEEM=cReRisk[,seq(6*9+1,7*9,1)]/cReRisk[,seq(1,9,1)]
rRiskI_MMA=cReRisk[,seq(7*9+1,8*9,1)]/cReRisk[,seq(1,9,1)]
rRiskTE_MMA=cReRisk[,seq(8*9+1,9*9,1)]/cReRisk[,seq(1,9,1)]
rRiskEW=cReRisk[,seq(9*9+1,10*9,1)]/cReRisk[,seq(1,9,1)]
rRiskCATE_JMA=cReRisk[,seq(10*9+1,11*9,1)]/cReRisk[,seq(1,9,1)]

Re=rbind(rRiskAIC,rRiskBIC,rRiskSAIC,rRiskSBIC,rRiskTECV,rRiskTEEM,rRiskI_MMA,rRiskTE_MMA,rRiskEW,rRiskCATE_JMA)
print(Re)
