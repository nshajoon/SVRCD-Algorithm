
library(bnlearn)
library("sparsebn")
library(igraph)
library(discretecdAlgorithm)
library(pcalg)

nSamples=50
nDatasets=20
lambda1=1.2
lambda2=0.9
etha=0.001
innerLoop=50
outerLoop=70
nNodes= 200


gr1=sample_pa(nNodes, p=-3, m=1)

gr1<-edgeList(gr1)
names(gr1)=c(1:nNodes)
func=function(n){runif(n,1,3)}
params=coef_gen(gr1, flip=FALSE)
dat5=generate_discrete_data(graph = gr1, n = nSamples, params = params)
while(1%in%auto_count_levels(dat5))
  dat5=generate_discrete_data(graph = gr1, n = nSamples, params = params)

dat5=list(dat5)
#gr1=list(gr1)

for (i in 2:nDatasets){
  dat5[[i]]=generate_discrete_data(graph = gr1, n = nSamples, params = params)
  while(1%in%auto_count_levels(dat5[[i]]))
    dat5[[i]]=generate_discrete_data(graph = gr1, n = nSamples, params = params)
}



trueDAG=gr1     #bipartite_edgeList2#bipartite_edgeList2#cytometryDiscrete$dag
bngr=to_bn(trueDAG)   #convert edgeList of gr to class bn
cgr=cpdag(bngr)       #cpdag of gr
matcgr=amat(cgr)      #convert cpdag of gr to adjacency matrix


dat=0
data=0
dags=0
error=0.4

dgs=0;

r=rep(0, nDatasets)
e=rep(0, nDatasets)
p=rep(0, nDatasets)
fp=rep(0, nDatasets)
m=rep(0, nDatasets)
tp=rep(0,nDatasets)
fdr=rep(0, nDatasets)
shd=rep(0, nDatasets)
ji=rep(0, nDatasets)

for (pat in 1:nDatasets) {
  
  dat=dat5[[pat]]
  data=sparsebnData(dat, type="d", ivn = NULL)
  #for(i in (1:50)){
  
 dgs=cd.run(data,adaptive=F,lambdas.length =2, error.tol= 0.8, convLb =0.5, upperbound = 50, alpha = 30, weight.scale = 0.001)
 #dgs=cd.run(data)
  
  
  dags=dgs
  
  reverse=rep(0, length(dags))
  expected=rep(0, length(dags))
  s0=nNodes-1
  overall=rep(0, length(dags))
  P=rep(0, length(dags))
  FP=rep(0, length(dags))
  M=rep(0, length(dags))
  TP=rep(0, length(dags))
  FDR=rep(0, length(dags))
  SHD=rep(0, length(dags))
  JI=rep(0, length(dags))
  
  for(k in 1:length(dags))
    P[[k]]=dags[[k]]$nedge
  
  for(k in 1:length(dags)){
    
    bndagTemp=0
    bncpdagTemp=0
    matTemp=0
    
    bndagTemp=to_bn(dags[[k]]$edges)
    bncpdagTemp=cpdag(bndagTemp)
    matTemp=amat(bncpdagTemp)
    
    for(h in 1:length(dags[[k]]$edges))
      for(j in 1:length(dags[[k]]$edges[[h]]))
        if(length(dags[[k]]$edges[[h]])!=0){
          if(!is.null(trueDAG[[dags[[k]]$edges[[h]][[j]]]]))
            if(h %in% trueDAG[[dags[[k]]$edges[[h]][[j]]]]){
              # 
              #           reverse[[k]]=reverse[[k]]+1
              #         }
              #         else{
              
              i_node=dags[[k]]$edges[[h]][[j]]
              j_node=h
              
              
              if((matTemp[i_node,j_node]==1 && matTemp[j_node,i_node]==0 && matcgr[i_node,j_node]==0 && matcgr[j_node,i_node]==1)||
                 (matTemp[i_node,j_node]==1 && matTemp[j_node,i_node]==0 && matcgr[i_node,j_node]==1 && matcgr[j_node,i_node]==1)||
                 #(matTemp[i_node,j_node]==1 && matTemp[j_node,i_node]==1 && matcgr[i_node,j_node]==1 && matcgr[j_node,i_node]==0)||
                 (matTemp[i_node,j_node]==1 && matTemp[j_node,i_node]==1 && matcgr[i_node,j_node]==0 && matcgr[j_node,i_node]==1)
                 # (matTemp[i_node,j_node]==0 && matTemp[j_node,i_node]==1 && matcgr[i_node,j_node]==1 && matcgr[j_node,i_node]==0)||
                 # (matTemp[i_node,j_node]==0 && matTemp[j_node,i_node]==1 && matcgr[i_node,j_node]==1 && matcgr[j_node,i_node]==1))
              )
                reverse[[k]]=reverse[[k]]+1
            }
          
          if(!is.null(trueDAG[[h]]))
            if( (dags[[k]]$edges[[h]][[j]] %in% trueDAG[[h]])||(h%in%trueDAG[[dags[[k]]$edges[[h]][[j]]]]))
              expected[[k]]=expected[[k]]+1
        }
    expected[[k]]=expected[[k]]-reverse[[k]]
    
    
    FP[[k]]=P[[k]]-expected[[k]]-reverse[[k]]    
    M[[k]]=s0-expected[[k]]-reverse[[k]]
    TP[[k]]=expected[[k]]/s0
    FDR[[k]]=(reverse[[k]]+FP[[k]])/P[[k]]
    SHD[[k]]=reverse[[k]]+M[[k]]+FP[[k]]
    JI[[k]]=expected[[k]]/(P[[k]]+s0-expected[[k]])
    
  }
  
  
  
  
  # for (index in 1:length(SHD)){
  #    if(SHD[index]<=min){
  #     min=SHD[index]
  #    ind=index
  # }
  #}*/
  
  #for(ind in (1:5)){
  ind=which.min(SHD)
  p[[pat]]=P[[ind]]
  r[[pat]]=reverse[[ind]]
  e[[pat]]=expected[[ind]]
  fp[[pat]]=FP[[ind]]
  m[[pat]]=M[[ind]]
  tp[[pat]]=TP[[ind]]
  fdr[[pat]]=FDR[[ind]]
  shd[[pat]]=SHD[[ind]]
  
  ji[[pat]]=JI[[ind]]

# cat("P=number of predicted edges= ", P[[ind]],"\n")
# cat("R=number of reverse edges= ", reverse[[ind]],"\n")
# cat("E=number of expected edges= ", expected[[ind]],"\n")
# cat("s0=number of edges in true dag= ", s0,"\n")
# cat("FP=false positive edges= ", FP[[ind]],"\n")
# cat("M=number of missing edges= ", M[[ind]],"\n")
# cat("TP=True positive rate= ", TP[[ind]],"\n")
# cat("FDR=false discovery rate= ", FDR[[ind]],"\n")
# cat("SHD=structural hamming distance= ", SHD[[ind]],"\n")
# cat("JI=jaccard index= ", JI[[ind]],"\n")

 }

mean(p)
mean(r)
mean(e)
mean(fp)
mean(m)
mean(tp)
mean(fdr)
mean(ji)
mean(shd)

var(p)
var(r)
var(e)
var(fp)
var(m)
var(tp)
var(fdr)
var(ji)
var(shd)

