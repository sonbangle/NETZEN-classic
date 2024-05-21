#!/usr/bin/env Rscript

#####################################################################################################




#score at each target divided by degree
# add option to include absolute expression level as a weight. (CellNet stat)
#use MI FDR as a weight (Viper stat)
# Leverage the Transcription Factor consensus sequence. Find the consensus Transcription Factors by deep learning . Bayesien framework.
#Leverage the patient survival.
#Let the iteration go to convergence.Objective function: sum of rank differencies between round
#Leverage the RNA seq isoforms.
#Deep learning on input: mutation (WGS), CNV, methylation  , network , transcription target consensus sequence, output:  survival, drug response
#Deep learning: mutation in non conding region -- effect on TF binding, histon mark, , expression, DNAse I hypersensitive sites
# out degree of the direct parent in propagation direction. Important.


#master_genes=c("SOX8","KLF9"),step=c(1,2),stat=c("logFC","pvalue","fdr","LR"),nround_iterative=c(1,2),weight=c("rho","MI","unw"),weight_power=c(1,2)
#library(IDPmisc)

# # print some progress messages to stderr if "quietly" wasn't requested
# if ( opt$verbose ) {
#   write("writing some verbose output to standard error...\n", stderr())
# }
# # do some operations based on user input
# if( opt$generator == "rnorm") {
#   cat(paste(rnorm(opt$count, mean=opt$mean, sd=opt$sd), collapse="\n"))
# } else {
#   cat(paste(do.call(opt$generator, list(opt$count)), collapse="\n"))
# }

print_v<-function(x)
{
  if (opt$verbose)print(paste(x, "at", Sys.time()))
}





#######################################################################################################
#initiation++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(igraph)

#library(profvis)
library(data.table)
#require(compiler)
#enableJIT(3)
#library(Rmpi)



#function to trim white space of string. x: string, can be a list, or vector
trim <- function (x) 
{
  result<-sapply(x, function(t) gsub("^\\s+|\\s+$", "", t))
  return(as.vector(result))              
}


pvalue_zscore<-function(pvalue)#convert pvalue to z score assuming normal distribution, pvalue is right tail
{
  #print("line 160, pvalue zscore")
  #print(paste("pvalue before z score",pvalue))
  pvalue[pvalue < 1e-100] = 1e-100

  pvalue_zscore<-qnorm(pvalue, lower.tail=F)
  #print(paste("pvalue zscore:",pvalue_zscore))
  pvalue_zscore[pvalue_zscore<0] <- 0.0001
  pvalue_zscore[is.infinite(pvalue_zscore)]<-0.0001
  print("done with pvalue z score")
  return(pvalue_zscore)
}

fdrvalue_confidence<-function(fdrvalue)
{
  
  return(1-fdrvalue)
}



#nodes<-data.frame(name=V(cn)$name)
#nodes$ID<-rownames(nodes)
#ego(cn,order=3,sample(1:17528,1))# test the number of steps needed to visit all nodes.
#Calculate shortest path from a source node to every other node, using inferred correlation

MI_to_rho<-function(MI)#convert MI to correlation, assuming the distribution is normal.
{
  rho<-sapply(MI, function(x)abs((1-exp(-2*x))^0.5))
  return(rho)
}


rho_to_MI<-function(rho)#convert correlation to MI, assuming the distribution is normal.
{
  MI<-sapply(rho, function(x)-log((1-x^2))/2)
  return(MI)
}





#calculate the estimated correlation between a node and all other nodes .
# distance_rho= r1*r2*r3*..., where r1, r2, r3.. are the rho of each step in the shortest path from source node to all other nodes in network



#return product of correlations based on correlation and the path as a list
#Path: shortest path between two nodes
distance_rho<-function(net=cn,paths) 
{
  distance_rho<-sapply(paths,function(p)
  {
    
    rhos<-E(net,path=p[[1]][[1]] ,directed = F)$rho
    prod_rho<-prod(rhos)
    return(prod_rho)
  }
  )
  return(distance_rho)
}


# return vertexes in n steps from the source node, shortest path from 
#neighborhood

#given source node, number of step, return :
#a) nbh_step list of target genes at each step, not including previous step
#b) nbh_path<-list of shortest path from target nodes at each step to the source node
#c) nbh_distance_rho -list of product of correlations of shortest path from source node to target nodes at each step
#output:n-neigbors nodes name only
#       np -nodes and path to nodes
#-------a  all: nodes, path to nodes, distance rho
distance_rho_step<-function(v1="KLF9", net=cn,step=2,output="n") 
{
  #print_v(v1)
  nbh<-list()#list of target genes that canbe reach from   the source in n step or less
  
  nbh[[1]]<-v1
  nbh_path<-list() #list of shortest path from target nodes at each step to the source node
  
  #find neighborhood 
  #nbh[[i+1]]<-names(unlist(ego(net,order=i, nodes = v1, mode="all")))
  
  if (step>1) 
  {
    subnet_previous<-make_ego_graph(net,order = step-1, node=v1, mode = "all")[[1]]
    nbh[[step]]<-V(subnet_previous)$name
  }
  subnet<-make_ego_graph(net,order = step, node=v1, mode = "all")[[1]]
  nbh[[step+1]]<-V(subnet)$name
  nbh_step<-setdiff(nbh[[step+1]],nbh[[step]])
  
  if (length(nbh_step)==0){return(NULL)} else
  {    
    if (output=="a")
    {
      if (step>1)
      {
        
        nbh_path<-lapply(nbh_step,function(target){shortest_paths(subnet, from=v1, to=target, mode="all", weights = E(subnet)$distance,output="vpath")})  
        nbh_distance_rho<-distance_rho(subnet,paths=nbh_path)
      }else
      {
        P<-vector() # making list of paths from source to target node when step is 1 (direct target) 
        for (j in 1:length(nbh_step))
        {
          P<-c(P,v1,nbh_step[[j]])
        }
        #nbh_distance_rho[[i]]<-E(net,P)$rho
        
        nbh_distance_rho<-E(subnet,P)$rho
        
      }
      
      names(nbh_distance_rho)<-nbh_step
    }
    
    
    
    
    if (output=="n")return(list(nbh_step))  
    # if (output=="np")return(list(nbh_step,nbh_path))  
    if (output=="a") return(list(nbh_step,nbh_path,nbh_distance_rho))
  }
  
  
}


# q<-profvis({
# nodes<-V(cn)$name
# newnodes<-sample(nodes, size = 1,replace=T)
# t<-lapply(newnodes,function(node)distance_rho_step(v1=node,step=2, net=cn,output = "a"))  
# 
# })






#################################################################################################################################
#rank the genes by sum of input gene importance score and step.
#step: number of step for calculating of local neigborhood
#weight: rho, MI, rho2,MI2 unw: unweighted
#weight_power: weight^power, for example weight square
#input_scores: data frame with format (gene,weight1, weight2,...weightn) new methods allowing multiple weights for a gene.
#option: divide by degree or not, influence of other nodes
#neighbor_aggregation_method: use average or sum of target scores  to calculate master score. If use sume, the node with high connection will always scores high.
master_score_step<-function(step=opt$step,input_scores=ges,weight=opt$edge_weight,weight_power=opt$weight_power,source_node_inclusion=opt$source_node_inclusion,neighbor_aggregation_method=opt$neighbor_aggregation_method,step_power=opt$step_power,steps_combined=opt$steps_combined,step_normalization=opt$step_normalization,input_normalization=opt$input_normalization,mpi=F,nslaves=10, cluster=F)
{
  
  colnames(input_scores)[1]<-c("gene")
  # for (i in 2:ncol(input_scores))
  # {input_scores[,i]<-scale_log(input_scores[,i])} #normalize all input scores.
  
  input_scores$source_node_score<-apply(data.frame(input_scores[,2:ncol(input_scores)]),1,prod)# aggreagate input  gene node statistics by multiplication
  if (input_normalization) input_scores$source_node_score<-scale_log(input_scores$source_node_score)  # normalize input score
  step_weight<-data.frame(nstep=c(1:step))
  step_weight$weight<-1/((step_weight$nstep)^step_power)
  
  #_________________________________________________________________________________________________________________________
  scores_df<-input_scores
  
  st<-1
  if (!steps_combined)  st<-step
  while(st<= step)
  {  
    
    if (!steps_combined)  st<-step # if steps scores not combined together, then only calculate score for the last step
    
    #running mpi
    if (mpi)
    {
      if (!is.loaded("mpi_initialize")) {library(Rmpi)}
      print_v("beginning mpi job")
      
      if (cluster)
      { nslaves <- mpi.universe.size() - 1 # only for cluster
      
      }
      
      mpi.spawn.Rslaves(nslaves=nslaves)# for multicore desktop
      mpi.bcast.cmd(scores_df<-mpi.bcast.Robj())
      mpi.bcast.Robj(scores_df) 
      mpi.bcast.cmd(st<-mpi.bcast.Robj())
      mpi.bcast.Robj(st)
      mpi.bcast.cmd(input_scores<-mpi.bcast.Robj())
      mpi.bcast.Robj(input_scores)
      mpi.bcast.cmd(weight<-mpi.bcast.Robj())
      mpi.bcast.Robj(weight)
      mpi.bcast.cmd(weight_power<-mpi.bcast.Robj())
      mpi.bcast.Robj(weight_power)
      mpi.bcast.cmd(neighbor_aggregation_method<-mpi.bcast.Robj())
      mpi.bcast.Robj(neighbor_aggregation_method)
      
      mpi.bcast.Robj2slave( comm = 1, all = T)
      mpi.bcast.Robj2slave(make_ego_graph)
      mpi.bcast.Robj2slave(V)
      input_scores<-input_scores[,c("gene","source_node_score")]
      mpi.bcast.cmd(input_scores<-mpi.bcast.Robj())
      mpi.bcast.Robj(input_scores)
      scores<-mpi.iapplyLB(scores_df$gene,function(gene) 
      {
        
        result<-neigborhood_step_score(gene=gene,st=st,input_scores=input_scores,weight=weight,weight_power=weight_power,neighbor_aggregation_method=neighbor_aggregation_method)
        
        return(result)
      }
      )
      mpi.close.Rslaves()
      
    }else
    {scores<- lapply(scores_df$gene,function(gene) 
    {
      #               t<-system.time(
      #                 { 
      #gc()
      result<-neigborhood_step_score(gene=gene,st=st,input_scores=input_scores[,c("gene","source_node_score")],weight=weight,weight_power=weight_power,neighbor_aggregation_method=neighbor_aggregation_method)
      #                }
      #                )
      #               print_v(paste("Step:",st))  
      #               print_v("time:")
      #                 print_v(t)
      #   cat(paste("time:",t,"\n"),file=logfile,append = T)
      return(result)
    }
    )
    }
    
    
    
    #scores_df<-data.frame(gene=V(cn)$name)
    #scores_df$score<-unlist(scores)
    scores_df[paste("step_",st,sep="")]<-unlist(scores)
    #       median_score<-median(unlist(scores_df[paste("step_",st,sep="")]),na.rm=T)
    #       mean_score<-mean(unlist(scores_df[paste("step_",st,sep="")]),na.rm=T)
    #       if (median_score>0)scores_df[paste("step_",st,sep="")]<- scores_df[paste("step_",st,sep="")]/median_score #normalized score
    #       else
    #       {
    #         print_v("median score is zero, normalized by mean")
    #         scores_df[paste("step_",st,sep="")]<- scores_df[paste("step_",st,sep="")]/mean_score
    #       }
    #scores_df[paste("step_",st,sep="")]<- scale(scores_df[paste("step_",st,sep="")])
    
    scores_df_st<-unlist(scores_df[paste("step_",st,sep="")])
    
    
    if (sum(is.infinite(scores_df_st))) scores_df_st[is.infinite(scores_df_st)]<-0 # convert score infinite to zero
    #        scores_df_st[scores_df_st==0]<-NA # convert where score = 0 to NA as log(0) is infininte.
    #        scores_df_log2<-log2(scores_df_st) # convert score to log
    #        score_sd<-sd(scores_df_log2, na.rm = TRUE)
    #        scores_df_log2_scale<-scale(scores_df_log2,center=TRUE, scale=score_sd )
    #        scores_df_st<-2^(scores_df_log2_scale)
    #        
    #        scores_df_st[is.na(scores_df_st)]<-0 # convert where score is NA to zero.
    #        
    if (step_normalization) scores_df[paste("step_",st,sep="")]<-scale_log(scores_df_st)
    
    
    # apply weight for different step
    scores_df[paste("step_",st,sep="")]<-scores_df[paste("step_",st,sep="")]*step_weight[st,"weight"]#weight of step   
    st<-st+1
    
  }
  
  #case: combine scores from different steps, equal sum of different steps.
  if (step==1|!steps_combined) scores_df$nbh_score<-scores_df[,paste("step_",step,sep="")] else 
  {
    step_columns<-vector()
    for (st in 1:step)
    {
      step_columns<-c(step_columns,paste("step_",st,sep=""))
    }
    scores_df$nbh_score<-apply(scores_df[,step_columns],1,function(x) sum(x)) 
    scores_df$nbh_score<-scale_log(scores_df$nbh_score)
  }
  
  #==============================================================================combine nbh score and source node score here:   
  
  
  switch(source_node_inclusion,
         
         n= {
           scores_df$score<-scores_df$nbh_score #do nothing
         },# no inculsion, not change any thing
         s={
           
           
           scores_df$score<-scores_df$nbh_score+scores_df$source_node_score
           
         },#sum of node and target score
         p={
           
           scores_df$score<-scores_df$nbh_score*scores_df$source_node_score
           scores_df$score<-scale_log(scores_df$score)
         }
         #product of node score and target score
         
         
         
  )    
  
  
  
  
  
  
  
  
  
  
  #=========================================================================================
  scores_df<-na.omit(scores_df)
  
  
  ngene<-nrow(scores_df)
  score_df_rank<-data.frame(gene=scores_df$gene)
  score_df_rank$rank<-  ngene -rank(scores_df$score)+1 
  return(list(scores=scores_df,ranks=score_df_rank))
}


#=====================================================================================================================================================
# the score for each gene with given step neighborhood (not combined yet)
#gene: gene name
#st: number of step from source node
#input_scores: format: gene, input_score
#neighbor_aggregation_method="Method to aggregate neighborhood score: sum, average or median (s,a,m)"
neigborhood_step_score<-function(gene="SOX2",st=opt$step,input_scores=ges[,c("gene","fdr")],weight=opt$edge_weight,weight_power=opt$weight_power ,neighbor_aggregation_method=opt$neighbor_aggregation_method) 
{  
  #print_v(paste("processing ",gene))
  
  if(weight=="unw")drs<-distance_rho_step(v1=gene,step=st,output="n") else drs<-distance_rho_step(v1=gene,step=st,output = "a")
  nbh<-unlist(drs[[1]]) # list of neiborhood genes at step i
  
  ges_gene_stat<-input_scores[input_scores$gene %in% nbh,]
  colnames(ges_gene_stat)[2]<-"input_score"
  ges_gene_stat<-ges_gene_stat[is.finite(ges_gene_stat$input_score),] # remove inf value in the vector
  ngene<-nrow(ges_gene_stat) # number of genes in neigborhood
  #node_stat<-input_scores[input_scores$gene==gene,"input_score"]
  #m: master node statistics only, equally products of all input statistics
  # if(source_node_inclusion=="m")score<-node_stat #only use source node statisctics to calculate score
  #else #calculate target statistics and combine
  # {
  if(weight=="unw") 
  {
    switch(neighbor_aggregation_method,
           m= {score<-median(ges_gene_stat$input_score,na.rm = T)},
           a = {score<-mean(ges_gene_stat$input_score,na.rm = T)},
           s =  {score<-sum(ges_gene_stat$input_score,na.rm=T)}
           
    )
  }else
  {
    #how to weight, absolute or relative. Sum of weight should be constant or not.
    
    nbh_distance_rho<-drs[[3]]
    ges_gene_stat$rho<-sapply(as.character(ges_gene_stat$gene),function(x)nbh_distance_rho[x])
    if(ngene>0)# at least one gene 
    {
      if(weight=="rho") 
      {
        target_scores<-ges_gene_stat$rho^weight_power*ges_gene_stat$input_score
        source_weight<-sum(ges_gene_stat$rho^weight_power,na.rm = T)
      }
      
      if(weight=="MI")
      {
        target_scores<-(rho_to_MI(ges_gene_stat$rho))^weight_power*ges_gene_stat$input_score
        source_weight<-sum((rho_to_MI(ges_gene_stat$rho))^weight_power,na.rm=T)
      }
      
      score<-sum(target_scores,na.rm = T)
      if (neighbor_aggregation_method=="m"|neighbor_aggregation_method=="a") score<-score/source_weight
      
    }else score<-0
  }
  
  #         if(neighbor_aggregation_method=="a" ) 
  #                 {if (ngene>0) score<-score/ngene
  #                   else score<-0}
  #                  print_v(paste("score=",score))
  #                      if (gene=="UNC80") 
  #                      {print_v("check here")}    
  
  #include source node statistics in the score or not.
  # n: no
  #s: sum of node and target score 
  #p: product of node score and target score
  #                 switch(source_node_inclusion,
  #                        
  #                      n= {
  #                         #do nothing
  #                           },# no inculsion, not change any thing
  #                      s={
  #                       
  #                         source_weight<-0
  #                         if(weight=="unw") source_weight<-ngene
  #                         if(weight=="rho") source_weight<-sum(ges_gene_stat$rho^weight_power)
  #                         if(weight=="MI") source_weight<-sum((rho_to_MI(ges_gene_stat$rho))^weight_power)
  #                         score<-score+source_weight*node_stat
  #                   
  #                       },#sum of node and target score
  #                      p={
  #                      
  #                         score<-node_stat*score
  #                       }
  #                      #product of node score and target score
  #                    
  #             
  #                   
  #                   )    
  
  
  if(is.nan(score))
  {
    print_v("Score is not a number")
    print_v(gene)
    print_v(ngene)
    
  }
  # } 
  return(score) #return score of function neigborhood_step_score
}







###############################################################################################################################################
#Calculate scores several round, the output scores rank of the previous round serve as the input of the next round
#nround = number of round
#input_scores: data frame with format (gene,scores) 
#return value: ranks: a list of rank and score of each iterative round. Each rank and score contain score at each step.
#use_rank: use rank as input instead of the score value
iterative_master_scores<-function(input_scores=ges[,c("gene",gene_statistics_list)],nround=opt$nround,step=opt$step,weight=opt$edge_weight,weight_power=opt$weight_power,use_rank=opt$use_rank,source_node_inclusion=opt$source_node_inclusion,neighbor_aggregation_method=opt$neighbor_aggregation_method,step_power=opt$step_power,steps_combined=opt$steps_combined,step_normalization=opt$step_normalization,input_normalization=opt$input_normalization,mpi=F,nslaves=10, cluster=F)
{
  i_scores<-input_scores  
  
  
  ranks<-list()
  for (round in 1:nround)
  {
    if(use_rank) # convert  value to rank
      for (column in 2:ncol(i_scores))
      { i_scores[column]<-rank(i_scores[column])}
    scores<-master_score_step(step = step,input_scores=i_scores,weight= weight,weight_power=weight_power,source_node_inclusion=source_node_inclusion,neighbor_aggregation_method=neighbor_aggregation_method,step_power=step_power,steps_combined = steps_combined,step_normalization=step_normalization,input_normalization=input_normalization,mpi=mpi,nslaves=nslaves, cluster=cluster)
    #ranks[[round]]<-scores[[2]]
    ranks[[round]]<-scores
    #i_scores<-scores[[2]][,1:2]
    i_scores<-scores[[1]][,c("gene","score")] #case use score from the first round as input for second round
    colnames(i_scores)[2]<-paste("round_",round,"_score",sep="")
    #i_scores$betweenness<-ges$betweenness
    ngenes<-nrow(i_scores)
    #i_scores[2]<-ngenes-(i_scores[2])+1 #update input score as rank of the resulting score
    
  }
  return(ranks[[nround]])
}


#check the rank of master driver after iterative calculation
#ranks : the return of the iterative_master_scores function
check_iterative_rank<-function(ranks,genes)
{
  
  rank_table<-ranks[[2]]
  return(rank_table[rank_table$gene %in% genes,][,])
  
  
}

iterative_master_scores_convergence<-function(input_scores=ges[,c("gene",gene_statistics_list)],step=opt$step,weight= opt$edge_weight,weight_power=opt$weight_power,use_rank=opt$use_rank,source_node_inclusion=opt$source_node_inclusion,neighbor_aggregation_method=opt$neighbor_aggregation_method,steps_combined=opt$steps_combined,convergence=5, mpi=F)
{
  i_scores<-input_scores  
  
  
  ranks<-list()
  round<-1
  out<-"ITERATIVE_RESULT"
  dir.create(out,showWarnings = F)
  dir.create(paste(out,"/scores",sep=""),showWarnings = F)
  dir.create(paste(out,"/ranks",sep=""),showWarnings = F)
  logfile<-paste(out,"/iterative_log.csv",sep="")
  file.create(logfile)
  convergence_score<-10
  while (convergence_score>convergence)
  {
    if(use_rank) # convert  value to rank
      for (column in 2:ncol(i_scores))
      { i_scores[column]<-rank(i_scores[column])}
    scores<-master_score_step(step = step,input_scores=i_scores,weight= weight,weight_power=weight_power,source_node_inclusion=source_node_inclusion,neighbor_aggregation_method=neighbor_aggregation_method,steps_combined=opt$steps_combined, mpi=mpi)
    scores_round<-scores[[1]]
    rank_round<-scores[[2]]
    write.table(scores_round,file =paste(out,"/scores/","round",round,sep = ""),quote=F )
    write.table(scores_round,file =paste(out,"/ranks/","round",round,sep = ""),quote=F )
    ranks[[round]]<-scores
    i_scores<-scores[[1]] #case use score from the first round as input for second round
    ngenes<-nrow(i_scores)
    master_genes_ranks<-rank_round[rank_round$gene %in% master_genes,]
    master_genes_ranks$round<-round
    write.table(master_genes_ranks,file=paste(out,"/master_gene_ranks.csv",sep=""),append = T)
    
    if (round>1)
    {
      rank_diff<-abs(ranks[[round]][[2]]$rank-ranks[[round-1]][[2]]$rank)
      mean_rank_diff<-mean(rank_diff,na.rm = T)
      print_v(paste("round:",round,"\n"))
      cat(paste("round:",round,"\n"), file=logfile,append=T)
      print_v(paste("mean_rank_diff:",mean_rank_diff,"\n"))
      cat(paste("mean_rank_diff:",mean_rank_diff,"\n"), file=logfile,append=T)
      convergence_score<- mean_rank_diff
    }
    
    round<-round+1
    
  }
  return(ranks)
}









#experiment indexing: use for job array when doing optimization condition of nScore.
run_experiment<-function(experiment_indexing=F)
  
{
  outdir<-opt$outdir 
  dir.create(outdir)
  sys_time<-gsub(" ", "_", format(Sys.time(), "%a %b %d %X %Y"))
  old_dir<-getwd()
  gep_filename<-tools::file_path_sans_ext(basename(gep_file))
  
  outdir<-paste0(outdir,"/",gep_filename)
  dir.create(outdir,showWarnings = F,recursive=T)
  if (experiment_indexing) newdir<-paste(outdir,"/experiment_run_",opt$index,sep = "") else newdir<-outdir
  
  dir.create(newdir)# create new directory and change working directry to this directory
  setwd(newdir)
  logfile<-paste("experiment",opt$index,".log", sep = "")
  file.create(logfile)
  print_v(logfile)
  cat(paste("Beginning calculation: ",format(Sys.time(), "%a %b %d %X %Y"),"\n"),file=logfile,append = TRUE)
  cat(paste("Master_genes:",master_genes,"\n"),file=logfile,append = TRUE)
  
  
  write.table(data.frame(opt),file="experiment_condition.csv",append = F,quote = F,row.names = F)
  runtime<-system.time({
    cat(paste("=================================================================================================\n"),file=logfile,append = T)
    experiment_rank<-iterative_master_scores(ges[,c("gene",gene_statistics_list)],nround=opt$nround,step=opt$step,weight= opt$edge_weight,weight_power=opt$weight_power,use_rank=opt$use_rank,source_node_inclusion=opt$source_node_inclusion,neighbor_aggregation_method=opt$neighbor_aggregation_method,step_power = opt$step_power)
    
    if (length(experiment_rank)==0) {
      cat(paste("no rank in round", round, "\n"),file=logfile,append=T)
    }else
    {
      
      round_scores<-as.data.frame(experiment_rank[[1]])
      
      #scores_sum = sum(round_scores$score)
      scores_median = median(round_scores$score)
      if (scores_median == 0)
      {
        scores_median = mean(round_scores$score)
      }
      #round_scores$score =  (round_scores$score/scores_sum) * nrow(round_scores)  # normalized the score so total sum would equal number of genes, average score would be 1
      round_scores$score =  round_scores$score/scores_median # normalized the score so median score would be 1
      round_scores<-round_scores[order(round_scores$score, decreasing = TRUE),]
      round_ranks<-as.data.frame(experiment_rank[[2]])
      round_ranks<-round_ranks[order(round_ranks$rank),]
      print("top 20 master genes:")
      print(round_ranks[1:20,])
      print(round_scores[1:20,])
      write.table(round_scores,file="scores.csv",quote=F,sep="\t",row.names =F)
      write.table(round_ranks,file="ranks.csv",quote=F,sep="\t",row.names  =F)
    }
    
    
    print_v("finish experiment")
    t<-check_iterative_rank(experiment_rank,master_genes)
    
    print_v(t)
    master_ranks_sum<-sum(t$rank)
    print_v(master_ranks_sum)
    if (nrow(t)==0) cat(paste( "zero master genes found in the ranking score. Master genes are out of top genes. \n"), file = logfile, append = T) else
    {
      
      print_v(t)
      write.table(t,file="master_genes_ranks.csv",append = T,quote = F,row.names = F)
      write.table(t,file=logfile,append = T,quote = F,row.names = F)
      
    }
    
    result<-list(scores=round_scores, ranks= round_ranks, master_ranks=t,
                 master_ranks_sum = master_ranks_sum,  de_data=ds)
    outfile<-paste0(gep_filename,"_master_score_result.RDS")
    saveRDS(result, file = outfile)
    
    
    
  }) # end of runtime
  
  cat(paste("Run time:", runtime,"\n"),file=logfile,append = T)
  cat(paste("Experiment finished, result saved in ","result.Rdata","\n"),file=logfile,append = T)
  cat(paste("=================================================================================================\n"),file=logfile,append = T)
  setwd(old_dir)   
  return(list(all_genes_rank=experiment_rank, master_genes_ranks=t))
  
  
  
}


#combine a list of data frames with the same columns into one big data frame
df_combine<-function(df_list)
{
  df_nrow<-unlist(sapply(df_list,function(x) nrow(x)))
  df_nrow<-sum(df_nrow)
  df <- data.frame(matrix(NA, ncol = ncol(df_list[[1]]), nrow = df_nrow))
  colnames(df)<-colnames(df_list[[1]])
  begining_row<-1
  for (i in 1:length(df_list))
  {
    df_l<-df_list[[i]]
    ending_row<-(begining_row+nrow(df_l)-1)
    if(length(begining_row:ending_row)==0)
    {
      
      print_v("need to stop")
    }
    rows<-c(begining_row:ending_row)
    if (length(rows)==0)
    {
      print_v("stop")
      
    }
    df[begining_row:ending_row,]<-df_l
    begining_row<-ending_row+1
    
  }
  rm(df_nrow)
  #gc()
  return(df)
  
}


#combine a list of data frames with the same columns into one big data frame, then sum of each column, out put a data frame with sum rum
df_combine_sum<-function(df_list,nodes)
{
  df<-df_combine(df_list)
  
  print_v("calculating node weight")
  print_v(paste("at",format(Sys.time(), "%a_%b_%d_%Y , %H : %M :%S")))
  df_colnames<-colnames(df)
  #sum each row:
  df<-data.frame(t(colSums(df,na.rm = T)))
  colnames(df)<-df_colnames# need to reassign back column names as the colSums operation change the col names for some character, for example from "SOX2-OT" to "SOX2.OT"
  print_v(" ending combine sum")
  print_v(paste("at",format(Sys.time(), "%a_%b_%d_%Y , %H : %M :%S")))
  return(df)
}

df_return<-function(asp,node_weight,ngene,df_colnames,nodes)
{
  
  #   node_weight<-data.table(node_weight)
  #   setkey(node_weight,node)
  df <- data.frame(matrix(NA, ncol = 1+2*ngene, nrow = length(asp)))
  colnames(df)<-df_colnames
  
  
  print_v("calculating node weight product")
  print_v(paste("at",format(Sys.time(), "%a_%b_%d_%Y , %H : %M :%S")))
  
  asp<-lapply(asp, function(x) names(unlist(x)))
  nsp<-length(asp)
  for (i in 1:nsp)
  {
    #print_v(i)
    sp<-asp[[i]]
    
    df[i,sp]<-1
    
    
    
    node_weights_product<-prod(node_weight[sp,"weight"])
    
    
    df[i,1]<-node_weights_product
    
  }
  
  
  
  
  
  #  for (node in nodes)
  #  {
  #    #print_v(node)
  #    # print_v(df[1,1:10])
  # 
  #    # print_v(paste("node",node)) #colnames(node_weight)<-c("node","weight")
  #    df[,paste(node,"_weight",sep="") ]<-df[,"node_weights_product"]* df[,node]
  #  }
  
  print_v("multiplying node weight")
  print_v(paste("at",  format(Sys.time(), "%a_%b_%d_%Y , %H : %M :%S")))
  df[,paste(nodes,"_weight",sep="") ]<-df[,"node_weights_product"]* df[,nodes] 
  
  df<- list(df)
  df<-df_combine_sum(df,nodes)
  print_v("finish df_return")
  print_v(paste("at",  format(Sys.time(), "%a_%b_%d_%Y , %H : %M :%S")))
  rm(node_weight,asp)
  return(df)
  
}


# get df from one gene , output is one summary row
get_df_list<-function(x,node_weight,ngene,df_colnames,nodes) 
{
  print_v("calculating all shortest paths...")
  print_v(paste("at", format(Sys.time(), "%a_%b_%d_%Y , %H : %M :%S")))
  asp<-all_shortest_paths(cn,from=x)[[1]]
  print_v("calculating df..")
  print_v(paste("at", "at", format(Sys.time(), "%a_%b_%d_%Y , %H : %M :%S")))
  df<-df_return(asp,node_weight,ngene,df_colnames,nodes)
  print_v("finish get_df_list")
  print_v(paste("at",  format(Sys.time(), "%a_%b_%d_%Y , %H : %M :%S")))
  rm(asp)
  return(df)
}

get_df_list_pair<-function(x,node_weight,ngene,df_colnames) 
{
  sampling_size<-1000
  asp<-list()
  g2<-sample(ges$gene,size=1)
  asp<-lapply(1:sampling_size,function(i) shortest_paths(cn, from = x, to =g2)$vpath)
  df<-df_return(asp,node_weight,ngene,df_colnames,nodes)
  rm(asp)
  return(df)
}


# node weighted betweennes
#convergence_option: ind change for each node comparing with previous bath, or m : mean of betweenness.  mean of betweenness would be quicker to converge, and contributed by more important genes.
#sampling_option: ind - randomly pick individual node and calculate all shortest paths from this node. p - pair: randomly pick a pair of nodes and calculate the shortest path from one node to another node
#mpi: use mpi in calculation

nw_betweenness<-function(node_weight= ges[,c("gene","logFC","fdr")],edge_weight="unw",convergence_criteria =2,convergence_option="ind",sampling_option="ind",mpi=T,nslaves=10,cluster=T, gene_index=1,batch_size=2) #total
{
  
  if (!is.loaded("mpi_initialize")) {library(Rmpi)}
  
  if (cluster)
  { 
    print_v("running on cluster...")
    nslaves <- mpi.universe.size() - 1 # only for cluster
    
  }
  
  
  
  
  
  
  print_v("beginning calculation betweenness")  
  #approximate calculation by sampling  two random nodes and calculate the shortest path, the weight. Use for all genes, the number of sampling increases incrementatally. Objective function
  # is the change of average betweenees. Will be data frame: gene1 gene 2 statistics1, statistics2  gene1 presenence, gene 2 presence....
  # another dataframe: gene1 gene 2 prod(statistics1, statistics2)  gene1 product , gene 2 presence....
  #batch size: number of random shortest path samples before recalculating betweenness
  #edge weight: use for finding shortest path.
  #convergence_criteria: the cut_off of change of average betweenneess to consider convergence
  
  
  
  node_weight$weight<-apply(data.frame(node_weight[,2:ncol(node_weight)]),1,prod)
  #colnames(node_weight)<-c("node","weight")
  colnames(node_weight)[1]<-"node"
  
  ngene<- nrow(node_weight)
  
  #sorting by node weight
  node_weight<-node_weight[order(node_weight$weight,decreasing = T),]
  
  nodes<-node_weight$node
  rownames(node_weight)<-nodes
  
  #df_colnames<-c("node1","node2","node1_weight", "node2_weight","node_weights_product",nodes, paste(nodes,"_weight",sep="") )
  df_colnames<-c("node_weights_product",nodes, paste(nodes,"_weight",sep="") )
  
  
  
  dirname<-paste("BETWEENNESS_",format(Sys.time(), "%a_%b_%d_%Y_%H_%M_%S"),sep="")
  dir.create(dirname)
  batch_size<-nslaves
  file.create(paste(dirname,"/mean_change.csv",sep=""))  
  betweenness<-vector()
  df <- data.frame(matrix(NA, ncol = 1+2*ngene, nrow = batch_size))
  colnames(df)<-df_colnames
  if (convergence_option=="m") criteria<-1 else criteria<-10^10
  batch=1
  df_sum_cumulative<-data.frame(Sum=colSums(df[,1:ncol(df)],na.rm = T))
  node_betweenness<-data.frame(node=nodes)
  node_betweenness_unw<-data.frame(node=nodes)# unweighted betweenness
  asp<-list()#list of shortest paths
  
  
  print_v("beginning checking for mpi")   
  try(mpi.close.Rslaves(),silent = T) 
  mpi.spawn.Rslaves(nslaves=nslaves)   
  mpi.bcast.cmd(cn<-mpi.bcast.Robj())
  mpi.bcast.Robj(cn) 
  mpi.bcast.cmd(ngene<-mpi.bcast.Robj())
  mpi.bcast.Robj(ngene)
  mpi.bcast.cmd(df_colnames<-mpi.bcast.Robj())
  mpi.bcast.Robj(df_colnames)
  mpi.bcast.cmd(node_weight<-mpi.bcast.Robj())
  mpi.bcast.Robj(node_weight)
  mpi.bcast.cmd(nodes<-mpi.bcast.Robj())
  mpi.bcast.Robj(nodes)
  mpi.bcast.Robj2slave(data.table)
  mpi.bcast.Robj2slave(setkey)
  mpi.bcast.Robj2slave(df_return)  
  mpi.bcast.Robj2slave(all_shortest_paths)
  mpi.bcast.Robj2slave(shortest_paths)
  mpi.bcast.Robj2slave(get_df_list)
  mpi.bcast.Robj2slave(get_df_list_pair)
  mpi.bcast.Robj2slave(df_combine)
  mpi.bcast.Robj2slave(df_combine_sum)
  mpi.bcast.cmd(ges<-mpi.bcast.Robj())
  mpi.bcast.Robj(ges) 
  mpi.remote.exec("ls") 
  gene_index<-1
  
  
  
  
  
  while (((criteria > convergence_criteria) | criteria==0) & gene_index<=ngene)
  {
    
    print_v(paste("beginning  batch", batch,  "at", format(Sys.time(), "%a_%b_%d_%Y , %H : %M :%S")))
    
    # do something to calculate the change
    print_v(paste("batch size: ", batch_size))
    #do mpi staff
    
    # the index of gene as a source node
    
    df <- data.frame(matrix(NA, ncol = 1+2*ngene, nrow = 0))
    colnames(df)<-df_colnames
    df_run<-list()
    for (run in 1:(batch_size/nslaves))
    {
      if (gene_index<=ngene)
        # if (!is.loaded("mpi_initialize")) {library(Rmpi)}
        
        # return a data frame from list of shortest paths      
      {
        print_v(paste("gene index:", gene_index))
        print_v(paste("run:", run,   "at", format(Sys.time(), "%a_%b_%d_%Y , %H : %M :%S")))   
        g1<-nodes[gene_index:min(gene_index+nslaves-1,ngene)]
        print_v(paste("beginning mpi jobs at", format(Sys.time(), "%a_%b_%d_%Y , %H : %M :%S")))
        
        print_v(g1)
        if (sampling_option=="ind")
        {
          
          df_list<-mpi.iapplyLB(g1, function(x)get_df_list(x,node_weight,ngene,df_colnames,nodes))
        }
        
        if (sampling_option=="p")   
          
        {
          df_list<-mpi.iapplyLB(g1,function(x) get_df_list_pair(x,node_weight,ngene,df_colnames,nodes))
        }
        
        mpi.bcast.cmd(gc())
        print_v(paste("ending mpi jobs at", format(Sys.time(), "%a_%b_%d_%Y , %H : %M :%S")))
        
        
        print_v("combine result of one run from different slaves")
        print_v(paste("at",format(Sys.time(), "%a_%b_%d_%Y , %H : %M :%S")))
        df_run[[run]]<-df_combine(df_list)
        gene_index<-gene_index + nslaves
        
        
        
      }
      
    }
    
    print_v(paste("combining different runs",   "at", format(Sys.time(), "%a_%b_%d_%Y , %H : %M :%S")))
    print_v(paste("at",format(Sys.time(), "%a_%b_%d_%Y , %H : %M :%S")))
    df<-df_combine(df_run)
    
    
    
    
    
    
    
    
    df_sum<-data.frame(Sum=colSums(df[,1:ncol(df)],na.rm = T))
    df_sum_cumulative<-df_sum_cumulative + df_sum
    Grand_weight<- df_sum["node_weights_product","Sum"]
    
    if(Grand_weight>0)node_betweenness[,paste("batch_",batch, sep="")]<-df_sum[paste(nodes,"_weight",sep=""),"Sum"]/Grand_weight else node_betweenness[,paste("batch_",batch, sep="")]<-0
    node_betweenness_unw[,paste("batch_",batch, sep="")]<-df_sum[nodes,"Sum"]/batch_size
    node_betweenness[node_betweenness==0]<-10^(-10)
    betweenness[batch]<-mean(node_betweenness[,paste("batch_",batch, sep="")], na.rm=T) #node_betweenness: vector of betweenness for each gene. 
    
    print_v(paste("batch:", batch,"; average betweenness:", betweenness[batch]))
    if (batch==1) {plot(x=batch, y=criteria, xlim=range(0,20), ylim=range(0,2) )} else
    {
      
      change<-abs(node_betweenness[,paste("batch_",batch, sep="")]-node_betweenness[,paste("batch_",batch-1, sep="")])/node_betweenness[,paste("batch_",batch-1, sep="")]# betweenness[batch] - the average of betweenness for this batch
      mean_change<-mean(change, na.rm = T)
      cat(paste("batch:", batch,". Mean change %:", mean_change *100,"\n"),file = paste(dirname,"/mean_change.csv",sep=""),append=T)
      #if (batch==2) plot(x=batch,y= log10(mean_change), xlim=range(0,50),ylim=range(0,15))
      #else points(x=batch,y= log10(mean_change))
      print_v(paste("mean change:", mean_change *100, "%"))
      switch(convergence_option,
             ind= {criteria<-mean_change},
             m={
               criteria<-abs(betweenness[batch]-betweenness[batch-1])/betweenness[batch-1]# % change of  betweenness
               
               
             }
             
      )
      
      points(x=batch, y=criteria)
      node_betweenness[,paste("batch_",batch, sep="")]<-df_sum_cumulative[paste(nodes,"_weight",sep=""),"Sum"]/df_sum_cumulative["node_weights_product","Sum"] # new node_betweenness updated, combined with the old one
      node_betweenness[node_betweenness==0]<-10^(-10)
      node_betweenness_unw[,paste("batch_",batch, sep="")]<-df_sum_cumulative[nodes,"Sum"]/(2*batch_size)
      betweenness[batch]<-mean(node_betweenness[,paste("batch_",batch, sep="")], na.rm=T)
      print_v(paste("ending calculation for batch:", batch, "at", format(Sys.time(), "%a_%b_%d_%Y , %H : %M :%S"),"; average betweenness after combination with old batch:", betweenness[batch]))
      print_v(paste("convergence value:",criteria))
      
      batch_size<-2*batch_size
    }
    print_v("Writing tables...)")
    write.table(node_betweenness,file=paste(dirname,"/node_betweenness_csv",sep=""),quote = F, sep="\t")
    write.table(node_betweenness_unw,file=paste(dirname,"/node_betweenness_unw.csv",sep=""),quote = F, sep="\t")
    write.table(betweenness,file=paste(dirname,"/betweenness.csv",sep=""),quote = F, sep="\t")
    print_v("Updating batch")
    batch<-batch+1
    print_v("Removing df_sum")
    rm(df,df_sum,Grand_weight)
    # gc()
    
    print_v("Closing mpi...")
    mpi.close.Rslaves() 
    mpi.finalize()
    
    
    
    #Node weigted betweennes = sum(logFC(node1)*logFC(node2)*PresenceofNode(v))/sum(logFC(node1)*logFC(node2))
    return(node_betweenness)
  }
  
  
  
  
  
  
  
}



#For each set of genes, selected by gene_index and batch size, calculate the df (intermediate tables)
nw_betweenness_nompi_df<-function(node_weight= ges[,c("gene","logFC","fdr")],edge_weight="unw", gene_index=1,batch_size=10,outdir="DFs")
{
  
  
  node_weight$weight<-apply(data.frame(node_weight[,2:ncol(node_weight)]),1,prod)
  #colnames(node_weight)<-c("node","weight")
  colnames(node_weight)[1]<-"node"
  
  ngene<- nrow(node_weight)
  
  #sorting by node weight
  node_weight<-node_weight[order(node_weight$weight,decreasing = T),]
  
  nodes<-node_weight$node
  rownames(node_weight)<-nodes
  
  #df_colnames<-c("node1","node2","node1_weight", "node2_weight","node_weights_product",nodes, paste(nodes,"_weight",sep="") )
  df_colnames<-c("node_weights_product",nodes, paste(nodes,"_weight",sep="") )
  
  
  if (gene_index>ngene)
  {
    return("gene_index > number of genes")
  }
  if (gene_index+batch_size-1>ngene) batch_size<-ngene-gene_index+1
  
  df_list<-lapply(1:batch_size,function(i)
  {
    current_index<-gene_index + i-1
    print_v(paste("run:",i))
    if (current_index<=ngene)
    {
      g1<-nodes[min(gene_index + i-1,ngene) ]
      print_v(g1)
      
      df<-get_df_list(g1,node_weight,ngene,df_colnames,nodes)
      
      return(df)
    }
  })
  
  df<-df_combine_sum(df_list,nodes)
  #gene_index<-gene_index + batch_size 
  
  write.table(df, file=paste(outdir,"/df_",gene_index,"_",gene_index+batch_size-1,".df",sep=""),quote=F, sep="\t")
  
  
}




#Reads the dfs files, combine and calculate the betweenness
nw_betweenness_nompi_finalize<-function(indir="DFs")
{
  files<-dir(indir,full.names = T)
  files<-grep(".df",files, value = T)
  
  
  
  df_rownames<-readLines(files[1])[1]
  df_rownames<-strsplit(df_rownames,"\t")[[1]]
  
  
  sum_df<-rep(0, length(df_rownames)+1) 
  i<-1
  for (f in files)
  {
    print_v(i)
    df<-readLines(f)[[2]]
    df<-as.numeric(strsplit(df,"\t")[[1]])
    i<- i+1
    sum_df<- sum_df + df
  }
  
  df<-data.frame(Item= df_rownames, value = sum_df[2:length(sum_df)])
  rownames(df)<-df_rownames
  
  ngene<-(length(df_rownames)-1)/2
  nodes<-df_rownames[2:(1+ngene)]
  
  node_weights_products<-df["node_weights_product","value"]
  
  node_betweenness<-data.frame(betweennness_nw=df[paste(nodes,"_weight",sep=""),2]/node_weights_products)
  
  
  rownames(node_betweenness)<-nodes
  colnames(node_betweenness)<-"betweenness_nw"
  node_betweenness$node<-nodes
  node_betweenness<-node_betweenness[,c("node","betweenness_nw")]
  node_betweenness<-node_betweenness[order(node_betweenness$betweenness_nw,decreasing = T),]
  node_betweenness$rank<-ngene +1 -rank(node_betweenness$betweenness_nw)
  total_paths<-sum(df[2:(1+ngene), "value"])
  node_betweenness_unw<-data.frame(df[nodes, "value"]/total_paths )
  colnames(node_betweenness_unw)<-"betweenness"
  node_betweenness_unw$node<-nodes
  node_betweenness_unw<-node_betweenness_unw[,c("node","betweenness")]
  node_betweenness_unw<-node_betweenness_unw[order(node_betweenness_unw$betweenness,decreasing = T),]
  node_betweenness_unw$rank<-ngene +1 -rank(node_betweenness_unw$betweenness)
  outdir<-paste(indir,"/BTW",sep="")
  dir.create(outdir)
  write.table( node_betweenness, file=paste(outdir,"/node_betweenness.csv",sep=""),quote=F, sep="\t",row.names = F)  
  write.table( node_betweenness_unw, file=paste(outdir,"/node_betweenness_unw.csv",sep=""),quote=F, sep="\t",row.names = F)
  return(list(node_betweenness,node_betweenness_unw))
}

#create gene index , write script for submitting to cluster
nw_betweenness_nompi_manager<-function(node_weight= ges[,c("gene","logFC","fdr")],batch_size=5,outdir="DFs",sbatch_templatefile="../master_regulator_template.sh", source_script="Master_Regulator_Scoring_v8.R",ncore=300)
{
  manager_script="Experiment_Managerv6.R"
  source(manager_script)
  
  old_dir<-getwd() 
  newdir<-paste("BETWEENNESS_",format(Sys.time(), "%a_%b_%d_%Y_%H_%M_%S"),sep="")
  
  dir.create(newdir)# create new directory and change working directry to this directory
  setwd(paste(old_dir,"/",newdir, sep = ""))
  logfile<-paste("btw.log",sep = "")
  file.create(logfile)
  logfile<-paste(getwd(),"/",logfile,sep="")#absolute path to logfile.
  
  dir.create("JOBS")
  ngene<- nrow(node_weight)
  gene_index<-1
  i<-1
  while (gene_index<=ngene)
  {
    jobfile<-paste("job",i,".sh",sep="")
    file.copy(from=sbatch_templatefile, to=paste("JOBS/",jobfile,sep=""))
    
    
    cmd_str<-paste("Rscript ../", source_script,  " --cmd btw --gene_index ", gene_index, " --network ../", opt$network, " --gene_ex_stat ../",opt$gene_ex_stat, " --batch_size ", batch_size, sep = "")
    cat(cmd_str,file=paste("JOBS/",jobfile,sep=""),append = T) 
    
    gene_index<-gene_index +batch_size
    i<-i+1
  }
  
  # write scripts
  
  
  print_c(paste("Total", i-1, "job files create. \n"),file=logfile,append = TRUE)
  DF_dir<-"DFs"
  dir.create(DF_dir)
  jobID<- array_submit(logfile=logfile, njob=ncore,nexp = i-1,template = "../btw_array_template.sh")
  print_c(jobID,logfile)
  
  print_c(getwd(),logfile)
  
  all_jobs_completed(logfile = logfile, jobID = jobID)
  
  nw_betweenness_nompi_finalize(DF_dir)
  
  setwd(old_dir)
  
  
}


#scale positive vector by converting to log scale, scale with mean and sd, and exponential back
scale_log<-function(x)
{       
  x[x==0]<-NA # convert where score = 0 to NA as log(0) is infininte.
  x_log2<-log2(x) # convert score to log
  x_sd<-sd(x_log2, na.rm = TRUE)
  x_log2_scale<-scale(x_log2,center=TRUE, scale=x_sd )[,1]
  x<-2^(x_log2_scale)
  
  x[is.na(x)]<-0 # convert where score is NA to zero.
  
  return(x)
  
}





##############################################################################################################################################################

library(optparse)
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stats"))

option_list <- list(
  
  make_option(c("-n", "--network"), action="store", default="/ufrc/dtran/son.le/Klf9_Oc/mouseIPSC_GSE19023/matrix.cyto.ncol.csv", #"h_esc_GSE73211.cyto.csv",#"GBM-rnaseq.final.tr.ncol",
              #"./humanIPSC_GSE14897/GSE73211.cyto.csv",#"./humanIPSC_GSE14897/matrix.cyto.ncol.csv",#"GBM-rnaseq.final.tr.ncol", ,
              help="The gene regulatory network in the tab delimited format gene1 gene 2 MI from generep pipeline", type = "character"),
  make_option(c("-x", "--gene_ex_stat"), action="store", default="/ufrc/dtran/son.le/Klf9_Oc/pipelines/SRA_RAHMAN_combined/DEA/ESC/Fibroblast2ESC_dea.RDS", #./humanIPSC_GSE14897/ges.csv",#"", #"ges_GBM_reprogramming.csv",#"./humanIPSC_GSE14897/ges.csv",#"./Paul_analysis/gene_and_probe_results_files/gene_ttable_nha_nSCORE.csv",#"ges_GBM_reprogramming.csv",# "NSC_vs_GSC.csv"  ,# , ,
              help="The gene expression statistics table in the tab delimited format gene   in exact sequence: gene,logFC,pvalue,fdr,LR ", type = "character"),
  make_option(c("-u","--consider_positive_values_only"), action="store_true", default=F,
              help="Only consider positive values of input statistic (like logFC) only", type = "logical"),
  make_option(c("-f","--fdr_to_confidence"), action="store", default=F,
              help="Convert fdr value to confidence score", type = "logical"),
  
  make_option(c("-s", "--step"), action="store", default=2,
              help="The number of steps from the source node to calculate neigborhood genes", type = "integer"),
  make_option(c("-g", "--top_gene_statistics"), action="store", default="pvalue",
              help="The gene expression statistics used to substract top genes subnetwork. Choose from: logFC,pvalue,fdr,LR ", type = "character"),
  make_option(c("-l", "--gene_statistics_list"), action="store", default="logFC,pvalue,pagerank",#"logFC,fdr,betweenness",
              help="The comma separated list of gene statistics used to calculate master score. Choose from: logFC,pvalue,fdr,LR,degree,betweenness,coreness,pagerank,eigen",
              type = "character"),
  
  
  make_option(c("-r", "--nround"), action="store", default=1,
              help="number of rounds in iterative calculation, the score from previous run serves as the input for the next run", type = "integer"),
  make_option(c("-e", "--edge_weight"), action="store", default="unw",
              help="The edge weight statistics. Choose from rho, MI, unw. rho: normalized MI likes correlation. MI: Mutual Information. unw: unweighted", type = "character"),
  make_option(c("-p", "--weight_power"), action="store", default=1,
              help="weight^power used in calculation, for example weight square", type = "integer"),
  make_option(c("-o", "--step_power"), action="store", default=0,
              help="step power used in  combination of step scores", type = "integer"),
  
  make_option(c("-z", "--steps_combined"), action="store", default=F,
              help="combine scores at different steps or only use the specified step score", type = "logical"),
  make_option(c("-y", "--step_normalization"), action="store", default=F,
              help="Normalize score at each step before step scores combination or not", type = "logical"),
  make_option(c("-w", "--input_normalization"), action="store", default=F,
              help="Normalize inputs before calculation", type = "logical"),
  
  make_option(c("-t", "--top_genes"), action="store", default=0.2,
              help="the proportion of top genes to extract subnetwork from whole network. Using top nodes subnetwork will decrease the computing time and may increase robustness but loosing sensitivity",
              type = "double"),
  
  make_option(c("-k","--use_rank"), action="store", default=FALSE,
              help="Use rank instead of gene differential expression value (LR, fdr, logFC) value as input", type = "logical"),
  make_option(c("-m","--master_genes"), action="store", default=c("POU5F1","SOX2","MYCN","NANOG","LIN28A"),
              help="comma separated list of master genes to validate", type = "character"),
  make_option(c("-d","--source_node_inclusion"), action="store", default="p",
              help="Include source node statistics in the calculation or not. Choose from n (no), s(sum) or p (product), m: use source node statistics only as inputs", type = "character"),
  
  make_option(c("-a","--neighbor_aggregation_method"), action="store", default="a",
              help="Method to aggregate neighborhood score: sum, average or median (s,a,m)", type = "character"),
  
  
  make_option(c("-v","--verbose"), action="store", default=TRUE,
              help="verbose printing messages", type = "logical"),
  make_option(c("-i","--index"), action="store", default=1,
              help="index of the experiment", type = "character"),
  make_option("--cmd", action="store", default="help", # various command btw: beetwenneess, btw_jobs: generate jobs
              help="command to execute, accepted values are run, btw, parameter_optimization, btw_jobs ", type = "character"),
  make_option("--gene_index", action="store", default="1", # 
              help="command to execute", type = "integer"),
  make_option("--batch_size", action="store", default="5", # 
              help="the number of genes that 1 core proceess to calculate betweenness per one batch run", type = "integer"),
  make_option("--outdir", action="store", default="nSCORE", 
              help="output directory", type = "character")  
  
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
program_description = " This program calculate master regulater scores from statistics file (for example differential expression statistics such as logFC, pvalue, fdr, LR) and networks"
parser = OptionParser(option_list=option_list, description = program_description)
opt <- parse_args(parser)

usage = function()
{
  print_help(parser)
  
}


command<-opt$cmd
if (command !="help")
{
  
  print(paste(1305,command))
  #input files: cell_network_file:  network table with 3 column in format: gene1, gene2, MI. gep_file:gene expresion profile file: table with columns in exact sequence: "gene","logFC","pvalue","fdr","LR"
  #setwd("/media/son/Data/Klf9_Oc")
  #gep_file<-"KLF9_ges.csv"
  gep_file<-opt$gene_ex_stat
  #print(gep_file)
  #cell_network_file<-"GBM-rnaseq.final.tr.ncol"
  cell_network_file<-opt$network
  #print_v(cell_network_file)
  top_genes<-opt$top_genes
  master_genes<-trim(strsplit(opt$master_genes,","))
  gene_statistics_list<-trim(strsplit(opt$gene_statistics_list,","))
  print_v(master_genes)
  
  # 
  # print_v(opt$network)
  # print_v(opt$gene_ex_stat)
  # print_v(opt$step)
  # print_v(opt$top_gene_statistics)
  # print_v(opt$nround)
  # print_v(opt$weight)
  # print_v(opt$weight_power)
  
  
  #input: cell_network, gene_exp_stat gene expression statistics table. cell_network: igraph network with weight equall to MI. gene_exp_stat: data frame with column: "gene","logFC","pvalue","fdr","LR"
  print_v("beginning..")
  ext<-substr(gep_file,nchar(gep_file)-3, nchar(gep_file))
  print_v(ext)
  if (ext == ".RDS") 
  {
    print_v(gep_file)
    indata<-readRDS(gep_file)
    if (!is.null(indata$lr_table)) ds<-indata$lr_table else ds<-indata
    print(paste("Total number of genes in input data:", nrow(ds)))
    print(ds[1:5,])
  } else ds <- read.delim(gep_file, stringsAsFactors=FALSE)
  print_v(summary(ds))
  print_v("done reading gep")
  cn<-read_graph(file = cell_network_file,format="ncol") #cell_network
  V(cn)$name<-toupper(V(cn)$name)
  #ges<-data.frame(gene=ds[,1],logFC=ds[,2],pvalue=ds[,3],fdr=ds[,4],LR=ds[,5])# gene expression statistic table
  ges<-data.frame(gene=ds[,1],logFC=ds[,2],pvalue=ds[,3],fdr=ds[,4])# gene expression statistic table
  ges$gene<-toupper(as.character(ges$gene))
  #print(ges[1:20,])
  #ges<-NaRV.omit(ges)# remove gene with missing value, infinite value
  ges<-na.omit(ges)# remove gene with missing value, infinite value
  for (i in 1:ncol(ges))
  {
    for (j in 1:nrow(ges))
    {
      if (is.infinite(ges[j,i]))ges[j,i]<-NA
    }
  }
  ges<-na.omit(ges)# remove gene with missing value, infinite value
  ges<-ges[!duplicated(ges$gene), ] #removed rows with duplicated gene name
  ges<-ges[!duplicated(ges$gene), ] #removed rows with duplicated gene name
  if(opt$consider_positive_values_only) ges[ges<0]<-0  #very important. Need to change again, strongly influence the result
  ges$logFC<-abs(ges$logFC) #convert to positive value
  #ges$LR<-abs(ges$LR) #convert to positive value
  ges$pvalue<-pvalue_zscore(ges$pvalue) #convert to zscore value
  if (opt$fdr_to_confidence==T) 
  {ges$fdr<-fdrvalue_confidence(ges$fdr)} else 
    
  {ges$fdr<-pvalue_zscore(ges$fdr)} #convert to zscore value
  #print(ges$fdr)
  
  for (i in 2:ncol(ges)){ges[,i]<-abs(ges[,i])}# convert all statistics to positive value
  #sort ges by the statistics of interest to get top genes
  ges<-ges[order(ges[,opt$top_gene_statistics],decreasing = T),]
  print("1389 order")
  #extract ges from top genes
  rows_to_extract<-round(nrow(ges)*top_genes)
  ges<-ges[1:rows_to_extract,]
  #print("line 228")
  #print(ges[1:20,])
  #extract subnetwork from top genes
  all_genes<-data.frame(gene=V(cn)$name)
  #print(all_genes)
  all_genes$ID<-rownames(all_genes)
  top_genes_ID<-all_genes[all_genes$gene %in% ges$gene,"ID"]
  cn<-induced_subgraph(cn,as.numeric(top_genes_ID))
  
  
  #remove genes that are lonely (not connected to other genes)
  
  cn_degree<-data.frame(degree=centr_degree(cn, mode = "all", loops = F, normalized = F)$res)
  cn_degree$ID<-rownames(cn_degree)
  not_lonely_genes<-as.numeric(cn_degree[cn_degree$degree>0,"ID"])
  cn<-induced_subgraph(cn,not_lonely_genes)
  rm(cn_degree,all_genes,not_lonely_genes)
  cn_genes<-data.frame(gene=V(cn)$name)
  ges_genes<-intersect(ges$gene, as.character(unlist(cn_genes)))
  ges<-ges[ges$gene %in% ges_genes,]
  #print(ges[1:20,])
  #print_v("calculate centrality measures")
  centralities_list<-c("degree","betweenness","coreness","pagerank","eigen")# list of all implemented topological statistics 
  centralities<-intersect(centralities_list ,gene_statistics_list)# list of centralities to calculate statistics
  #degree
  for (centr in centralities)
  {
    centr_value<-data.frame(gene=V(cn)$name)
    switch (centr,
            
            degree={
              
              centr_value[centr]<-centr_degree(cn, mode = "all", loops = F, normalized = F)$res
              
            } ,
            betweenness=
            {
              centr_value[centr]<-centr_betw(cn, directed = F,normalized = T)$res
              
            },
            coreness=
            {
              centr_value[centr]<-coreness(cn, mode = "all")
              
            } , 
            
            pagerank=
            {
              centr_value[centr]<-page_rank(cn, directed = F, damping = 0.85, personalized = NULL, weights = NA,options = NULL)$vector
            },
            
            eigen=
            {
              centr_value[centr]<-centr_eigen(cn, directed = FALSE, scale = TRUE, options = arpack_defaults, normalized = TRUE)$vector
            }
            
    )
    
    
    ges<-merge(ges,centr_value,by="gene") 
  }
  
  #print_v("done calculation centralities measures")
  
  
  E(cn)$MI<-E(cn)$weight
  E(cn)$rho<-MI_to_rho(E(cn)$MI)
  E(cn)$distance<-1/E(cn)$MI
  
  
  
  
  
switch(command,
       
       btw={
         outdir<-"DFs"
         dir.create(outdir)
         nw_betweenness_nompi_df(gene_index=opt$gene_index, outdir = outdir, batch_size = opt$batch_size)
       },
       parameter_optimization={},
       btw_jobs={nw_betweenness_nompi_manger()},
       run ={
         out = run_experiment()
         }
)
  
  
  
} else
{
  usage()
}
  
#node weighted betweenes, input is node weight data frame, from ges, in format gene;weight

#nw_betweenness(sampling_option = "ind",convergence_option = "ind",convergence_criteria = 0.05,mpi=F,cluster=T)
# library(profvis)
#  t<-profvis(
# {test<-nw_betweenness(sampling_option = "ind",convergence_option = "ind",convergence_criteria = 0.5,mpi=F,nslaves=4,cluster=F)})
#   print_v(t)
