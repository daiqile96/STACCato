library(tensorregress)
library(dplyr)
library(grid)
library(pracma)
library(parallel)
library(reticulate)
library(reshape2)
library(CJIVE)

####. Chordal Distance #####
compare_norm = function(res1, res2){
  
  w2.norm = chord.norm.diff(res1$W$W2, res2$W$W2)
  w3.norm = chord.norm.diff(res1$W$W3, res2$W$W3)
  w4.norm = chord.norm.diff(res1$W$W4, res2$W$W4)
  
  norm = c(w2.norm, w3.norm, w4.norm)
  
  return(norm)
}

#### GSEA ############

GSEA_by_cell_pair = function(sender, receiver, pathways,
                             minSize,
                             maxSize,
                             nPerm){
  
  eff_by_cell_pair = disease_eff %>% 
    filter(Sender == sender & Receiver == receiver) %>% 
    select(LR, effect)
  
  eff = eff_by_cell_pair$effect
  names(eff) = eff_by_cell_pair$LR
  
  fgseaRes <- fgsea(pathways = pathways, 
                    stats    = eff,
                    minSize  = minSize,
                    maxSize  = maxSize,
                    nPermSimple = nPerm)
  
  fgseaRes = fgseaRes[, c('pathway', 'pval', 'NES')] %>% 
    mutate(sender = sender,
           receiver = receiver)
  
  fgseaRes = fgseaRes[complete.cases(fgseaRes)]
  
  return(fgseaRes)
  
}


##### composition of cell factors #########


# pdf(file.path(wk.dir, 'cell_comp.pdf'), width = 12, height = 14)
plot_sender_receiver = function(dcomp_res, cell_level, cell_label){
  
  n = max(ncol(dcomp_res$W$W3), ncol(dcomp_res$W$W4))
  p2 = plot_cell_score('W3', dcomp_res, n, cell_level, cell_label)
  p3 = plot_cell_score('W4', dcomp_res, n, cell_level, cell_label)
  
  my_plot <- data.frame(score = 1:length(cell_level),
                        cells = factor(cell_level, 
                                       levels = cell_level,
                                       labels = cell_label)) %>% 
    ggplot(aes(x=1:length(score), 
               y=score,
               fill=cells)) +
    geom_bar(position="dodge", stat="identity") +
    guides(fill=guide_legend(title="Cell Type"))
  
  legend <- cowplot::get_legend(my_plot)
  
  p_all = gridExtra::arrangeGrob(grobs = list(p2, p3, legend), ncol = 3, widths=c(4.5,4.5,1))
  return(p_all)
}


plot_cell_score = function(mode, dcomp_res, n, cell_level, cell_label){
  
  score.mat = dcomp_res$W[[mode]]
  cells = row.names(score.mat)
  
  if (mode == 'W3'){
    name = 'Sender Cells'
  } else if (mode == 'W4'){
    name = 'Receiver Cells'
  }
  
  p_list = list()
  
  for (i in 1:n){
    
    if (i <= ncol(score.mat)){
      tmp.sender.score = score.mat[,i]
      p_list[[i]] = data.frame(score = unlist(tmp.sender.score),
                               cells = factor(cells, 
                                              levels = cell_level,
                                              labels = cell_label)) %>% 
        arrange(cells) %>% 
        ggplot(aes(x=1:length(tmp.sender.score), 
                   y=score,
                   fill=cells)) +
        geom_bar(position="dodge", stat="identity") +
        ylab(paste('Factor', i)) +
        xlab('') +
        theme_bw() +
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              legend.position = 'none',
              axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank())
      
    } else {
      p_list[[i]] = rectGrob(gp=gpar(fill="white", col="white"))
    }
    
  }
  
  g <- gridExtra::arrangeGrob(grobs = p_list, ncol = 1, 
                              top = name, newpage = T)
  
  return(g)
  
}

##### get disease effects and pvalues ######


get_eff_sig = function(eff_pval_df, eff_cutoff, pval_cutoff){
  
  
  eff_pval_df = eff_pval_df %>% 
    filter(pval_cutoff <= pval_cutoff & abs(effect) >= eff_cutoff) %>% 
    select(LR, Sender, Receiver, effect) %>% 
    mutate(LR = as.character(LR),
           Sender = as.character(Sender),
           Receiver = as.character(Receiver))
  
  return(eff_pval_df)
  
}

plot_eff_pval_sig = function(eff_pval_df, eff_cutoff, pval_cutoff){
  
  eff_pval_df = eff_pval_df %>% 
    mutate(Direction = 
             case_when(effect >= eff_cutoff & pval <= pval_cutoff ~ 'Up',
                       effect <= -1 * eff_cutoff & pval <= pval_cutoff ~ 'Down',
                       .default = 'Not Significant')) %>% 
    mutate(Direction = factor(Direction,
                              levels = c('Down', 'Up', 'Not Significant'),
                              labels = c('Down', 'Up', 'Not Significant')))
  
  ggplot(eff_pval_df, aes(x=effect, y=-1 * log(pval, 10), 
                         col=Direction)) +
    geom_point() + 
    ylab('-log10(p-value)') +
    xlab('Estimated effects') +
    theme(text = element_text(size = 10),
          legend.position = 'bottom') +
    geom_hline(yintercept=-1*log(pval_cutoff,10), linetype="dashed", color = "red") +
    scale_color_manual(values=c("blue", "red", "black")) +
    geom_vline(xintercept=c(-1 * eff_cutoff, eff_cutoff), 
               col="red", linetype="dashed") + 
    scale_x_continuous(breaks = c(round(min(eff_pval_df$effect), 2), 
                                  -1 * eff_cutoff, 
                                  0, 
                                  eff_cutoff,
                                  round(max(eff_pval_df$effect),2)),
                       limits = c(round(min(eff_pval_df$effect), 2) - 0.01,
                                  round(max(eff_pval_df$effect),2) + 0.01)) 
  
  
}

get_full_effect = function(B, var, val){
  
  eff = B[var, , ,]
  lrs = dimnames(B)[[2]]
  senders = dimnames(B)[[3]]
  receivers = dimnames(B)[[4]]
  
  new.names = expand.grid(senders, receivers) %>% 
    mutate(new_name = paste(Var1, Var2, sep = ' - ')) %>% 
    dplyr::select(new_name) %>% 
    unlist() %>% 
    unname()
  
  dim_eff = dim(eff)
  B.arr.reshape = array_reshape(eff, 
                                c(dim_eff[1], dim_eff[2]*dim_eff[3]), 
                                order = 'F')
  B.arr.reshape = array(B.arr.reshape, 
                        c(dim_eff[1], dim_eff[2]*dim_eff[3]),
                        dimnames = list(lrs, new.names))
  eff_long = melt(B.arr.reshape) 
  colnames(eff_long) = c('LR', 'CellPair', val)
  
  eff_long = eff_long %>% 
    mutate(Sender = sapply(CellPair, 
                           function(x) strsplit(as.character(x), split = ' - ')[[1]][1]),
           Receiver = sapply(CellPair, 
                             function(x) strsplit(as.character(x), split = ' - ')[[1]][2])) 
    
    
    return(eff_long[, c('Sender', 'Receiver', 'LR', val)])
  
}

plot_disease_eff_by_receiver_cell = function(df){
  
  max = max(abs(df$effect))
  
  p = ggplot(df, aes(Sender, LR)) +    # Create default ggplot2 heatmap
    geom_tile(aes(fill = effect)) +
    scale_fill_gradient2(name = 'Disease Effect',
                         limits = c(-1 * max, max),
                         low = "blue", 
                         mid = "white", 
                         high = "red") +
    theme_classic() +
    theme(text = element_text(size = 14),
          axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
          legend.position = 'right') +
    facet_wrap(~Receiver, ncol = length(unique(df$Receiver))) +
    xlab('Sender Cell Types') +
    ylab('Ligand - Receptor Pairs')
  
  return(p)
  
}

##### contribution of cell factors ######
get_contribution = function(comp, mode, cov, dcomp_res){
  
  W1 = dcomp_res$W[['W1']]
  W2 = dcomp_res$W[['W2']]
  W3 = dcomp_res$W[['W3']]
  W4 = dcomp_res$W[['W4']]
  G = dcomp_res$G
  r = dim(G)
  
  if (mode == 'cov'){
    
    W1_sub = W1[-comp, ]
    lizt = list(W1_sub, W2, W3, W4)
    C_ts_sub = ttl(as.tensor(G), lizt, ms = c(1,2,3,4))
    U_sub = ttm(C_ts_sub, cov[, -comp], m = 1)
    
  } else if (mode == 'W2'){
    
    G_sub = array(G[,-comp,,], dim = c(r[1], r[2] - 1, r[3], r[4]))
    W2_sub = matrix(W2[, -comp], ncol = r[2] - 1)
    
    lizt = list(W1, W2_sub, W3, W4)
    C_ts_sub = ttl(as.tensor(G_sub), lizt, ms = c(1,2,3,4))
    U_sub = ttm(C_ts_sub, cov, m = 1)
    
  } else if (mode == 'W3'){
    
    G_sub = array(G[,,-comp,], dim = c(r[1], r[2], r[3] - 1, r[4]))
    W3_sub = matrix(W3[, -comp], ncol = r[3] - 1)
    
    lizt = list(W1, W2, W3_sub, W4)
    C_ts_sub = ttl(as.tensor(G_sub), lizt, ms = c(1,2,3,4))
    U_sub = ttm(C_ts_sub, cov, m = 1)
    
  } else if (mode == 'W4'){
    
    G_sub = array(G[,,,-comp], dim = c(r[1], r[2], r[3], r[4] - 1))
    W4_sub = matrix(W4[, -comp], ncol = r[4] - 1)
    
    lizt = list(W1, W2, W3, W4_sub)
    C_ts_sub = ttl(as.tensor(G_sub), lizt, ms = c(1,2,3,4))
    U_sub = ttm(C_ts_sub, cov, m = 1)
    
  }
  
  mean((dcomp_res$U - U_sub@data)^2)
  
  
}

get_contribution_by_mode = function(mode, dcomp_res){
  
  cov = dcomp_res$X_covar1
  
  if (mode == 'cov'){
    r_mode = 2:ncol(cov)
    names = colnames(cov)[r_mode]
  } else {
    r_mode = 1:(ncol(dcomp_res$W[[mode]]))
    names = factor(paste0('Factor', r_mode))
  }

  contr = sapply(r_mode, 
                 get_contribution, 
                 mode = mode, 
                 cov = cov, 
                 dcomp_res = dcomp_res)
  
  contr_df = data.frame(var = names,
                        mse = contr)
  
  if (mode == 'cov'){
    title = 'Sample-level Variables'
  } else if (mode == 'W3'){
    title = 'Sender Cell Types'
  } else if (mode == 'W4'){
    title = 'Receiver Cell Types'
  }
  
  p = ggplot(data=contr_df, aes(x=var, y=mse)) +
    geom_bar(stat="identity", alpha = 0.4) + xlab('') +
    ylab('Contribution') +
    scale_y_continuous(trans = 'sqrt') +
    ggtitle(title) +
    theme(text = element_text(size=10),
          plot.title = element_text(size=10)) 
  
  return(p)
  
}

##### decomposition ######
staccato = function(tsr,
                    X_covar1 = NULL, 
                    lr.names = NULL, 
                    sender.names = NULL, 
                    receiver.names = NULL,
                    core_shape){
  
  # initial: "random" for random initialization; "QR_tucker" for QR-based tucker initialization
  Y_1 = unfold(tsr, row_idx = 1, col_idx = c(2,3,4))@data
  Y_2 = unfold(tsr, row_idx = 2, col_idx = c(1,3,4))@data
  Y_3 = unfold(tsr, row_idx = 3, col_idx = c(1,2,4))@data
  Y_4 = unfold(tsr, row_idx = 4, col_idx = c(1,2,3))@data
  d1 = dim(tsr)[1] ; d2 = dim(tsr)[2] ; d3 = dim(tsr)[3]; d4 = dim(tsr)[4]
  r1 = core_shape[1] ; r2 = core_shape[2] ; r3 = core_shape[3]; r4 = core_shape[4]
  
  # check if sample-level covariates are provided 
  un_m1 = FALSE ; 
  if(is.null(X_covar1)|(identical(X_covar1,diag(d1)))) {X_covar1 = diag(d1) ; un_m1 = TRUE}
  X_covar2 = diag(d2) ; un_m2 = TRUE
  X_covar3 = diag(d3) ; un_m3 = TRUE
  X_covar4 = diag(d4) ; un_m4 = TRUE
  p1 = dim(X_covar1)[2] ; p2 = dim(X_covar2)[2] ; p3 = dim(X_covar3)[2]; p4 = dim(X_covar4)[2]

  # QR based tucker initialization
  # use new QR generalization
  qr1 = qr(X_covar1); qr2 = qr(X_covar2); qr3 = qr(X_covar3); qr4 = qr(X_covar4)
  Q1 = qr.Q(qr1); Q2 = qr.Q(qr2); Q3 = qr.Q(qr3); Q4 = qr.Q(qr4)
  R1 = qr.R(qr1); R2 = qr.R(qr2); R3 = qr.R(qr3); R4 = qr.R(qr4)
  
  new_y = ttl(tsr, list_mat = list(t(Q1), t(Q2),t(Q3), t(Q4)), c(1,2,3,4)) # Y \times Q^T = B \times R
  res_un = tucker(new_y,ranks = core_shape) # HOOI, not random
  C_ts = ttl(res_un$est, list(solve(R1),solve(R2), solve(R3), solve(R4)), c(1,2,3,4))
  
  # output factors need to be orthogonal
  ortho_decomp = tucker(C_ts, ranks = core_shape)
  
  W1 = ortho_decomp$U[[1]]
  W2 = ortho_decomp$U[[2]]
  W3 = ortho_decomp$U[[3]]
  W4 = ortho_decomp$U[[4]]
  
  G = ortho_decomp$Z
  
  U = ttl(C_ts, list(X_covar1, X_covar2, X_covar3, X_covar4),c(1,2,3,4))
  sigma_est=mean((tsr@data-U@data)^2)
  lglk=dnorm(tsr@data,U@data,sqrt(sigma_est),log=TRUE)
  dimnames(C_ts@data) <- list(colnames(X_covar1),
                         lr.names,
                         sender.names,
                         receiver.names)
  rownames(W2) = lr.names
  rownames(W3) = sender.names
  rownames(W4) = receiver.names
  
  return(list(W = list(W1 = W1,W2 = W2,W3 = W3, W4 = W4),
              X_covar1 = X_covar1,
              G = G@data,U=U@data, 
              C_ts = C_ts@data,
              lglk = lglk, 
              sigma=sigma_est))

  
  
}


cal_eigen_varexp = function(tensor, mode0, modes, varexp){
  
  arr = unfold(tensor, row_idx=mode0, col_idx=modes)@data
  mat = arr %*% t(arr) 
  
  # calculate eigenvalue
  e = eigen(mat)
  e.value = e$values
  e.prop = e.value / sum(e.value)
  
  prop.df = data.frame(k = 1:length(e.prop), value = e.prop)
  
  if (nrow(prop.df) > 20){
    prop.df = prop.df[1:20, ]
  }
  
  p = ggplot(data=prop.df, aes(x=factor(k), y=value*100)) + 
    geom_bar(stat="identity", alpha = 0.4) + xlab('') +
    geom_hline(yintercept=varexp, linetype="dashed", color = "red") +
    scale_y_continuous(trans = 'sqrt',
                       breaks = c(0, varexp, 20, 40, 60, 80, 100)) +
    ylab('explained variance (%)')
  
  return(list(p = p, rank = max(which(e.prop > varexp/100))))
  
}

##### bootstrap ########
boot = function(idx, dcomp_res){
  
  set.seed(idx)
  U = dcomp_res$U
  d = dim(U)
  B = dcomp_res$C_ts
  r = dim(dcomp_res$G)
  sd = sqrt(dcomp_res$sigma)
  
  # parametric bootstrap 
  new.residuals =  array(rnorm(U, sd = sd), dim = d)
  new.Y = as.tensor(U + new.residuals)
  new.res = staccato(new.Y, dcomp_res$X_covar1, 
                     core_shape=r)
  
  return(abs(new.res$C_ts - B) > abs(B))
  
}

boot_p = function(n_boot, dcomp_res, n_thread,
                  lr.names, sender.names, receiver.names){
  
  boot.res = mclapply(1:n_boot, 
                      boot, 
                      dcomp_res = dcomp_res, 
                      mc.cores = n_thread)
  
  boot.p = (Reduce('+', boot.res) + 1) / (n_boot + 1)
  
  dimnames(boot.p) = list(colnames(dcomp_res$X_covar1),
                          lr.names,
                          sender.names,
                          receiver.names)
  
  return(boot.p)
  
}


sim_comm_tsr = function(seed=NA, 
                        whole_shape = c(20,20,20,20), 
                        core_shape = c(3,3,3,3),
                        p1,
                        dup, 
                        signal){
  
  d1 = whole_shape[1] ; d2 = whole_shape[2] ; d3 = whole_shape[3]  ; d4 = whole_shape[4]
  r1 = core_shape[1] ; r2 = core_shape[2] ; r3 = core_shape[3]; r4 = core_shape[4]
  p1 = p1; 
  
  #### warning for r should larger than 0
  if(r1<=0 | r2 <= 0|r3<= 0|r4<= 0){
    warning("the rank of coefficient tensor should be larger than 0",immediate. = T)
    return()
  }
  
  if(p1 > d1){
    warning("the number of covariates at each mode should be no larger than the dimension of the tensor",immediate. = T)
  }
  
  #### warning for p should larger than r
  if(p1<r1 & p1>0){
    warning("the rank of coefficient tensor should be no larger than the number of covariates",immediate. = T)
  }
  
  ####-------- generate data
  if(is.na(seed)==FALSE) set.seed(seed)
  X_covar1 = X_covar2 = X_covar3 = X_covar4 =  NULL
  
  if(p1<=0){
    X_covar1=diag(1,d1)
    p1=d1
  }else{
    X_covar1 = matrix(rnorm(d1*p1,mean = 0, sd = 1/sqrt(d1)),d1,p1)
  }

  #### if use block, p and r should larger than 1 of smaller than 0
  W1 =as.matrix(randortho(p1)[,1:r1])
  A =  ## factor matrix
  
  W2 = as.matrix(randortho(d2)[,1:r2]) ## factor matrix
  
  W3 = as.matrix(randortho(d3)[,1:r3]) ## factor matrix
  
  W4 = as.matrix(randortho(d4)[,1:r4]) ## factor matrix
  
  ### G: core tensor
  G = as.tensor(array(runif(r1*r2*r3*r4,min=-1,max=1),dim = c(r1,r2,r3,r4)))
  
  ### U: linear predictor
  U = ttl(G,list(X_covar1%*%W1,W2,W3,W4), ms = c(1,2,3,4))@data
  
  ## rescale to entrywise constraint
  G=G/max(abs(U))*signal 
  U=U/max(abs(U))*signal
  
  ## C_ts: coefficient
  C_ts=ttl(G,list(W1,W2,W3,W4),ms = c(1,2,3,4))@data 
  
  ### tsr:
  tsr = lapply(seq(dup), function(x) array(rnorm(d1*d2*d3*d4,U),dim = c(d1,d2,d3,d4)))

  return(list(tsr = tsr,
              X_covar1 = X_covar1, 
              W = list(W1 = W1, W2 = W2, W3 = W3, W4 = W4), 
              G=G@data, 
              U=U,
              C_ts=C_ts))
}

