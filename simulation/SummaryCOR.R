setwd("/home/xingjies/mammot/simulation/mammotPartVScommCOR") 
source("/home/xingjies/R/ggplot_theme_publication.R")

nz_ti <- c(1, 2)[2]
Ti = 3
method <- c("MAMMO", "CoMM", "PrediXcan", "TWAS")
h_z <- 0.01
ymax <- 1#ifelse(rhoX < 0.8, 0.5, 0.75)
gg_fdr <- NULL
gg_pow <- NULL
for(i in 1:5) {
  for (j in 1:3) {
    for (k in 1:3)  {

      h_c <- c(0.025, 0.05, 0.1, 0.2, 0.4)[i] 
      rhoW = c(0.2, 0.5, 0.8)[j] 
      rhoX <- c(0.2, 0.5, 0.8)[k]

      path <- paste0("part_hc", h_c*10, "_rhoX", rhoX*10, "_rhoW", rhoW*10, "nz_ti", nz_ti, "_batch-*")

      batch <- list.files(pattern=path)
      
      pvalue_mammot <- NULL
      pvalue_mammotpart <- NULL
      pvalue_comm <- NULL
      pvalue_predixcan <- NULL
      pvalue_TWAS <- NULL

      for (l in batch)  {
        res <- readRDS(l)
        pvalue_mammot <- rbind(pvalue_mammot, res$pvalue_mammot)
        pvalue_mammotpart <- rbind(pvalue_mammotpart, res$pvalue_mammotpart)
        pvalue_comm <- rbind(pvalue_comm, res$pvalue_comm)
        pvalue_predixcan <- rbind(pvalue_predixcan, res$pvalue_predixcan)
        pvalue_TWAS <- rbind(pvalue_TWAS, res$pvalue_TWAS)
      }

      Rep <- nrow(pvalue_TWAS)
      cutoff <- 0.05
      power <- c(#sum(pvalue_mammotpart[pvalue_mammot < cutoff/Rep, (1:nz_ti)] < 0.05/Ti)/Rep, 
             mean(pvalue_mammotpart[, (1:nz_ti)] < 0.05/Ti),
             mean(pvalue_comm[, (1:nz_ti)] < 0.05/Ti),  
             mean(pvalue_predixcan[, (1:nz_ti)] < 0.05/Ti), 
             mean(pvalue_TWAS[, (1:nz_ti)] < 0.05/Ti)
             )

      fpr <- c(#sum(pvalue_mammotpart[pvalue_mammot < cutoff/Rep, -(1:nz_ti)] < 0.05/Ti)/Rep, 
              mean(pvalue_mammotpart[, -(1:nz_ti)] < 0.05/Ti),
              mean(pvalue_comm[, -(1:nz_ti)] < 0.05/Ti), 
             mean(pvalue_predixcan[, -(1:nz_ti)] < 0.05/Ti),
             mean(pvalue_TWAS[, -(1:nz_ti)] < 0.05/Ti))

      gg_pow <- rbind(gg_pow, data.frame(method, power, h_c, rhoX, rhoW))

      gg_fdr <- rbind(gg_fdr, data.frame(method, fpr, h_c, rhoX, rhoW))

    }
  }
}

gg_fdr$method <- ordered(gg_fdr$method, levels = c("MAMMO", "CoMM", "PrediXcan", "TWAS"))
gg_fdr$h_c <- factor(gg_fdr$h_c)
#gg_fdr$nz_ti <- factor(gg_fdr$nz_ti)
#gg_fdr$s <- factor(gg_fdr$s)
    
gg_pow$method <- ordered(gg_pow$method, levels = c("MAMMO", "CoMM", "PrediXcan", "TWAS"))
gg_pow$h_c <- factor(gg_pow$h_c)
#gg_pow$nz_ti <- factor(gg_pow$nz_ti)
#gg_pow$s <- factor(gg_pow$s)


library(ggplot2)

outfile = paste0("../fpr_PartTestCOR_nzTis", nz_ti, ".png")

pl_fdr <- ggplot(gg_fdr, aes(x = h_c, y = fpr, fill = method)) + labs(x = expression(h[c]^2), y="FPR") + 
geom_bar(stat="identity", position=position_dodge())+
coord_cartesian(ylim=c(0, ymax)) + 
#facet_grid(s~nz_ti, labeller = label_both, scales = "free_y") + 
facet_grid(rhoW~rhoX, labeller = label_bquote(rows =  rho[W]: .(rhoW), cols= rho[X]: .(rhoX))) + 

scale_fill_manual(values=c("#FDBF6F", "#cccccc", "#969696", "#525252"),##FDBF6F
                  labels=c("TisCoMM", "CoMM", "PrediXcan", "TWAS")) + theme_Publication(base_size=28, border_col=NA)
ggsave(outfile, plot=pl_fdr, width = 360, height = 360, units = "mm", dpi=300, type = "cairo")


outfile = paste0("../pow_PartTestCOR_nzTis", nz_ti, ".png")

pl_pow <- ggplot(gg_pow, aes(x = h_c, y = power, fill = method)) + labs(x = expression(h[c]^2), y="Power") + 
geom_bar(stat="identity", position=position_dodge())+
coord_cartesian(ylim=c(0, 1)) + 
#facet_grid(s~nz_ti, labeller = label_both, scales = "free_y") + 
facet_grid(rhoW~rhoX, labeller = label_bquote(rows =  rho[W]: .(rhoW), cols= rho[X]: .(rhoX))) + 
 
scale_fill_manual(values=c("#FDBF6F", "#cccccc", "#969696", "#525252"),
                  labels=c("TisCoMM", "CoMM", "PrediXcan", "TWAS")) +theme_Publication(base_size=28, border_col=NA)
ggsave(outfile, plot=pl_pow, width = 360, height = 360, units = "mm", dpi=300, type = "cairo")

