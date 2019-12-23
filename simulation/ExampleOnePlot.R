source("/home/xingjies/mammot/qqplot5.R")
source("/home/xingjies/R/ggplot_theme_publication.R")
############################################################################################################
############################################################################################################
################################################ QQ plot (type 1 error) ####################################
############################################################################################################
############################################################################################################

for(k in 1:5){

case <- 1
plot_list <- vector("list", 9)

	for (j in 1:3) {
		for (i in 1:3)	{

			h_z <- c(0, 0.01)[1]
  			h_c <- c(0.025, 0.05, 0.1, 0.2, 0.4)[k] 
			rhoX = c(0.2, 0.5, 0.8)[i]
			s = c(0.1, 0.5, 0.9)[j] 


			batch <- list.files(pattern = paste0("pvalue_hz", h_z*10, "_hc", h_c*10, "_rhoX", rhoX*10, "_s", s*10, "_batch-*"))
			pvalue <- NULL
			for (l in batch)	{
				pvalue <- rbind(pvalue, read.table(l, header=F))
			}

			pvalue[is.na(pvalue)] <- 1

 			pv <- list(MAMMOT = pvalue[, 1], 
 				      MAMMOT_S = pvalue[, 2], 
 				      MultiXcan = pvalue[, 3],
			 		  MultiXcan_S = pvalue[, 4],
			  		  UTMOST = pvalue[, 5])
	
 			plot_list[[case]] <- qqunif.plot(pv, xlab="Expected: -log10(p-value)", ylab = "Observed: -log10(p-value)",
 											mytitle = bquote((s~","~rho[X])==(.(s)~","~.(rhoX))))
 			case <- case + 1
		}
	}
 

library(gridExtra) # also loads grid
outfile = paste0("../example1_qqplot_", "h_c_", k,".png")
png(filename = outfile, width = 480, height = 480, units = "mm", res=200, type = "cairo")
grid.arrange(grobs=plot_list, nrow = 3)
#grid_arrange_shared_legend(grobs=plot_list, nrow = 3, position='bottom')

dev.off()
}

############################################################################################################
############################################################################################################
#################################################### power plot ############################################
############################################################################################################
############################################################################################################
method <- c("MAMMO", "S-MAMMO", "MultiXcan", "S-MultiXcan", "UTMOST")

h_z <- c(0, 0.01)[2]
ymax <- ifelse(h_z==0, 0.2, 1)
gg_df <- NULL
for(i in 1:5)	{
	for (j in 1:3) {
		for (k in 1:3)	{
			h_c <- c(0.025, 0.05, 0.1, 0.2, 0.4)[i] 
			rhoX = c(0.2, 0.5, 0.8)[j]
			s = c(0.1, 0.5, 0.9)[k] 
			batch <- list.files(pattern = paste0("pvalue_hz", h_z*10, "_hc", h_c*10, "_rhoX", rhoX*10, "_s", s*10, "_batch-*"))
			pvalue <- NULL
			for (l in batch)	{
				pvalue <- rbind(pvalue, read.table(l, header=F))
			}

			pvalue[is.na(pvalue)] <- 1
			power <- apply(pvalue < 0.05, 2, mean, na.rm=T) 
			gg_df <- rbind(gg_df, data.frame(method, power, h_c, rhoX, s))

		}
	}
}
  	

gg_df$method <- ordered(gg_df$method, levels = c("MAMMO", "S-MAMMO", "MultiXcan", "S-MultiXcan", "UTMOST"))
gg_df$h_c <- factor(gg_df$h_c)
library(ggplot2)

outfile = paste0("../power_hz", h_z*100, ".png")

pl <- ggplot(data=gg_df, aes(x = h_c, y = power, fill = method)) + labs(x = expression(h[c]^2), y="Power") + 
geom_bar(stat="identity", position=position_dodge())+
coord_cartesian(ylim=c(0, ymax)) + 
facet_grid(s~rhoX, labeller = label_bquote(rows =  s: .(s), cols= rho[X]: .(rhoX))) + 
scale_fill_manual(values=c("#FDBF6F", "#FF7F00", "#A6CEE3", "#1F78B4","#6A3D9A"),
				  labels=c("TisCoMM", expression(TisCoMM-S^2), "MultiXcan", "S-MultiXcan", "UTMOST")) +
#scale_colour_manual(labels=c("TisCoMM", "S-MAMMO", "MultiXcan", "S-MultiXcan", "UTMOST")) +
theme_Publication(base_size=30, border_col=NA)

ggsave(outfile, plot=pl, width = 360, height = 360, units = "mm", dpi=300, type = "cairo")

