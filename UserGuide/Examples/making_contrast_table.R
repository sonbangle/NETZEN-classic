contrast=paste0(c(rep(0,9),1), collapse = ",")
network ="GBM_net.csv"
contrast_network= data.frame(contrast= contrast, network=network)
write.table(contrast_network, file="contrast_network_table.csv", quote=FALSE, row.names = FALSE, sep="\t")