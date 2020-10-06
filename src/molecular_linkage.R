#!/usr/bin/env Rscript

# load packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(furrr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tictoc))
suppressPackageStartupMessages(library(dbscan))

# parse arguments
args <- commandArgs(trailingOnly = T)

# set up parallel sessions
plan(multisession, workers = args[4])

# read and pre-process dataframes
read <- function(filename, kvalue) {
	
	df <- fread(filename, sep = "\t", header = F)
	colnames(df) <- c("reference", "query", "dist", "pval", "jaccard")
	
	df <- df %>% 
		separate(jaccard, into = c("shared_k", "total_k"), sep = "/") %>% 
		mutate(kvalue = kvalue,
					 shared_k = as.numeric(shared_k),
					 total_k = as.numeric(total_k),
					 pmatch = shared_k/total_k) %>% 
		select(reference, kvalue, pmatch, dist)
		
	
	return(df)
}

# merge mash results with k=13,15,17,19,21,23
load <- function(dir, sample_name) {
	
	k <- c(13,15,17,19,21,23)
	df <- future_map_dfr(k, function(x) read(paste(paste(dir, sample_name, sep = "/"), x, "tsv", sep = "."), x))
	
	# identify irrelevant strains
	ref_ids <- df %>% 
		filter(kvalue == 21, 
					 dist <= 0.01) %>% 
		pull(reference)
	# filter irrelevant strains and create nested dataframes
	df <- df %>% 
		filter(reference %in% ref_ids) %>% 
		select(-dist) %>% 
		group_by(reference) %>% 
		nest() %>% 
		mutate(reference = gsub(".*/|.fasta|.fna", "", reference)) %>%  # remove paths in reference ids
		ungroup()
	
	return(df)
}

# construct linear fit
linear_fit <- function(df) {
	
	# log transform pmatch
	df <- df %>% 
		mutate(pmatch = log(pmatch))
	# design matrix
	dmatrix <- matrix(df$kvalue, length(df$kvalue), 1)
	# linear fit
	linear.fit <- lm.fit(cbind(1, dmatrix), df$pmatch)
	intercept <- linear.fit$coefficients[[1]]
	slope <- linear.fit$coefficients[[2]]
	
	return(paste(intercept, slope, sep = ","))
}

# compute accessory and core genome distances
compute <- function(df, linear_regression) {
	
	distances <- df %>% 
		cbind(linear_fit = linear_regression) %>% 
		separate(linear_fit, into = c("intercept", "slope"), sep = ",") %>% 
		mutate(intercept = as.numeric(intercept),
					 slope = as.numeric(slope),
					 accessory = 1-exp(1)^intercept,
					 core = 1-exp(1)^slope,
					 core = case_when(core < 0 ~ 0,
					 								 T ~ core),
					 accessory = case_when(accessory < 0 ~ 0,
					 											T ~ accessory))
	
	return(distances)
}

# optimize hdbscan hyperparameters
optimize <- function(df) {
	
	optimize_tbl <- tibble(minPts = seq(5,105,10),
												 low_proportion = 0)
	optimize_tbl$low_proportion <- map_dbl(seq(5,105,10), function(x) {
		
		# cluster data
		cluster_res <- hdbscan(df, minPts = x)
		# compute proportion of low confidence assignments (p < 0.05)
		low_count <- length(cluster_res$membership_prob[which(cluster_res$membership_prob < 0.05)])
		total_count <- length(cluster_res$membership_prob)
		low_proportion <- low_count/total_count
		
		return(low_proportion)
	})
	
	# identify optimal minPts
	opt_minPts <- optimize_tbl %>% 
		filter(low_proportion == min(low_proportion)) %>% 
		arrange(-desc(minPts)) %>% 
		slice(1) %>% 
		pull(minPts)
	
	p <- optimize_tbl %>% 
		ggplot(aes(x = minPts, y = low_proportion)) +
		geom_line()+
		ylab("Fraction of Low Probability Assignments")+
		xlab("minPts")
	#ggsave(paste0(args[1], "/", args[2], "_optimization.png"), p)
	
	return(opt_minPts)
}

# plot cluster assignments 
plot_cluster <- function(df, cutoff) {
	
	plot <- ggplot() +
		geom_point(data = df %>% 
							 	filter(cluster == 0),
							 mapping = aes(x = core, y = accessory), 
							 alpha = 0.1, 
							 color = "black") +
		geom_point(data = df %>% 
							 	filter(cluster != "0"),
							 mapping = aes(x = core, y = accessory, color = cluster),
							 alpha = 0.1) +
		stat_ellipse(data = df %>% 
							 	filter(cluster != "0"),
							 geom = "polygon",
							 alpha = 0.2,
							 mapping = aes(x = core, y = accessory, color = cluster, fill = cluster)) +
		geom_label(data = df %>% 
							 	filter(cluster != "0") %>% 
							 	group_by(cluster) %>% 
							 	summarize(avg_core = mean(core), # average X and Y to find cluster centroid
							 						avg_accessory = mean(accessory)),
							 mapping = aes(x = avg_core, y= avg_accessory, label = cluster, color = cluster))+
		#geom_abline(slope = cutoff[[1]], intercept = cutoff[[2]], color = "red") +
		#geom_abline(slope = -0.029, intercept = 0.016, color = "red") +
		geom_hline(yintercept = cutoff[[2]], color = "red") +
		geom_vline(xintercept = cutoff[[1]], color = "red") +
		xlab("Core genome distance") +
		ylab("Accessory genome distance") +
		guides(color = F,
					 fill = F) +
		theme_minimal() #+
	##scale_x_continuous(limits = c(0,0.001)) +
	#scale_y_continuous(limits = c(0,0.05))
	
	return(plot)
}

# identify threshold cut-offs defined by the cluster closest to origin
threshold <- function(df) {
	
	# determine the two cluster ids closest to origin
	cluster_dist <- df %>% 
		filter(cluster != 0) %>% # remove noise
		group_by(cluster) %>% 
		summarize(avg_core = mean(core), # average X and Y to find cluster centroid
							avg_accessory = mean(accessory)) %>% 
		mutate(distance = sqrt(avg_core^2 + avg_accessory^2))
	
	cluster_id <- cluster_dist %>% 
		arrange(-desc(distance)) %>% 
		slice(1:2) %>% 
		pull(cluster)
	
	# calculate slope, midpoint and intercept between the two closest clusters
	# cluster1 <- as.numeric((cluster_dist %>% filter(cluster == cluster_id[1]) %>% select(avg_core, avg_accessory))[1,])
	# cluster2 <- as.numeric((cluster_dist %>% filter(cluster == cluster_id[2]) %>% select(avg_core, avg_accessory))[1,])
	# 
	# slope <- ((cluster2[2]-cluster1[2])/(cluster2[1]-cluster1[1]))
	# midpoint <- c((cluster2[1]+cluster1[1])/2, (cluster2[2]+cluster1[2])/2)
	# print(midpoint)
	# intercept <- midpoint[2]-slope*midpoint[1]
	
	#return(list(slope, intercept))
	
	# determine max core and accessory boundaries of cluster closest to origin
	max_core <- df %>% 
		filter(cluster == cluster_id[1]) %>% 
		pull(core) %>% 
		max()
	
	max_accessory <- df %>% 
		filter(cluster == cluster_id[1]) %>% 
		pull(accessory) %>% 
		max()
	
	return(list(max_core, max_accessory))
}

# filter entire database given cutoff
database_filter <- function(df, cutoff) {
	
	neighbours <- df %>% 
		filter(core <= cutoff[[1]],
					 accessory <= cutoff[[2]])
	
	return(neighbours)
}

# write output
output <- function(df, output_path) {
	
	neighbours_id <- df %>% 
		pull(reference)
	
	write.table(neighbours_id, file = output_path, quote = F, col.names = F, row.names = F)
}

# main
main <- function(dir, sample_name, out_dir) {
	
	tic()
	# load data
	message(paste("Loading mash results for", sample_name))
	mash_res <- load(dir, sample_name)
	# fit linear model
	message("Fitting linear model")
	linear.fit <- future_map_chr(mash_res$data, ~linear_fit(.))
	# compute accessory and core genome distances
	mash_res <- compute(mash_res, linear.fit)
	# identify optimal hyperparameters for hdbscan
	message("Optimizing hyperparameters for clustering")
	opt_minpts <- optimize(mash_res[c("core", "accessory")])
	# cluster data using optimal parameters
	cluster_res <- hdbscan(mash_res[c("core", "accessory")], minPts = opt_minpts)
	mash_res <- cbind(mash_res, cluster = as.character(cluster_res$cluster))
	# identify threshold cut-off
	cutoff <- threshold(mash_res)
	# plot cluster assignments
	p <- plot_cluster(mash_res, cutoff)
	ggsave(paste0(out_dir, "/", sample_name, "_cluster_res.png"), p, width = 10, height = 10)
	# filter database
	mash_res_filter <- database_filter(mash_res, cutoff)
	# write out neighbour ids
	output(mash_res_filter, paste0(out_dir, "/", sample_name, "_neighbours.list"))
	toc()
}

# call main
main(args[1], args[2], args[3])
