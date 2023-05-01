### Combining Cerberus replicates and mean-normalizing the expression signal

library(dplyr)
library(stringr)

data_type <- "Cerberus"

data_dir <- "data_dir/"
bed_dir <- paste0(data_dir, "raw_beds/Cerberus_replicates/")
output_dir <- paste0(data_dir, "raw_beds/Cerberus/")


repl_conv_file <- paste0(bed_dir, "replicate_to_experiment.txt")
repl_conv <- read.delim(repl_conv_file, header = F)
colnames(repl_conv) <- c("Experiment", "Accession")
rownames(repl_conv) <- repl_conv$Experiment


beds <- list.files(bed_dir, "*.bed", full.names = T)
peaks_lst <- lapply(beds, read.delim2, header = T)
names(peaks_lst) <- gsub(".bed", "", basename(beds))
peaks_lst <- lapply(peaks_lst, function (df) {
	df %>% mutate(across("tpm", as.numeric)) %>% mutate(dataset = toupper(dataset))
})

all_repls <- names(peaks_lst)

get_exper_name <- function(repl_name) {
	gsub("_\\d+$", "", repl_name)
}

get_exper_accession <- function(exper) {
	repl_conv[exper, "Accession"]	
}

get_exper_repls <- function(exper) {
	all_repls[which(get_exper_name(all_repls) == exper)]
}

get_cell_line_name <- function(sample) {
	gsub("_.*", "", sample)
}

all_experiments <- unique(get_exper_name(names(peaks_lst)))

for (exper in all_experiments) {
	print(exper)
	repls <- get_exper_repls(exper)
	repls_dfs <- peaks_lst[repls]
	all_regions <- data.table::rbindlist(repls_dfs) %>% select(-tpm, -dataset)
	all_regions <- all_regions[!duplicated(all_regions), ]
	
	agg0 <- data.table::rbindlist(repls_dfs) %>% select(Name, tpm, dataset)
	repl_means <- aggregate(tpm ~ dataset, agg0, mean) %>% pull(tpm, dataset)
	all_means <- mean(agg0$tpm)
	agg0 <- agg0 %>% mutate(normalized_tpm = agg0$tpm * (all_means / repl_means[agg0$dataset]))  # mean-normalizing every replicate
	agg <- aggregate(normalized_tpm ~ Name, agg0 %>% select(Name, normalized_tpm), mean)

	exper_df <- merge(all_regions, agg, by = "Name") %>% mutate(dataset = exper)
	colnames(exper_df)[which(colnames(exper_df) == "normalized_tpm")] <- "tpm"
	exper_df <- exper_df %>% relocate(tpm, .after = gene_id) %>% relocate(Name, .after = End) %>%
			arrange(Name) %>% mutate(tpm = round(tpm, 2))
	exper_full_name <- paste0(get_cell_line_name(exper), "_", get_exper_accession(exper))
	output_file <- paste0(output_dir, exper_full_name, ".bed")
	write.table(exper_df, output_file, row.names = F, col.names = T, sep = "\t", quote = F)
}


