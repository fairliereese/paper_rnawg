### Combining LAPA replicates and mean-normalizing the expression signal

library(dplyr)
library(stringr)

data_type <- "LAPA"

data_dir <- "data_dir/"
bed_dir <- paste0(data_dir, "raw_beds/LAPA_replicates/")
output_dir <- paste0(data_dir, "raw_beds/LAPA/")


repl_conv_file <- paste0(bed_dir, "replicate_to_experiment.txt")
repl_conv <- read.delim(repl_conv_file, header = F)
colnames(repl_conv) <- c("Experiment", "Accession")
rownames(repl_conv) <- repl_conv$Experiment


beds <- list.files(bed_dir, "*.bed", full.names = T)
peaks_lst0 <- lapply(beds, read.delim, header = F)
names(peaks_lst0) <- gsub(".bed", "", basename(beds))
# peaks_lst <- lapply(peaks_lst, function (df) {
peaks_lst <- lapply(names(peaks_lst0), function (repl_name) {
	df <- peaks_lst0[[repl_name]]
	colnames(df) <- c("chr", "start", "end", "polyA_site", "count", "strand", "feature",
                                         "gene_id", "TPM", "gene_count", "fracA", "annotated_site")
    df %>% mutate(across(c("TPM", "fracA"), as.numeric)) %>%
		mutate(dataset = gsub("_\\d$", "", repl_name)) %>%
		mutate(Name = paste(chr, start, end, gene_id, sep = "_")) %>%
		filter(feature == "five_prime_utr")
})

names(peaks_lst) <- names(peaks_lst0)
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
	all_regions <- data.table::rbindlist(repls_dfs) %>% select(-TPM, -dataset)
	all_regions <- all_regions[!duplicated(all_regions), ]
	
	agg0 <- data.table::rbindlist(repls_dfs) %>% select(Name, TPM, dataset)
	repl_means <- aggregate(TPM ~ dataset, agg0, mean) %>% pull(TPM, dataset)
	all_means <- mean(agg0$TPM)
	agg0 <- agg0 %>% mutate(normalized_TPM = agg0$TPM * (all_means / repl_means[agg0$dataset]))  # mean-normalizing every replicate
	agg <- aggregate(normalized_TPM ~ Name, agg0 %>% select(Name, normalized_TPM), mean)

	exper_df <- merge(all_regions, agg, by = "Name") %>% mutate(dataset = exper)
	colnames(exper_df)[which(colnames(exper_df) == "normalized_TPM")] <- "TPM"
	exper_df <- exper_df %>% relocate(TPM, .after = gene_id) %>% relocate(Name, .after = end) %>%
			arrange(Name) %>% mutate(TPM = round(TPM, 2))
	exper_full_name <- paste0(get_cell_line_name(exper), "_", get_exper_accession(exper))
	output_file <- paste0(output_dir, exper_full_name, ".bed")
	write.table(exper_df, output_file, row.names = F, col.names = T, sep = "\t", quote = F)
}



