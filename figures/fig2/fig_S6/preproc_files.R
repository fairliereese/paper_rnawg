### change column names based on the input type and
# filter out peaks based on chromosomal location and
# sort .bed files based on chr location

library(dplyr)
library(stringr)

# data_type <- "LAPA"
# data_type <- "Cerberus"
data_type <- commandArgs(trailingOnly = T)[1]

data_dir <- "data_dir/"
bed_dir <- paste0(data_dir, "raw_beds/", data_type, "/") 
output_dir <- paste0(data_dir, "processed_beds/", data_type, "/")

hd <- (data_type == "Cerberus")

beds <- list.files(bed_dir, "*.bed", full.names = T)
peaks_lst <- lapply(beds, read.delim2, header = hd)
names(peaks_lst) <- gsub(".bed", "", basename(beds))
if (data_type == "LAPA") {
        peaks_lst <- lapply(peaks_lst, function (df) {
        colnames(df) <- c("chr", "start", "end", "polyA_site", "count", "strand", "feature",
                                         "gene_id", "TPM", "gene_count", "fracA", "annotated_site")
        df %>% mutate(across(c("TPM", "fracA"), as.numeric))
	})
} else {
        peaks_lst <- lapply(peaks_lst, function (df) {
        colnames(df) <- tolower(colnames(df))
		colnames(df)[1] <- "chr"
        colnames(df)[which(colnames(df) == "tpm")] <- "TPM"
	    df %>% mutate(across("TPM", as.numeric))
                #mutate(length = end - start, .after = end)
#				mutate(peak_id = paste(chr, start, end, sep = "_"))
	})
}

valid_chroms <- paste0("chr", c(1:22, "X", "Y"))

peaks_lst <- lapply(peaks_lst, function (df) {
	df <- df %>% mutate(length = end - start, .after = end) %>%
			mutate(peak_id = paste(chr, start, end, sep = "_"), .after = end) %>%
			filter(chr %in% valid_chroms)
	df <- df[!duplicated(df$peak_id), ]
	df <- df %>% select(chr, start, end, peak_id, TPM)
})

for (f in names(peaks_lst)) {
	fname <- paste0(output_dir, "/", f, ".bed")
	tmp_fname <- paste0(fname, ".tmp")
	write.table(peaks_lst[[f]], tmp_fname, sep = "\t", quote = F, col.names = F, row.names = F)
	cmd1 <- sprintf("sort -k1,1 -k2,2n %s > %s", tmp_fname, fname)
	system(cmd1)
	cmd2 <- paste("rm", tmp_fname)
	system(cmd2)
}

