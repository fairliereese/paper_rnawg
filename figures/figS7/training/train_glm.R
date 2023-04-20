### After training_diff_models.R and finding that logit(label~DHS+TPM+length) is the best model
### Training logit(label~DHS+TPM+peak_length)


library(dplyr)
library(reshape2)
library(stringr)

# LR_type <- "LAPA"
# LR_type <- "Cerberus"
LR_type <- commandArgs(trailingOnly = T)[1]

train_dir <- file.path("../data_dir/labeled_beds/split_chr_with_DHS", LR_type, "train")
model_dir <- file.path("models/individual_exper_same_cell_line", LR_type)

train_beds <- list.files(train_dir, pattern = "*.bed", full.name = T)

for (train_bed in train_beds) {
	df <- read.delim2(train_bed, header = F)
	colnames(df) <- c("chr", "start", "end", "peak_id",
					  "TPM", "DHS", "full_label", "label")
	df <- df %>% mutate(across(c("TPM", "DHS"), as.numeric)) %>%
			mutate(length = end - start, .after = peak_id)
	to_log <- c("TPM", "DHS", "length")
	for (col in to_log) {
		df[[col]] <- log2(df[[col]] + 1)
	}
	model <- glm(label ~ TPM + DHS + length, data=df, family = "binomial")
	experiment <- gsub("\\..*", "", basename(train_bed))
	model_file <- file.path(model_dir, paste0(experiment, ".RDS"))
	saveRDS(model, model_file)
}



