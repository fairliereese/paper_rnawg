### Cross-cell-type prediction using the logit model trained on the other cell type
### Ex: test set=chr2 & chr3 in K562, train set=chr1 & chr4-22 of GM12878, model = logit(label~DHS+TPM+peak_length)

library(dplyr)
library(reshape2)
library(stringr)
library(PRROC)

# LR_type <- "Cerberus"
# LR_type <- "LAPA"
LR_type <- commandArgs(trailingOnly = T)[1]

data_dir <- "../data_dir/labeled_beds/split_chr_with_DHS"
test_dir <- file.path(data_dir, LR_type, "test")
eval_dir <- file.path(data_dir, LR_type, "test_with_labels")

training_dir <- "training"
model_dir <- file.path("models/individual_exper_same_cell_line", LR_type)
pred_dir <- file.path("predicts/individual_exper_diff_cell_line", LR_type)

plot_dir <- file.path("plots/individual_exper_diff_cell_line/", LR_type)

test_beds <- list.files(test_dir, pattern = "*.bed", full.name = T)

experiments <- gsub("\\..*", "", basename(test_beds))
models <- setNames(lapply(file.path(model_dir, paste0(experiments, ".RDS")), readRDS), experiments)

get_other_cell_line <- function (cell_line) {
	if ("K562" %in% cell_line) {
		other_cell_line <- "GM12878_ENCSR962BVU"
	} else {
		other_cell_line <- "K562_ENCSR589FUJ"
	}
	other_cell_line
}

auc_df <- data.frame()
for (test_bed in test_beds) {
		print(test_bed)
        df <- read.delim2(test_bed, header = F)
        colnames(df) <- c("chr", "start", "end", "peak_id", "TPM", "DHS")
        df <- df %>% mutate(across(c("TPM", "DHS"), as.numeric)) %>%
                        mutate(length = end - start, .after = peak_id)
        to_log <- c("TPM", "DHS", "length")
        for (col in to_log) {
                df[[col]] <- log2(df[[col]] + 1)
        }
        experiment <- gsub("\\..*", "", basename(test_bed))
		cell_line <- str_extract(experiment, "GM12878|K562")
		############
		other_cell_line <- get_other_cell_line(cell_line)
		other_cell_line_model <- models[[other_cell_line]]
#		preds <- as.data.frame(lapply(other_cell_line_models, predict, df, type = "response"))
#		avg_pred <- rowMeans(preds)
		pred <- predict(other_cell_line_model, df, type="response")
		pred_file <- file.path(pred_dir, paste0(experiment, ".out"))
		pred_df <- df %>% select(chr, start, end, peak_id) %>% mutate(prediction = pred)
		write.table(pred_df, pred_file, row.names = F, col.names = F, quote = F, sep = "\t")
		############################# eval:
		eval_file <- file.path(eval_dir, paste0(experiment, ".testlabels.bed"))
		binary_labels <- read.table(eval_file, header = F) %>% pull(tail(names(.), 1))

		PRROC_obj <- roc.curve(scores.class0 = pred, weights.class0 = binary_labels, curve = T)
		auroc_plot_file <- file.path(plot_dir, paste0(experiment, "_AUROC_diff_cell_line.pdf"))
		pdf(auroc_plot_file, 7, 7)
		#auroc_plot_title <- paste("Train on", other_cell_line, ", prediction on", experiment, "-", LR_type, "- AUROC:", round(PRROC_obj$auc, 2))
		other_cell_line_name <- str_extract(other_cell_line, "GM12878|K562")
		auroc_plot_title <- paste("Train on", other_cell_line_name, ", prediction on", cell_line, "-", LR_type,
								  "\nAUROC:", round(PRROC_obj$auc, 3))
		plot(PRROC_obj, color=F, main = auroc_plot_title,
			 auc.main = F, cex.lab=1.5, cex.axis=1.5, cex.main=1.3, cex.sub=1.3)
		dev.off()
		auc_df <- rbind(auc_df, c(experiment, PRROC_obj$auc))
}
colnames(auc_df) <- c("Experiment", "AUROC")
auc_df <- auc_df %>% mutate(AUROC = as.numeric(AUROC))
write.table(auc_df, paste0(pred_dir, "/auc.txt"), col.names = T, row.names = F, quote = F, sep = "\t") 
