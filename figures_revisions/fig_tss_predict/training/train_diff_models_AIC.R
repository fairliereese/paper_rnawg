### Training GLM (logit) models based on all 2^3-1 = 7 combinations of DHS, TPM expr and peak length
### Calculating AIC and finding the best model overall

library(dplyr)
library(reshape2)
library(stringr)
library(ggplot2)

theme_set(theme_classic(base_size = 16))
theme_update(axis.text = element_text(color = "black"),
                axis.ticks.length = unit(0.3, "cm"),
                plot.title = element_text(hjust = 0.3),
                panel.border = element_rect(color = "black", fill = NA, size = 1))


# LR_type <- "Cerberus"
# LR_type <- "LAPA"
LR_type <- commandArgs(trailingOnly = T)[1]

train_dir <- file.path("../data_dir/labeled_beds/split_chr_with_DHS", LR_type, "train")
# model_dir <- file.path("training", "models/diff_models_AIC_individual_exper_same_cell_line", LR_type)
model_dir <- file.path("models/diff_models_AIC_individual_exper_same_cell_line")
print(model_dir)  ###
train_beds <- list.files(train_dir, pattern = "*.bed", full.name = T)
print(train_beds) ###
model_params_list <- list(c("TPM", "DHS", "length"),
					 c("TPM", "DHS"),
					 c("TPM", "length"),
					 c("TPM"),
					 c("DHS", "length"),
					 c("DHS"),
					 c("length"))
fm_to_str <- function(fm) {
	ch <- as.character(fm)  # "~", "label", "TPM + DHS + length"
	fm_str <- paste(unlist(str_split(ch[3], " \\+ ")), collapse = "_")
	return(fm_str)
}

aic_df <- data.frame()
for (model_params in model_params_list) {
	fm <- formula(paste("label ~", paste(model_params, collapse = " + ")))
	fm_str <- fm_to_str(fm)
	
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
		model <- glm(formula = fm, data=df, family = "binomial")
		aic <- model$aic
		experiment <- gsub("\\..*", "", basename(train_bed))
		aic_df <- rbind(aic_df, c(experiment, fm_str, aic))
	}
}
colnames(aic_df) <- c("Experiment", "Model", "AIC")
print(aic_df)
aic_df <- aic_df %>% mutate(AIC = as.numeric(AIC))
write.table(aic_df, file.path(model_dir, paste0(LR_type, "_AIC.out")), row.names = F , col.names = T, quote = F, sep = "\t")

aic_rank_df <- data.frame()
for (train_bed in train_beds) {
	experiment <- gsub("\\..*", "", basename(train_bed))
	exper_df <- aic_df %>% filter(Experiment == experiment) %>%
			mutate(AIC = as.numeric(AIC))
	exper_df <- exper_df %>% mutate(Order = rank(AIC))
	aic_rank_df <- rbind(aic_rank_df, exper_df)
}
agg <- aggregate(Order ~ Model, aic_rank_df, mean) %>% arrange(desc(Order)) %>%
		mutate(Model = factor(Model, levels = Model))

plot_png <- file.path("plots", paste0(LR_type, "_AIC.png"))
plot_pdf <- file.path("plots", paste0(LR_type, "_AIC.pdf"))

ggplot(agg) + geom_bar(aes(x = Model, y = Order), size = 0.6, fill="chocolate3", color="black", stat = "identity") + 
	scale_y_continuous(breaks = scales::pretty_breaks(n=7)) +
	labs(x = "Logit Model Parameters",
		 y = paste0("Avg. AIC Rank (", length(train_beds), " Samples)"),
		 title = paste("AIC of Logit Models Trained on", LR_type)) +
	theme(axis.text.x=element_text(angle = +30, vjust = 0.5))
ggsave(plot_png, width=7, height=7)
ggsave(plot_pdf, width=7, height=7)



