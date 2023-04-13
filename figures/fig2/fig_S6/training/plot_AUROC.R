# Boxplot for AUROC of cross-cell-type predictions  (SFig 6f)

library(dplyr)
library(reshape2)
library(stringr)
library(ggplot2)
library(RColorBrewer)

theme_set(theme_classic(base_size = 16))
theme_update(axis.text = element_text(color = "black"),
                axis.ticks.length = unit(0.3, "cm"),
                plot.title = element_text(hjust = 0.3),
                panel.border = element_rect(color = "black", fill = NA, size = 1))


data_type_dir <- "predicts/individual_exper_diff_cell_line"
pdf_file <- paste0("plots/", "AUROC_diff_cell_line_boxplot.pdf")

auroc_lapa <- read.delim2(file.path(data_type_dir, "LAPA/auc.txt"))
auroc_cerb <- read.delim2(file.path(data_type_dir, "Cerberus/auc.txt"))

auroc_df <- rbind(auroc_lapa %>% mutate(Data_Type = "LAPA"), auroc_cerb %>% mutate(Data_Type = "Cerberus")) %>%
		mutate(AUROC = as.numeric(AUROC))

plot_df <- auroc_df

png_file <- gsub(".pdf", ".png", pdf_file)
pal <- brewer.pal(n = 6, "Set2")
cols <- setNames(pal[5:6], c("GM12878", "K562"))
ggplot(plot_df) +
		geom_boxplot(aes(x = Data_Type, y = AUROC, fill = str_extract(Experiment, "GM12878|K562")), width = 0.4, color = "black") +
		scale_fill_manual(values = cols, name = "Prediction Cell Line") +
		scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
		labs(x = "", y = "AUROC", title = "Cross-cell type prediction")
		# labs(x = "", y = "AUROC", title = "Same-cell type prediction")
ggsave(pdf_file, width = 7, height = 7)
ggsave(png_file, width = 7, height = 7)


