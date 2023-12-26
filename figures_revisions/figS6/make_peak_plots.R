library(stringr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

theme_set(theme_classic(base_size = 16))
theme_update(axis.text.x = element_text(vjust = -0.15),
                         axis.text = element_text(color = "black"),
                         axis.ticks.length = unit(0.5, "cm"),
                         plot.title = element_text(hjust = 0.5),
                         panel.border = element_rect(color = "black", fill = NA, size = 1))

peaks <- c()
read_peaks_lst <- function (LR_type) {
bed_dir <- paste0("data_dir/labeled_beds/all_chr/", LR_type, "/")

beds <- list.files(bed_dir, "*.bed", full.names = T)
peaks_lst <- lapply(beds, read.delim2, header = F)
names(peaks_lst) <- gsub(".bed", "", basename(beds))
peaks_lst <- lapply(peaks_lst, function (df) {
	colnames(df) <- c("chr", "start", "end", "peak_id", "TPM", "full_label", "label")
	df <- df %>% mutate(TPM = as.numeric(TPM)) %>%
		mutate(length = end - start, .after = peak_id)
})
return(peaks_lst)
}
#################################
cerb_peaks_lst <- read_peaks_lst("Cerberus")
names(cerb_peaks_lst) <- paste0(gsub("\\..*", "", gsub("Allchr.", "", names(cerb_peaks_lst))), ".Cerberus")
lapa_peaks_lst <- read_peaks_lst("LAPA")
names(lapa_peaks_lst) <- paste0(gsub("\\..*", "", gsub("Allchr.", "", names(lapa_peaks_lst))), ".LAPA")
peaks_lst <- c(cerb_peaks_lst, lapa_peaks_lst)

cerb_samps <- c("GM12878_ENCSR962BVU.Cerberus", "K562_ENCSR589FUJ.Cerberus")
lapa_samps <- c("GM12878_ENCSR962BVU.LAPA", "K562_ENCSR589FUJ.LAPA")
samps <- c(cerb_samps, lapa_samps)
pal <- brewer.pal(n = 6, "Set2")

######
overlap_status <- as.data.frame(lapply(peaks_lst, function (df) {setNames(as.numeric(table(df$label)), c("No Overlap", "CAGE or RAMPAGE"))})) %>%
		t() %>% as.data.frame() %>% mutate(sample = rownames(.))
m_overlap_status <- melt(overlap_status)
cols <- c("CAGE or RAMPAGE" = pal[1], "No Overlap" = pal[2])
plot_df <- m_overlap_status %>% mutate(cell_line = gsub("_.*", "", sample))
ggplot(plot_df %>% filter(sample %in% samps)) +
	   geom_bar(aes(x = sample, y = value, fill = variable), position = "stack", stat = "identity") +
	   scale_fill_manual(values = cols, name = "TSS Peak Overlap") +
	   scale_y_continuous(breaks = scales::pretty_breaks(n=7)) +
	   labs(x = "", y = "No. TSS Peaks", title = "Long-read overlap with RAMPAGE or CAGE") +
	   theme(axis.text.x=element_text(angle = +60, size = 12, vjust = +1, hjust=+1))
overlap_plot_pdf <- "plots/overlap.pdf"
overlap_plot_png <- gsub(".pdf", ".png", overlap_plot_pdf)
ggsave(overlap_plot_pdf, width = 7, height = 6)
ggsave(overlap_plot_png, width = 7, height = 6)
#######
test_train_size <- as.data.frame(lapply(peaks_lst, function(df) {c(df %>% filter(chr %in% c("chr2", "chr3")) %>% nrow, df %>% filter(!(chr %in% c("chr2", "chr3"))) %>% nrow)})) %>% t %>% as.data.frame %>% mutate(sample = rownames(.))
colnames(test_train_size)[1:2] <- c("Test Set Peaks", "Training Set Peaks")
m_test_train <- melt(test_train_size)
cols <- setNames(pal[3:4], colnames(test_train_size)[1:2])
plot_df <- m_test_train
ggplot(plot_df %>% filter(sample %in% samps)) +
       geom_bar(aes(x = sample, y = value, fill = variable), position = "fill", stat = "identity") +
	   scale_y_continuous(breaks = scales::pretty_breaks(n=6)) +
       scale_fill_manual(values = cols, name = "Set Size") +
       labs(x = "", y = "%TSS Peaks", title = "Test set (chr2&3) and training set size") +
       theme(axis.text.x=element_text(angle = +60, size = 12, vjust = +1, hjust=+1))
test_train_plot_pdf <- "plots/test_train_set_peaks.pdf"
test_train_plot_png <- gsub(".pdf", ".png", test_train_plot_pdf)
ggsave(test_train_plot_pdf, width = 7, height = 6)
ggsave(test_train_plot_png, width = 7, height = 6)

########################
m_encode <- data.frame(RAMPAGE = c(14647, 14619), CAGE = c(33178, 27681), row.names = c("GM12878", "K562")) %>% mutate(sample = rownames(.)) %>% melt
colnames(m_encode) <- c("Sample", "Assay", "TSS_Peaks")
cols <- setNames(pal[5:6], c("GM12878", "K562"))
ggplot(m_encode) + geom_bar(aes(x = Assay, y = TSS_Peaks, fill = Sample), position="dodge", stat="identity") + 
	scale_y_continuous(breaks = scales::pretty_breaks(n=6)) +
	scale_fill_manual(values = cols, name = "Cell Line") +
	labs(x = "", y = "No. TSS Peaks", title = "RAMPAGE & CAGE TSS peaks") +
	theme(axis.text.x=element_text(angle = +45, size = 12, vjust = +1, hjust=+1))
m_encode_pdf <- "plots/encode_rampage_cage.pdf"
m_encode_png <- gsub(".pdf", ".png", m_encode_pdf)
ggsave(m_encode_pdf, width = 7, height = 6)
ggsave(m_encode_png, width = 7, height = 6)


