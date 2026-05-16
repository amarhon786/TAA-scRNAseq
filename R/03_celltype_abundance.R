# ct fractions by disease x sex, cells are already typed so no deconvolution

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})
source("R/theme_nature.R")

obj <- readRDS("data/processed/seurat_annotated.rds")

abund <- obj@meta.data %>%
  count(sample, sex, disease, celltype, name = "n") %>%
  group_by(sample) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

write.csv(abund, "results/tables/celltype_abundance.csv", row.names = FALSE)

p_stack <- ggplot(abund, aes(sample, frac, fill = celltype)) +
  geom_col(width = 0.85) +
  facet_grid(. ~ disease + sex, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent_format(1), expand = c(0, 0)) +
  scale_fill_manual(values = OKABE_ITO) +
  labs(x = NULL, y = "Cell-type fraction", title = "Composition per sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        panel.spacing.x = unit(1, "mm"))

p_box <- ggplot(abund, aes(disease, frac, fill = sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, linewidth = 0.3) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15), size = 0.8,
             color = "black", stroke = 0) +
  scale_fill_manual(values = SEX_COLORS) +
  facet_wrap(~ celltype, scales = "free_y") +
  labs(y = "Fraction of sample", x = NULL,
       title = "Cell-type abundance by disease × sex")

save_nature(p_stack, "results/figures/composition_stack.pdf",
            width_mm = NATURE_DOUBLE_COL, height_mm = 95)
save_nature(p_box, "results/figures/composition_box.pdf",
            width_mm = NATURE_DOUBLE_COL, height_mm = 110)
pdf("results/figures/composition.pdf", width = 11, height = 8)
print(p_stack); print(p_box)
dev.off()

# plain lm per ct as a quick read, beta regression would suit fractions better
fit_tab <- abund %>%
  group_by(celltype) %>%
  group_modify(~ {
    fit <- lm(frac ~ disease * sex, data = .x)
    broom::tidy(fit)
  }) %>%
  filter(term != "(Intercept)") %>%
  ungroup()

write.csv(fit_tab, "results/tables/celltype_abundance_tests.csv", row.names = FALSE)

message("Composition tables + figures written.")
message("Next: source('R/04_pseudobulk_de.R')")
