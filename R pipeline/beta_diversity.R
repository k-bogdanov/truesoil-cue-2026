# Prokaryotes
ps_nz_16s <- readRDS("~/NZ_CUE_16S_phyloseq.rds") %>%
  ps_mutate(Location = factor(Location, levels = c("Lincoln", "Kumara")))
ps_nz_16s

# Recode values
sample_data(ps_nz_16s)$Urine <- ifelse(
  sample_data(ps_nz_16s)$Urine == "Yes", "+",
  ifelse(sample_data(ps_nz_16s)$Urine == "No", "-", sample_data(ps_nz_16s)$Urine))

# Sample order
sample_data(ps_nz_16s)$Treatment = factor(
  sample_data(ps_nz_16s)$Treatment,
  levels = c("Control",
             "Roots",
             "Urine",
             "Both"))

sample_data(ps_nz_16s)$Roots = case_when(sample_data(ps_nz_16s)$Treatment %in% c("Roots", "Both") ~ "Yes", TRUE ~ "No")
# rarefy at 5k
ps_nz_16S_r <- rarefy_even_depth(ps_nz_16s, 5000, rngseed = TRUE)

# Lincoln only as an example, same for Kumara
bray_16s <- phyloseq::distance(NZUT_16S_merge_clean %>% ps_filter(Location == "Lincoln"), method = "bray")
meta_16s <- data.frame(sample_data(NZUT_16S_merge_clean %>% ps_filter(Location == "Lincoln")))
meta_16s$Urine <- factor(meta_16s$Urine)
meta_16s$Time_d <- factor(meta_16s$Time_d)

# PERMANOVA with the interaction term
permanova_16s <- adonis2(bray_16s ~ Urine * Time_d, data = meta_16s,
                         permutations = 999, by = "terms")
permanova_16s

# PERMDISP - dispersion by Urine
disp_urine <- betadisper(bray_16s, meta_16s$Urine)
anova(disp_urine)
plot(disp_urine)          
boxplot(disp_urine)

# PERMDISP - dispersion by Day
disp_Day <- betadisper(bray_16s, meta_16s$Time_d)
anova(disp_Day)
plot(disp_Day)          
boxplot(disp_Day)

# Boxplot of distances to centroid

df_box <- data.frame(
  group    = disp_Day$group,
  distance = disp_Day$distances
)

# order days numerically rather than alphabetically
day_levels <- sort(unique(as.numeric(as.character(df_box$group))))
df_box$group <- factor(df_box$group, levels = day_levels)

centroid_boxplot_Lincoln <- ggplot(df_box, aes(x = group, y = distance)) +
  geom_boxplot(outlier.shape = 1, width = 0.6) +
  #geom_jitter(width = 0.08, alpha = 0.6, size = 1.5) +
  labs(x = "Time (d)", y = "Bray distance to centroid") +
  revised_theme +            # swap for your own theme
  theme(legend.position = "none")
centroid_boxplot_Lincoln

# Ordination plot with hulls

eig <- disp_Day$eig
pct <- round(100 * eig / sum(eig[eig > 0]), 1)

df_pts <- data.frame(
  PCoA1 = disp_Day$vectors[, 1],
  PCoA2 = disp_Day$vectors[, 2],
  group = factor(disp_Day$group, levels = day_levels)
)

df_cent <- data.frame(
  group  = factor(rownames(disp_Day$centroids), levels = day_levels),
  cPCoA1 = disp_Day$centroids[, 1],
  cPCoA2 = disp_Day$centroids[, 2]
)

df_spider <- df_pts %>%
  left_join(df_cent, by = "group")

df_hull <- df_pts %>%
  group_by(group) %>%
  slice(chull(PCoA1, PCoA2)) %>%
  ungroup()

pcoa_xlim <- c(-0.3, 0.55)
pcoa_ylim <- c(-0.35, 0.45)

centroid_PCoA_Lincoln <- ggplot() +
  #geom_polygon(data = df_hull, aes(PCoA1, PCoA2),
  #alpha = 0.15, linewidth = 0.4) +
  geom_segment(data = df_spider,
               aes(x = cPCoA1, y = cPCoA2, xend = PCoA1, yend = PCoA2, alpha = group)) +
  scale_alpha_discrete(breaks = c("1", "7", "28", "72")) +
  geom_point(data = df_pts, aes(PCoA1, PCoA2, alpha = group), size = 2) +
  geom_point(data = df_cent, aes(cPCoA1, cPCoA2),
             shape = 4, size = 4, stroke = 1.2, color = "#FFC107") +
  geom_label_repel(data = df_cent, aes(cPCoA1, cPCoA2, label = group),
                   fill = "white", fontface = "bold", size = 4.5,
                   label.size = 0.5, segment.color = "grey50",
                   show.legend = FALSE) +
  labs(x = paste0("PCoA1 (", pct[1], "%)"),
       y = paste0("PCoA2 (", pct[2], "%)")) +
  coord_equal() +
  revised_theme +
  theme(legend.position = "none") +
  coord_equal(xlim = pcoa_xlim, ylim = pcoa_ylim)
centroid_PCoA_Lincoln