files <- list.files(path = "~/FBA_results/tables", pattern = "\\.xls$", full.names = TRUE)

# Read and combine all files
FBA_results <- purrr::map_df(files, function(file) {
  df <- read_excel(file, sheet = 1)
  sample <- gsub("\\.xls$", "", basename(file))
  df$bin_id <- sample
  return(df)
}) %>%
  dplyr::relocate(bin_id, .before = id) %>%
  separate(
    bin_id,
    into = c("seqcoast_id", "bin_id"),
    sep = "_bin\\.",
    remove = FALSE
  ) %>%
  mutate(bin_id = paste0("bin.", bin_id))

FBA_results

# calculate uptake
C_uptake_flux <- FBA_results %>%
  filter(uptake > 0, str_detect(formula, "^C\\d+") | str_detect(formula, "^CH")) %>%
  mutate(C_atoms = case_when(
    str_detect(formula, "^CH") ~ 1,
    TRUE ~ as.numeric(str_extract(formula, "(?<=^C)\\d+"))
  )) %>%
  mutate(
    C_uptake_rate = C_atoms * uptake
  ) %>%
  dplyr::relocate(C_uptake_rate, C_atoms, .before = charge) %>%
  group_by(seqcoast_id, bin_id) %>%
  dplyr::summarise(sum(C_uptake_rate))

# calculate secretion
C_secretion_flux <- FBA_results %>%
  filter(uptake < 0, formula == "CO2" | str_detect(formula, "^C\\d+") | str_detect(formula, "^CH"), name != "Biomass [c0]") %>%
  mutate(C_atoms = case_when(
    str_detect(formula, "^CH") ~ 1,
    str_detect(formula, "CO2") ~ 1,
    TRUE ~ as.numeric(str_extract(formula, "(?<=^C)\\d+"))
  )) %>%
  mutate(
    C_secretion_rate = C_atoms * uptake * -1
  ) %>%
  dplyr::relocate(C_secretion_rate, C_atoms, .before = charge) %>%
  group_by(seqcoast_id, bin_id) %>%
  dplyr::summarise(sum(C_secretion_rate))

# combine them and calculate CUE
CUE_results <- C_uptake_flux %>%
  full_join(C_secretion_flux, by = c("seqcoast_id", "bin_id")) %>%
  dplyr::rename("uptake" = "sum(C_uptake_rate)", "secretion" = "sum(C_secretion_rate)") %>%
  mutate(cue = (uptake - secretion)/uptake)

# this is already combined with measured CUE
combined_results <- readRDS("~/combined_results.rds")

combined_results$Treatment_combined <- case_when(combined_results$Treatment %in% c("Control", "Roots") ~ "No urine", TRUE ~ "Urine applied")

combined_results_sum <- combined_results %>%
  filter(!is.na(cue)) %>%
  group_by(Location, Treatment_combined, Time_d) %>%
  dplyr::summarise(
    across(
      cue,
      list(median = ~median(.x, na.rm = TRUE), 
           IQR = ~IQR(.x, na.rm = TRUE),
           min = ~min(.x, na.rm = TRUE), 
           max = ~max(.x, na.rm = TRUE),
           mean = ~mean(.x, na.rm = TRUE))
    ),
    .groups = "drop"
  ) %>%
  mutate(Time_d = factor(Time_d, levels = c("1", "2", "7", "14", "19", "20", "21", "22", "24", "28", "50", "72", "82")))


# combined
combined_results_sum$color <- sapply(combined_results_sum$Treatment_combined, color_mapping_combo)
combined_results_sum$Treatment_combined <- factor(combined_results_sum$Treatment_combined, levels = c("No urine", "Urine applied"))

CUE_MAGs_plot_L <- ggplot(combined_results_sum %>% filter(Location == "Lincoln"), aes(y = cue_mean, x = Time_d, group = Treatment_combined, color = Treatment_combined, shape = Treatment_combined)) +
  geom_vline(xintercept = 3.8, colour = lines_color, linetype = "dashed") +
  geom_vline(xintercept = 1.2, colour = lines_color, linetype = "dotted") +
  geom_line(aes(color = color), linewidth = 0.7) +
  geom_point(aes(color = color), size = 3, position = position_dodge(width = 0.15)) +
  geom_errorbar(aes(ymin = cue_min, ymax = cue_max, color = color), width =.1, position = position_dodge(width = 0.15)) +
  scale_color_identity(name = "Treatment", breaks = c(col2, col4), labels = c("No urine", "Urine applied"), guide = "legend") +
  scale_shape_manual(name = "Treatment", values = c(17, 16), labels = c("No urine", "Urine applied"), guide = "legend") +
  labs(x = "Time (d)", y = "Predicted CUE") +
  scale_x_discrete(labels = scales::label_wrap(2)) +
  revised_theme +
  theme(plot.margin = margin(l = 10, b = 10, r = 10, t = 10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("Low-MAP soil (Lincoln)")

CUE_MAGs_plot_L


CUE_MAGs_plot_K <- ggplot(combined_results_sum %>% filter(Location == "Kumara"), aes(y = cue_mean, x = Time_d, group = Treatment_combined, color = Treatment_combined, shape = Treatment_combined)) +
  geom_vline(xintercept = 2.8, colour = lines_color, linetype = "dashed") +
  geom_vline(xintercept = 1.2, colour = lines_color, linetype = "dotted") +
  geom_line(aes(color = color), linewidth = 0.7, show.legend = FALSE) +
  geom_point(aes(color = color), size = 3, position = position_dodge(width = 0.15)) +
  geom_errorbar(aes(ymin = cue_min, ymax = cue_max, color = color), width =.1, position = position_dodge(width = 0.15), show.legend = FALSE) +
  scale_color_identity(name = "Treatment", breaks = c(col2, col4), labels = c("No urine", "Urine applied"), guide = "legend") +
  scale_shape_manual(name = "Treatment", values = c(17, 16), labels = c("No urine", "Urine applied"), guide = "legend") +
  labs(x = "Time (d)", y = "Predicted CUE") +
  scale_x_discrete(labels = scales::label_wrap(2)) +
  revised_theme +
  theme(plot.margin = margin(l = 10, b = 10, r = 10, t = 10),
        legend.position = c(0.75, 0.75),
        plot.title = element_text(hjust = 0.5),
        legend.title.align = 0.5,
        legend.box = "horizontal",        # stack legends in rows instead of side-by-side
        legend.box.just = "center",
        legend.justification = "center") +
  #axis.text.y = element_blank()) +
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("High-MAP soil (Kumara)")

CUE_MAGs_plot_K
