#### 13C from roots only ####
# background 13C flux (from unlabeled sources)
background_flux <- metadata %>%
  filter(Treatment %in% c("Control", "Urine"),  # No labeled roots
         !is.na(X13C_CO2_delta_VPDB)) %>%
  mutate(
    CO2_13C_flux_ug_g_h = (CO2_ppm * 0.85 * 44 * 12) / (22.4 * 44 * Soil_g * 1.5) *
      (((((X13C_CO2_delta_VPDB * 0.0111132)/1000) + 0.0111132) /
          (((X13C_CO2_delta_VPDB * 0.0112372)/1000) + 1.0111132)) * 100) / 100
  ) %>%
  group_by(Time_d, Location, Treatment) %>%
  dplyr::summarise(
    background_13C_flux = mean(CO2_13C_flux_ug_g_h, na.rm = TRUE),
    .groups = "drop"
  )

background_for_roots <- background_flux %>%
  filter(Treatment == "Control") %>%
  select(Time_d, Location, background_13C_flux) %>%
  dplyr::rename(control_background = background_13C_flux)

background_for_both <- background_flux %>%
  filter(Treatment == "Urine") %>%
  select(Time_d, Location, background_13C_flux) %>%
  dplyr::rename(urine_background = background_13C_flux)

# calculate label-derived 13C (only from roots)
label_derived <- metadata %>%
  filter(Treatment %in% c("Roots", "Both"),
         !is.na(X13C_CO2_delta_VPDB)) %>%
  mutate(
    CO2_13C_flux_ug_g_h = (CO2_ppm * 0.85 * 44 * 12) / (22.4 * 44 * Soil_g * 1.5) *
      (((((X13C_CO2_delta_VPDB * 0.0111132)/1000) + 0.0111132) /
          (((X13C_CO2_delta_VPDB * 0.0112372)/1000) + 1.0111132)) * 100) / 100
  ) %>%
  left_join(background_for_roots, by = c("Time_d", "Location")) %>%
  left_join(background_for_both, by = c("Time_d", "Location")) %>%
  mutate(
    # subtract appropriate background
    background_13C_flux = case_when(
      Treatment == "Roots" ~ control_background,
      Treatment == "Both" ~ urine_background,
      TRUE ~ NA_real_
    ),
    
    # label-derived 13C flux (only from roots)
    label_13C_flux = pmax(0, CO2_13C_flux_ug_g_h - background_13C_flux)
  ) %>%
  # now calculate cumulative label-derived flux
  group_by(Location, Treatment, Tech_rep) %>%
  arrange(Time_h, Gas_d) %>%
  mutate(interval_h = Time_h - lag(Time_h, default = first(Time_h))) %>%
  mutate(interval_h = case_when(soil_id %in% c("K7", "K8", "K9") ~ 164, TRUE ~ interval_h)) %>%
  mutate(avg_flux = (label_13C_flux + lag(label_13C_flux, default = first(label_13C_flux))) / 2,
         interval_cumulative = avg_flux * interval_h,
         cumulative_13C_ug_per_g = cumsum(interval_cumulative)
  ) %>%
  ungroup()


# summarise
label_derived_sum <- label_derived %>%
  mutate(Gas_d = factor(Gas_d, levels = c("1", "2", "7", "14", "19", "19 (2h)", "20", "21 (2h)", "21 (9h)", "22", "24", "28", "50", "72", "82"))) %>%
  group_by(Location, Treatment, Time_h, Gas_d) %>%
  dplyr::summarise_at(vars(label_13C_flux), funs(median(.), IQR(.)))

# map colours
label_derived_sum$color <- sapply(label_derived_sum$Treatment, color_mapping)
label_derived_sum$Treatment <- factor(label_derived_sum$Treatment, levels = c("Control", "Roots", "Urine", "Both"))

roots_initial <- label_derived_sum %>%
  filter(Location == "Lincoln", 
         Treatment == "Roots",
         Time_h < 506)

both_21_d <- label_derived_sum %>%
  filter(Location == "Lincoln", 
         Treatment == "Both",
         Time_h == 506)

dashed_line_df <- bind_rows(roots_initial, both_21_d)

# plot with ribbon
label_derived_flux_per_h <- ggplot(label_derived_sum %>% filter(Location == "Lincoln"), aes(y = median, x = Gas_d, group = Treatment, color = Treatment, shape = Treatment)) +
  geom_vline(xintercept = 4.8, colour = lines_color, linetype = "dashed") +
  geom_vline(xintercept = 1.2, colour = lines_color, linetype = "dotted") +
  geom_line(aes(color = color), linewidth = 0.7) +
  geom_point(aes(color = color), size = 3) +
  geom_ribbon(aes(ymin = 0, ymax = median, fill = color, color = color), alpha = 0.2, show.legend = FALSE) +
  geom_errorbar(aes(ymin = median - IQR, ymax = median + IQR, color = color), width =.1) +
  scale_color_identity(name = "Treatment", breaks = c(col1, col2, col3, col4), labels = c("C", "R", "U", "U + R"), guide = "legend") +
  scale_fill_identity(guide = "none") +
  scale_shape_manual(name = "Treatment", values = c(17, 16, 15, 18), labels = c("C", "R", "U", "U + R"), guide = "legend") +
  labs(x = "Time (d)", y = expression(bold(CO[2]*' '^{"13"}~C*' flux '*"(µg "*g^{-1}*' '*h^{-1}*")"))) +
  scale_x_discrete(labels = scales::label_wrap(2)) +
  truesoil_theme +
  theme(plot.margin = margin(l = 10, b = 10, r = 10, t = 10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  geom_ribbon(data = dashed_line_df, 
              aes(ymin = 0, ymax = median, xmax = Gas_d, group = 1), 
              color = col4,
              fill = col4,
              alpha = 0.2,
              show.legend = FALSE,
              linetype = "dashed") +
  ggtitle("Lincoln")

label_derived_flux_per_h

# Kumara
roots_initial <- label_derived_sum %>%
  filter(Location == "Kumara", 
         Treatment == "Roots",
         Time_h < 480)

both_19_d <- label_derived_sum %>%
  filter(Location == "Kumara", 
         Treatment == "Both",
         Time_h == 480)

dashed_line_df <- bind_rows(roots_initial, both_19_d)

Kumara_label_derived_flux_per_h <- ggplot(label_derived_sum %>% filter(Location == "Kumara"), aes(y = median, x = Gas_d, group = Treatment, color = Treatment, shape = Treatment)) +
  geom_vline(xintercept = 3.8, colour = lines_color, linetype = "dashed") +
  geom_vline(xintercept = 1.2, colour = lines_color, linetype = "dotted") +
  geom_line(aes(color = color), linewidth = 0.7) +
  geom_point(aes(color = color), size = 3) +
  geom_ribbon(aes(ymin = 0, ymax = median, fill = color, color = color), alpha = 0.2, show.legend = FALSE) +
  geom_errorbar(aes(ymin = median - IQR, ymax = median + IQR, color = color), width =.1) +
  scale_color_identity(name = "Treatment", breaks = c(col2, col4), labels = c("R", "U + R"), guide = "legend") +
  scale_fill_identity(guide = "none") +
  scale_shape_manual(name = "Treatment", values = c(17, 16), labels = c("R", "U + R"), guide = "legend") +
  labs(x = "Time (d)", y = expression(bold(CO[2]*' '^{"13"}~C*' flux '*"(µg "*g^{-1}*' '*h^{-1}*")"))) +
  scale_x_discrete(labels = scales::label_wrap(2)) +
  truesoil_theme +
  theme(plot.margin = margin(l = 10, b = 10, r = 10, t = 10),
        legend.position = c(0.8, 0.75),
        plot.title = element_text(hjust = 0.5)) +
  geom_ribbon(data = dashed_line_df, 
              aes(ymin = 0, ymax = median, xmax = Gas_d, group = 1), 
              color = col4,
              fill = col4,
              alpha = 0.2,
              show.legend = FALSE,
              linetype = "dashed") +
  ggtitle("Kumara Junction") +
  coord_cartesian(ylim = c(0, 0.25))

Kumara_label_derived_flux_per_h
#### Cumulative flux ####
# Get the Roots cumulative value at day 20 (or closest time point)
roots_day20_cumulative <- label_derived %>%
  filter(Location == "Lincoln", 
         Treatment == "Roots",
         Time_h == 506) %>%
  dplyr::summarise_at(vars(cumulative_13C_ug_per_g), funs(median(.))) %>%
  pull(cumulative_13C_ug_per_g)

# Create adjusted cumulative for Both treatment
label_derived_adjusted <- label_derived %>%
  mutate(
    cumulative_13C_ug_per_g_adj = case_when(
      Treatment == "Both" ~ cumulative_13C_ug_per_g + roots_day20_cumulative,
      TRUE ~ cumulative_13C_ug_per_g
    )
  )

# summarise
cumulative_label_derived_sum <- label_derived_adjusted %>%
  mutate(Gas_d = factor(Gas_d, levels = c("1", "2", "7", "14", "19", "19 (2h)", "20", "21 (2h)", "21 (9h)", "22", "24", "28", "50", "72", "82"))) %>%
  group_by(Location, Treatment, Time_h, Gas_d) %>%
  dplyr::summarise_at(vars(cumulative_13C_ug_per_g_adj), funs(median(.), IQR(.)))

# map colours
cumulative_label_derived_sum$color <- sapply(cumulative_label_derived_sum$Treatment, color_mapping)
cumulative_label_derived_sum$Treatment <- factor(cumulative_label_derived_sum$Treatment, levels = c("Control", "Roots", "Urine", "Both"))

# plots
cumulative_flux <- ggplot(cumulative_label_derived_sum %>% filter(Location == "Lincoln"), aes(y = median, x = Gas_d, group = Treatment, color = Treatment, shape = Treatment)) +
  geom_vline(xintercept = 4.8, colour = lines_color, linetype = "dashed") +
  geom_vline(xintercept = 1.2, colour = lines_color, linetype = "dotted") +
  geom_line(aes(color = color), linewidth = 0.7) +
  geom_point(aes(color = color), size = 3) +
  geom_errorbar(aes(ymin = median - IQR, ymax = median + IQR, color = color), width =.1) +
  scale_color_identity(name = "Treatment", breaks = c(col1, col2, col3, col4), labels = c("C", "R", "U", "U + R"), guide = "legend") +
  scale_shape_manual(name = "Treatment", values = c(17, 16, 15, 18), labels = c("C", "R", "U", "U + R"), guide = "legend") +
  labs(x = "Time (d)", y = expression(bold("Cumulative"*' '^{"13"}~C*' '*"(µg "*g^{-1}*")"))) +
  scale_x_discrete(labels = scales::label_wrap(2)) +
  truesoil_theme +
  theme(plot.margin = margin(l = 10, b = 10, r = 10, t = 10),
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 40))

cumulative_flux

# Kumara
# Get the Roots cumulative value at day 20 (or closest time point)
roots_day20_cumulative <- label_derived %>%
  filter(Location == "Kumara", 
         Treatment == "Roots",
         Time_h == 480) %>%
  dplyr::summarise_at(vars(cumulative_13C_ug_per_g), funs(median(.))) %>%
  pull(cumulative_13C_ug_per_g)

# Create adjusted cumulative for Both treatment
label_derived_adjusted <- label_derived %>%
  mutate(
    cumulative_13C_ug_per_g_adj = case_when(
      Treatment == "Both" ~ cumulative_13C_ug_per_g + roots_day20_cumulative,
      TRUE ~ cumulative_13C_ug_per_g
    )
  )

# summarise
cumulative_label_derived_sum <- label_derived_adjusted %>%
  mutate(Gas_d = factor(Gas_d, levels = c("1", "2", "7", "14", "19", "19 (2h)", "20", "21 (2h)", "21 (9h)", "22", "24", "28", "50", "72", "82"))) %>%
  group_by(Location, Treatment, Time_h, Gas_d) %>%
  dplyr::summarise_at(vars(cumulative_13C_ug_per_g_adj), funs(median(.), IQR(.)))


# map colours
cumulative_label_derived_sum$color <- sapply(cumulative_label_derived_sum$Treatment, color_mapping)
cumulative_label_derived_sum$Treatment <- factor(cumulative_label_derived_sum$Treatment, levels = c("Control", "Roots", "Urine", "Both"))



Kumara_cumulative_flux <- ggplot(cumulative_label_derived_sum %>% filter(Location == "Kumara"), aes(y = median, x = Gas_d, group = Treatment, color = Treatment, shape = Treatment)) +
  geom_vline(xintercept = 3.8, colour = lines_color, linetype = "dashed") +
  geom_vline(xintercept = 1.2, colour = lines_color, linetype = "dotted") +
  geom_line(aes(color = color), linewidth = 0.7) +
  geom_point(aes(color = color), size = 3) +
  geom_errorbar(aes(ymin = median - IQR, ymax = median + IQR, color = color), width =.1) +
  scale_color_identity(name = "Treatment", breaks = c(col1, col2, col3, col4), labels = c("C", "R", "U", "U + R"), guide = "legend") +
  scale_shape_manual(name = "Treatment", values = c(17, 16, 15, 18), labels = c("C", "R", "U", "U + R"), guide = "legend") +
  
  labs(x = "Time (d)", y = expression(bold("Cumulative"*' '^{"13"}~C*' '*"(µg "*g^{-1}*")"))) +
  scale_x_discrete(labels = scales::label_wrap(2)) +
  truesoil_theme +
  theme(plot.margin = margin(l = 10, b = 10, r = 10, t = 10),
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 40))

Kumara_cumulative_flux

#### Label-derived 13C HWEC and plot ####
# calculate background 13C flux (from unlabeled sources)
background_13C_HWEC <- metadata %>%
  select(Location, Treatment, Time_d, HWC_TC_mg_L, HW_13C_DeltaPDB) %>%
  filter(Treatment %in% c("Control", "Urine"),
         !is.na(HW_13C_DeltaPDB)) %>%
  mutate(
    HWEC_13C_ug_g = (HWC_TC_mg_L * 30 / 3) *
      (((((HW_13C_DeltaPDB * 0.0111132)/1000) + 0.0111132) /
          (((HW_13C_DeltaPDB * 0.0112372)/1000) + 1.0111132)) * 100) / 100
  ) %>%
  group_by(Time_d, Location, Treatment) %>%
  dplyr::summarise(
    background_HWEC_13C_ug_g = mean(HWEC_13C_ug_g, na.rm = TRUE),
    .groups = "drop"
  )

background_for_roots_HWEC <- background_13C_HWEC %>%
  filter(Treatment == "Control") %>%
  select(Time_d, Location, background_HWEC_13C_ug_g) %>%
  dplyr::rename(control_background = background_HWEC_13C_ug_g)

background_for_both_HWEC <- background_13C_HWEC %>%
  filter(Treatment == "Urine") %>%
  select(Time_d, Location, background_HWEC_13C_ug_g) %>%
  dplyr::rename(urine_background = background_HWEC_13C_ug_g)

# calculate label-derived 13C (only from roots)
label_derived_13C_HWEC <- metadata %>%
  select(soil_id, Location, Treatment, Time_d, HWC_TC_mg_L, HW_13C_DeltaPDB) %>%
  filter(Treatment %in% c("Roots", "Both"),
         !is.na(HW_13C_DeltaPDB)) %>%
  mutate(
    HWEC_13C_ug_g = (HWC_TC_mg_L * 30 / 3) *
      (((((HW_13C_DeltaPDB * 0.0111132)/1000) + 0.0111132) /
          (((HW_13C_DeltaPDB * 0.0112372)/1000) + 1.0111132)) * 100) / 100
  ) %>%
  left_join(background_for_roots_HWEC, by = c("Time_d", "Location")) %>%
  left_join(background_for_both_HWEC, by = c("Time_d", "Location")) %>%
  mutate(
    # Subtract appropriate background
    background_13C_HWEC = case_when(
      Treatment == "Roots" ~ control_background,
      Treatment == "Both" ~ urine_background,
      TRUE ~ NA_real_
    ),
    
    # Label-derived 13C (only from roots)
    excess_13C_HWEC = pmax(0, HWEC_13C_ug_g - background_13C_HWEC)
  )

# summarize
label_derived_13C_HWEC_sum <- label_derived_13C_HWEC %>%
  filter(!is.na(excess_13C_HWEC)) %>%
  mutate(Time_d = factor(Time_d, levels = c("1","7","14", "19", "20", "21", "22", "24", "28", "50", "72", "82"))) %>%
  group_by(Location, Time_d, Treatment) %>%
  dplyr::summarize_at(vars(excess_13C_HWEC), funs(median(.), IQR(.), mean(.), min(.), max(.)))

label_derived_13C_HWEC_sum$color <- sapply(label_derived_13C_HWEC_sum$Treatment, color_mapping)
label_derived_13C_HWEC_sum$Treatment <- factor(label_derived_13C_HWEC_sum$Treatment, levels = c("Control", "Roots", "Urine", "Both"))

#label_derived_13C_HWEC_sum$Treatment_Bag <- paste(label_derived_13C_HWEC_sum$Treatment, "_", label_derived_13C_HWEC_sum$Bag)

# plot
HWEC_13C <- ggplot(label_derived_13C_HWEC_sum %>% filter(Location == "Lincoln"), aes(y = mean, x = Time_d, group = Treatment, color = Treatment, shape = Treatment)) +
  geom_vline(xintercept = 2.8, colour = lines_color, linetype = "dashed") +
  geom_vline(xintercept = 1.2, colour = lines_color, linetype = "dotted") +
  geom_line(aes(color = color), linewidth = 0.7) +
  geom_point(aes(color = color), size = 3) +
  scale_x_discrete(labels = scales::label_wrap(2)) +
  geom_errorbar(aes(ymin = min, ymax = max, color = color), width =.1) +
  scale_color_identity(name = "Treatment", breaks = c(col1, col2, col3, col4), labels = c("C", "R", "U", "U + R"), guide = "legend") +
  scale_shape_manual(name = "Treatment", values = c(17, 16, 15, 18), labels = c("C", "R", "U", "U + R"), guide = "legend") +
  labs(x = "Time (d)", y = expression(bold("HWEC"*' '^{"13"}~C*' '*"(µg "*g^{-1}*")"))) +
  truesoil_theme +
  theme(plot.margin = margin(l = 10, b = 10, r = 10, t = 10),
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 20))

HWEC_13C


# Kumara
Kumara_HWEC_13C <- ggplot(label_derived_13C_HWEC_sum %>% filter(Location == "Kumara"), aes(y = mean, x = Time_d, group = Treatment, color = Treatment, shape = Treatment)) +
  geom_vline(xintercept = 1.8, colour = lines_color, linetype = "dashed") +
  geom_vline(xintercept = 0.6, colour = lines_color, linetype = "dotted") +
  geom_line(aes(color = color), linewidth = 0.7) +
  geom_point(aes(color = color), size = 3) +
  scale_x_discrete(labels = scales::label_wrap(2)) +
  geom_errorbar(aes(ymin = min, ymax = max, color = color), width =.1) +
  scale_color_identity(name = "Treatment", breaks = c(col1, col2, col3, col4), labels = c("C", "R", "U", "U + R"), guide = "legend") +
  scale_shape_manual(name = "Treatment", values = c(17, 16, 15, 18), labels = c("C", "R", "U", "U + R"), guide = "legend") +
  labs(x = "Time (d)", y = expression(bold("HWEC"*' '^{"13"}~C*' '*"(µg "*g^{-1}*")"))) +
  truesoil_theme +
  theme(plot.margin = margin(l = 10, b = 10, r = 10, t = 10),
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 20))

Kumara_HWEC_13C

#### CUE ####
CUE_df <- label_derived_13C_HWEC %>%
  select(soil_id, Location, Treatment, Time_d, excess_13C_HWEC, HWC_TC_mg_L) %>%
  filter(Time_d != "1") %>%
  left_join(label_derived_adjusted %>%
              select(soil_id, Location, Treatment, Time_d, cumulative_13C_ug_per_g_adj, HWC_TC_mg_L), by = c("soil_id", "HWC_TC_mg_L", "Treatment", "Location", "Time_d")) %>%
  select(-HWC_TC_mg_L) %>%
  mutate(CUE = excess_13C_HWEC / (excess_13C_HWEC + cumulative_13C_ug_per_g_adj)) %>%
  group_by(Location, Treatment, Time_d) %>%
  dplyr::summarise(
    across(
      CUE,
      list(median = ~median(.x, na.rm = TRUE), 
           IQR = ~IQR(.x, na.rm = TRUE),
           mean = ~mean(.x, na.rm = TRUE), 
           min = ~min(.x, na.rm = TRUE),
           max = ~max(.x, na.rm = TRUE))
    ),
    .groups = "drop"
  ) %>%
  mutate(Time_d = factor(Time_d, levels = c("1", "7", "14", "19", "20", "21", "22", "24", "28", "50", "72", "82")))

CUE_df$color <- sapply(CUE_df$Treatment, color_mapping)
CUE_df$Treatment <- factor(CUE_df$Treatment, levels = c("Control", "Roots", "Urine", "Both"))

#CUE_df$Treatment_Bag <- paste(CUE_df$Treatment, "_", CUE_df$Bag)

write.csv(label_derived_13C_HWEC, "C:/Users/bogda/Documents/PhD/TRUESOIL/urine_experiment/label_derived_13C_HWEC.csv")

CUE_plot_Lincoln <- ggplot(CUE_df %>%filter(Location == "Lincoln"), aes(y = CUE_mean, x = Time_d, group = Treatment, color = Treatment, shape = Treatment)) +
  geom_vline(xintercept = 1.8, colour = lines_color, linetype = "dashed") +
  geom_vline(xintercept = 0.8, colour = lines_color, linetype = "dotted") +
  geom_line(aes(color = color), linewidth = 0.7) +
  geom_point(aes(color = color), size = 3, position = position_dodge(width = 0.1)) +
  geom_errorbar(aes(ymin = CUE_min, ymax = CUE_max, color = color), width =.1, position = position_dodge(width = 0.1)) +
  scale_color_identity(name = "Treatment", breaks = c(col2, col4), labels = c("R", "U + R"), guide = "legend") +
  scale_shape_manual(name = "Treatment", values = c(17, 16, 15, 18), labels = c("R", "U + R"), guide = "legend") +
  labs(x = "Time (d)", y = "Measured CUE") +
  scale_x_discrete(labels = scales::label_wrap(2)) +
  revised_theme +
  theme(plot.margin = margin(l = 10, b = 10, r = 10, t = 10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(0, 1))
CUE_plot_Lincoln

legend <- get_legend(CUE_plot_Lincoln)

CUE_plot_Lincoln <- CUE_plot_Lincoln + theme(legend.position = "none")


# Kumara
CUE_plot_Kumara <- ggplot(CUE_df %>% filter(Location == "Kumara"), aes(y = CUE_mean, x = Time_d, group = Treatment, color = Treatment, shape = Treatment)) +
  geom_vline(xintercept = 1.8, colour = lines_color, linetype = "dashed") +
  geom_vline(xintercept = 0.6, colour = lines_color, linetype = "dotted") +
  geom_line(aes(color = color), linewidth = 0.7) +
  geom_point(aes(color = color), size = 3, position = position_dodge(width = 0.1)) +
  geom_errorbar(aes(ymin = CUE_min, ymax = CUE_max, color = color), width =.1, position = position_dodge(width = 0.1)) +
  scale_color_identity(name = "Treatment", breaks = c(col2, col4), labels = c("R", "U + R"), guide = "legend") +
  scale_shape_manual(name = "Treatment", values = c(17, 16, 15, 18), labels = c("R", "U + R"), guide = "legend") +
  labs(x = "Time (d)", y = "Measured CUE") +
  scale_x_discrete(labels = scales::label_wrap(2)) +
  revised_theme +
  theme(plot.margin = margin(l = 10, b = 10, r = 10, t = 10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  #axis.text.y = element_blank())+
  coord_cartesian(ylim = c(0, 1))
CUE_plot_Kumara