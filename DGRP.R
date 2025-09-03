library(readxl)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(rstatix)
library(dplyr)
library(binom)


#Data from DGRP-517+ crossed with 9 DGRP lines and w1118
DGRP<-read_excel("./1DGRP_cross.xlsx")
#Data for maternal or paternally infected embryos
eggs<-read_excel("./DGRP405_DGRP530_embryos.xlsx")
#Data for maternal and paternally infected embryos at 2 days and 8 days
egg_28152<-read_excel("./DGRP189_embryos_over_time.xlsx")


##Shared Y range to standaardize the Y-axis on subsequent plots
shared_y_range <- c(1e-9, 1e5)

#Make pseudocounts for 0 infections so they are farther away from the postiive but lowly infected samples
DGRP <- DGRP %>%
  mutate(
    dCt_modified = ifelse(dCt == 0, 1e-10, dCt),
    is_pseudocount = ifelse(dCt == 0, "pseudocount", "real")
  )

# Calculate p-values and convert to significance stars; Wilcoxon test excluding 0 values

p_values <- DGRP %>%
  group_by(Generation, Cross) %>%
   filter((Infected)=="Yes")%>%
  summarise(
    p_value = tryCatch(
      wilcox.test(dCt ~ Parent_Infected)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    signif_label = case_when(
      is.na(p_value) ~ NA_character_,
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )


####Medians for each Generation, Cross, Parent Infected Type
medians <- DGRP %>%
  filter(Infected == "Yes") %>%
  group_by(Generation, Cross, Parent_Infected) %>%
  summarise(
    median_value = median(dCt_modified, na.rm = TRUE),
    .groups = "drop"
  )


###### number of infected flies per generation and parent infected
summary_table <- DGRP %>%
  group_by(Parent_Infected, Infected, Generation) %>%
  summarise(count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = Infected,
    values_from = count,
    values_fill = 0
  )

summary_table




# Figure 1A
plot <- ggplot(data = DGRP, aes(x = Parent_Infected, y = dCt)) +
  geom_jitter(
    height = 0, width = 0.2,
    aes(color = as.factor(Replicate), shape = Offspring_Sex), size=1.5
  , show.legend=FALSE) +
  scale_alpha_manual(values = c("real" = 1, "pseudocount" = 0.3), guide = FALSE) +
  facet_grid(Generation ~ Cross, drop = TRUE) +
  theme(strip.text.x.top = element_text(size = 6)) +
  scale_y_log10(limits= shared_y_range)+
  ylab("Galbut virus RNA relative to RpL32 mRNA") +
  xlab("Original Transmission Type") +
  labs(color = "Replicate Number", shape = "Offspring Sex") +
  
  # Add asterisks instead of p-values
  geom_text(
    data = filter(p_values, !is.na(signif_label)),
    aes(x = 1.5, y = 1e-5, label = signif_label),
    inherit.aes = FALSE,
    size = 4,
    color = "black"
  ) +
  coord_cartesian(clip = "off")+theme_bw()+
  scale_color_brewer(palette= "Dark2")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# Show the plot
plot


###Adds the X to DGRP-405
missing_points <- tibble(
  Cross = rep("DGRP-405", 3),                    
  Generation = c("F1", "F2", "F3"),              
  Parent_Infected = "Paternal", levels = levels(DGRP$Parent_Infected),
  dCt = 1e-3                                     
)

# Add black X's to the plot
plot <- plot +
  geom_point(
    data = missing_points,
    aes(x = Parent_Infected, y = dCt),
    inherit.aes = FALSE,
    color = "black",
    shape = 4,
    size = 3
  )

plot
  

#Figure 1C
stats_asterisks_gen <- DGRP %>%
  filter(Infected %in% c("Yes", "No"),
         Generation %in% c("F1", "F2", "F3")) %>%
  group_by(Cross, Parent_Infected, Generation) %>%
  summarise(
    infected_count   = sum(Infected == "Yes", na.rm = TRUE),
    uninfected_count = sum(Infected == "No",  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Cross, Parent_Infected) %>%
  summarise(
    comparisons = list({
      gen_pairs <- combn(unique(Generation), 2, simplify = FALSE)
      
      map_dfr(gen_pairs, function(pair) {
        subdata <- filter(cur_data(), Generation %in% pair)
        
        mat <- matrix(
          c(
            subdata$infected_count[subdata$Generation == pair[1]],
            subdata$uninfected_count[subdata$Generation == pair[1]],
            subdata$infected_count[subdata$Generation == pair[2]],
            subdata$uninfected_count[subdata$Generation == pair[2]]
          ),
          nrow = 2,
          byrow = TRUE
        )
        
        tibble(
          comparisons = list(pair),  
          p_value = fisher.test(mat)$p.value
        )
      })
    }),
    .groups = "drop"
  ) %>%
  unnest(comparisons) %>%
  mutate(
    p.adj = p.adjust(p_value, method = "bonferroni"),
    signif = case_when(
      p.adj <= 0.0001 ~ "****",
      p.adj <= 0.001  ~ "***",
      p.adj <= 0.01   ~ "**",
      p.adj <= 0.05   ~ "*",
      TRUE            ~ ""
    )
  )


# Plot
  gen_plot <- ggplot(data = DGRP, aes(x = Generation, y = dCt, group = Cross, Parent_Infected)) +
    geom_jitter(
      height = 0, width = 0.2,
      aes(color = as.factor(Replicate), shape = Offspring_Sex), size=1.5
    ) +
    scale_alpha_manual(values = c("real" = 1, "pseudocount" = 0.3), guide = FALSE) +
    facet_grid(Parent_Infected ~ Cross, drop = TRUE) +
    theme(strip.text.x.top = element_text(size = 6)) +
    scale_y_log10(limits= shared_y_range)+
    ylab("Galbut virus RNA relative to RpL32 mRNA") +
    xlab("Generation") +
    labs(color = "Replicate Number", shape = "Offspring Sex") +
    coord_cartesian(clip = "off")+theme_bw()+
    scale_color_brewer(palette= "Dark2")+
    geom_point(
      data = missing_points,
      aes(x = "F2", y = dCt),
      inherit.aes = FALSE,
      color = "black",
      shape = 4,
      size = 3
    )

  

# Show the plot
gen_plot



####Figure 1B
mvf_p_values <- DGRP %>%
  group_by(Generation) %>%
  filter((Infected)=="Yes")%>%
  summarise(
    p_value = tryCatch(
      wilcox.test(dCt_modified ~ Parent_Infected)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    signif_label = case_when(
      is.na(p_value) ~ NA_character_,
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )



mvf<- ggplot(DGRP, aes(x = Parent_Infected, y = dCt, group = Parent_Infected, shape=Offspring_Sex, color=as.factor(Replicate))) +
  geom_jitter(width = 0.2, show.legend = FALSE) +
  ggtitle("") +
  facet_grid(Generation~All) +
  ylab("Galbut virus Ct relative to RpL32") +
  xlab("Original Transmission Type") +
  scale_y_log10(limits= shared_y_range)+
  coord_cartesian(clip = "off")+theme_bw()+
  geom_text(
    data = filter(mvf_p_values, !is.na(signif_label)),
    aes(x = 1.5, y = 1e-5, label = signif_label),
    inherit.aes = FALSE,
    size = 4,
    color = "black"
  )+ scale_color_brewer(palette= "Dark2")+ theme(axis.title.y = element_blank(),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

mvf


#Figure 1D

gen_p_values <- DGRP %>%
  group_by(Parent_Infected) %>%
  filter((Infected)=="Yes")%>%
  summarise(
    p_value = tryCatch(
      wilcox.test(dCt_modified ~ Parent_Infected)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    signif_label = case_when(
      is.na(p_value) ~ NA_character_,
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )



gen_all<- ggplot(DGRP, aes(x = Generation, y = dCt, shape=Offspring_Sex, color=as.factor(Replicate))) +
  geom_jitter(width = 0.2, show.legend = FALSE) +
  ggtitle("") +
  facet_grid(Parent_Infected~All) +
  ylab("Galbut virus Ct relative to RpL32") +
  xlab("Generation") +
  scale_y_log10(limits= shared_y_range)+
  coord_cartesian(clip = "off")+theme_bw()+
 scale_color_brewer(palette= "Dark2")+ theme(axis.title.y = element_blank(),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

gen_all




# Make the top row first with widths
top_row <- plot + mvf + 
  plot_layout(widths = c(3.5, .5), guides = "collect", axes = "collect")
bottom_row <- gen_plot + gen_all + 
  plot_layout(widths = c(3.5, .5), guides = "collect", axes = "collect")

# Stack top row over gen_plot
final <- top_row / bottom_row +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

#Final Figure 1
final



##Supplemental Figure 1
ggqqplot(DGRP$dCt)
shapiro<-shapiro.test(DGRP$dCt)
shapiro_p <- shapiro$p.value
pval_text <- paste0("Shapiro-Wilk p = ", signif(shapiro_p, 3))

normal<- ggplot(data= DGRP, aes(x=dCt))+
  geom_histogram()+scale_x_log10()+  annotate("text", 
                                              x = Inf, y = Inf, 
                                              label = pval_text, 
                                              hjust = 1.1, vjust = 1.5, 
                                              size = 5, 
                                              fontface = "italic") +xlab("Galbut virus RNA compared to RpL32 mRNA")+
  ylab("Count")

normal


#Figure 2
infection_data <- DGRP %>%
  filter(Infected %in% c("Yes", "No")) %>%
  group_by(Generation, Cross, Parent_Infected) %>%
  summarise(
    infected_count = sum(Infected == "Yes", na.rm = TRUE),
    total_count = n(),
    .groups = "drop"
  ) %>%
  filter(total_count > 0) %>%
  mutate(
    infection_freq = infected_count / total_count,
    percent_label = paste0(round(infection_freq * 100, 1), "%")
  ) %>%
  bind_cols(
    binom.confint(x = .$infected_count, n = .$total_count, methods = "wilson")[, c("lower", "upper")]
  )


# Prepare asterisks table
stats_asterisks <- DGRP %>%
  filter(Infected %in% c("Yes", "No")) %>%
  group_by(Cross, Generation) %>%
  filter(all(c("Maternal", "Paternal") %in% Parent_Infected)) %>%
  group_by(Cross, Generation, Parent_Infected) %>%
  summarise(
    infected_count = sum(Infected == "Yes", na.rm = TRUE),
    uninfected_count = sum(Infected == "No", na.rm = TRUE),
    .groups = "drop_last"
  ) %>%
  summarise(
    p_value = {
      mat <- matrix(
        c(
          infected_count[Parent_Infected == "Maternal"],
          uninfected_count[Parent_Infected == "Maternal"],
          infected_count[Parent_Infected == "Paternal"],
          uninfected_count[Parent_Infected == "Paternal"]
        ),
        nrow = 2,
        byrow = TRUE
      )
      fisher.test(mat)$p.value
    },
    .groups = "drop"
  ) %>%
  mutate(
    p.adj = p.adjust(p_value, method = "bonferroni"),
    signif = case_when(
      p.adj <= 0.0001 ~ "****",
      p.adj <= 0.001  ~ "***",
      p.adj <= 0.01   ~ "**",
      p.adj <= 0.05   ~ "*",
      TRUE            ~ ""
    )
  )

# Add y.position = max of both Maternal/Paternal points in infection_data
stats_asterisks <- stats_asterisks %>%
  left_join(
    infection_data %>%
      group_by(Cross, Generation) %>%
      summarise(y.position = max(infection_freq) + 0.05, .groups = "drop"),
    by = c("Cross", "Generation")
  )

# Plot
ggplot(infection_data, aes(x = Generation, y = infection_freq, color = Parent_Infected, group = Parent_Infected)) +
  geom_point(size = 2) +
  geom_line() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  facet_wrap(~ Cross, ncol = 2) +
  labs(
    x = "Generation",
    y = "Frequency of Infection",
    color = "Original Transmission Type"
  ) +
  scale_color_manual(values = c("Maternal" = "hotpink", "Paternal" = "skyblue3")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1.2)) +
  geom_text(
    data = stats_asterisks,
    aes(x = Generation, y = y.position, label = signif),
    inherit.aes = FALSE,
    size = 5
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 8),
    legend.position = "bottom"
  )




#########Frequency plot separated by replicate
infection_data_reps <- DGRP %>%
  filter(Infected %in% c("Yes", "No")) %>%
  group_by(Generation, Cross, Parent_Infected, Replicate) %>% 
  summarise(
    infected_count = sum(Infected == "Yes", na.rm = TRUE),
    total_count = n(),
    .groups = "drop"
  ) %>%
  filter(total_count > 0) %>%
  mutate(
    infection_freq = infected_count / total_count
  ) %>%
  bind_cols(
    binom.confint(x = .$infected_count, n = .$total_count, methods = "wilson")[, c("lower", "upper")]
  )


ggplot(infection_data_reps, aes(
  x = Generation, 
  y = infection_freq, 
  color = as.factor(Replicate))) +
  geom_point(size = 2) +
  geom_line() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  facet_grid(Parent_Infected ~ Cross) +
  labs(
    x = "Generation",
    y = "Percent Infected",
    color = "Replicate"
  ) +
  scale_color_brewer(palette = "Dark2") +  # Dark2 palette
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 8),
    legend.position = "bottom"
  )

#Supplemental Figure 2
maternal_p <-
  ggplot(filter(infection_data_reps, Parent_Infected == "Maternal"),
         aes(
           x = Generation,
           y = infection_freq * 100,
           color = as.factor(Replicate),
           group = Replicate)) +
  geom_point(size = 1.5) +
  geom_line(linewidth = 0.5, linetype = "dotted") +
  facet_grid(Replicate ~ Cross) +
  labs(
    x = "",
    y = "Percent Infected",
    color = "Replicate"
  ) +
  scale_color_brewer(palette = "Dark2") +  
  scale_y_continuous(breaks = c(0,50,100)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 8),
    legend.position = "none"
  )


paternal_p <-
  ggplot(filter(infection_data_reps, Parent_Infected == "Paternal"),
         aes(
           x = Generation,
           y = infection_freq * 100,
           color = as.factor(Replicate),
           group = Replicate)) +
  geom_point(size = 1.5) +
  geom_line(linewidth = 0.5, linetype = "dotted") +
  facet_grid(Replicate ~ Cross) +
  labs(
    x = "Generation",
    y = "Percent Infected",
    color = "Replicate"
  ) +
  scale_color_brewer(palette = "Dark2") +  
  scale_y_continuous(breaks = c(0,50,100)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 8),
    legend.position = "bottom"
  )


maternal_p / paternal_p+plot_annotation(tag_levels = "A")

ggsave("prevalence_by_replicate.pdf", width = 7.5, height=7, units="in")


#Median for Offspring Sex

medians_offspring <- DGRP %>%
  filter(Infected == "Yes") %>%
  group_by(Generation, Offspring_Sex) %>%
  summarise(
    median_value = median(dCt_modified, na.rm = TRUE),
    .groups = "drop"
  )

#Figure 4

sex_plot<-ggplot(DGRP, aes(x=Offspring_Sex, y=dCt_modified, color=Offspring_Sex, group=Offspring_Sex))+
  geom_jitter(width=.2, aes(alpha = is_pseudocount), show.legend = FALSE)+
  scale_alpha_manual(values = c("real" = 1, "pseudocount" = 0.15), guide="none")+
  ggtitle("")+ scale_y_log10()+
  ylab("Galbut virus RNA \nrelative to RpL32 mRNA")+
  xlab("Offspring Sex")+coord_cartesian(clip="off")+
  labs(color="", alpha="")+theme_bw()+ 
  facet_wrap(~Generation)+scale_color_manual(values= c("Female"= "mediumorchid", "Male"="palegreen2"))+
  stat_compare_means(paired=FALSE, method="wilcox.test", label="p.signif", label.x.npc = "center", show.legend = FALSE)


sex_plot


infection_data_sex <- DGRP %>%
  filter(Infected %in% c("Yes", "No")) %>%
  group_by( Offspring_Sex, Generation) %>%
  summarise(
    infected_count = sum(Infected == "Yes", na.rm = TRUE),
    total_count = n(),
    .groups = "drop"
  ) %>%
  filter(total_count > 0) %>%
  mutate(
    infection_freq = infected_count / total_count,
    percent_label = paste0(round(infection_freq * 100, 1), "%")
  ) %>%
  bind_cols(
    binom.confint(x = .$infected_count, n = .$total_count, methods = "wilson")[, c("lower", "upper")]
  )



# Calculate Fisher's exact test for each Generation
fisher_results_sex <- infection_data_sex %>%
  group_by(Generation) %>%
  summarise(
    fisher_p = {
      sexes <- sort(unique(Offspring_Sex))
      if (length(sexes) == 2) {
        tab <- matrix(
          c(
            infected_count[Offspring_Sex == sexes[1]],
            total_count[Offspring_Sex == sexes[1]] - infected_count[Offspring_Sex == sexes[1]],
            infected_count[Offspring_Sex == sexes[2]],
            total_count[Offspring_Sex == sexes[2]] - infected_count[Offspring_Sex == sexes[2]]
          ),
          nrow = 2,
          byrow = TRUE,
          dimnames = list(
            Sex = sexes,
            Infection = c("Infected", "Not Infected")
          )
        )
        fisher.test(tab)$p.value
      } else {
        NA_real_
      }
    },
    .groups = "drop"
  )

# Make label data frame for plotting (one label per Generation)
pval_labels_sex <- infection_data_sex %>%
  distinct(Generation) %>%
  left_join(fisher_results_sex, by = "Generation") %>%
  mutate(
    label = paste0("p-value = ", signif(fisher_p, 2)),
    y_pos = 0  
  )




sex_freq<-ggplot(infection_data_sex, aes(x = Offspring_Sex, y = infection_freq, color=Offspring_Sex)) +
  geom_point(show.legend = FALSE) +
  ylab("Prevalence of \nGalbut virus Infection") +
  xlab("Offspring Sex") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 9),
    panel.border = element_rect(color = "black", size = 0.9)
    
  )+
  facet_wrap(~Generation, nrow=1)+
  scale_color_manual(values=c("mediumorchid", "palegreen2"))+
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.2, show.legend=FALSE)+   geom_text(
                  data = pval_labels_sex,
                  aes(x = "Male", y = y_pos, label = label),
                  inherit.aes = FALSE
                )
  

sex_freq

(sex_freq/sex_plot)+plot_layout(axes="collect")+plot_annotation(tag_levels = "A")


#Figure 5
eggs <- eggs %>%
  mutate(
    dCt_modified = ifelse(dCt == 0, 1e-6, dCt),
    is_pseudocount = ifelse(dCt == 0, "pseudocount", "real")
  )

eggs$Infected <- factor(eggs$Infected, levels=c("Yes", "No"))


egg_plot<-ggplot(eggs, aes(x=Parent_Infected, y=dCt_modified, color=Parent_Infected, group=Parent_Infected))+
  geom_jitter(width=.2, aes(alpha = is_pseudocount), show.legend = FALSE)+
  scale_alpha_manual(values = c("real" = 1, "pseudocount" = 0.15), guide="none")+
  ggtitle("")+ scale_y_log10()+
  ylab("Galbut virus RNA \nrelative to RpL32 mRNA")+
  xlab("Transmission Type")+coord_cartesian(clip="off")+
  stat_compare_means(paired=FALSE, method="wilcox", show.legend = FALSE)+
  labs(color="", alpha="")+theme_bw()+facet_grid(~Cross)+  scale_color_manual(values= c("Maternal"= "orange2", "Paternal"="darkred"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


egg_plot


infection_data_eggs <- eggs %>%
  filter(Infected %in% c("Yes", "No")) %>%
  group_by(Cross, Parent_Infected) %>%
  summarise(
    infected_count = sum(Infected == "Yes", na.rm = TRUE),
    total_count = n(),
    .groups = "drop"
  ) %>%
  filter(total_count > 0) %>%
  mutate(
    infection_freq = infected_count / total_count,
    percent_label = paste0(round(infection_freq * 100, 1), "%")
  ) %>%
  bind_cols(
    binom.confint(x = .$infected_count, n = .$total_count, methods = "wilson")[, c("lower", "upper")]
  )

# Fisher’s exact test for each Cross
fisher_results <- infection_data_eggs %>%
  group_by(Cross) %>%
  summarise(
    fisher_p = {
      if (all(c("Maternal", "Paternal") %in% Parent_Infected)) {
        tab <- matrix(
          c(
            infected_count[Parent_Infected == "Maternal"],
            total_count[Parent_Infected == "Maternal"] - infected_count[Parent_Infected == "Maternal"],
            infected_count[Parent_Infected == "Paternal"],
            total_count[Parent_Infected == "Paternal"] - infected_count[Parent_Infected == "Paternal"]
          ),
          nrow = 2,
          byrow = TRUE,
          dimnames = list(
            Parent = c("Maternal", "Paternal"),
            Infection = c("Infected", "Not Infected")
          )
        )
        fisher.test(tab)$p.value
      } else {
        NA
      }
    },
    .groups = "drop"
  )

# Merge p-values into main table
infection_data_eggs <- left_join(infection_data_eggs, fisher_results, by = "Cross")


#make a p-value label 
pval_labels <- infection_data_eggs %>%
  distinct(Cross, fisher_p) %>%
  mutate(
    label = paste0("p-value = ", signif(fisher_p, 2)),
    Parent_Infected = "Paternal",
    y_pos = 1.05
  )



####plot
egg_freq<-ggplot(infection_data_eggs, aes(x = Parent_Infected, y = infection_freq, color = Parent_Infected)) +
  geom_point(show.legend = FALSE) +
  facet_grid(~Cross) +
  ylab("Prevalence of \nGalbut virus Infection") +
  xlab("Transmission Type") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 9),
    panel.border = element_rect(color = "black", fill = NA, linewidth = .85)
    
    
  )+
  scale_color_manual(values=c("orange2", "darkred"))+
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 0.9), width = 0.2, show.legend=FALSE)+
  geom_text(
    data = pval_labels,
    aes(x = Parent_Infected, y = y_pos, label = label),
    inherit.aes = FALSE
  )

egg_freq

(egg_freq/ egg_plot)+plot_layout(axes="collect")+plot_annotation(tag_levels = "A")



#Figure 6
egg_28152 <- egg_28152 %>%
  mutate(
    dCt_modified = ifelse(dCt == 0, 1e-7, dCt),
    is_pseudocount = ifelse(dCt == 0, "pseudocount", "real")
  )
egg_28152$Parent_Infected <- factor(egg_28152$Parent_Infected, levels=c("Maternal", "Paternal", "Neither"))

egg_28152$Infected <- factor(egg_28152$Infected, levels=c("Yes", "No"))

two<- ggplot(egg_28152, aes(x=Parent_Infected, y=dCt_modified, group=Parent_Infected, color=dpi))+
  geom_jitter(width=.2, aes(alpha = is_pseudocount), show.legend = FALSE)+
  scale_alpha_manual(values = c("real" = 1, "pseudocount" = 0.15), guide="none")+
  ggtitle("")+ scale_y_log10()+
  ylab("Galbut virus RNA \nrelative to RpL32 mRNA")+
  xlab("")+coord_cartesian(clip="off")+
  labs(color="Tested \nPositive", alpha="")+theme_bw()+
  geom_signif(
    comparisons = list(c("Maternal", "Paternal")), 
    y_position = c(1e0), 
    tip_length = 0.01,
    textsize = 4, show.legend = FALSE, color="black")+
  facet_grid(~dpi)+
  scale_color_brewer(palette= "Dark2")+  xlab("Transmission Type")

two


parent_data_eggs <- egg_28152 %>%
  filter(Infected %in% c("Yes", "No")) %>%
  group_by(dpi, Parent_Infected) %>%
  summarise(
    infected_count = sum(Infected == "Yes", na.rm = TRUE),
    total_count = n(),
    .groups = "drop"
  ) %>%
  filter(total_count > 0) %>%
  mutate(
    infection_freq = infected_count / total_count,
    percent_label = paste0(round(infection_freq * 100, 1), "%")
  ) %>%
  bind_cols(
    binom.confint(x = .$infected_count, n = .$total_count, methods = "wilson")[, c("lower", "upper")]
  )

# Reorder factor levels
parent_data_eggs <- parent_data_eggs %>%
  mutate(Parent_Infected = factor(Parent_Infected, levels = c("Maternal", "Paternal", "Neither")))


# Calculate Fisher's exact test per dpi (compare Maternal vs Paternal only)
fisher_results_parent <- parent_data_eggs %>%
  filter(Parent_Infected %in% c("Maternal", "Paternal")) %>%  
  group_by(dpi) %>%
  summarise(
    fisher_p = {
      parents <- c("Maternal", "Paternal")
      if (all(parents %in% Parent_Infected)) {
        tab <- matrix(
          c(
            infected_count[Parent_Infected == "Maternal"],
            total_count[Parent_Infected == "Maternal"] - infected_count[Parent_Infected == "Maternal"],
            infected_count[Parent_Infected == "Paternal"],
            total_count[Parent_Infected == "Paternal"] - infected_count[Parent_Infected == "Paternal"]
          ),
          nrow = 2,
          byrow = TRUE,
          dimnames = list(
            Parent = parents,
            Infection = c("Infected", "Not Infected")
          )
        )
        fisher.test(tab)$p.value
      } else {
        NA_real_
      }
    },
    .groups = "drop"
  )

# Create label data frame for plotting (one label per dpi)
pval_labels_parent <- dpi_data_eggs %>%
  distinct(dpi) %>%                     
  left_join(fisher_results_parent, by = "dpi") %>%
  mutate(
    label = paste0("p-value = ", signif(fisher_p, 2)),
    y_pos = 1.05                          
  )




parent_freq<-ggplot(parent_data_eggs, aes(x = Parent_Infected, y = infection_freq, color = dpi)) +
  geom_point(show.legend = FALSE) +
  facet_grid(~dpi) +
  ylab("Prevalence of \nGalbut virus Infection") +

  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 9),
    panel.border = element_rect(color = "black", size = 0.9))+
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 0.9), width = 0.2, show.legend=FALSE)+
  scale_color_brewer(palette="Dark2")+
  geom_text(
    data = pval_labels_parent,
    aes(x = "Maternal", y = y_pos, label = label),
    inherit.aes = FALSE,
    hjust = -0.1 
  )+xlab("")


parent_freq

one<- ggplot(egg_28152, aes(x=dpi, y=dCt_modified, color=dpi))+
  geom_jitter(width=.2, aes(alpha = is_pseudocount), show.legend = FALSE)+
  scale_alpha_manual(values = c("real" = 1, "pseudocount" = 0.15), guide="none")+
  ggtitle("")+ scale_y_log10()+
  ylab("Galbut virus RNA \nrelative to RpL32 mRNA")+
  xlab("")+coord_cartesian(clip="off")+
  labs(color="Tested \nPositive", alpha="")+theme_bw()+
  geom_signif(
    comparisons = list(c("2 days post mating", "8 days post mating")), 
    map_signif_level = TRUE,
    y_position = c(1e0), 
    tip_length = 0.01,
    textsize = 4,
    color= "black")+
  facet_grid(~Parent_Infected)+
  scale_color_brewer(palette= "Dark2")+  xlab("Days Post Mating")+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


one

dpi_data_eggs <- egg_28152 %>%
  filter(Infected %in% c("Yes", "No")) %>%
  group_by(dpi, Parent_Infected) %>%
  summarise(
    infected_count = sum(Infected == "Yes", na.rm = TRUE),
    total_count = n(),
    .groups = "drop"
  ) %>%
  filter(total_count > 0) %>%
  mutate(
    infection_freq = infected_count / total_count,
    percent_label = paste0(round(infection_freq * 100, 1), "%")
  ) %>%
  bind_cols(
    binom.confint(x = .$infected_count, n = .$total_count, methods = "wilson")[, c("lower", "upper")]
  )


fisher_results_eggs <- dpi_data_eggs %>%
  group_by(Parent_Infected) %>%
  summarise(
    fisher_p = {
      dpis <- sort(unique(dpi))
      if (length(dpis) == 2) {
        tab <- matrix(
          c(
            infected_count[dpi == dpis[1]],
            total_count[dpi == dpis[1]] - infected_count[dpi == dpis[1]],
            infected_count[dpi == dpis[2]],
            total_count[dpi == dpis[2]] - infected_count[dpi == dpis[2]]
          ),
          nrow = 2,
          byrow = TRUE,
          dimnames = list(
            DPI = dpis,
            Infection = c("Infected", "Not Infected")
          )
        )
        fisher.test(tab)$p.value
      } else {
        NA_real_
      }
    },
    .groups = "drop"
  )

# Make label data frame for plotting (one label per facet)
pval_labels_dpi <- dpi_data_eggs %>%
  distinct(Parent_Infected) %>%     
  left_join(fisher_results_eggs, by = "Parent_Infected") %>%
  mutate(
    label = paste0("p-value = ", signif(fisher_p, 2)),
    y_pos = 1.05 
  )




dpi_freq<-ggplot(dpi_data_eggs, aes(x = dpi, y = infection_freq, color = dpi)) +
  geom_point(show.legend = FALSE) +
  facet_grid(~Parent_Infected) +
  ylab("Prevalence of \nGalbut virus Infection") +

  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "black", size = 0.9)
    
  )+
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 0.9), width = 0.2, show.legend=FALSE)+
  scale_color_brewer(palette="Dark2")+
  geom_text(
    data = pval_labels_dpi,
    aes(x = "2 days post mating", y = y_pos, label = label),
    inherit.aes = FALSE,
    hjust = -0.05 
  )+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+xlab("")

dpi_freq


dpi_freq+parent_freq+one+two+plot_layout(ncol=2) +plot_annotation(tag_levels = list(c('A','C','B','D')))+plot_layout(axes="collect")


compare_means(dCt~dpi, data=egg_28152, method="wilcox", paired = FALSE)
compare_means(dCt~Parent_Infected, data=egg_28152, method="wilcox", paired = FALSE)
compare_means(dCt~Treatment, data=egg_28152, method="wilcox", paired = FALSE)

summary_table_eggs <- egg_28152 %>%
  group_by(Parent_Infected, dpi, Infected) %>%
  summarise(count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = Infected,
    values_from = count,
    values_fill = 0
  )

summary_table_eggs




#Supplemental Figure 3
pct <- egg_28152 %>%
  group_by(Parent_Infected, Treatment, dpi) %>%
  summarize(
    total = n(),
    infected = sum(Infected == "Yes"),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    prop_test = list(binom.test(infected, total)),
    percent_infected = (infected / total) * 100,
    lower = prop_test$conf.int[1] * 100,
    upper = prop_test$conf.int[2] * 100
  ) %>%
  ungroup() %>%
  select(-prop_test)

ggplot(pct, aes(x = Parent_Infected, y = percent_infected, color = Treatment, group = dpi)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.2,
    position = position_dodge(width = 0.5)
  ) +
  labs(
    x = "Parent Infected",
    y = "Percent Infected",
    color = "Treatment"
  ) +
  ylim(0, 100) +
  theme_bw()+facet_grid(~dpi)

pct %>%
  group_by(Parent_Infected, dpi) %>%
  get_summary_stats(percent_infected, type = "mean_sd")
res.aov <- pct %>% anova_test(percent_infected ~ Parent_Infected * dpi)
res.aov


egg_28152 %>%
  group_by(Parent_Infected, Treatment, dpi) %>%
  summarize(
    total = n(),
    infected = sum(Infected == "Yes"),
    percent_infected = (infected / total) * 100,
    .groups = "drop"
  )



