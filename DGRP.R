library(readxl)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(rstatix)
library(dplyr)
library(binom)
library(patchwork)


DGRP<-read_excel("./1DGRP_cross.xlsx")

eggs<-read_excel("./eggs.xlsx")
egg_28152<-read_excel("./28152_eggs.xlsx")


###################################
###################################
###################################LARGE DGRP PLOT
###################################


shared_y_range <- c(1e-9, 1e5)

#Make pseudocounts for 0 infections
DGRP <- DGRP %>%
  mutate(
    dCt_modified = ifelse(dCt == 0, 1e-10, dCt),
    is_pseudocount = ifelse(dCt == 0, "pseudocount", "real")
  )

# Calculate p-values and convert to significance stars

p_values <- DGRP %>%
  group_by(Generation, Cross) %>%
  filter(n_distinct(Parent_Infected) == 2) %>%   # filter here so the test runs only when valid
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


# Plot
plot <- ggplot(data = DGRP, aes(x = Parent_Infected, y = dCt, group = Cross, Generation)) +
  geom_jitter(
    height = 0, width = 0.2,
    aes(color = as.factor(Replicate), shape = Offspring_Sex, alpha = is_pseudocount), size=1.5
  ) +
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


  

########## GEN ON THE X
p_values_by_generation <- DGRP %>%
  group_by(Cross, Parent_Infected) %>%
  summarise(
    p_value = tryCatch(
      kruskal.test(dCt_modified ~ Generation)$p.value,
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

  # Plot
  gen_plot <- ggplot(data = DGRP, aes(x = Generation, y = dCt, group = Cross, Parent_Infected)) +
    geom_jitter(
      height = 0, width = 0.2,
      aes(color = as.factor(Replicate), shape = Offspring_Sex, alpha = is_pseudocount), size=1.5
    ) +
    scale_alpha_manual(values = c("real" = 1, "pseudocount" = 0.3), guide = FALSE) +
    facet_grid(Parent_Infected ~ Cross, drop = TRUE) +
    theme(strip.text.x.top = element_text(size = 6)) +
    scale_y_log10(limits= shared_y_range)+
    ylab("Galbut virus RNA relative to RpL32 mRNA") +
    xlab("Original Transmission Type") +
    labs(color = "Replicate Number", shape = "Offspring Sex") +
    coord_cartesian(clip = "off")+theme_bw()+
    scale_color_brewer(palette= "Dark2")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    stat_compare_means(method="kruskal.test")

# Show the plot
gen_plot



####male v female parent infected comparison
mvf_p_values <- DGRP %>%
  group_by(Generation) %>%
  filter(n_distinct(Parent_Infected) == 2) %>%   # filter here so the test runs only when valid
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



mvf<- ggplot(DGRP, aes(x = Parent_Infected, y = dCt, group = Parent_Infected,  alpha = is_pseudocount, shape=Offspring_Sex, color=as.factor(Replicate))) +
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

(plot + mvf) + 
  plot_layout(widths=c(1,.1), guides= "collect", axes="collect")+plot_annotation(tag_levels = "A")&theme(legend.position = "bottom")












#####b ig plot arranged differnetly

# Plot
plot <- ggplot(data = DGRP, aes(x = Cross, y = dCt, group = Cross)) +
  geom_jitter(
    height = 0, width = 0.2,
    aes(color = as.factor(Replicate), shape = Offspring_Sex, alpha = is_pseudocount)
  ) +
  scale_alpha_manual(values = c("real" = 1, "pseudocount" = 0.3), guide = FALSE) +
  facet_grid(Parent_Infected ~ Generation, drop = TRUE) +
  theme(strip.text.x.top = element_text(size = 6)) +
  scale_y_log10() +
  ylab("Galbut virus relative to RpL32") +
  xlab("Infected parent") +
  labs(color = "Replicate", shape = "Offspring Sex") +
theme_bw()

# Show the plot
plot












######Frequency plot
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




ggplot(infection_data, aes(x = Cross, y = infection_freq, fill = factor(Cross))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), show.legend = TRUE) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 0.9), width = 0.2) +

  facet_grid(Parent_Infected~Generation) +
  labs(
    x = "Parent Originally Infected",
    y = "Frequency of Infection",
    title = "Frequency of Infection by Parent Infected, Generation, and Cross"
  ) +
  scale_fill_brewer(palette= "Dark2")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylim(0,1.2)+ labs(fill = "Cross")




####line version
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
  scale_color_manual(values= c("Maternal"= "hotpink", "Paternal"="skyblue3")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 8)
  ) +
  ylim(0, 1.2)+
  theme(legend.position="bottom")






#####function to plot individual crosses
#############SEE IF STATS ARE ACCURATE, MIGHT NEED TO CHANGE THE WAY THEY ARE CALCULATED!!!!!!!!!!!

plot_cross_data <- function(data, cross_value) {
  # Filter for each cross
  filtered_data <- data[data$Cross == cross_value, ]
  
  # Add a pseudocount to the dCt values where they are 0 (replace 0 with 1e-10)
  filtered_data$dCt_modified <- ifelse(filtered_data$dCt == 0, 1e-10, filtered_data$dCt)
  
  # Plot
  p <- ggplot(filtered_data, aes(
    x = Parent_Infected,
    y = dCt_modified,
    color = as.factor(Replicate),
    shape = Offspring_Sex,
    group = Parent_Infected
  )) +
    geom_jitter(height = 0, width = 0.1) +
    facet_wrap(~Generation) +
    scale_y_log10() +
    ylab("Galbut virus relative to RpL32 (log scale, pseudocounts applied)") +
    xlab("Initial infected parent") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    theme(strip.text.x.top = element_text(size = 6)) +
    

    stat_compare_means(
      aes(x = Parent_Infected, y = dCt_modified),
      data = filtered_data,  # Use the modified data for stats
      method = "wilcox", paired = FALSE
    )+
  theme_bw()+
    labs(color= "Replicate", shape = "Offspring Sex")
  
  return(p)
}

plot_cross_data(DGRP, "W1118")









#####offspring sex
ggplot(DGRP, aes(x=Offspring_Sex, y=dCt_modified, color=Cross, group=Offspring_Sex))+
  geom_jitter(width=.2)+stat_compare_means(method= "wilcox", paired= FALSE, label= "p.signif", label.x.npc = "center")+
  ylab("Galbut virus Ct relative to RpL32")+xlab("Offspring Sex")+labs(color="Fly Line")+
  facet_grid(~Generation)+scale_y_log10()






sex_plot<-ggplot(DGRP, aes(x=Offspring_Sex, y=dCt_modified, color=Offspring_Sex, group=Offspring_Sex))+
  geom_jitter(width=.2, aes(alpha = is_pseudocount), show.legend = FALSE)+
  scale_alpha_manual(values = c("real" = 1, "pseudocount" = 0.15), guide="none")+
  ggtitle("")+ scale_y_log10()+
  ylab("Galbut virus RNA \nrelative to RpL32 mRNA")+
  xlab("Offspring Sex")+coord_cartesian(clip="off")+
  labs(color="", alpha="")+theme_bw()+ 
  facet_wrap(~Generation)+scale_color_manual(values= c("Female"= "mediumorchid", "Male"="palegreen2"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  stat_compare_means(paired=FALSE, method="wilcox", label="p.signif", label.x.npc = "center", show.legend = FALSE)



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




sex_freq<-ggplot(infection_data_sex, aes(x = Offspring_Sex, y = infection_freq, color=Offspring_Sex)) +
  geom_point(show.legend = FALSE) +
  ylab("Prevalence of \nGalbut virus Infection") +
  xlab("Offspring Sex") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 9),
    panel.border = element_rect(color = "black", size = 0.9)
    
  )+
  facet_wrap(~Generation, nrow=1)+
  scale_color_manual(values=c("mediumorchid", "palegreen2"))+
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.2, show.legend=FALSE)

sex_freq

(sex_plot/ sex_freq)+plot_layout(axes="collect")+plot_annotation(tag_levels = "A")







compare_means(dCt~Parent_Infected, data=DGRP, method="wilcox", paired= FALSE)

compare_means(dCt~Offspring_Sex, data=DGRP, method="wilcox", paired = FALSE)



######EMBRYOS
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
  stat_compare_means(paired=FALSE, method="wilcox", label="p.signif", label.x.npc = "center", show.legend = FALSE)+
  labs(color="", alpha="")+theme_bw()+facet_grid(~Cross)+  scale_color_manual(values= c("Maternal"= "orange2", "Paternal"="darkred"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

compare_means(dCt~Parent_Infected, data=eggs, method="wilcox", paired = FALSE)




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
    panel.border = element_rect(color = "black", size = 0.9)
    
  )+
  scale_color_manual(values=c("orange2", "darkred"))+
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 0.9), width = 0.2, show.legend=FALSE)

egg_freq

(egg_plot/ egg_freq)+plot_layout(axes="collect")+plot_annotation(tag_levels = "A")


#######eggs over time
egg_time<-read_excel("./egg_trial.xlsx")


ggplot(data=egg_time, aes(x=Combo, y=Count, color=Combo, fill=Combo))+
  geom_col()+facet_wrap(~Hours, nrow=1)+ theme_bw()+scale_color_brewer(palette="Spectral")+scale_fill_brewer(palette="Spectral")+
  labs(fill="Transmission", color= "Transmission")+xlab("Transmission Type")+ylab("Total Egg Count Per Day")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggplot(data=egg_time, aes(x=Hours, y=Count, fill = as.factor(Hours)))+
  geom_col()+facet_wrap(~Combo, nrow=1)+ theme_bw()+labs(fill= "Hours")


compare_means(Count~Combo, data=egg_time, method="wilcox", paired=FALSE)


#########25192
egg_time_25192<-read_excel("./Egg_Trial_25192.xlsx")


ggplot(data=egg_time_25192, aes(x=Parent_Infected, y=Count, color=Parent_Infected, fill=Parent_Infected))+
  geom_col()+facet_wrap(~Hours, nrow=1)+ theme_bw()+scale_color_brewer(palette="Spectral")+scale_fill_brewer(palette="Spectral")+
  labs(fill="Transmission", color= "Transmission")+xlab("Transmission Type")+ylab("Total Egg Count Per Day")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggplot(data=egg_time_25192, aes(x=Hours, y=Count, fill = as.factor(Hours)))+
  geom_col()+facet_wrap(~Parent_Infected, nrow=1)+ theme_bw()+labs(fill= "Hours")+
  scale_fill_brewer(palette = "Spectral")


compare_means(Count~Combo, data=egg_time, method="wilcox", paired=FALSE)




##########28152
egg_28152 <- egg_28152 %>%
  mutate(
    dCt_modified = ifelse(dCt == 0, 1e-7, dCt),
    is_pseudocount = ifelse(dCt == 0, "pseudocount", "real")
  )
egg_28152$Parent_Infected <- factor(egg_28152$Parent_Infected, levels=c("Maternal", "Paternal", "Neither"))

egg_28152$Infected <- factor(egg_28152$Infected, levels=c("Yes", "No"))

two<- ggplot(egg_28152, aes(x=Parent_Infected, y=dCt_modified, group=Parent_Infected, color=Infected))+
  geom_jitter(width=.2, aes(alpha = is_pseudocount))+
  scale_alpha_manual(values = c("real" = 1, "pseudocount" = 0.15), guide="none")+
  ggtitle("")+ scale_y_log10()+
  ylab("Galbut virus Ct relative to RpL32 mRNA")+
  xlab("Transmission Type")+coord_cartesian(clip="off")+
  labs(color="Tested \nPositive", alpha="")+theme_bw()+
  geom_signif(
    comparisons = list(c("Maternal", "Paternal")), # replace with actual levels
    map_signif_level = TRUE,
    y_position = c(1e0), # Adjust to fit your y scale (log10)
    tip_length = 0.01,
    textsize = 4)+
  facet_grid(~dpi)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_color_brewer(palette= "Dark2")


one<- ggplot(egg_28152, aes(x=dpi, y=dCt_modified, color=Infected))+
  geom_jitter(width=.2, aes(alpha = is_pseudocount))+
  scale_alpha_manual(values = c("real" = 1, "pseudocount" = 0.15), guide="none")+
  ggtitle("")+ scale_y_log10()+
  ylab("Galbut virus Ct relative to RpL32 mRNA")+
  xlab("Transmission Type")+coord_cartesian(clip="off")+
  labs(color="Tested \nPositive", alpha="")+theme_bw()+
  geom_signif(
    comparisons = list(c("2 days post mating", "8 days post mating")), # replace with actual levels
    map_signif_level = TRUE,
    y_position = c(1e0), # Adjust to fit your y scale (log10)
    tip_length = 0.01,
    textsize = 4,
    color= "black")+
  facet_grid(~Parent_Infected)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_color_brewer(palette= "Dark2")+xlab("Days Post Mating")


one/two+ plot_annotation(tag_levels = "A")


compare_means(dCt~dpi, data=egg_28152, method="wilcox", paired = FALSE)
compare_means(dCt~Parent_Infected, data=egg_28152, method="wilcox", paired = FALSE)
compare_means(dCt~Treatment, data=egg_28152, method="wilcox", paired = FALSE)

ggplot(egg_28152, aes(x=Parent_Infected, y=dCt_modified, group=Treatment, color= Treatment))+
  geom_jitter(width=.2, aes(alpha = is_pseudocount))+
  scale_alpha_manual(values = c("real" = 1, "pseudocount" = 0.15), guide="none")+
  ggtitle("")+ scale_y_log10()+
  ylab("Galbut virus Ct relative to RpL32")+
  xlab("")+coord_cartesian(clip="off")+
  labs(color="Treatment", alpha="")+theme_bw()+
  geom_signif(
    comparisons = list(c("Bleach", "EtOH")), # replace with actual levels
    map_signif_level = TRUE,
    y_position = c(1e0), # Adjust to fit your y scale (log10)
    tip_length = 0.01,
    textsize = 4)+
scale_color_brewer(palette="Dark2")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))





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



