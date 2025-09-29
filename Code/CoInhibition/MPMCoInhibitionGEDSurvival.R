#' Analyze Panel: GED vs Survival
#'
#' This script examines the relationship between Graph Edit Distance (GED) of
#' single-cell phenotypes and patient survival outcomes for the Co-inhibition panel.
#' It loads precomputed GED values and patient response data, merges them with
#' overall and progression-free survival information, and visualizes the results.
#' The script produces line plots of survival versus GED for responders and
#' non-responders, as well as boxplots of survival and GED distributions by response group.
#' All outputs are saved as SVG figures.
#' 
setwd(dirname(this.path::this.path()))
source("../MPMUtilities.R")
setwd(dirname(this.path::this.path()))


coinhib_ged_path = "../../Data/Raw Data/Co-inhibition/Final Cell Assigned Data/Co-inhibition_GED.xlsx"
coinhib_ged = read.xlsx(coinhib_ged_path, sheet = "Sheet1")
coinhib_ged = na.omit(coinhib_ged)
coinhib_ged$Response = c("Responder","Responder","Responder","Responder","Non-responder","Non-responder",
                         "Non-responder","Responder","Responder","Non-responder","Non-responder",
                         "Non-responder","Responder")
coinhib_ged$Response = factor(coinhib_ged$Response, levels = c("Non-responder", "Responder"))


os_pfs = read.xlsx("../../Data/Raw Data/OS_PFS.xlsx", sheet = "Sheet1")
os_pfs = os_pfs[-3, ]
os_pfs$Response = coinhib_ged$Response
os_pfs$GED = coinhib_ged$Baseline_After3Vac


# Reshape the data to long format for easier plotting
os_pfs_long <- os_pfs %>%
  pivot_longer(cols = c(OS, PFS),
               names_to = "SurvivalType",
               values_to = "SurvivalValue")

# Plot for Responders
ggplot(os_pfs_long[os_pfs_long$Response == "Responder", ],
       aes(x = GED, y = SurvivalValue, color = SurvivalType)) +
  #geom_jitter(width = 0.2, height = 0.2, size = 2) +
  geom_line()+
  #geom_smooth(method = "loess", se = FALSE) +
  labs(title = "Responder: Survival vs GED",
       x = "GED",
       y = "Survival Time",
       color = "Survival Type") +
  theme_minimal()

# Plot for Non-Responders
ggplot(os_pfs_long[os_pfs_long$Response == "Non-responder", ],
       aes(x = GED, y = SurvivalValue, color = SurvivalType)) +
  #geom_jitter(width = 0.2, height = 0.2, size = 2) +
  geom_line()+
  #geom_smooth(method = "loess", se = FALSE) +
  labs(title = "Non-responder: Survival vs GED",
       x = "GED",
       y = "Survival Time",
       color = "Survival Type") +
  theme_minimal()


# Create boxplot
g1=ggplot(os_pfs_long, aes(x = Response, y = SurvivalValue, fill = SurvivalType)) +
  geom_boxplot(position = position_dodge(0.5), width = 0.4, alpha=0.7, size=0.3, outlier.size = 0.5) +
  scale_fill_manual(values = mycolor[c(10, 12)]) +
  theme(legend.position="right", legend.text=element_text(size=8), legend.title = element_text(size=9),
        legend.key = element_blank(),
        #legend.box.margin=margin(-10,-10,-10,-10),
        plot.title = element_text(hjust = 0.5, size=10),
        plot.subtitle = element_text(hjust = 0.5, size=10),
        axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        panel.background = element_rect(fill='white', colour='grey'),
        panel.grid.major = element_line(color = "grey", linewidth=0.25, linetype=2)) +
  labs(title = "Survival distributions by Response", x = "Response", y = "Months from CRS-HIPEC")
svglite("../../Output/ged/CoInhib_Survival_Response.svg", width = 3.5, height = 3)
plot(g1)
dev.off()

# Boxplot of GED by Response group
g2 = ggplot(os_pfs, aes(x = Response, y = GED, fill = Response)) +
  geom_boxplot(width = 0.2, alpha=0.8, size=0.3, outlier.size = 0.5) +
  scale_fill_manual(values = mycolor[c(2, 1)]) +
  theme(legend.position="right", legend.text=element_text(size=8), legend.title = element_text(size=9),
        legend.key = element_blank(),
        #legend.box.margin=margin(-10,-10,-10,-10),
        plot.title = element_text(hjust = 0.5, size=10),
        plot.subtitle = element_text(hjust = 0.5, size=10),
        axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        panel.background = element_rect(fill='white', colour='grey'),
        panel.grid.major = element_line(color = "grey", linewidth=0.25, linetype=2)) +
  labs(title = "GED distribution by Response", x = "Response", y = "GED ( Baseline, After 3 vac )")
svglite("../../Output/ged/CoInhib_GED_Response.svg", width = 3.9, height = 3)
plot(g2)
dev.off()
