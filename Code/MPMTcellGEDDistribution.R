#' Analyze and Visualize Graph Edit Distances (GED) Across Panels and Patients
#'
#' This script processes Graph Edit Distance (GED) data for Co-inhibition, Co-stimulation, 
#' and Cytokine panels across multiple patients and time points. Specifically, it:
#' 
#' 1. Loads GED data from Excel files for each panel.
#' 2. Cleans the data by removing missing values and orders patients based on response status.
#' 3. Transforms the data into a long format suitable for plotting, including creating 
#'    factors for time points.
#' 4. Generates horizontal scatter and bar plots of GED values for individual patients, 
#'    highlighting responders vs non-responders and time-dependent changes.
#' 5. Performs statistical analysis:
#'    - Paired t-tests to determine whether GED changes from baseline are significant.
#'    - Repeated measures ANOVA to assess differences across time points within each panel.
#'    - Tukey's HSD post-hoc tests for multiple comparisons.
#' 6. Combines GED data from all panels to perform overall statistical comparisons.
#' 7. Saves plots as high-resolution SVG files and visualizes GED distributions using the 
#'    `compareGEDDistribution()` function.
#'
#' @note Patients with missing data (rows 3, 11, 12) are excluded from analysis. 
#'       The script assumes the presence of response labels for Co-inhibition patients 
#'       and a predefined color palette `mycolor` for plotting.

setwd(dirname(this.path::this.path()))
source("MPMUtilities.R")
setwd(dirname(this.path::this.path()))


coinhib_ged_path = "../Data/Raw Data/Co-inhibition/Final Cell Assigned Data/Co-inhibition_GED.xlsx"
costim_ged_path = "../Data/Raw Data/Co-stimulation/Final Cell Assigned Data/Co-stimulation_GED.xlsx"
cytokine_ged_path = "../Data/Raw Data/Cytokine/Final Cell Assigned Data/Cytokine_GED.xlsx"


coinhib_ged = read.xlsx(coinhib_ged_path, sheet = "Sheet1")
costim_ged = read.xlsx(costim_ged_path, sheet = "Sheet1")
cytokine_ged = read.xlsx(cytokine_ged_path, sheet = "Sheet1")


coinhib_ged = na.omit(coinhib_ged)
coinhib_ged$Response = c("Responder","Responder","Responder","Responder","Non-responder","Non-responder",
                         "Non-responder","Responder","Responder","Non-responder","Non-responder",
                         "Non-responder","Responder")
coinhib_ged$Response = factor(coinhib_ged$Response, levels = c("Non-responder", "Responder"))

coinhib_ged <- coinhib_ged[order(coinhib_ged$Response), ]
coinhib_ged <- coinhib_ged %>%
  arrange(Response, (Baseline_After1Vac)) %>%
  mutate(Baseline.ID = factor(Baseline.ID, levels = unique(Baseline.ID)))

coinhib_ged_2 = data.frame(Patient=rep(coinhib_ged$Baseline.ID, 3), Response=rep(coinhib_ged$Response, 3),
                           GED=c(rep(0, nrow(coinhib_ged)), coinhib_ged$Baseline_After1Vac, coinhib_ged$Baseline_After3Vac),
                           Time = c(rep("Baseline", nrow(coinhib_ged)), rep("After 1 vac", nrow(coinhib_ged)), rep("After 3 vac", nrow(coinhib_ged))) )
coinhib_ged_2$Time <- factor(coinhib_ged_2$Time, levels = c("Baseline", "After 1 vac", "After 3 vac"))

g = ggplot() +
  geom_point(data=coinhib_ged_2, aes(x = Patient, y = GED, color = Response, shape = Time), cex=2, alpha=0.8) +
  geom_line(data=coinhib_ged_2,aes(x = Patient, y = GED, color = Response, group = Patient), linetype="dotted", show.legend = F) +
  #geom_point(data=coinhib_ged, aes(x = Baseline.ID, y = Baseline_After3Vac, color = Baseline.ID), pch=7, cex=3) +
  coord_flip() + 
  theme(legend.position="right", legend.text=element_text(size=8), legend.title = element_text(size=10),
        legend.key=element_blank(),
        #legend.box.margin=margin(-10,-10,-10,-10),
        plot.title = element_text(hjust = 0.5, size=10),
        plot.subtitle = element_text(hjust = 0.5, size=10),
        axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        panel.background = element_rect(fill='white', colour='grey'),
        panel.grid.major = element_line(color = "white", linewidth=0.25, linetype=2)) +
  scale_colour_manual(values = mycolor[c(2, 1)]) +
  #guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 3)) +
  labs(title="GED between FC samples", x = "Patient ID", y = "GED")
svglite("../Output/ged/CoInhib_GED_Horz_No13_legend.svg", width = 6, height = 3)
plot(g)
dev.off()

#################################################################
coinhib_ged_3 = data.frame(Patient=rep(coinhib_ged$Baseline.ID, 3), 
                           GED=c(rep(0, nrow(coinhib_ged)), coinhib_ged$Baseline_After1Vac, coinhib_ged$After1Vac_After3Vac),
                           Time = c(rep("Baseline", nrow(coinhib_ged)), rep("After 1 vac", nrow(coinhib_ged)), rep("After 3 vac", nrow(coinhib_ged))) )
coinhib_ged_3$Time <- factor(coinhib_ged_3$Time, levels = c("After 3 vac", "After 1 vac", "Baseline"))
ggplot() +
  geom_bar(data=coinhib_ged_3, aes(x = Patient, y = GED), 
           stat = "identity", alpha=0.9, color="black", fill="white", width=0.5) +
  coord_flip() +
  theme(legend.position="right", legend.text=element_text(size=10), legend.title = element_text(size=12),
        #legend.box.margin=margin(-10,-10,-10,-10),
        plot.title = element_text(hjust = 0.5, size=10),
        plot.subtitle = element_text(hjust = 0.5, size=10),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        panel.background = element_rect(fill='white', colour='black'),
        panel.grid.major = element_line(color = "grey", linewidth=0.25, linetype=2)) +
  labs(title = "", x = "Patient", y = "Graph Edit Distance")


#### paired t-test and ANOVA
coinhib_patient = na.omit(coinhib_ged[-c(3, 11, 12), 1])
coinhib_ged = na.omit(coinhib_ged[-c(3, 11, 12), c(4,6)])
summary(coinhib_ged)

# whether the mean of these differences is significantly different from zero
# if there is a significant change in patient status from baseline to after 1 vaccination.
t.test(coinhib_ged$Baseline_After1Vac, mu=0)
#t.test(coinhib_ged$After1Vac_After3Vac, mu=0)
t.test(coinhib_ged$Baseline_After3Vac, mu=0)

stacked_coinhib_ged = utils::stack(coinhib_ged)
stacked_coinhib_ged$patient = rep(coinhib_patient, 2)
coinhib_anova = aov(values~ind+Error(patient/ind), data = stacked_coinhib_ged)
summary(coinhib_anova)


######################################################
costim_patient = na.omit(costim_ged[-c(3, 11, 12), 1])
costim_ged = na.omit(costim_ged[-c(3, 11, 12), c(4,6)])
summary(costim_ged)

t.test(costim_ged$Baseline_After1Vac, mu=0)
# t.test(costim_ged$After1Vac_After3Vac, mu=0)
t.test(costim_ged$Baseline_After3Vac, mu=0)

stacked_costim_ged = utils::stack(costim_ged)
stacked_costim_ged$patient = rep(costim_patient, 2)
costim_anova = aov(values~ind+Error(patient/ind), data = stacked_costim_ged)
summary(costim_anova)

# Perform Tukey's HSD Test
posthoc <- TukeyHSD(aov(values~ind, data = stacked_costim_ged))
print(posthoc)
plot(posthoc)


######################################################
cytokine_patient = na.omit(cytokine_ged[-c(3, 11, 12), 1])
cytokine_ged = na.omit(cytokine_ged[-c(3, 11, 12), c(4,6)])
summary(cytokine_ged)

t.test(cytokine_ged$Baseline_After1Vac, mu=0)
# t.test(cytokine_ged$After1Vac_After3Vac, mu=0)
t.test(cytokine_ged$Baseline_After3Vac, mu=0)

stacked_cytokine_ged = utils::stack(cytokine_ged)
stacked_cytokine_ged$patient = rep(cytokine_patient, 2)
cytokine_anova = aov(values~ind+Error(patient/ind), data = stacked_cytokine_ged)
summary(cytokine_anova)


######################################################
all_ged = rbind(coinhib_ged, costim_ged, cytokine_ged)
t.test(all_ged$Baseline_After1Vac, mu=0)
#t.test(all_ged$After1Vac_After3Vac, mu=0)
t.test(all_ged$Baseline_After3Vac, mu=0)


stacked_all_ged = rbind(stacked_coinhib_ged, stacked_costim_ged, stacked_cytokine_ged)
stacked_all_ged$panel = c(rep("coinhib", nrow(stacked_coinhib_ged)),
                          rep("costim", nrow(stacked_costim_ged)),
                          rep("cytokine", nrow(stacked_cytokine_ged)))
all_anova = aov(values~(ind*panel)+Error(patient/(ind*panel)), data = stacked_all_ged)
summary(all_anova)
# Perform Tukey's HSD Test
posthoc <- TukeyHSD(aov(values~ind, data = stacked_all_ged))
print(posthoc)

svglite(paste("../Output/ged/", "all_ged.svg", sep=""), width = 4, height = 4)
compareGEDDistribution(all_ged)
dev.off()
