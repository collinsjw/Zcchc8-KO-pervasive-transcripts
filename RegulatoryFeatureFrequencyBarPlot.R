library(tidyverse)

#Genomic Regulatory features from NIH Sims cells
#Bar plots in R are weird. I found it easier to simply enter the data manually than try to import and manipulate a dataset.
#These data come from the SharedRegFeaturesVennDiagrams.R file
df <- data.frame(
  "Gene" = rep(c("Zcchc8 KO"), 12), 
  "Feature" = c(rep("Enhancers", 2),
                rep("Promoters", 2), 
                rep("Promoter Flanking \n Regions", 2),
                rep("CTCF Sites", 2),
                rep("Transcription Factor \n Binding Sites", 2),
                rep("Open Chromatin \n Regions", 2)),
  "Regulation" = rep(c("Up", "Down"), 6), 
  "Number" = c(243, 28, 105, 22, 471, 51, 247, 19, 45, 6, 142, 9)
)

g <- ggplot() +
  geom_bar(data = df, aes(x = Feature, y = Number, fill = Regulation), stat = "identity", position = position_dodge(0.4), width = 0.5)

g + 
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c("dodgerblue", "firebrick")) + 
  labs(fill = "Expression") +
  xlab("") + ylab("") +
  guides(fill = guide_legend(reverse=TRUE)) +
  theme(
    #aspect.ratio = 1/1,
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.position = c(0.9, 0.92),
    legend.key.size = unit(0.1, "in"),
    legend.direction = "vertical",
    legend.title.align = 0,
    legend.background = element_rect(fill="transparent"),
    legend.key = element_blank(),
    axis.text.x=element_text(angle=45, hjust = 1, color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 8),
    panel.border = element_blank(),
    panel.grid.major.y = element_line(color = "black", linetype = "dashed", size = 0.1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = 'white'),
    # Change axis line
    axis.line = element_line(color = "black")
  )

#Genomic Regulatory features from Hopkins KO mice
df <- data.frame(
  "Gene" = rep(c("Zcchc8 KO"), 12), 
  "Feature" = c(rep("Enhancers", 2),
                rep("Promoters", 2), 
                rep("Promoter Flanking \n Regions", 2),
                rep("CTCF Sites", 2),
                rep("Transcription Factor \n Binding Sites", 2),
                rep("Open Chromatin \n Regions", 2)),
  "Regulation" = rep(c("Up", "Down"), 6), 
  "Number" = c(775, 12, 438, 25, 1342, 40, 1159, 14, 198, 1, 638, 7)
)

g <- ggplot() +
  geom_bar(data = df, aes(x = Feature, y = Number, fill = Regulation), stat = "identity", position = position_dodge(0.4), width = 0.5)

g + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1400, 500), minor_breaks = seq(0, 1300, 100), limits = c(0, 1400)) +
  scale_fill_manual(values = c("dodgerblue", "firebrick")) + 
  labs(fill = "Expression") +
  xlab("") + ylab("") +
  guides(fill = guide_legend(reverse=TRUE)) +
  theme(legend.position = "none",
    #aspect.ratio = 1/1,
    #legend.title = element_text(size = 8),
    #legend.text = element_text(size = 8),
    #legend.position = c(0.9, 0.89),
    #legend.key.size = unit(0.1, "in"),
    #legend.direction = "horizontal",
    #legend.title.align = 0,
    #legend.background = element_rect(fill="transparent"),
    legend.key = element_blank(),
    axis.text.x=element_text(angle=45, hjust = 1, color = "black", size = 7),
    axis.text.y = element_text(color = "black", size = 7),
    panel.border = element_blank(),
    panel.grid.major.y = element_line(color = "black", size = 0.25),
    panel.grid.minor.y = element_line(color = "black", linetype = "dashed", size = 0.1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = 'white'),
    # Change axis line
    axis.line = element_line(color = "black")
  )


#Genomic Regulatory features from Shanghai KO mice
df <- data.frame(
  "Gene" = rep(c("Zcchc8 KO"), 12), 
  "Feature" = c(rep("Enhancers", 2),
                rep("Promoters", 2), 
                rep("Promoter Flanking \n Regions", 2),
                rep("CTCF Sites", 2),
                rep("Transcription Factor \n Binding Sites", 2),
                rep("Open Chromatin \n Regions", 2)),
  "Regulation" = rep(c("Up", "Down"), 6), 
  "Number" = c(152, 47, 143, 68, 235, 86, 101, 27, 21, 9, 61, 17)
)

g <- ggplot() +
  geom_bar(data = df, aes(x = Feature, y = Number, fill = Regulation), stat = "identity", position = position_dodge(0.4), width = 0.5)

g + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1400, 500), minor_breaks = seq(0, 1300, 100), limits = c(0, 1400)) +
  scale_fill_manual(values = c("dodgerblue", "firebrick")) + 
  labs(fill = "Expression") +
  xlab("") + ylab("") +
  guides(fill = guide_legend(reverse=TRUE)) +
  theme(
    #aspect.ratio = 1/1,
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.position = c(0.75, 0.97),
    legend.key.size = unit(0.1, "in"),
    legend.direction = "horizontal",
    legend.title.align = 0,
    legend.background = element_rect(fill="transparent"),
    legend.key = element_blank(),
    axis.text.x = element_text(angle=45, hjust = 1, color = "black", size = 7),
    axis.text.y = element_text(color = "transparent", size = 7),
    panel.border = element_blank(),
    panel.grid.major.y = element_line(color = "black", size = 0.25),
    panel.grid.minor.y = element_line(color = "black", linetype = "dashed", size = 0.1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = 'white'),
    # Change axis line
    axis.line = element_line(color = "black")
  )
