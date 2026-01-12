##packages 
install.packages("readxl")
install.packages("RColorBrewer")
install.packages("ggplot2")
install.packages("reshape2")
install.packages("extrafont")
install.packages("patchwork")
install.packages("gridExtra")
install.packages("dplyr")

##libraries
library(readxl)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(extrafont)
library(patchwork)
library(gridExtra)
library(dplyr)

#clean the environment
rm(list = ls())


## import data from excel file
data <- read_excel("Growth_reduction.xlsx", 
                   sheet = "R", col_types = c("text", "numeric", 
                                              "numeric","numeric","numeric", "numeric"))


## prepare data matrix
# set effector names as row names
effector_names <- data[[1]]         # extracts the effector names (first column)
data_matrix <- as.matrix(data[,-1]) # converts the rest of the data to a matrix
rownames(data_matrix) <- effector_names # sets the row names 

# Replace values below 0 with 0
data_matrix <- pmax(pmin(data_matrix, 1), 0)

print(data_matrix) # check matrix format

# Melt the data matrix to a long format for ggplot2
data_melted <- melt(data_matrix)
colnames(data_melted) <- c("Row", "Col", "Value")



## Change order of bacterial strains according to their sensitivity
# Summarize data: calculate the mean value for each row (effector)
row_means <- data_melted %>%
  group_by(Row) %>%
  summarize(mean_value = mean(Value, na.rm = TRUE)) # Calculate the mean inhibition percentage per row

# Order rows by the calculated mean_value
ordered_rows <- row_means %>%
  arrange(mean_value) %>%   # Change to arrange(desc(mean_value)) for descending order
  pull(Row)

# Reverse the order of the rows
reversed_rows <- rev(ordered_rows)

# Set the order of the 'Row' factor in the melted data according to the new order
data_melted$Row <- factor(data_melted$Row, levels = reversed_rows)



## Change order of effector proteins
# Reverse the order of the columns as well
ordered_cols <- unique(data_melted$Col)  
reversed_cols <- rev(ordered_cols)       
data_melted$Col <- factor(data_melted$Col, levels = reversed_cols)



# heatmap
ggplot(data_melted, aes(x = Row, y = Col, fill = Value)) +
  geom_tile() +
  scale_fill_gradientn(colors = brewer.pal(9, "OrRd"), na.value = "lightgrey") +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 12, color = "#333A3E", face = "italic"),
        axis.text.y = element_text(size = 12, color = "#333A3E"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",          
        legend.direction = "horizontal",      
        legend.box = "horizontal",
        legend.title = element_blank(),  
        legend.text = element_text(size = 12, color = "#333A3E"),  
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(1, "cm")
  )+
  coord_fixed()


# Save the plot as a PDF file
ggsave("BacHeatmap.pdf", device = "pdf", width = 8.27, height = 4.4)



##palettes of color brewer
display.brewer.all()
