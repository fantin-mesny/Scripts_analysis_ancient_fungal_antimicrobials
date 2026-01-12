## install packages and load libraries

install.packages("readxl")
install.packages("RColorBrewer")
install.packages("ggplot2")
install.packages("reshape2")
install.packages("extrafont")
install.packages("patchwork")
install.packages("gridExtra")
install.packages("dplyr")

library(readxl)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(extrafont)
library(patchwork)
library(gridExtra)
library(dplyr)



# Import and load fonts
font_import()
loadfonts(device = "win")

# List available fonts
fonts()


## import data from excel file
data <- read_excel("Biomass_quan_ImageJ.xlsx", 
                   sheet = "Fil", col_types = c("text", "numeric", "numeric", 
                                              "numeric"))
data2 <- read_excel("Biomass_quan_ImageJ.xlsx", 
                   sheet = "yea", col_types = c("text", "numeric", "numeric", 
                                                "numeric"))

## prepare data matrix
# set effector names as row names
effector_names <- data[[1]]         # extracts the effector names (first column)
data_matrix <- as.matrix(data[,-1]) # converts the rest of the data to a matrix
rownames(data_matrix) <- effector_names # sets the row names 

effector_names2 <- data2[[1]] 
data_matrix2 <- as.matrix(data2[,-1])
rownames(data_matrix2) <- effector_names2

# Replace values below 0 with 0
data_matrix <- pmax(pmin(data_matrix, 1), 0)
data_matrix2 <- pmax(pmin(data_matrix2, 1), 0)

print(data_matrix) # check matrix format
print(data_matrix2)

# Melt the data matrix to a long format for ggplot2
data_melted <- melt(data_matrix)
colnames(data_melted) <- c("Row", "Col", "Value")

data_melted2 <- melt(data_matrix2)
colnames(data_melted2) <- c("Row", "Col", "Value")

# Summarize data: calculate the mean value for each row (effector)
row_means <- data_melted %>%
  group_by(Row) %>%
  summarize(mean_value = mean(Value, na.rm = TRUE)) # Calculate the mean inhibition percentage per row

row_means2 <- data_melted2 %>%
  group_by(Row) %>%
  summarize(mean_value2 = mean(Value, na.rm = TRUE))

# Order rows by the calculated mean_value
ordered_rows <- row_means %>%
  arrange(mean_value) %>%   # Change to arrange(desc(mean_value)) for descending order
  pull(Row)

ordered_rows2 <- row_means2 %>%
  arrange(mean_value2) %>%   
  pull(Row)

# Reverse the order of the rows
reversed_rows <- rev(ordered_rows)
reversed_rows2 <- rev(ordered_rows2)

# Set the order of the 'Row' factor in the melted data according to the new order
data_melted$Row <- factor(data_melted$Row, levels = reversed_rows)
data_melted2$Row <- factor(data_melted2$Row, levels = reversed_rows2)

# Reverse the order of the columns as well
ordered_cols <- unique(data_melted$Col)  
reversed_cols <- rev(ordered_cols)       

# Set the order of the 'Col' factor in the melted data according to the new order
data_melted$Col <- factor(data_melted$Col, levels = reversed_cols)
data_melted2$Col <- factor(data_melted2$Col, levels = reversed_cols)

# heatmap
heatmap1 <- ggplot(data_melted, aes(x = Row, y = Col, fill = Value)) +
  geom_tile() +
  scale_fill_gradientn(colors = brewer.pal(9, "OrRd"), na.value = "lightgrey", limits = c(0, 1)) +
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

print(heatmap1)

heatmap2 <- ggplot(data_melted2, aes(x = Row, y = Col, fill = Value)) +
  geom_tile() +
  scale_fill_gradientn(colors = brewer.pal(9, "OrRd"), na.value = "lightgrey", limits = c(0, 1)) +
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

print(heatmap2)

display.brewer.all()

# Save the plots as a PDF file
ggsave("FilamentousFungi.pdf", device = "pdf", plot = heatmap1, width = 4, height = 3)
ggsave("Yeasts.pdf", device = "pdf", plot = heatmap2, width = 4, height = 3)


