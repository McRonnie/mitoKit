alldata <- read.delim( file = "/home/jingyu/Desktop/mito-chisq", stringsAsFactors = FALSE, header = FALSE )
library( ggplot2 )
colnames( alldata ) <- c( "pos", "pvalue", "or" )

library(readxl)
data <- read_excel("~/Documents/data.xlsx")
colnames( data ) = c( "pos", "key" )
View(data)


g <- ggplot( data, aes( x = pos, y = key, color = key ) )
g + geom_point() + 
  xlim( 0,16570 ) +
  ylim( -0.5, 3 ) +
  theme(
    panel.spacing = unit (0.1, "lines"),
    panel.background = element_rect(fill = "grey96")) +
  labs(x = "", y = "-log(P-value)" ) +
  coord_polar( theta = "x", direction = -1 )
#scale_y_continuous(expand = c(0,0)) +
#scale_x_continuous(expand = c(0,0))
