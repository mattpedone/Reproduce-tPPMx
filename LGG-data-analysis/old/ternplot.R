rm(list=ls())
set.seed(121)
library(ggplot2)
library(ggtern)
library(compositions)
library(robCompositions)
library(dirmult)
library(scales)

hex2 <- hue_pal()(2)
hex3 <- hue_pal()(3)

# Treatment 1
eta <- c(0.1253, 0.4178, 0.4569)
reta <- data.frame(rdirichlet(n = 1600, exp(eta)*10))
colnames(reta) <- c("PD", "PS", "CR")
deta11 <- ggtern(data=reta, aes(x=PD, y=PS, z=CR)) +              
  geom_density_tern(colour = "#F8766D", bdl = 0.00000001) + 
  ggtitle("Group 1")
  
eta <- c(0.3241, 0.2252, 0.4508)
reta <- data.frame(rdirichlet(n = 1600, exp(eta)*10))
colnames(reta) <- c("PD", "PS", "CR")
deta12 <- ggtern(data=reta, aes(x=PD, y=PS, z=CR)) +              
  geom_density_tern(colour = "#00BFC4", bdl = 0.00000001) + 
  ggtitle("Group 2")

tern1 <- ggtern::grid.arrange(deta11, deta12, ncol = 2)#, common.legend = TRUE, legend="bottom")

# Treatment 2
eta <- c(-0.5378, -0.2301, -0.7604)
reta <- data.frame(rdirichlet(n = 1600, exp(eta)*10))
colnames(reta) <- c("PD", "PS", "CR")
deta21 <- ggtern(data=reta, aes(x=PD, y=PS, z=CR)) +              
  geom_density_tern(colour = "#F8766D", bdl = 0.00000001) + 
  ggtitle("Group 1")

eta <- c(-1.3823, -0.1178, -0.5438)
reta <- data.frame(rdirichlet(n = 1600, exp(eta)*10))
colnames(reta) <- c("PD", "PS", "CR")
deta22 <- ggtern(data=reta, aes(x=PD, y=PS, z=CR)) +              
  geom_density_tern(colour = "#00BA38", bdl = 0.00000001) + 
  ggtitle("Group 2")

eta <- c(-0.9260, 0.0485, -1.4872)
reta <- data.frame(rdirichlet(n = 1600, exp(eta)*10))
colnames(reta) <- c("PD", "PS", "CR")
deta25 <- ggtern(data=reta, aes(x=PD, y=PS, z=CR)) +              
  geom_density_tern(colour = "#619CFF", bdl = 0.00000001) + 
  ggtitle("Group 5")

tern2 <- ggtern::grid.arrange(deta21, deta22, deta25, ncol = 3)#, common.legend = TRUE, legend="bottom")

ggsave(tern1, device = "pdf", path = "figs", filename = "tern1.pdf")
ggsave(tern2, device = "pdf", path = "figs", filename = "tern2.pdf")
