library(dplyr)
library(tidyr)
library(ggplot2)
acc <- read.csv("mean_kcount_accuracy.csv")%>%
  select(contains("k"))

#sd() berekenen
#------------------------------------------
row_stdev <- apply(acc, 2, sd, na.rm=TRUE)
k7_stdev <- sd(k7_acc)
stdev <- data.frame(c(row_stdev[],k7_stdev))
colnames(stdev)<- "standaarddeviatie"
rownames(stdev) <- c("k1", "k2", "k3", "k4","k5", "k6","k7")
write.table(stdev, "standaarddeviatie")
#-------------------------------------------

k7_acc <- c(0.6811404,0.6697368,0.665309,0.6671053,0.6846491)
k7_gem <- mean(k7_acc)
gemiddelde_accuraatheid <- data.frame(gem_acc = c(colMeans(acc),k7_gem),
                                      kcount = c(1, 2,3,4,5,6,7))
rownames(gemiddelde_accuraatheid) <- c("k1", "k2", "k3", "k4","k5", "k6","k7")
#plot met accuracies
ggplot(gemiddelde_accuraatheid, aes(x = kcount, y = gem_acc)) +
  geom_point(colour = "darkblue") +
  geom_line(colour = "darkblue") +
  scale_x_continuous(breaks = seq(1, 7, by = 1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
  labs(title = "Gemiddelde accuraatheid ifv. lengte k-meren",
       x = "Lengte k-meer",
       y = "Gemiddelde accuraatheid")
#------------------------------------------------------------------------
#confidence interval
n <- 20
s <- 0.0119
xbar <- gemiddelde_accuraatheid$gem_acc[6]
margin <- qt(0.975,df=n-1)*s/sqrt(n)

lowerinterval <- xbar - margin
upperinterval <- xbar + margin
CI <- c(lowerinterval,upperinterval)
#-----------------------------------------------------------------------------
minimum <- min(test3$sequentielengte)
maximum <- max(test3$sequentielengte)

temp <- merge(x = juist,
              y = test3, 
              by = "Genus")
temp <- group_by(temp, Genus)%>%
  summarise(avg_lengte = mean(sequentielengte),
            avg_match = mean(matching))
cor(temp$avg_lengte,temp$avg_match)

pdata <- group_by(test3,Genus)%>%
  summarise(avg = mean(sequentielengte),
            sequentielengte = mean)
# Calculate overall mean sequentielengte
overall_mean <- mean(test3$sequentielengte)

minimum <- min(pdata$avg)
maximum <- max(pdata$avg)

# Create plot
sequentieplot <- ggplot(data = pdata, mapping = aes(x=Genus, y = avg)) +
  geom_col() +
  scale_colour_viridis_c() +
  theme(axis.text.x = element_blank()) +
  ylab("Gemiddelde sequentielengte") +
  labs(title = "Gemiddelde sequentielengte per genus") +
  # Add horizontal red line for overall mean sequentielengte
  geom_hline(yintercept = overall_mean, color = "red")

library(gridExtra)
grid.arrange(sequentieplot, plot_juist, ncol = 1)
#---------------------------------------------------------------------
pltdata <- mutate(juist, sequentielengte = pdata$avg)
my_colors <- c("#FDFFD0", "#FEB24C", "#FC4E2A", "#BD0026")

my_breaks <- c(0, 500, 1000, 2000, 3000)

ggplot(data = pltdata, aes(x = sequentielengte, y = matching)) +
  geom_point(colour = "darkblue") +
  labs(title = "Spreidingsplot van sequentielengte vs accuraatheid",
       x = "Gemiddelde sequentielengte van genus",
       y = "Procent juist geïdentificeerd")

ggplot(data = pltdata, aes(x = Frequentie, y = matching)) +
  geom_point(colour = "darkblue") +
  labs(title = "Spreidingsplot van aantal sequenties genus vs accuraatheid",
       x = "Aantal sequenties genus in dataset",
       y = "Procent juist geïdentificeerd")

#-------------------------------------------------------------------------------

test22 <- read.csv("sampled_data/test2.csv")
test22 <- test22[c(2, 3, 4, 5)]
plotmatrix <- mutate(plotmatrix, Genus = Prediction)
plotmatrix2 <- merge(x = test22,
                     y = plotmatrix, 
                     by = "Genus")

plotmatrix3 <- plotmatrix2%>%
  group_by(Familie, Genus, Reference)%>%
  summarise(across(c(normalized), mean))

colnames(plotmatrix3) <- c("Familie", "Prediction", "Reference", "normalized")

families <- data.frame(Familie = unique(plotmatrix3$Familie))

families <- families%>%
  mutate(volgorde = order(families$Familie))

plotmatrix4 <- merge(x = plotmatrix3, 
                     y = families,
                     by = "Familie")
colnames(plotmatrix4) <- c("Familie", "Prediction",
                           "Reference", "normalized", "volgorde_prediction")

familie_Prediction <- plotmatrix4[, 1:2]

colnames(familie_Prediction) <- c("Familie_reference", "Reference")

plotmatrix5 <- merge(x = plotmatrix4, 
                     y = familie_Prediction, 
                     by = "Reference")


plotmatrix6 <- plotmatrix5%>%
  group_by(Familie, Prediction, Reference, Familie_reference)%>%
  summarise(across(c(normalized, volgorde_prediction), mean))

families_reference <- data.frame(Familie_reference =
                                   unique(plotmatrix6$Familie_reference))
families_reference <- families_reference%>%
  mutate(volgorde_reference = order(Familie_reference))

plotmatrix6 <- plotmatrix6%>%
  arrange(Familie_reference)

plotmatrix6 <- merge(x = plotmatrix6, 
                     y = families_reference, 
                     by = "Familie_reference")

#plotten van confusionmatrix
ggplot(plotmatrix6, aes(x = reorder(Reference, volgorde_reference),
                        y = reorder(Prediction, volgorde_prediction),
                        fill = normalized)) + 
  geom_tile() + theme_bw() + coord_equal() +
  scale_fill_distiller(direction=1) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  ggtitle("Verwarringsmatrix")+
  xlab("Referentie") + ylab("Predictie")+
  labs(fill = "overeenkomst") +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank()  #remove y axis ticks
  )
#--------------------------------------------------------------
enterobacter <- filter(plotmatrix, Reference == "Enterobacter")%>%
  filter(normalized != 0)
salmonella <- filter(plotmatrix, Reference == "Salmonella")%>%
  filter(normalized != 0)
listeria <- filter(plotmatrix, Reference == "Listeria")%>%
  filter(normalized != 0)
merged_df <- as.data.frame(rbind(enterobacter, blautia, salmonella))
ggplot(salmonella, mapping = aes(x = Reference, y = Freq))+
  geom_col(fill = Prediction)+
  facet_wrap(~Reference)

#--------------------------------------------------------
#verdeling lengtes
ggplot(temp, aes(x = avg_lengte)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30,
                 colour = "darkblue",
                 fill = "lightblue")+
  labs(title = "Dichtheidshistogram van de gemiddelde sequentielengte per genus")+
  xlab("Gemiddelde sequentielengte")+
  ylab("Dichtheid")
#statistische analyse sequentielengte ~ match
temp_sorted <- temp[order(temp$avg_lengte),]

# Bepaal de lengte van de dataframe
n <- nrow(temp_sorted)

# Bereken de helft van de lengte
half <- floor(n/2)

# Maak de twee dataframes
temp_half1 <- temp_sorted[1:half,]
temp_half2 <- temp_sorted[(half+1):n,]




# Bekijk de resultaten
t_result

#controleren voorwaarde normaliteit
ggplot(temp_half1, aes(sample=avg_match)) +
  geom_qq() +
  geom_qq_line()
ggplot(temp_half2, aes(sample=avg_match)) +
  geom_qq() +
  geom_qq_line()

#de data is niet normaal verdeeld

#controleren voorwaarde homoskedasticiteit
temp_half1 <- mutate(temp_half1,lengte = "hoog")
temp_half2 <- mutate(temp_half2,lengte = "laag")

tempo <- rbind(temp_half1, temp_half2)
tempo$lengte <- factor(tempo$lengte)
var.test(avg_match ~ lengte, data = tempo)
#de data is heteroskedastisch
ggplot(tempo, mapping = aes(x = lengte, y = avg_match))+
  geom_boxplot()
# Voer de t-test uit
t_result <- wilcox.test(avg_match ~ lengte, data = tempo, alternative = "less")
t_result
#-------------------------------------------------------------------------
#--------------------------------------------------------
#statistische analyse frequentie ~ match
top30 <- juist %>% 
  top_n(30, Frequentie)

rest <- juist %>% 
  anti_join(top30, by = "Genus")

#controleren voorwaarde normaliteit
ggplot(top30, aes(sample = matching)) +
  geom_qq() +
  geom_qq_line()
ggplot(rest, aes(sample = matching)) +
  geom_qq() +
  geom_qq_line()


#testen homoskedasticiteit
top30 <- mutate(top30,lengte = "hoog")
rest <- mutate(rest,lengte = "laag")

tempo_freq <- rbind(top30, rest)
tempo_freq$lengte <- factor(tempo$lengte)

var.test(matching ~ lengte, data = tempo_freq)
#de data is heteroskedastisch
ggplot(tempo_freq, mapping = aes(x = lengte, y = matching))+
  geom_boxplot()
# Voer de t-test uit
t_result_freq <- wilcox.test(matching ~ lengte, data = tempo_freq, 
                             alternative = "less")
t_result_freq

# Bekijk de resultaten
t_result



#controleren voorwaarde homoskedasticiteit
temp_half1 <- mutate(temp_half1,lengte = "hoog")
temp_half2 <- mutate(temp_half2,lengte = "laag")

