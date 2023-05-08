library(dplyr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(naivebayes)
library(kmer)
library(ape)
library(adegenet)
library(tidyselect)
library(e1071)
library(caret)
library(R.utils)

#inlezen dataset

test2 <- read.csv("test2.csv")
test2 <- test2[c(2, 3, 4, 5)]

#volledig <- gunzip("SILVA_138.1_LSUParc_tax_silva.fasta.gz")
database <- read.fasta(file = "sample_test_ 66234(5procent).fasta")
#matrix prealloceren
test <- data.frame(matrix(ncol = 3, nrow = length(database)))
#database van fasta formaat omzetten in leesbare dataframe
for (i in 1:length(database)) {
  new_elem <- getAnnot(database[[i]])
  split <- strsplit(new_elem, split = ";")[[1]]
  genus <- split[length(split)-1]
  if (length(split) > 2){
    familie <- split[length(split)-2]
  }
  domein <- split[1]
  string_spatie <- paste(database[[i]], collapse = " ")
  string <- gsub(" ", "", string_spatie)
  if (grepl("Bacteria", domein)){
    test[i, 1] <- string
    test[i, 2] <- genus
    test[i, 3] <- familie
  }
}



test <- na.omit(test)


colnames(test) <- c("Sequentie", "Genus", "Familie")

frequentie <- data.frame(table(test[, 2]))
colnames(frequentie) <- c("Genus", "Frequentie")

test <- right_join(test, frequentie, by='Genus')
#filteren zodat enkel genera voorkomen die meer dan 14 keer voorkomen 
test <- test%>%  
  filter(Frequentie > 14)%>%
  filter(Genus != "uncultured")




tmp <- data.frame()


library(stringr)
#sequenties met random letters eruit filteren
test2 <- test[str_detect(test$Sequentie, regex("^[AGCU]*$", ignore_case = TRUE)),]

write.csv(test2, "test2.csv")


uniek <- unique(test2$Genus)


#dataframe met accuraatheden van verschillende aantal sequenties maken

aantal_seq <- length(test2$Sequentie)
k <- c(6)
herhalingen <- 1
accuraatheid <- data.frame(accuraatheid = matrix(nrow = herhalingen,
                                                 ncol = length(k)))

#--------------------------------------------------------------------------
start_time <- Sys.time()
previous_time <- start_time

for (u in 1:length(k)){
  
  # initialiseer timer variabelen voor deze iteratie van 'u'
  iter_start_time <- Sys.time()
  iter_previous_time <- iter_start_time
  
  for(y in 1:herhalingen){
    set.seed(y)
    mat <- matrix(nrow = aantal_seq,ncol=4^k) #preallocatie matrix
    for ( i in 1:aantal_seq){   #matrix maken
      oef <- test2[i,1]
      var1 <- strsplit(oef, "")
      mat[i,] <- kcount(var1, k = k)/nchar(test2$Sequentie[i]) #genormaliseerde matrix
    }
    mat <- data.frame(mat)
    mattrain <- mat%>%
      mutate(genus = factor(test2$Genus[1:aantal_seq]))
    index <- createDataPartition(mattrain$genus, p = .9, list = FALSE) #stratisfier zodat soorten in verhoudingen worden gesampled
    
    train <- mattrain[index, ] #train dataset
    
    testdata <- mattrain[-index, ] #test dataset
    
    test_mat <- testdata[, 1:ncol(testdata)-1]
    test.label <- testdata[, ncol(testdata)] #labels ter controle van predicties
    
    model <- naiveBayes(train[, 1:ncol(train)-1],train[, ncol(train)]) #naive bayes model
    
    uitkomst <- predict(model, test_mat) #predictions
    
    
    accuraatheid[y, u] <- sum(uitkomst == test.label)/length(uitkomst)
    
    # bereken tijd voor deze iteratie van 'y'
    iter_current_time <- Sys.time()
    iter_elapsed_time <- iter_current_time - iter_previous_time
    iter_previous_time <- iter_current_time
    
    # print timer output voor deze iteratie van 'y'
    cat(sprintf("  Iteratie %d van y uitgevoerd in %s\n", y, format(iter_elapsed_time)))
    
  }
  
  # bereken tijd voor deze iteratie van 'u'
  current_time <- Sys.time()
  elapsed_time <- current_time - previous_time
  previous_time <- current_time
  
  # print timer output voor deze iteratie van 'u'
  cat(sprintf("Iteratie %d van u uitgevoerd in %s\n", u, format(elapsed_time)))
  
}

# bereken totale tijd
end_time <- Sys.time()
total_time <- end_time - start_time

# print timer output voor totale uitvoering van de code
cat(sprintf("\nTotaal uitgevoerd in %s\n", format(total_time)))
#------------------------------------------------------------------------

colnames(accuraatheid) <- c("k6")

gemiddelde_accuraatheid <- data.frame(gem_acc = colMeans(na.omit(accuraatheid)))

#write.csv(accuraatheid, "accuraatheid_k7.csv")

#write.csv(gemiddelde_accuraatheid, "mean_kcount7_accuracy.csv")
#-------------------------------------------------------------------------
#data confusionmatrix

k5_accuraatheid <- data.frame(k5_accuraatheid = matrix(nrow = herhalingen,
                                                    ncol = 1))
  for(y in 1:herhalingen){
    set.seed(herhalingen + y)
  
      mat <- matrix(nrow = aantal_seq,ncol=4^5) #preallocatie matrix
      for ( i in 1:aantal_seq){   #matrix maken
        oef <- test2[i,1]
        var1 <- strsplit(oef, "")
        mat[i,] <- kcount(var1, k = 5)/nchar(test2$Sequentie[i]) #genormaliseerde matrix
      }
      mat <- data.frame(mat)
      mattrain <- mat%>%
        mutate(genus = factor(test2$Genus[1:aantal_seq]))
      index <- createDataPartition(mattrain$genus, p = .9, list = FALSE) #stratisfier zodat soorten in verhoudingen worden gesampled
      
      train <- mattrain[index, ] #train dataset
      
      testdata <- mattrain[-index, ] #test dataset
      
      test_mat <- testdata[, 1:ncol(testdata)-1]
      test.label <- testdata[, ncol(testdata)] #labels ter controle van predicties
      
      model <- naiveBayes(train[, 1:ncol(train)-1],train[, ncol(train)]) #naive bayes model
      
      uitkomst <- predict(model, test_mat) #predictions
      
      
      k5_accuraatheid[y, 1] <- sum(uitkomst == test.label)/length(uitkomst)

  }

pos_mediaan <- which(k5_accuraatheid$k5_accuraatheid == median(k5_accuraatheid$k5_accuraatheid))[1]

write.csv(k5_accuraatheid, "k5_accuraatheid.csv")




set.seed(herhalingen + pos_mediaan)

mat <- matrix(nrow = aantal_seq,ncol=4^5) #preallocatie matrix
for ( i in 1:aantal_seq){   #matrix maken
  oef <- test2[i,1]
  var1 <- strsplit(oef, "")
  mat[i,] <- kcount(var1, k = u)/nchar(test2$Sequentie[i]) #genormaliseerde matrix
}
mat <- data.frame(mat)
mattrain <- mat%>%
  mutate(genus = factor(test2$Genus[1:aantal_seq]))
index <- createDataPartition(mattrain$genus, p = .9, list = FALSE) #stratisfier zodat soorten in verhoudingen worden gesampled

train <- mattrain[index, ] #train dataset

testdata <- mattrain[-index, ] #test dataset

test_mat <- testdata[, 1:ncol(testdata)-1]
test.label <- testdata[, ncol(testdata)] #labels ter controle van predicties

model <- naiveBayes(train[, 1:ncol(train)-1],train[, ncol(train)]) #naive bayes model

uitkomst <- predict(model, test_mat) #predictions
#-------------------------------------------------------------------------
#plotten van accuraatheid tov kcount
ggplot(accuraatheid, aes(x = kcount, y = accuraatheid)) +
  geom_line() + ylim(0, 1)


#set.seed(1234)

test3 <- test2%>%
  mutate(sequentielengte = nchar(test2$Sequentie))


#confusion matrix

cm <- confusionMatrix(data= uitkomst,
                      reference = test.label)


a <- test2[1:3000,]
uniek <- unique(a$Genus)

labels <- data.frame(labels = matrix(nrow = length(test2$Genus), ncol = 1))
#genera selecteren die meer dan 250 keer voorkomen om ze dan te highlighten
for (i in 1:length(test2$Genus)){
  if (test2$Frequentie[i] >= 250){
    labels[i, ] <- test2$Genus[i]
  }
}

#frequentieplot maken 

gekleurd <- unique(na.omit(labels$labels))

plot <- cbind(test2, labels)

plot <- plot%>%
  arrange(desc(Frequentie))

ggplot(plot, aes(x = Genus, fill = labels)) +
  geom_bar(stat = "count")+
  ggtitle("Aantal per genus (train -en testdata)") + 
  theme(axis.text.x=element_blank(), #remove x axis labels
  ) + 
  ylab("Aantal")+
  scale_fill_discrete(labels = c(gekleurd, "Andere genera"))









plotmatrix <- data.frame(cm$table)





colnames(genus_familie) <- c("Genus", "Familie")


genus_familie <- genus_familie%>%
  group_by(Genus, Familie)%>%
  summarise(.groups = "drop")

colnames(genus_familie) <- c("Reference", "Familie")



norm <- plotmatrix %>% group_by(Reference) %>% 
  summarise(som=sum(Freq),
            .groups = 'drop')

norm <- right_join(norm, genus_familie, by = "Reference")



plotmatrix <- right_join(plotmatrix, norm, by = "Reference")

plotmatrix <- right_join(plotmatrix, aantal_families, by = "Familie")

plotmatrix <- plotmatrix%>%    #normalizeren van frequenties zodat confusion matrix duidelijker af te lezen is
  mutate(normalized = Freq/som)

aantal_families <- data.frame(table(plotmatrix$Familie)/length(unique(test2$Genus[1:aantal_seq])))
colnames(aantal_families) <- c("Familie", "Aantal")


plotmatrix$normalized[is.nan(plotmatrix$normalized)] <- 0

colnames(plotmatrix) <- c("Genus", "Reference", "Freq", "som", "normalized")

plotmatrix2 <- merge(x = test2,
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
families_reference <- data.frame(Familie_reference = arrange(families_reference, Familie_reference))

families_reference <- families_reference%>%
  mutate(volgorde_reference = order(Familie_reference))

plotmatrix6 <- plotmatrix6%>%
  arrange(Familie_reference)

plotmatrix6 <- merge(x = plotmatrix6, 
                     y = families_reference, 
                     by = "Familie_reference")

write.csv(plotmatrix6, "plotmatrix6.csv")


#plotten van confusionmatrix
ggplot(plotmatrix6, aes(x = reorder(Reference, volgorde_reference),
                        y = reorder(Prediction, volgorde_prediction),
                        fill = normalized)) + 
  geom_tile() + theme_bw() + coord_equal() +
  scale_fill_distiller(direction=1) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  ggtitle("Verwarringsmatrix")+
  xlab("Genera in testdataset") + ylab("Voorspelde genera")+
  labs(fill = "overeenkomst") +
   theme(axis.text.x=element_blank(), #remove x axis labels
         axis.ticks.x=element_blank(), #remove x axis ticks
         axis.text.y=element_blank(),  #remove y axis labels
         axis.ticks.y=element_blank()  #remove y axis ticks
 )


jpeg("C:/Users/bogae/OneDrive/Documenten/3e bach/Bachelorproef/BAP/Confusionmatrix.jpg",
     quality = 100)


ggplot(plotmatrix, aes(x = Reference,
                        y = Genus,
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


# kader <- data.frame(matrix(nrow = length(plotmatrix$Prediction),
#                            ncol = length(plotmatrix)))
# 
# for (i in 1:length(plotmatrix$Prediction)){
#   if(plotmatrix$Aantal[i] >= 3){
#     kader[i, ] <- plotmatrix[i, ]
#   }
# }
# 
# colnames(kader) <- colnames(plotmatrix)
# posities <- which(!is.na(kader$Familie))



#distributieplot train -en testdata

genera_traindata <- train[ncol(train)]

genera_traindata <- data.frame(Genus = as.character(genera_traindata$genus))
genera_traindata <- data.frame(table(genera_traindata$Genus))
colnames(genera_traindata) <- c("Genus","Frequentie")
highlight <- data.frame(highlight = matrix(ncol = 1,
                                           nrow = length(genera_traindata$Genus)))
for (i in 1:length(genera_traindata$Genus)){
  if (genera_traindata$Frequentie[i] >= 100){
    highlight[i, ] <- as.character(genera_traindata$Genus[i])
  }
}

genera_traindata <- cbind(genera_traindata, highlight)

gekleurd_train <- unique(na.omit(genera_traindata$highlight))

#plot traindata

ggplot(genera_traindata, aes(x = Genus, y = Frequentie, fill = highlight)) +
  geom_bar(stat = "identity")+ 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  ggtitle("Frequentieverdeling trainingsdata") + 
  theme(axis.text.x=element_blank()) + 
  labs(fill = "Meest voorkomende genera") + 
  scale_fill_discrete(labels = c(gekleurd_train, "Andere genera"))


genera_testdata <- testdata[ncol(testdata)]

genera_testdata <- data.frame(Genus = as.character(genera_testdata$genus))
genera_testdata <- data.frame(table(genera_testdata$Genus))
colnames(genera_testdata) <- c("Genus","Frequentie")
highlight_test <- data.frame(highlight = matrix(ncol = 1,
                                                nrow = length(genera_testdata$Genus)))
for (i in 1:length(genera_testdata$Genus)){
  if (genera_testdata$Frequentie[i] >= 10){
    highlight_test[i, ] <- as.character(genera_testdata$Genus[i])
  }
}

genera_testdata <- cbind(genera_testdata, highlight_test)

gekleurd_test <- unique(na.omit(genera_testdata$highlight))

#plot testdata

ggplot(genera_testdata, aes(x = Genus, y = Frequentie, fill = highlight)) +
  geom_bar(stat = "identity")+ 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ggtitle("Frequentieverdeling testdata") + 
  labs(fill = "Meest voorkomende genera") + 
  scale_fill_discrete(labels = c(gekleurd_test, "Andere genera"))

#tabel van meest voorkomende misparen

foutief <- data.frame(fout = matrix(nrow = length(plotmatrix$Prediction), 
                                    ncol = 3))

for (i in 1:length(plotmatrix$Prediction)){
  if (plotmatrix$Prediction[i] != plotmatrix$Reference[i] && plotmatrix$normalized[i] != 0){
    foutief[i, 1] <- plotmatrix$normalized[i]
    foutief[i, 2] <- as.character(plotmatrix$Prediction[i])
    foutief[i, 3] <- as.character(plotmatrix$Reference[i])
  }
  
}

foutief <- na.omit(foutief)
colnames(foutief) <- c("mismatching","Prediction", "Reference")

foutief <- foutief%>%
  arrange(desc(mismatching))

library(gt)

gt(foutief)


frequentieplot <- ggplot(plot, aes(x = reorder(Genus, -Frequentie), fill = labels)) +
  geom_bar(stat = "count")+
  ggtitle("Aantal per genus (train -en testdata)") + 
  theme(axis.text.x=element_blank(), #remove x axis labels
  ) + 
  ylab("Aantal")+
  scale_fill_discrete(labels = c(gekleurd, "Andere genera"))






#plot van juiste matching
juist <- data.frame(juist = matrix(nrow = length(plotmatrix$Prediction), 
                                   ncol = 2))

for (i in 1:length(plotmatrix$Prediction)){
  if (plotmatrix$Prediction[i] == plotmatrix$Reference[i]){
    juist[i, 1] <- plotmatrix$normalized[i]
    juist[i, 2] <- as.character(plotmatrix$Prediction[i])
  }
  
}
juist <- na.omit(juist)
colnames(juist) <- c("matching","Genus")

juist <- merge(x = test2,
                    y = juist, 
                    by = "Genus")



juist <- juist%>%
  group_by(Genus)%>%
  summarise(across(c(Frequentie, matching), mean))

plot_juist <- ggplot(juist, aes(x = reorder(Genus, -Frequentie), y = matching)) + 
  geom_bar(stat = "identity") +
  theme_light()+
  scale_fill_distiller(direction=1) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ggtitle("Percentage van juiste toekenningen")

# library(gridExtra)
# 
# grid.arrange(frequentieplot, plot_juist, ncol=1)



write.csv(foutief, "foutief.csv")


k <- c("1","2","3","4","5","6","7")

tijd <- c(4.9, 5.1, 5.4, 5.9, 8.6, 19.6, 193)

k_tijd <- data.frame(k = k, minuten = tijd)

colnames(k_tijd) <- c("lengte k-meer", "tijd (min.)")

library(gt)

k_tijd_plot <- gt(k_tijd)

tab_options(k_tijd_plot, 
  table.width = px(300),
  table.margin.left = px(0)
  
)
