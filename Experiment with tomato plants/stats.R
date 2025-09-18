library(multcompView)
#Parse data
csv<-read.csv('data/WPFMAllReps.csv')
csv$Treatment<- paste(csv$Condition, csv$Genotype, sep = '_')

print(csv)

#Linear model for ANOVA testing

print('Testing for normal distribution of the data...')

print('-- non-transformed data --')
LM<-lm(Freshweight~Treatment,data=csv)
SW<-shapiro.test(residuals(LM))
print(SW)

print('-- log-transformed data --')
LM<-lm(log(Freshweight)~Treatment,data=csv)
SW<-shapiro.test(residuals(LM))
print(SW)

print('-- sqrt-transformed data --')
LM<-lm(sqrt(Freshweight)~Treatment,data=csv)
SW<-shapiro.test(residuals(LM))
print(SW)


#ANOVA
A<-aov(LM)
print(summary(A))

#TukeyHSD
T<-TukeyHSD(A)
print(T)

#Letters
L<-multcompLetters4(A, T)
print(L)
