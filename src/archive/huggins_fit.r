data <- read.csv("TC_Data_Charles.csv");
head(data);
Tab <- table(data[,7])
T.D <- table(data$surveydate)

D <- as.Date(data$surveydate,format="%d/%m/%Y")
d.1 <- min(D)
D.1 <- as.numeric(D-d.1)


F.T <- table(data$idfish)


m <- which.max(F.T)
n.m <- names(F.T)[m]
m.1 <- which(data$idfish==n.m)
TC.1 <- data[m.1,]

plot(D.1[m.1],data$weight[m.1],type="l");

plot(D.1[m.1],data$totallength[m.1],type="l");

plot(D.1[m.1],data$idsite[m.1]);

cbind(D[m.1],data$idsite[m.1]);

table(data$idsite[m.1])

plot(D.1[m.1],1:7,type="s")
