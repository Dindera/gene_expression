#import data
genes = read.csv("C:/Users/chide/OneDrive/Desktop/Statistical Course/gene_dataa.csv")

#----------------Task 1-------------------------

#------------ check for summary of data ---------#
summary(genes)



#------ create data matrix ------#

df = data.matrix(genes)

#--------- select each columns ---------#


df
x0 = df[, 1]
x1 = df[, 2]
x2 = df[, 3]
x3 = df[, 4]
x4 = df[, 5]
x5 = df[, 6]

#sumDf = sum(x1,x2,x3,x4,x5)
#sampleMean = sumDf/301
#sampleMean

#------- Time series plots ----------#

plot.ts(genes, main = "Time series plot")

#------ Distribution for each gene --------#


distribution = function(gene, gname = ""){
  gene_mean = mean(gene)
  gene_sd = sd(gene)
  x_lab = gname
  
  print(gene_mean)
  print(gene_sd)
  
#histogram and normal curve plot
  h = hist(gene, breaks = 10, col = "darkblue", xlab = x_lab , ylab = "Density", main = "Histogram with Normal Curve") 
  xfit = seq(min(gene), max(gene), length = 70) 
  yfit = dnorm(xfit, mean = gene_mean, sd = gene_sd) 
  yfit = yfit * diff(h$mids[1:2]) * length(gene) 
  
  lines(xfit, yfit, col = "red", lwd = 2)
}

#----- distribution plots --------#
distribution(genes$x1, "Gene x1")
distribution(genes$x2, "Gene x2")
distribution(genes$x3, "Gene x3")
distribution(genes$x4, "Gene x4")
distribution(genes$x5, "Gene x5")


#Correlation of genes
# library(corrplot)

cor(df,  method = "pearson")
cor.test(x3,x5, method = "pearson")
cor.test(x3,x4, method = "pearson")
#corrplot(cor(genes), type="upper", method="color", addCoef.col = "black",tl.col = "black", number.cex = 0.6,mar=c(0,0,0,0))

#Scatter Plots of genes
plot(genes[2:6], col="darkred", main="Genes Scatter Plot")

#----------------- Linear regression of genes -------------------# 

dff = as.matrix(data.frame(x1, x2, x4, x5))
y = x3
yH = list()

  for(c in 1:length(dff)){
        n = length(dff[, c])
        ones = matrix(1 , n , 1)
        
        x = dff[, c]
        
        X = cbind(ones, x)
        thetaHat = solve(t(X) %*% X) %*% t(X) %*% y
          
        yhat = X %*% thetaHat
        yHatm = data.frame(yhat)
       
        
        errr = y - yhat
        yH[[c]] = yHatm
        
        sse = norm(errr , type = "2")^2
        print(sse)
  }

yHat = do.call(cbind, yH)

par(mfrow=c(1,3))

plot(dff[,1],y, xlab = "x1", ylab = "x3")
lines(dff[,1],yHat[[1]], col="red")

plot(dff[,2],y, xlab = "x2", ylab = "x3")
lines(dff[,2],yHat[[2]], col="red")

plot(dff[,3],y, xlab = "x4", ylab = "x3")
lines(dff[,3],yHat[[3]], col="red")

plot(dff[,4],y, xlab = "x5", ylab = "x3")
lines(dff[,4],yHat[[4]], col="red")



#---------------- Task 2 ----------------------------#


X = t(as.matrix(genes[2:6])) 

par(mfrow=c(2,3)) # plot the subsequent graphs in matrix of 2x3. Try par(mfcol=c(2,3)) and see the difference

plot( X[2,] , X[1,] , main="1,2")
plot( X[3,] , X[1,] , main="1,3")
plot( X[4,] , X[1,] , main="1,4")
plot( X[5,] , X[1,] , main="1,5")

plot( X[3,] , X[2,] , main="2,3")
plot( X[4,] , X[2,] , main="2,4")
plot( X[5,] , X[2,] , main="2,4")

plot( X[4,] , X[3,] , main="4,3")
plot( X[5,] , X[3,] , main="5,3")

plot( X[5,] , X[4,] , main="5,4")

#----- Eigen decomposition of genes --------------#
X_mean = matrix(rowMeans( X ) , 5 , 301 ) # Mean matrix of rowise means
X_without_mean = X - X_mean # remove the means of genes
A = X_without_mean %*% t(X_without_mean)

Cx = (1/301) * (A)

E = eigen(Cx) 

#------------- PCA ----------------#

P = t(E$vectors)

Y = P %*% X_without_mean # Y is the new reduced data
P
sdev = sqrt(diag((1/(dim(X_without_mean)[2]-1)* P %*% A  %*% t(P)))) # standard dev of X without mean

variance = sdev^2


prop_of_var = variance/sum(variance)

prop_of_var
print(sum(prop_of_var))
colo = c("red", "purple", "green", "blue", "yellow")

par(mfrow=c(1,3))
plot(Y[1,], Y[2,], xlab = "PC1", ylab = "PC2", col=colo, main = "PCA plot", pch=19)
plot(prop_of_var,  type = "b")
plot(cumsum(prop_of_var),  type = "b")


#---------- Task 3 --------------

#------------ split data ---------------

df = as.matrix(genes[4:6])

sample_size = floor(0.80 * nrow(df))

set.seed(123)
trainTest = sample(seq_len(nrow(df)), size = sample_size)
trainData = df[trainTest, ]
testData = df[-trainTest, ]

xtrain = trainData[,2:3]
xtest = testData[,2:3]
ytrain = trainData[,1]
ytest = testData[,1]

#-------------- model -----------------#

x4_tr =  matrix(xtrain[,1])
x4_te =  matrix(xtest[,1])
x5_tr =  matrix(xtrain[,2])
x5_te =  matrix(xtest[,2])
ones = matrix(1 , length(x4_tr), 1)
onesT = matrix(1, length(x4_te), 1)

#------------ selection 1 ---------------#

#model terms = {ones, x4, x4^2,x4^3, x5, x5^2,  x5^3}
model_terms = data.frame(one = ones, x4 = x4_tr,  x4_2 = x4_tr^2, x4_3 = x4_tr^3, x5 = x5_tr, x5_2 = x5_tr^2, x5_3 = x5_tr^3)
model_terms2 = data.frame(ones = onesT, x4 = x4_te, x4_2 = x4_te^2, x4_3 = x4_te^3, x5 = x5_te, x5_2 = x5_te^2,  x5_3 = x5_te^3)

for (c in 1:length(model_terms)){
  X = cbind(model_terms[,c])
  y = ytrain
  
  thetaHat = solve(t(X) %*% X) %*% t(X) %*% y
}

for (c in 1:length(model_terms2)){
  y = ytest
  x = cbind(model_terms2[,c])
    
  yHat = x %*% thetaHat
  err = y - yHat
  MSE = t(err)%*%(err)
  print(MSE)
}

# smallest mse for model selection 1 = x5^3

#---------------- selection 2 ----------------#

model_terms = data.frame(one = ones, x4 = x4_tr, x4_2 = x4_tr^2,  x4_3 = x4_tr^3, x5 = x5_tr, x5_2 = x5_tr^2)
model_terms2 = data.frame(ones = onesT, x4 = x4_te, x4_2 = x4_te^2, x4_3 = x4_te^3, x5 = x5_te, x5_2 = x5_te^2)

selectTerm = x5_tr^3
selectTerm2 = x5_te^3

for (c in 1:length(model_terms)){
  X = cbind(selectTerm, model_terms[,c])
  y = ytrain
  
  thetaHat = solve(t(X) %*% X) %*% t(X) %*% y
}

for (c in 1:length(model_terms2)){
  y = ytest
  x = cbind(selectTerm2 ,model_terms2[,c])
  
  yHat = x %*% thetaHat
  err = y - yHat
  MSE = t(err)%*%(err)
  print(MSE)
}

# smallest mse for model selection 2 = x5^2

#---------------- selection 3 ----------------#

model_terms = data.frame(one = ones, x4 = x4_tr,  x4_2 = x4_tr^2, x4_3 = x4_tr^3, x5 = x5_tr)
model_terms2 = data.frame(ones = onesT, x4 = x4_te, x4_2 = x4_te^2, x4_3 = x4_te^3, x5 = x5_te)

x5_2 = x5_tr^2
x5_2t = x5_te^2

x5_3 = x5_tr^3
x5_3t = x5_te^3

for (c in 1:length(model_terms)){
  X = cbind(x5_3 , x5_2, model_terms[,c])
  y = ytrain
  
  thetaHat = solve(t(X) %*% X) %*% t(X) %*% y
}

for (c in 1:length(model_terms2)){
  y = ytest
  x = cbind(x5_3t, x5_2t, model_terms2[,c])
  
  yHat = x %*% thetaHat
  err = y - yHat
  MSE = t(err)%*%(err)
  print(MSE)
}

# smallest mse for model selection 3 = x5

#---------------- selection 4 ----------------#

model_terms = data.frame(one = ones, x4 = x4_tr,  x4_2 = x4_tr^2, x4_3 = x4_tr^3)
model_terms2 = data.frame(ones = onesT, x4 = x4_te, x4_2 = x4_te^2, x4_3 = x4_te^3)

x5_3 = x5_tr^3
x5_3t = x5_te^3

x5_2 = x5_tr^2
x5_2t = x5_te^2

x5 = x5_tr
x5t = x5_te

for (c in 1:length(model_terms)){
  X = cbind(x5_3, x5_2 , x5, model_terms[,c])
  y = ytrain
  
  thetaHat = solve(t(X) %*% X) %*% t(X) %*% y
}

for (c in 1:length(model_terms2)){
  y = ytest
  x = cbind(x5_3t, x5_2t, x5t, model_terms2[,c])
  
  yHat = x %*% thetaHat
  err = y - yHat
  MSE = t(err)%*%(err)
  print(MSE)
}

# smallest mse for model selection 4 = x4^3


# final model x3 =  x5 + x5^2 + x4^3 + x5^3

#--------------Pairwise Plots --------------------------------#
x5 = x5_tr
x5t = x5_te

x4_3 = x4_tr^3
x4_3t = x4_te^3

x5_2 = x5_tr^2
x5_2t = x5_te^2

X5_3 = x5_tr^3
X5_3t = x5_te^3

X1 = cbind(x5, x4_3)
X2 = cbind(x5, x5_2)
X3 = cbind(x5, x5_3)
X4 = cbind(x5_2, x4_3)
X5 = cbind(x4_3, x5_3)

X1t = cbind(x5t, x4_3t)
X2t = cbind(x5t, x5_2t)
X3t = cbind(x5t, x5_3t)
X4t = cbind(x5_2t, x4_3t)
X5t = cbind(x4_3t, x5_3t)



pairwise = function(X, Xt, theta_1_range, theta_2_range) {
  yy = ytrain
  y = ytest
  n = length(X[,1])
  A = t(X) %*% X
  
  thetaHat = solve(A) %*% t(X) %*% yy
  print(thetaHat)
  
  yHat = Xt %*% thetaHat # model predicted
  err = y - yHat
  sse = norm(err , type = "2")^2
  
 
  
  sigma_2 = sse/( n - 1 ) # error variance sigma^2
  
  cov_thetaHat = 1 * (solve(A)) 
  cov_thetaHat_inv = (A) * (1/sigma_2) # inverse of cov_thetaHat
  det_cov_thetaHat = det(cov_thetaHat) # determinent of cov_thetaHat

  
  no_points = 20 # no point on the plot
  number_of_parameters = 2
  
  p_thetaHat_D = matrix(0 , no_points , no_points)
  
  for(r in 1:20){
    for(c in  1:20){
      
      theta_12 = matrix( c( theta_1_range[r] , theta_2_range[c] ) , number_of_parameters , 1)
      thetaHat_theta = theta_12 - thetaHat
      
      p_thetaHat_D[r,c] = ( 1/sqrt( ( (2*pi)^number_of_parameters ) * det_cov_thetaHat) ) * 
        exp( -0.5 * t(-thetaHat_theta) %*% cov_thetaHat_inv %*% -thetaHat_theta )
      
        
    }
  }
  par(mfrow = c(1,1))
  #contour(theta_1_range, theta_2_range, p_thetaHat_D, main="Contour of pairwise parameters", xlab = "X5", ylab="X4^3") #Uncomment to see contour
  persp(theta_1_range, theta_2_range, p_thetaHat_D , theta = 30 , phi = 15, xlab = "Theta 1 range", ylab = "Theta 2 range", zlab = "P.d.f of parameter", main = "3D of pairwise parameters", col = "springgreen", shade = 0.5) # theta changes the rotation (left-right), phi changes the rotation (up-down)
  
}

theta_X1 = seq(0.24 , 0.27 , length=20)
theta_X1_2 = seq(0.21 , 0.22 , length=20)

theta_X2 = seq(-0.45 , -0.35 , length=20)
theta_X2_2 = seq(0.55 , 0.7 , length=20)

theta_X3 = seq(0.09 , 0.13 , length=20)
theta_X3_2 = seq(0.15 , 0.25 , length=20)

theta_X4 = seq(0.30 , 0.45 , length=20)
theta_X4_2 = seq(0.01 , 0.025 , length=20)

theta_X5 = seq(-0.15 , -0.08  , length=20)
theta_X5_2 = seq(0.25 , 0.35 , length=20)

help("persp")

pairwise(X1, X1t, theta_X1, theta_X1_2)
pairwise(X2, X2t, theta_X2, theta_X2_2)
pairwise(X3, X3t, theta_X3, theta_X3_2)
pairwise(X4, X4t, theta_X4, theta_X4_2)
pairwise(X5, X5t, theta_X5, theta_X5_2)

#---------------- Model Validation ------------------------------#

df = as.matrix(genes[4:6])

sample_size = floor(0.85 * nrow(df))

set.seed(123)
trainTest = sample(seq_len(nrow(df)), size = sample_size)
trainData = df[trainTest, ]
testData = df[-trainTest, ]

xtrain = trainData[,2:3]
xtest = testData[,2:3]
ytrain = trainData[,1]
ytest = testData[,1]

x4_tr =  matrix(xtrain[,1])
x4_te =  matrix(xtest[,1])
x5_tr =  matrix(xtrain[,2])
x5_te =  matrix(xtest[,2])


model_terms = data.frame(x4_3 = x4_tr^3, x5 = x5_tr, x5_2 = x5_tr^2, x5_3 = x5_tr^3)
model_terms2 = data.frame(x4_3 = x4_te^3, x5 = x5_te, x5_2 = x5_te^2,  x5_3 = x5_te^3)

Xx = cbind(model_terms[,2], model_terms[,1], model_terms[,3], model_terms[,4])
yy = ytrain
y = ytest
A = t(Xx) %*% Xx

thetaHat = solve(A) %*% t(Xx) %*% yy
thetaHat

xx = cbind(model_terms2[,2], model_terms2[,1], model_terms2[,3], model_terms2[,4])

yHat = xx %*% thetaHat # model predicted
err = y - yHat
MSE = t(err)%*%(err)
print(MSE)

par(mfrow = c(1,1))
plot(yHat, err)


#------------------ CI ---------------------------#

cov_thetaHat = 1 * (solve(t(xx) %*% xx))
n = 46
number_of_parameters = 4

var_y_hat = matrix(0 , n , 1)

for( i in 1:n){
  X_i = matrix( xx[i,] , 1 , number_of_parameters ) # X[i,] creates a vector. Convert it to matrix

  var_y_hat[i,1] = X_i %*% cov_thetaHat %*% t(X_i) # same as sigma_2 * ( X_i %*% ( solve(t(X) %*% X)  ) %*% t(X_i) )

}

CI = 2 * sqrt(var_y_hat) # Confidance interval
yHat
plot(xx[,1], yHat , type = "l")
segments(xx, yHat - CI, xx, yHat+CI) # Adds error bars to the indivigual data points



par(mfrow = c(1, 1))
plot(model_terms[,2], yy , type = "p", main = "Model Result on Test Data", xlab = "X value", ylab = "Y value")
curve(thetaHat[1]*x + thetaHat[3]*x^2 + thetaHat[4]*x^2 + thetaHat[2]*x^3 + thetaHat[5]*x^3  , -4, 4, add = TRUE, col = "red", lw = 2)


