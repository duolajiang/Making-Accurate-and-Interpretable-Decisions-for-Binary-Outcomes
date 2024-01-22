Note that BART_LReg_*.RData in this directory is generated in the setting: 
1. f1 is BART fitted to {(x,z,y)}
2. f2 is logistic regression with main effect of X,Z and interaction effect between X*Z, fitted to {(x,z,p(z*=1|f1))}
3. f3 is logistic regression with the same functional form as f2, fitted to {(x,z,y)}
