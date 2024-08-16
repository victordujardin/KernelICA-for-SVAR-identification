## Create exogenous variables
## Dummies for tax event
time    <- seq(as.Date("1960/1/1"), by = "quarter", length.out = 152)
tD0     <- as.numeric(time == "1975-04-01")
tD1     <- as.numeric(time == "1975-07-01")
tD2     <- as.numeric(time == "1975-10-01")
tD3     <- as.numeric(time == "1976-01-01")
tD4     <- as.numeric(time == "1976-04-01")

## Time trends & const
constant <- rep(1,152)
t.lin  <- seq(from = 1, length.out = 152,by = 1)
t.sq   <- t.lin^2

## Seasonal Dimmies
q1     <- as.numeric(grepl(x = time,pattern = "-01-"))
q2     <- as.numeric(grepl(x = time,pattern = "-04-"))
q3     <- as.numeric(grepl(x = time,pattern = "-07-"))
q4     <- as.numeric(grepl(x = time,pattern = "-10-"))

## Lag is taken from Hmisc package
q1tax  <- q1*tsy[,"TT"]
q1tax1 <- q1*(Lag(x = tsy[,"TT"],shift = 1))
q1tax2 <- q1*(Lag(x = tsy[,"TT"],shift = 2))
q1tax3 <- q1*(Lag(x = tsy[,"TT"],shift = 3))
q1tax4 <- q1*(Lag(x = tsy[,"TT"],shift = 4))
q1gcn  <- q1*tsy[,"GG"]
q1gcn1 <- q1*(Lag(x = tsy[,"GG"],shift = 1))
q1gcn2 <- q1*(Lag(x = tsy[,"GG"],shift = 2))
q1gcn3 <- q1*(Lag(x = tsy[,"GG"],shift = 3))
q1gcn4 <- q1*(Lag(x = tsy[,"GG"],shift = 4))
q1gdp <-  q1*tsy[,"XX"]
q1gdp1 <- q1*(Lag(x = tsy[,"XX"],shift = 1))
q1gdp2 <- q1*(Lag(x = tsy[,"XX"],shift = 2))
q1gdp3 <- q1*(Lag(x = tsy[,"XX"],shift = 3))
q1gdp4 <- q1*(Lag(x = tsy[,"XX"],shift = 4))
q2tax  <- q2*tsy[,"TT"]
q2tax1 <- q2*(Lag(x = tsy[,"TT"],shift = 1))
q2tax2 <- q2*(Lag(x = tsy[,"TT"],shift = 2))
q2tax3 <- q2*(Lag(x = tsy[,"TT"],shift = 3))
q2tax4 <- q2*(Lag(x = tsy[,"TT"],shift = 4))
q2gcn  <- q2*tsy[,"GG"]
q2gcn1 <- q2*(Lag(x = tsy[,"GG"],shift = 1))
q2gcn2 <- q2*(Lag(x = tsy[,"GG"],shift = 2))
q2gcn3 <- q2*(Lag(x = tsy[,"GG"],shift = 3))
q2gcn4 <- q2*(Lag(x = tsy[,"GG"],shift = 4))
q2gdp  <- q2*tsy[,"XX"]
q2gdp1 <- q2*(Lag(x = tsy[,"XX"],shift = 1))
q2gdp2 <- q2*(Lag(x = tsy[,"XX"],shift = 2))
q2gdp3 <- q2*(Lag(x = tsy[,"XX"],shift = 3))
q2gdp4 <- q2*(Lag(x = tsy[,"XX"],shift = 4))
q3tax  <- q3*tsy[,"TT"]
q3tax1 <- q3*(Lag(x = tsy[,"TT"],shift = 1))
q3tax2 <- q3*(Lag(x = tsy[,"TT"],shift = 2))
q3tax3 <- q3*(Lag(x = tsy[,"TT"],shift = 3))
q3tax4 <- q3*(Lag(x = tsy[,"TT"],shift = 4))
q3gcn  <- q3*tsy[,"GG"]
q3gcn1 <- q3*(Lag(x = tsy[,"GG"],shift = 1))
q3gcn2 <- q3*(Lag(x = tsy[,"GG"],shift = 2))
q3gcn3 <- q3*(Lag(x = tsy[,"GG"],shift = 3))
q3gcn4 <- q3*(Lag(x = tsy[,"GG"],shift = 4))
q3gdp  <- q3*tsy[,"XX"]
q3gdp1 <- q3*(Lag(x = tsy[,"XX"],shift = 1))
q3gdp2 <- q3*(Lag(x = tsy[,"XX"],shift = 2))
q3gdp3 <- q3*(Lag(x = tsy[,"XX"],shift = 3))
q3gdp4 <- q3*(Lag(x = tsy[,"XX"],shift = 4))
q4tax  <- q4*tsy[,"TT"]
q4tax1 <- q4*(Lag(x = tsy[,"TT"],shift = 1))
q4tax2 <- q4*(Lag(x = tsy[,"TT"],shift = 2))
q4tax3 <- q4*(Lag(x = tsy[,"TT"],shift = 3))
q4tax4 <- q4*(Lag(x = tsy[,"TT"],shift = 4))
q4gcn  <- q4*tsy[,"GG"]
q4gcn1 <- q4*(Lag(x = tsy[,"GG"],shift = 1))
q4gcn2 <- q4*(Lag(x = tsy[,"GG"],shift = 2))
q4gcn3 <- q4*(Lag(x = tsy[,"GG"],shift = 3))
q4gcn4 <- q4*(Lag(x = tsy[,"GG"],shift = 4))
q4gdp  <- q4*tsy[,"XX"]
q4gdp1 <- q4*(Lag(x = tsy[,"XX"],shift = 1))
q4gdp2 <- q4*(Lag(x = tsy[,"XX"],shift = 2))
q4gdp3 <- q4*(Lag(x = tsy[,"XX"],shift = 3))
q4gdp4 <- q4*(Lag(x = tsy[,"XX"],shift = 4))

EXOG <- as.matrix(cbind(#constant,
#                        t.lin, 
                        t.sq,
                        q1, q2, q3,
                        tD0, tD1, tD2, tD3, tD4,
                        as.vector(q2gcn1),
                        as.vector(q2tax1),
                        as.vector(q2gdp1),
                        as.vector(q2gcn2),
                        as.vector(q2tax2),
                        as.vector(q2gdp2),
                        as.vector(q2gcn3),
                        as.vector(q2tax3),
                        as.vector(q2gdp3),
                        as.vector(q2gcn4),
                        as.vector(q2tax4),
                        as.vector(q2gdp4),
                        as.vector(q3gcn1),
                        as.vector(q3tax1),
                        as.vector(q3gdp1),
                        as.vector(q3gcn2),
                        as.vector(q3tax2),
                        as.vector(q3gdp2),
                        as.vector(q3gcn3),
                        as.vector(q3tax3),
                        as.vector(q3gdp3),
                        as.vector(q3gcn4),
                        as.vector(q3tax4),
                        as.vector(q3gdp4),
                        as.vector(q4gcn1),
                        as.vector(q4tax1),
                        as.vector(q4gdp1),
                        as.vector(q4gcn2),
                        as.vector(q4tax2),
                        as.vector(q4gdp2),
                        as.vector(q4gcn3),
                        as.vector(q4tax3),
                        as.vector(q4gdp3),
                        as.vector(q4gcn4),
                        as.vector(q4tax4),
                        as.vector(q4gdp4)))


colnames(EXOG) <- c(#"constant",
                    #"linear", 
                    "quadratic",
                    "q1", "q2", "q3", 
                    "dtax0","dtax1", "dtax2", "dtax3", "dtax4",
                    "q2gcn1",
                    "q2tax1",
                    "q2gdp1",
                    "q2gcn2",
                    "q2tax2",
                    "q2gdp2",
                    "q2gcn3",
                    "q2tax3",
                    "q2gdp3",
                    "q2gcn4",
                    "q2tax4",
                    "q2gdp4",
                    "q3gcn1",
                    "q3tax1",
                    "q3gdp1",
                    "q3gcn2",
                    "q3tax2",
                    "q3gdp2",
                    "q3gcn3",
                    "q3tax3",
                    "q3gdp3",
                    "q3gcn4",
                    "q3tax4",
                    "q3gdp4",
                    "q4gcn1",
                    "q4tax1",
                    "q4gdp1",
                    "q4gcn2",
                    "q4tax2",
                    "q4gdp2",
                    "q4gcn3",
                    "q4tax3",
                    "q4gdp3",
                    "q4gcn4",
                    "q4tax4",
                    "q4gdp4")
time.dummy <- cbind(tD0,tD1,tD2,tD3,tD4)

EXOG2    <- as.matrix(cbind(#constant, 
                            t.lin, 
                            #t.sq,
                            q1, q2, q3,
                            tD0, tD1, tD2, tD3, tD4))
colnames(EXOG2) <- c(#"constant",
                     "linear", 
                    #"quadratic",
                    "q1", "q2", "q3", 
                    "dtax0","dtax1", "dtax2", "dtax3", "dtax4")
