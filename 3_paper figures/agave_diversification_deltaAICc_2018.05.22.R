aicc <- c(196.2890, 189.9803, 189.7409, 194.0388, 202.4494, 199.4794, 220.229, 196.965, 194.0175, 192.3932, 193.7995, 189.7856, 200.314, 196.8731, 208.0882, 200.5535, 194.4207, 192.7458, 194.1521, 199.0309, 199.2042)
delta <- aicc-min(aicc)
delta
getwd()
write.table(delta, file = "/Users/luna/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2018.01/agave_diversification_deltaAICc_2018.05.22.txt", row.names = F, col.names = F, sep = "\t")
