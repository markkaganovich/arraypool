lines = 31
reportname = 'MK3.txt'

G <- as.matrix(read.table('chosenlinesSfile', sep=",", nrows = 1123483))
lcls.uniform <- rep(1/lines, lines)
S.array <- as.matrix(read.table(paste('arrayBfreq',reportname), sep=","))
