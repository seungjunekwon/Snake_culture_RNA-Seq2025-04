library(AnnotationForge)

gaf <- read.delim("snake.gaf", skip=9, header=FALSE) # gaf from ncbi ftp

fGO <- gaf[, c(2,5,7)]
colnames(fGO) <- c("GID", "GO", "EVIDENCE")
fGO <- fGO[!duplicated(fGO), ] # Remove duplicated rows

fSym <- read.csv("id.csv")

makeOrgPackage(gene_info = fSym, go=fGO,
               version = "1.1",
               maintainer = "Seung June Kwon <udnjf251@gmail.com>",
               author = "Seung June Kwon",
               outputDir = ".",
               tax_id = "94885",
               genus = "Pantherophis",
               species = "guttatus",
               goTable = "go")

install.packages("./org.Pguttatus.eg.db/", repos = NULL, type = "source")