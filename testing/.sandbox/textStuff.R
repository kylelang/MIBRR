

getwd()

lic <- readLines("../../source/MIBRR/LICENSE")


lic

writeLines(readLines("../../source/MIBRR/LICENSE"))

start <- grep("15. Disclaimer of Warranty", lic)
end   <- grep("END OF TERMS AND CONDITIONS", lic) - 1

writeLines(lic[start : end])

?read.dcf
?system.file

readLines(system.file("LICENSE", package = "MIBRR"))
