####slide 3####
#input 2 sequences
DNA <- bio3d::read.fasta(r'(P53DNA.fasta)')
DNA <- paste(DNA$ali, collapse = '')
DNA <- substr(DNA, 45*3+1, nchar(DNA))
DNA <- Biostrings::DNAString(DNA)
AA <- as.character(Biostrings::translate(DNA))
DNA <- as.character(DNA)
AA
DNA


#I want to get the 5th residue
getResidue <- function(seq, pos){
  return(substr(seq,pos,pos))
}

#lets try on the two sequences
getResidue(AA,5)
getResidue(DNA,5)
#getResidue() can't get the 5th residue for a DNA sequence, it only gives me the 5th base

#make the two sequences into two different class of objects
P53DNA <- list(name = 'P53', seq = DNA)
class(P53DNA) <- 'DNA'

P53AA <- list(name = 'P53', seq = AA)
class(P53AA) <- 'AA'

class(P53DNA)
class(P53AA)

#write a generic function
getResidue <- function(seq, pos) UseMethod("getResidue")

#write the specific methods for the two classes
getResidue.DNA <- function(seq, pos){
  codon <- substr(seq$seq, pos*3-2, pos*3)
  resi <- as.character(Biostrings::GENETIC_CODE[[codon]])
  out <- list(resi = resi)
  class(out) <- 'AminoAcid'
  return(out)
}

getResidue.AA <- function(seq, pos){
  resi <- substr(seq$seq, pos, pos)
  out <- list(resi = resi)
  class(out) <- 'AminoAcid'
  return(out)
}

#see how it works
getResidue(P53AA,5)
getResidue(P53DNA,5)

#what if we have multiple sequences in the data set?
dt <- list(P53DNA,P53AA,P53AA,P53DNA,P53DNA,P53AA,P53DNA)
unlist(lapply(dt, getResidue, pos = 5))



####slide 4####
#get training data of speed and distance
df <- datasets::cars
head(df)
#I want to predict distance according to speed
speeds <- data.frame(speed = c(10,11,12,13,14,15))

#building linear model
linear_model <- lm(dist~speed, data = df)
summary(linear_model)

#building random forest model
install.packages('randomForest')
random_forest_model <- randomForest::randomForest(dist~speed, data = df)
summary(random_forest_model)

#building SVM model
install.packages('e1071')
SVM_model <- e1071::svm(dist~speed, data = df)
summary(SVM_model)

#predicting with different models
predict(linear_model, speeds)
predict(random_forest_model, speeds)
predict(SVM_model, speeds)

#whats going on with predict()?
install.packages('sloop')
sloop::ftype(predict)

####slide 5####
#how if we write the AA and DNA classes in S4?
setClass('DNA',
         slots = list(name = 'character',
                      source = 'character',
                      seq = 'character'))

setClass('AA',
         slots = list(name = 'character',
                      source = 'character',
                      seq = 'character'))

setClass('AminoAcid',
         slots = list(residueType = 'character'))

setGeneric("getResidue", function(seq, pos){
  standardGeneric("getResidue")
})

setMethod("getResidue",
          c(seq = "DNA", pos = 'numeric'),
          function(seq, pos){
            codon <- substr(seq@seq, pos*3-2, pos*3)
            resi <- as.character(Biostrings::GENETIC_CODE[[codon]])
            out <- new('AminoAcid', residueType = resi)
            return(out)
          })

setMethod("getResidue",
          c(seq = "AA", pos = 'numeric'),
          function(seq, pos){
            resi <- substr(seq@seq, pos, pos)
            out <- new('AminoAcid', residueType = resi)
            return(out)
          })

#see how it works
P53DNA <- new('DNA', name = 'P53', source = 'Hs', seq = DNA)
P53AA <- new('AA', name = 'P53', source = 'Hs', seq = AA)
P53DNA
getResidue(P53AA,5)
getResidue(P53DNA,5)

#what if we have multiple sequences in the data set?
dt <- list(P53DNA,P53AA,P53AA,P53DNA,P53DNA,P53AA,P53DNA)
unlist(lapply(dt, getResidue, pos = 5))

#what to do if we want a vector of those residues?
setGeneric("getChar", function(x){
  standardGeneric("getChar")
})

setMethod("getChar",
          c(x = 'AminoAcid'),
          function(x){
            return(x@residueType)
          })

AminoAcids <- lapply(dt, getResidue, pos = 5)
unlist(lapply(AminoAcids, getChar))

#classes can inherit from other classes
#Lets create classes for certified p53 sequences
setClass("p53_DNA",
         slots = list(mutated = 'logical',
                      mut_pos = 'numeric'),
         contains = "DNA")

setClass("p53_AA",
         slots = list(mutated = 'logical',
                      mut_pos = 'numeric'),
         contains = "AA")

#why is this important? Let's say we only want to mutate P53 sequences that are not mutated
#(i.e. only care about SNPs)
#DNA and AA mutates in different manner, mutate residues for AA, mutate bases for DNA
setClass('AminoAcid',
         slots = list(residueType = 'character'))

setClass("Base",
         slots = list(baseType = 'character'))


setGeneric("mutateP53", function(seq, pos, mutTo){
  standardGeneric("mutateP53")
})

setMethod("mutateP53",
          c(seq = "p53_DNA", pos = 'numeric', mutTo = 'Base'),
          function(seq, pos, mutTo){
            if(seq@mutated){
              return(print(paste0('sequence already mutated at position: ', seq@mut_pos)))
            }
            substr(seq@seq, pos, pos) <- mutTo@baseType
            seq@mutated <- T
            seq@mut_pos <- pos
            return(seq)
          })

setMethod("mutateP53",
          c(seq = "p53_AA", pos = 'numeric', mutTo = 'AminoAcid'),
          function(seq, pos, mutTo){
            if(seq@mutated){
              return(print(paste0('sequence already mutated at position: ', seq@mut_pos)))

            }
            substr(seq@seq, pos, pos) <- mutTo@residueType
            seq@mutated <- T
            seq@mut_pos <- pos
            return(seq)
          })

b <- new('Base', baseType = 'A')
a <- new('AminoAcid', residueType = 'K')

mutateP53(P53DNA, 13, b)
P53DNA_new <- new('p53_DNA', name = P53DNA@name, source = 'Hs', seq = DNA,  mutated = F, mut_pos = 0)
mutP53DNA <- mutateP53(P53DNA_new, 13, b)

#Let's check the mutated residue
getResidue(P53DNA_new, 5)
getResidue(mutP53DNA,5)
#we did not define getResidue() for class "p53_DNA", why would this function work?


predict.lm

mutateP53(P53AA, 5, a)
P53AA_new <- new('p53_AA', name = 'P53', source = 'Hs', seq = AA,  mutated = F, mut_pos = 0)
mutP53AA <- mutateP53(P53AA_new, 5, a)
getResidue(P53AA_new, 5)
getResidue(mutP53AA,5)
#Can we mutate a mutated P53 sequence?
mutateP53(mutP53AA, 10, a)
