rm(list=ls())

setwd('D:/university/Amasters/firstsemester/PRG/exercise_07')
# Motif Search

## The Brute-Force Motif Search
### Task 1
* In R, create a function `Score()`, that calculates the score for a consensus string.

* Input:
    * An array of starting indexes.
    * `DNAStringSet` object of sequences (for example file `seq_score.fasta`).
    * Motif length.
   
* Output:
    * The score for the consensus string.

library(Biostrings)

Score <- function(start_positions, motif_len, sequences, l) {
  
# Výběr motivů z každé sekvence podle startovních pozic
  motifs <- character(length(sequences))
  for (i in seq_len(length(sequences))) {
    motifs[i] <- as.character(subseq(sequences[[i]], start = start_positions[i], width = l))
  }
  
# Rozdělení motivů na jednotlivé znaky
  motif_matrix <- do.call(rbind, strsplit(motifs, split = ""))
  
# Inicializace skóre
  total_score <- 0
  
# Pro každý sloupec spočítáme četnosti písmen a přičteme nejvyšší
  for (j in seq_len(l)) {
    column <- motif_matrix[, j]
    freqs <- table(column)
    total_score <- total_score + max(freqs)
  }
  
  return(total_score)
}

# Umělá data
seqs_fasta <- readDNAStringSet("seq_score.fasta")
# seqs <- DNAStringSet(c(
#   "ATGCACTGAC",
#   "ATGAACTTAC",
#   "ATGCGCTAAC"))

# Startovací pozice pro motiv délky 5
starts <- c(1, 1, 1, 1, 1)

Score(starts, seqs_fasta, 8)


### Task 2
* In R, create function `NextLeaf()` according to the following pseudocode.

* Input:
    * `s` An array of starting indexes `s = (s1 s2 … st)`, where *t* is the number of sequences.
    * `t` Number of sequences.
    * `k` `k = n - l + 1`, where `n` is length of sequences and `l` is motif length.

* Output:
    * `s` An array of starting indexes that corresponds to the next leaf in the tree.

```
NextLeaf(s, t, k)
1   for i ← t to 1
2     if s[i] < k
3       s[i] ← s[i] + 1
4       return s
5     s[i] ← 1
6   return s
```

NextLeaf <- function(s, t, k){
  for (i in t:1){
    if (s[i] < k){
      s[i] <- s[i] + 1
      return (s)
    }
    s[i] <- 1
  }
  return (s)
}

# Počet sekvencí
t <- 3

# k = n - l + 1 (např. n = 10, l = 6)
k <- 5

# Počáteční vektor
s <- c(1, 1, 1)

# Iteruj několikrát:
for (i in 1:10) {
  print(s)
  s <- NextLeaf(s, t, k)
}


### Task 3
* In R, create a function `BFMotifSearch()` according to the following pseudocode.

* Input:
    * `DNA` `DNAStringSet` object of sequences (for example file `seq_motif.fasta`).
    * `t` Number of sequences.
    * `n` Length of each sequence.
    * `l` Motif length.

* Output:
    * `bestMotif` An array of starting positions for each sequence with the best score for the consensus string.

```
BFMotifSearch(DNA, t, n, l)
1   s ← (1, 1, ... , 1)
2   bestScore ← Score(s, DNA, l)
3   while forever
4     s ← NextLeaf(s, t, n − l + 1)
5       if Score(s, DNA, l) > bestScore
6         bestScore ← Score(s, DNA, l)
7         bestMotif ← (s1, s2, . . . , st)
8       if s = (1, 1, . . . , 1)
9         return bestMotif
```

BFMotifSearch <- function(DNA, t, n, l){
  k <- n - l + 1  # počet možných startů v jedné sekvenci
  s <- rep(1, t)  # inicializace vektoru startů
  
  bestScore <- Score(s, DNA, l)
  bestMotif <- s
  
  repeat {
    s <- NextLeaf(s, t, k)
    if (Score(s, DNA, l) > bestScore){
      bestScore <- Score(s, DNA, l)
      bestMotif <- s
    }
    if (all(s == rep(1, t))) {
      break
    }
  }
  return(bestMotif)
}

# Načtení FASTA souboru
DNA <- readDNAStringSet("seq_motif.fasta")

t <- length(DNA)        # počet sekvencí
n <- width(DNA[1])      # délka jedné sekvence (všechny stejné)
l <- 6                  # délka motivu

bestMotif <- BFMotifSearch(DNA, t, n, l)

cat("Nejlepší startovací pozice motivu:", bestMotif, "\n")

## The Branch-and-Bound Motif Search
### Task 4
* In R, create a function `NextVertex()` according to the following pseudocode.

* Input:
    * `s` An array of starting indexes `s = (s1 s2 … st)`, where *t* is the number of sequences.
    * `i` Level of vertex.
    * `t` Number of sequences.
    * `k` `k = n - l + 1`, where `n` is length of sequences and `l` is motif length.

* Output:
    * `s` The next vertex in the tree.
    * Current level of vertex.

```
NextVertex(s, i, t, k)
1   if i < t
2     s[i + 1] ← 1
3     return (s, i + 1)
4   else
5     for j ← t to 1
6       if s[j] < k
7         s[j] ← s[j] + 1
8         return (s, j)
9   return (s, 0)
```

NextVertex <- function(s, i, t, k){
  if (i < t){
    s[i+1] <- 1
    return(list(s = s, i = i + 1))
  } else {
    for (j in t:1){
      if (s[j] < k){
      s[j] <- s[j] +1
      return(list(s = s,i = j))
      }
    }
  }
  return(list(s = s, i = 0))
}

# Parametry
t <- 3   # počet sekvencí
k <- 4   # počet možných startů (n - l + 1)

# Počáteční vrchol
s <- c(1, 1, 1)
i <- 1

# Posouváme se stromem
res <- NextVertex(s, i, t, k)
print(res)

# Další krok
res <- NextVertex(res$s, res$i, t, k)
print(res)

### Task 5
* In R, create a function `ByPass()` according to the following pseudocode.

* Input:
    * `s = (s1 s2 … st)`; an array of starting indexes, where *t* is the number of sequences
    * `i`; level of vertex
    * `t`; number of DNA sequences
    * `k = n - l + 1`, where `n` is length of DNA sequences and `l` is motif length

* Output:
    * the next leaf after a skip of a subtree
    * current level of vertex

```
ByPass(s, i, t, k)
1   for j ← i to 1
2     if s[j] < k
3       s[j] ← s[j] + 1
4       return (s, j)
5   return (s, 0)
```

ByPass <- function(s, i, t, k){
  for (j in i:1){
    if (s[j] < k){
      s[j] <- s[j] + 1
      return(list(s=s, i=j))
    }
  }
  return(list(s=s, i=0))
}

t <- 4    # počet sekvencí
k <- 5    # počet možných startovacích pozic
s <- c(1, 3, 5, 5)
i <- 4    # aktuální úroveň

res <- ByPass(s, i, t, k)
print(res)


### Task 6
* In R, create a function `BBMotifSearch()` according to the following pseudocode.

* Input:
    * `DNA` `DNAStringSet` object of sequences (for example file `seq_motif.fasta`).
    * `t` Number of sequences.
    * `n` Length of each sequence.
    * `l` Motif length.
* Output:
    * `bestMotif` An array of starting positions for each sequence with the best score for the consensus string.

* Modify function `Score()` to calculate score for the consensus string of the first `i` sequences of `DNA`.

```
BBMotifSearch(DNA, t, n, l)
1   s ← (1, ... , 1)
2   bestScore ← 0
3   i ← 1
4   while i > 0
5     if i < t
6       optimisticScore ← Score(s, i, DNA, l) + (t - i) * l
7       if optimisticScore < bestScore
8         (s, i) ← ByPass(s, i, t, n - l + 1)
9       else
10        (s, i) ← NextVertex(s, i, t, n − l + 1)
11    else
12      if Score(s, t, DNA, l) > bestScore
13        bestScore ← Score(s, t, DNA, l)
14        bestMotif ← (s1, s2, ... , st)
15      (s, i) ← NextVertex(s, i, t, n − l + 1)
16  return bestMotif
```
BBMotifSearch <- function(DNA, t, n, l){
  s <- rep(1, t)
  bestScore <- 0
  i <- 1
  while (i > 0){
    if (i < t){
      optimisticScore <- Score(s, i, DNA, l) + ((t - i) *1)
      if (optimisticScore < bestScore){
        temp <- ByPass(s, i, t, n - l + 1)
        s <- temp$s
        i <- temp$i 
      } else {
          temp <- NextVertex(s, i, t, n - l + 1)
          s <- temp$s
          i <- temp$i 
        }
    } else {
      if (Score(s, t, DNA, l) > bestScore){
          bestScore <- Score(s, t, DNA, l)
          bestMotif <- s
      }
      temp <- NextVertex(s, i, t, n - l + 1)
      s <- temp$s
      i <- temp$i 
    }
  }
  return(bestMotif)
}

DNA <- readDNAStringSet('seq_motif.fasta')
t <- length(DNA)   # počet sekvencí
n <- width(DNA)[1] # délka sekvencí
l <- 5             # délka motivu

bestMotif <- BBMotifSearch(DNA, t, n, l)
cat("Best motif start positions:", bestMotif, "\n")

<details>
<summary>Download files from GitHub</summary>
<details>
<summary>Basic Git settings</summary>

> * Configure the Git editor
>     ```bash
>     git config --global core.editor notepad
>     ```
> * Configure your name and email address
>     ```bash
>     git config --global user.name "Zuzana Nova"
>     git config --global user.email z.nova@vut.cz
>     ```
> * Check current settings
>     ```bash
>     git config --global --list
>     ```
>
</details>

* Create a fork on your GitHub account. 
  On the GitHub page of this repository find a <kbd>Fork</kbd> button in the upper right corner.
  
* Clone forked repository from your GitHub page to your computer:
    ```bash
    git clone <fork repository address>
    ```
* In a local repository, set new remote for a project repository:
```bash
git remote add upstream https://github.com/mpa-prg/exercise_07.git
```

#### Send files to GitHub
Create a new commit and send new changes to your remote repository.
* Add file to a new commit.
    ```bash
    git add <file_name>
    ```
* Create a new commit, enter commit message, save the file and close it.
    ```bash
    git commit
    ```
* Send a new commit to your GitHub repository.
    ```bash
    git push origin main
    ```

</details>
