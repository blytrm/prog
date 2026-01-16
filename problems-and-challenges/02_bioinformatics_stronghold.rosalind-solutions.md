# Rosalind Problems + Answers

Solutions to bioinformatics problems from [Rosalind](http://rosalind.info/), implemented in both Python and R.

# ROSALIND Solutions – Python & R

## 1. Counting DNA Nucleotides

**Problem**  
Given: A DNA string `s` of length at most 1000 bp.  
Return: Four integers (separated by spaces) counting the respective number of occurrences of `A`, `C`, `G`, and `T` in `s`.

**Python**

```
with open("Counting DNA Nucleotides (1).txt", "r") as file:
    S = file.readline().strip()

bases = ("A", "C", "G", "T")
chr_count = {}

for base in S:
    chr_count[base] = chr_count.get(base, 0) + 1

counts = [chr_count.get(base, 0) for base in bases]
print(" ".join(str(count) for count in counts))
```

**R**

```
txt_path <- "Counting DNA Nucleotides (1).txt"
S <- readLines(txt_path)
bases <- strsplit(S, "")[][1]
count <- table(bases)
print(count)
```

---

## 2. Transcribing DNA into RNA

**Problem**  
Given: A DNA string `t` corresponding to a coding strand.  
Return: The transcribed RNA string `u` formed by replacing all `T` with `U`.

**Python**

```
with open("Transcribing DNA to RNA.txt", "r") as file:
    t = file.readline().strip()

transc = t.replace("T", "U")
print(transc)
```

**R**

```
txt <- "Transcribing DNA to RNA.txt"
t <- readLines(txt)
transc <- gsub(pattern = "T", replacement = "U", x = t)
print(transc)

# or
transc2 <- chartr(old = "T", new = "U", x = t)
print(transc2)
```

---

## 3. Complementing a Strand of DNA

**Problem**  
Given: A DNA string `s`.  
Return: The reverse complement of `s`.

**Python**

```
from Bio.Seq import Seq

with open("DNA Complement Strand.txt", "r") as file:
    s = file.readline().strip()

seq = Seq(s)
rev_comp = seq.reverse_complement()
print(str(rev_comp))
```

**R**

```
txt <- "DNA Complement Strand.txt"
s <- readLines(txt)[1]

rc <- function(x) {
  ss    <- unlist(strsplit(x, NULL))
  ssrev <- rev(ss)

  for (i in seq_along(ssrev)) {
    if (ssrev[i] == "T") {
      ssrev[i] <- "A"
    } else if (ssrev[i] == "A") {
      ssrev[i] <- "T"
    } else if (ssrev[i] == "G") {
      ssrev[i] <- "C"
    } else if (ssrev[i] == "C") {
      ssrev[i] <- "G"
    }
  }

  paste(ssrev, collapse = "")
}

rc_seq <- rc(s)
cat(rc_seq, "\n")
```

---

## 4. Rabbits and Recurrence Relations

**Problem**  
Given: Positive integers `n` and `k`. Start with 1 pair of rabbits; each reproduction‑age pair produces `k` new pairs each month.  
Return: The total number of rabbit pairs after `n` months.

**R**

```
rabbit_pairs_month <- c()
n <- 28  # months
k <- 2   # new pairs per reproducing pair

for (i in 1:n) {
  if (i == 1) {
    rabbit_pairs_month[i] <- 1
  } else if (i == 2) {
    rabbit_pairs_month[i] <- 1
  } else {
    rabbit_pairs_month[i] <- rabbit_pairs_month[i - 1] + k * rabbit_pairs_month[i - 2]
  }
}

rabbit_pairs_month[length(rabbit_pairs_month)]
```

**Python**

```
n, k = 28, 2

def mortal_r(n, k):
    if n == 1 or n == 2:
        return 1

    F_prev2 = 1  # F1
    F_prev1 = 1  # F2

    for _ in range(3, n + 1):
        F_curr = F_prev1 + k * F_prev2
        F_prev2, F_prev1 = F_prev1, F_curr

    return F_curr

print(mortal_r(n, k))


def fib(n, k):
    previous1, previous2 = 1, 1
    for i in range(2, n):
        current = previous1 + k * previous2
        previous2 = previous1
        previous1 = current
    return current
```

---

## 5. GC Content

**Problem**  
Given: At most 10 DNA strings in FASTA format (length ≤ 1 kbp).  
Return: The ID of the string with the highest GC content, followed by its GC percentage (to within 0.001).

**Python**

```
def readFile(filePath):
    with open(filePath, "r") as f:
        return [l.strip() for l in f.readlines()]

ffile = readFile("GC Content Analysis.txt")

fdict = {}
flabel = ""

for line in ffile:
    if ">" in line:
        flabel = line
        fdict[flabel] = ""
    else:
        fdict[flabel] += line

def gc_content(seq):
    return round(((seq.count("C") + seq.count("G")) / len(seq) * 100), 6)

resultsDict = {key: gc_content(value) for (key, value) in fdict.items()}

maxGC_key = max(resultsDict, key=resultsDict.get)
print(maxGC_key)
print(f"{maxGC_key[1:]}\n{resultsDict[maxGC_key]}")
```

**R**

```
file_path <- "GC Content Analysis.txt"
lines <- readLines(file_path)

seq_list  <- list()
current_id <- NULL

for (line in lines) {
  if (startsWith(line, ">")) {
    current_id <- substring(line, 2)
    seq_list[[current_id]] <- ""
  } else {
    seq_list[[current_id]] <- paste0(seq_list[[current_id]], line)
  }
}

gc_content <- function(seq) {
  chars <- strsplit(seq, "")[][1]
  g <- sum(chars == "G")
  c <- sum(chars == "C")
  gc <- (g + c) / length(chars) * 100
  round(gc, 6)
}

gc_values <- sapply(seq_list, gc_content)

max_id <- names(gc_values)[which.max(gc_values)]
max_gc <- gc_values[max_id]

cat(max_id, "\n", max_gc, "\n", sep = "")
```

---

## 6. Counting Point Mutations

**Problem**  
Given: Two DNA strings `s` and `t` of equal length.  
Return: The Hamming distance between them (number of differing positions).

**Python**

```
def readFile(filePath):
    with open(filePath, "r") as file:
        return [l.strip() for l in file.readlines()]

ffile = readFile("Counting Point Mutations (1).txt")
s = ffile
t = ffile[1]

difrncs = 0

if len(s) == len(t):
    for char1, char2 in zip(s, t):
        if char1 != char2:
            difrncs += 1

print(difrncs)
```

**R**

```
file_path <- "Counting Point Mutations (1).txt"
lines <- readLines(file_path)
s <- lines[1]
t <- lines[2]

hamming_distance <- function(a, b) {
  count <- 0
  for (i in seq_len(nchar(a))) {
    if (substr(a, i, i) != substr(b, i, i)) {
      count <- count + 1
    }
  }
  count
}

d <- hamming_distance(s, t)
cat(d, "\n")
```

---

## 7. Mendel’s First Law

**Problem**  
Given: Three positive integers `k`, `m`, `n` representing a population of `k + m + n` organisms:  
- `k` homozygous dominant (`AA`),  
- `m` heterozygous (`Aa`),  
- `n` homozygous recessive (`aa`).  

Return: The probability that two randomly selected organisms will produce an offspring with at least one dominant allele.

**Python**

```
with open("Mendel's First Law (1).txt", "r") as f:
    line = f.readline().strip()

values = line.split()
k, m, n = map(int, values)
N = k + m + n

P_AA_AA = (k / N) * ((k - 1) / (N - 1))
P_AA_Aa = (k / N) * (m / (N - 1)) + (m / N) * (k / (N - 1))
P_AA_aa = (k / N) * (n / (N - 1)) + (n / N) * (k / (N - 1))
P_Aa_Aa = (m / N) * ((m - 1) / (N - 1))
P_Aa_aa = (m / N) * (n / (N - 1)) + (n / N) * (m / (N - 1))
P_aa_aa = (n / N) * ((n - 1) / (N - 1))

Pr_dom = (
    P_AA_AA * 1.0 +
    P_AA_Aa * 1.0 +
    P_AA_aa * 1.0 +
    P_Aa_Aa * 0.75 +
    P_Aa_aa * 0.5 +
    P_aa_aa * 0.0
)

print(f"{Pr_dom:.5f}")
```

**R**

```
line <- readLines("Mendel's First Law (1).txt")
line <- trimws(line)

values <- strsplit(line, " ")[][1]
nums <- as.numeric(values)
k <- nums[1]
m <- nums[2]
n <- nums[3]
N <- k + m + n

P_AA_AA <- (k / N) * ((k - 1) / (N - 1))
P_AA_Aa <- (k / N) * (m / (N - 1)) + (m / N) * (k / (N - 1))
P_AA_aa <- (k / N) * (n / (N - 1)) + (n / N) * (k / (N - 1))
P_Aa_Aa <- (m / N) * ((m - 1) / (N - 1))
P_Aa_aa <- (m / N) * (n / (N - 1)) + (n / N) * (m / (N - 1))
P_aa_aa <- (n / N) * ((n - 1) / (N - 1))

Pr_dom <-
  P_AA_AA * 1.0 +
  P_AA_Aa * 1.0 +
  P_AA_aa * 1.0 +
  P_Aa_Aa * 0.75 +
  P_Aa_aa * 0.5 +
  P_aa_aa * 0.0

cat(sprintf("%.5f\n", Pr_dom))
```

---

## 8. Translating RNA into Protein

**Problem**  
Given: An RNA string `s` corresponding to an mRNA strand.  
Return: The protein string encoded by `s`.

**Python**

```
from Bio.Seq import Seq

def readFile(filePath):
    with open(filePath, "r") as f:
        content = f.readline().strip()
    return content

s = readFile("RNA to Protein (2).txt")

rna = Seq(s)
protein = rna.translate(to_stop=True)
print(protein)
```

**R**

```
library(Biostrings)
library(bioseq)

file_path <- "RNA to Protein.txt"
s <- trimws(readLines(file_path))[1]

# Biostrings
rna_bs <- RNAString(s)
protein_bs <- translate(rna_bs)
print(protein_bs)

# bioseq
rna_bq <- rna(s)
protein_bq <- seq_translate(rna_bq)
print(protein_bq)
```

---

## 9. Finding a Motif in DNA

**Problem**  
Given: Two strings `s` and `t`, with `t` a substring of `s`.  
Return: All locations (1‑based) where `t` appears in `s`.

**Python – simple double loop**

```
def read_file(path):
    with open(path) as f:
        lines = [l.strip() for l in f]
    return lines, lines[1]

s, t = read_file("Finding a Motif in DNA.txt")

n = len(s)
m = len(t)
positions = []

for i in range(n - m + 1):
    has_match = True
    for j in range(m):
        if s[i + j] != t[j]:
            has_match = False
            break
    if has_match:
        positions.append(i + 1)

print(" ".join(map(str, positions)))
```

**R (Biostrings)**

```
library(Biostrings)

lines <- readLines("Finding a Motif in DNA.txt")
s <- lines[1]
t <- lines[2]

dna   <- DNAString(s)
motif <- DNAString(t)

hits <- matchPattern(motif, dna)
pos  <- start(hits)

cat(paste(pos, collapse = " "), "\n")
```

---

## 10. Consensus and Profile

**Problem**  
Given: A collection of DNA strings in FASTA format, all of the same length.  
Return: A consensus string and a profile matrix (counts of A, C, G, T at each position).

**Python**

```
from Bio import motifs
from Bio.Seq import Seq
from Bio import SeqIO
from io import StringIO

def readFile(filePath):
    with open(filePath, "r") as f:
        return [l.strip() for l in f.readlines()]

lines = readFile("Consensus Profile Rosalind.txt")
fasta_content = "\n".join(lines)

seq_strings = []
for record in SeqIO.parse(StringIO(fasta_content), "fasta"):
    seq_strings.append(str(record.seq))

def profile_matrix(sequences):
    seq_length = len(sequences)
    profile = {
        "A":  * seq_length,
        "C":  * seq_length,
        "G":  * seq_length,
        "T":  * seq_length,
    }
    for seq in sequences:
        for i, n in enumerate(seq):
            if n in profile:
                profile[n][i] += 1
    return profile

def consensus(profile):
    cols = len(next(iter(profile.values())))
    con = []
    for i in range(cols):
        col_count = {nt: profile[nt][i] for nt in "ACGT"}
        con.append(max(col_count, key=col_count.get))
    return "".join(con)

profile = profile_matrix(seq_strings)
con = consensus(profile)

print(con)
for nt in "ACGT":
    print(f"{nt}: {' '.join(map(str, profile[nt]))}")

seq_obj = [Seq(s) for s in seq_strings]
m = motifs.create(seq_obj)
print(m.consensus)
print(m.counts)
```

**R**

```
library(Biostrings)
library(seqinr)

fasta_file <- "Consensus Proile Rosalind (2).txt"
seq_list   <- read.fasta(file = fasta_file, seqtype = "DNA", as.string = FALSE)

seq_strings <- vapply(seq_list, paste0, collapse = "", FUN.VALUE = character(1))

prof_matrix <- function(sequences) {
  seq_length <- nchar(sequences)[1]
  profile <- list(
    A = integer(seq_length),
    C = integer(seq_length),
    G = integer(seq_length),
    T = integer(seq_length)
  )
  for (seq in sequences) {
    nts <- strsplit(seq, "")[][1]
    for (i in seq_along(nts)) {
      n <- nts[i]
      if (n %in% names(profile)) {
        profile[[n]][i] <- profile[[n]][i] + 1
      }
    }
  }
  profile
}

profile <- prof_matrix(seq_strings)

consensus_matrix_acgt <- function(profile) {
  cols <- length(profile[])[1]
  consensus <- character(cols)
  for (i in seq_len(cols)) {
    col_counts <- sapply(c("A", "C", "G", "T"), function(nt) profile[[nt]][i])
    consensus[i] <- names(which.max(col_counts))
  }
  paste0(consensus, collapse = "")
}

cons <- consensus_matrix_acgt(profile)
cat(cons, "\n")

for (nt in c("A", "C", "G", "T")) {
  cat(nt, ":", paste(profile[[nt]], collapse = " "), "\n")
}
```
---

## 11. Mortal Fibonacci Rabbits

**Problem**  
Given: Positive integers `n` and `m`.  
Return: The total number of rabbit pairs that will remain after the `n`-th month if all rabbits live for `m` months.

**Python**

```python
def read(path):
    with open(path, 'r') as f:
        return [l.strip() for l in f.readlines()]

fi = read('Fibonacci Rabbits (2).txt')
n, m = map(int, fi[0].split())
k = 1  # each pair produces 1 new pair

def mortal_rabbits(n, m, k):
    ages = [0] * m  # age groups: index 0 = newborns
    ages[0] = 1     # month 1: 1 newborn pair
    
    for month in range(1, n):
        breeders = sum(ages[1:])  # all except newborns
        newborns = breeders * k
        ages = [newborns] + ages[:-1]  # shift ages, oldest die
        
    return sum(ages)

print(mortal_rabbits(n, m, k))
```

**R**

```r
library(gmp)  # for large integers

fi <- readLines('Fibonacci Rabbits (2).txt')
n <- as.integer(strsplit(fi[1], " ")[[1]][1])
m <- as.integer(strsplit(fi[1], " ")[[1]][2])
k <- as.bigz(1)

mortal_rabbits <- function(n, m, k = 1) {
  ages <- c(as.bigz(1), rep(as.bigz(0), m - 1))
  
  for (i in 1:(n - 1)) {
    newborns <- k * sum(ages[-1])  # sum of all except newborns
    ages <- c(newborns, ages[-length(ages)])  # newborns + all except oldest
  }
  
  return(sum(ages))
}

cat(as.character(mortal_rabbits(n, m, k)))
```

---

## 12. Overlap Graphs

**Problem**  
Given: A collection of DNA strings in FASTA format (length ≤ 10 kbp).  
Return: The adjacency list corresponding to O₃ (k=3). You may return edges in any order.

**Python**

```python
def read_fasta(path):
    with open(path, 'r') as f:
        lines = [l.strip() for l in f.readlines()]
    
    sequences = {}
    current_label = ""
    
    for line in lines:
        if line.startswith(">"):
            current_label = line[1:]
            sequences[current_label] = ""
        else:
            sequences[current_label] += line
    
    return sequences

sequences = read_fasta('Overlap Graphs.txt')
k = 3

for s_label, s_dna in sequences.items():
    for t_label, t_dna in sequences.items():
        if s_label != t_label and s_dna[-k:] == t_dna[:k]:
            print(f"{s_label} {t_label}")
```

**R (base R)**

```r
read_fasta <- function(path) {
  lines <- readLines(path)
  sequences <- list()
  current_label <- ""
  
  for (line in lines) {
    if (startsWith(line, ">")) {
      current_label <- substring(line, 2)
      sequences[[current_label]] <- ""
    } else {
      sequences[[current_label]] <- paste0(sequences[[current_label]], line)
    }
  }
  
  return(sequences)
}

sequences <- read_fasta('Overlap Graphs (1).txt')
k <- 3
labels <- names(sequences)

for (s_name in labels) {
  for (t_name in labels) {
    if (s_name != t_name) {
      s_dna <- sequences[[s_name]]
      t_dna <- sequences[[t_name]]
      
      suffix <- substr(s_dna, nchar(s_dna) - k + 1, nchar(s_dna))
      prefix <- substr(t_dna, 1, k)
      
      if (suffix == prefix) {
        cat(s_name, t_name, "\n")
      }
    }
  }
}
```

**R (with stringr)**

```r
library(Biostrings)
library(stringr)

dna_obj <- readDNAStringSet("Overlap Graphs (1).txt")
dna_seqs <- as.character(dna_obj)
k <- 3
ids <- names(dna_seqs)

for (s in ids) {
  for (t in ids) {
    if (s != t && str_sub(dna_seqs[[s]], -k) == str_sub(dna_seqs[[t]], 1, k)) {
      cat(s, t, "\n")
    }
  }
}
```

---

## 13. Calculating Expected Offspring

**Problem**  
Given: Six nonnegative integers corresponding to the number of couples with genotypes:  
1. AA-AA  
2. AA-Aa  
3. AA-aa  
4. Aa-Aa  
5. Aa-aa  
6. aa-aa  

Return: The expected number of offspring displaying the dominant phenotype in the next generation, assuming each couple has exactly two offspring.

**Python (NumPy)**

```python
import numpy as np

def read_expected(path):
    with open(path, 'r') as f:
        data = list(map(int, f.read().split()))
    
    couples = np.array(data)
    probs = [1.0, 1.0, 1.0, 0.75, 0.50, 0.0]
    offspring_per_couple = 2
    
    expected = np.sum(couples * probs * offspring_per_couple)
    return expected

result = read_expected('Calculating Expected Offspring.txt')
print(result)
```

**Python (direct calculation)**

```python
def calculate_expected_offspring(path):
    with open(path, 'r') as f:
        a, b, c, d, e, f = map(int, f.read().split())
    
    expected = 2 * (a + b + c) + 1.5 * d + 1.0 * e
    return expected

print(calculate_expected_offspring('Calculating Expected Offspring.txt'))
```

**R**

```r
read_expected <- function(path) {
  couples <- scan(path, what = integer(), quiet = TRUE)[1:6]
  probs <- c(1.0, 1.0, 1.0, 0.75, 0.50, 0.0)
  offspring_per_couple <- 2
  
  expected <- sum(couples * probs * offspring_per_couple)
  return(expected)
}

result <- read_expected('Calculating Expected Offspring.txt')
cat(result, "\n")
```

---

## 14. Finding a Shared Motif

**Problem**  
Given: A collection of k DNA strings in FASTA format (length ≤ 1 kbp each).  
Return: A longest common substring of the collection.

**Python**

```python
def read_fasta(path):
    sequences = []
    current_dna = ""
    
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith(">"):
                if current_dna:
                    sequences.append(current_dna)
                current_dna = ""
            else:
                current_dna += line
        
        if current_dna:
            sequences.append(current_dna)
    
    return sequences

def longest_common_substring(sequences):
    if not sequences:
        return ""
    
    reference = min(sequences, key=len)
    ref_len = len(reference)
    
    for window_size in range(ref_len, 0, -1):
        for start in range(ref_len - window_size + 1):
            substring = reference[start:start + window_size]
            if all(substring in seq for seq in sequences):
                return substring
    
    return ""

sequences = read_fasta('Shared Motif LCSM.txt')
result = longest_common_substring(sequences)
print(result)
```

**R**

```r
read_fasta_manual <- function(path) {
  lines <- readLines(path)
  sequences <- character()
  current_dna <- ""
  
  for (line in lines) {
    line <- trimws(line)
    if (nchar(line) == 0) next
    
    if (startsWith(line, ">")) {
      if (nchar(current_dna) > 0) {
        sequences <- c(sequences, current_dna)
      }
      current_dna <- ""
    } else {
      current_dna <- paste0(current_dna, line)
    }
  }
  
  if (nchar(current_dna) > 0) {
    sequences <- c(sequences, current_dna)
  }
  
  return(sequences)
}

longest_common_substring <- function(seqs) {
  seqs <- toupper(seqs)
  idx_min <- which.min(nchar(seqs))
  ref <- seqs[idx_min]
  others <- seqs[-idx_min]
  ref_len <- nchar(ref)
  
  for (window_size in ref_len:1) {
    starts <- 1:(ref_len - window_size + 1)
    for (start in starts) {
      sub <- substr(ref, start, start + window_size - 1)
      if (all(vapply(others, grepl, logical(1), 
                     pattern = sub, fixed = TRUE))) {
        return(sub)
      }
    }
  }
  
  return("")
}

sequences <- read_fasta_manual('Shared Motif LCSM (1).txt')
result <- longest_common_substring(sequences)
cat(result, "\n")
```

---

## 15. Independent Alleles

**Problem**  
Given: Two positive integers k and N.  
Return: The probability that at least N AaBb organisms will belong to the k-th generation (starting with AaBb in generation 0). Assume each organism always mates with an AaBb organism.

**Python**

```python
import math

def read_alleles(path):
    with open(path, 'r') as f:
        k, n = map(int, f.readline().split())
    return k, n

k, n = read_alleles('Independent Alleles.txt')

def probability_at_least_n(k, n):
    m = 2 ** k  # total organisms in k-th generation
    p = 0.25    # probability of AaBb in offspring
    q = 1 - p
    
    total_prob = 0
    for i in range(n, m + 1):
        comb = math.comb(m, i)
        prob_i = comb * (p ** i) * (q ** (m - i))
        total_prob += prob_i
    
    return total_prob

print(probability_at_least_n(k, n))

# Alternative using scipy
# from scipy.stats import binom
# print(binom.sf(n - 1, 2**k, 0.25))
```

**R**

```r
data <- readLines('Independent Alleles.txt')
k <- as.integer(strsplit(data[1], " ")[[1]][1])
n <- as.integer(strsplit(data[1], " ")[[1]][2])

m <- 2^k
p <- 0.25

# P(X >= n) = 1 - P(X <= n-1)
prob <- pbinom(n - 1, size = m, prob = p, lower.tail = FALSE)
cat(prob, "\n")
```

---
