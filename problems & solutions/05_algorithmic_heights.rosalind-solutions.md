# Rosalind Solutions - _Algorithmic Heights_
### Solutions to bioinformatics problems from [Rosalind](https://rosalind.info/problems/list-view/?location=algorithmic-heights), implemented in both Python and R.

---

## 1. Fibonacci Numbers

**Problem**  
Given: A positive integer n < 26.  
Return: The value of F(n) (the n-th Fibonacci number).

**Python**

```python
with open("/content/Fibonacci Numbers.txt", 'r') as f:
    n = int(f.read().split()[0])

def fib(n):
    if n <= 0:
        return 0
    if n == 1:
        return 1
    a, b = 0, 1
    for i in range(2, n + 1):
        a, b = b, a + b
    return b 

print(fib(n))
```

**R**

```r
n <- as.integer(readLines("/content/Fibonacci Numbers.txt", n = 1))

fib <- function(n) {
    if (n <= 0) return(0)
    if (n == 1) return(1)
    a <- 0
    b <- 1
    for (i in 2:n) {
        temp <- b
        b <- a + b
        a <- temp
    }
    return(b)
}

cat(fib(n), "\n")
```

---

## 2. Binary Search

**Problem**  
Given: Two positive integers n ≤ 10^5 and m ≤ 10^5, a sorted array A[1...n] of integers, and a list of m integers K₁, K₂, ..., Kₘ.  
Return: For each Kᵢ, output an index 1 ≤ j ≤ n such that A[j] = Kᵢ, or -1 if there is no such index.

**Python**

```python
def read(file):
    with open(file, 'r') as f:
        lines = [l.strip() for l in f.readlines()]
    n = int(lines[0])
    m = int(lines[1])
    array = [int(x) for x in lines[2].split()]
    keys = [int(x) for x in lines[3].split()]
    return n, m, array, keys

def binary_search(list, search):
    first = 0
    last = len(list) - 1

    while first <= last:
        middle = (first + last) // 2
        if search == list[middle]:
            return middle + 1  # 1-based index
        elif search < list[middle]:
            last = middle - 1
        else:
            first = middle + 1
    return -1

def main():
    path = "/content/Binary Search Rosalind.txt"
    n, m, array, keys = read(path)
    results = [binary_search(array, k) for k in keys]
    print(*(results))

if __name__ == "__main__":
    main()
```

**R**

```r
read_input <- function(file) {
    lines <- readLines(file)
    n <- as.integer(lines[1])
    m <- as.integer(lines[2])
    array <- as.integer(strsplit(lines[3], " ")[[1]])
    keys <- as.integer(strsplit(lines[4], " ")[[1]])
    return(list(n = n, m = m, array = array, keys = keys))
}

binary_search <- function(arr, target) {
    left <- 1
    right <- length(arr)
    
    while (left <= right) {
        mid <- floor((left + right) / 2)
        if (arr[mid] == target) {
            return(mid)  # 1-based index in R
        } else if (arr[mid] < target) {
            left <- mid + 1
        } else {
            right <- mid - 1
        }
    }
    return(-1)
}

main <- function() {
    path <- "/content/Binary Search Rosalind.txt"
    input <- read_input(path)
    results <- sapply(input$keys, function(k) binary_search(input$array, k))
    cat(paste(results, collapse = " "), "\n")
}

main()
```

---
