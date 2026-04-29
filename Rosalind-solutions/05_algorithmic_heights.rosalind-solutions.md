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

## 3. Degree Array

**Problem**
Given: undirected/simple graph with n vertices in edge list format
Return: an array D[1..n] where D[i] is the degree of vertex i
- in simple graph, degree d(u) of vertex u = # of neighbours u has = # number of edges incident upon it

**Python**

```python

def degree_array(data):
    lines = data.strip().split('\n')
    n, m = map(int, lines[0].split())

    degrees = [0] * (n + 1)

    for line in lines[1:]:
        u, v = map(int, line.split())
        degrees[u] += 1
        degrees[v] += 1
    return " ".join(map(str, degrees[1:]))

file_path = "/Users/billytrim/Desktop/2026/learning/rosalind/Degree Array.txt"
with open(file_path, 'r') as f:
    data = f.read()
print(degree_array(data))

# or w OOP
class Graph:
    def __init__(self, file_path):
        self.n = 0
        self.m = 0
        self.edges = []
        self._load_data(file_path)

    def _load_data(self, file_path):
        """Reads the file and parses the number of vertices/edges."""
        with open(file_path, 'r') as f:
            lines = f.readlines()
            
        # Parse the first line for n and m
        self.n, self.m = map(int, lines[0].split())
        
        # Store edges as a list of tuples
        for line in lines[1:]:
            if line.strip():
                u, v = map(int, line.split())
                self.edges.append((u, v))

    def get_degree_array(self):
        """Calculates the degree for each vertex."""
        # Initialize degree array with zeros (n+1 for 1-based indexing)
        degrees = [0] * (self.n + 1)
        
        for u, v in self.edges:
            degrees[u] += 1
            degrees[v] += 1
            
        # Return as a space-separated string (skipping index 0)
        return " ".join(map(str, degrees[1:]))

path = "/Users/billytrim/Desktop/2026/learning/rosalind/Degree Array.txt"
my_graph = Graph(path)
print(my_graph.get_degree_array())
```
---

## 4. Insertion Sort

**Problem**
Given: A positive integer 'n' (<10^3) and an array A[1..n] of ints
Return: the number of swaps performed by insertion sort algorithm
**Python**
```python
def readin(path):
    with open(path, 'r') as f:
        lines = [line.strip() for line in f.readlines()]
        return lines
l = readin("./bin/rosalind_ins.txt")
a_len = int(l[0])
arry = [int(x) for x in l[1].split()]
# or
    # arry = list(map(int, l[1].split())
print(a_len)
def insertionsort(arr):
    swaps = 0
    for i in range(1, len(arr)):
        key = arr[i] # start @ 2nd ; because indx 0 already sorted
        j = i - 1 # j = indx of item left of key
        # compare key with elements to left
        # move elements > key right 1 position
        while j >= 0 and arr[j] > key:
            # 1) ensure index > 0
            # 2) left item > right item
            arr[j + 1] = arr[j]
            # move j right 1
            j -= 1
            swaps += 1
            # move pointer left
        arr[j + 1] = key
        # loop stops when no smaller number than j
    return swaps, arr
result = insertionsort(arry)
print(result)
```
---

